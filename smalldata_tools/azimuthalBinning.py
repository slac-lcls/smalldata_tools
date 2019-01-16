import numpy as np
import pylab as plt
import sys
from scipy import hypot,arcsin,arccos
import os
import time
import h5py
from scipy.interpolate import griddata
import utilities as util
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpiSize = comm.Get_size()

g_az = None

def displayimg(img,**kwargs):
    plt.imshow(img.transpose(),origin="lower",**kwargs)

class azimuthalBinning:
    def __init__(self, x, y, xcen, ycen, img=None, **kwargs):
    #new are: ADU/photon, gainImg. darkImg,phiBins
        """ 
        This function azumithally averages images into q & phi bins
        it applies geometrical & (X-ray) polarization corrections
        correctedImage = (Image-darkImg)/gainImg/geom_correction/pol_correction

        Parameters
        ----------
        x,y            = pixel coordinate (1D array each); note: they should be the center of the pixels
        xcen,ycen = center beam position
        tx,ty = angle of detector normal with respect to incoming beam (in deg)
                        zeros are for perpendicular configuration
        darkImg    = darkImage to subtract, by default input is pedestal subtracted.
        ADU_per_photon : used to estimate errors
        qbin = rebinning q (def 0.01)
        phiBins = bin in azimuthal angle (def: one bin)
        Pplane = Polarization (1 = horizontal, 0 = vertical)
        d         = distance of center of detector to sample (in m)
        lam     = wavelength in Ang
        img is used only for displaying corrections
        """
        # save parameters for later use
        self.mask = kwargs.pop("mask",None)
        self.gainImg = kwargs.pop("gainImg",None)
        self.darkImg = kwargs.pop("darkImg",None)
        self.debug = kwargs.pop("debug",False)
        self.ADU_per_photon = kwargs.pop("ADU_per_Photon",1.)
        self.d = kwargs.pop("d",100e-3)
        phiBins = kwargs.pop("phiBins",1)
        lam = kwargs.pop("lam",1.)
        Pplane = kwargs.pop("Pplane",0)
        tx = kwargs.pop("tx",0.)
        ty = kwargs.pop("ty",0.)
        qbin = kwargs.pop("qbin",5e-3)
        self.xcen = xcen
        self.ycen = ycen

        if self.mask is not None: self.mask = np.asarray(self.mask,dtype=np.bool)
        if rank==0:
            if self.mask is not None: 
                print('initialize azimuthal binning, mask %d pixel for azimuthal integration'%self.mask.sum())
            else:
                print('no mask has been passed, will return None')
                return None

        tx = np.deg2rad(tx)
        ty = np.deg2rad(ty)
        self.xcen = float(self.xcen)
        self.ycen = float(self.ycen)
        # equations based on J Chem Phys 113, 9140 (2000) [logbook D30580, pag 71]
        (A,B,C) = (-np.sin(ty)*np.cos(tx),-np.sin(tx),-np.cos(ty)*np.cos(tx))
        (a,b,c) = (self.xcen+self.d*np.tan(ty),float(self.ycen)-self.d*np.tan(tx),self.d)

        r = np.sqrt( (x-a)**2+(y-b)**2+c**2)
        self.r = r
        
        self.msg("calculating theta...",cr=0)
        matrix_theta = np.arccos( (A*(x-a)+B*(y-b)-C*c )/r )
        self.matrix_theta = matrix_theta
        self.msg("...done")
        
        self.msg("calculating phi...",cr=0)
        matrix_phi = np.arccos( ((A**2+C**2)*(y-b)-A*B*(x-a)+B*C*c )/ \
                np.sqrt((A**2+C**2)*(r**2-(A*(x-a)+B*(y-b)-C*c)**2)))
        idx = (y>self.ycen) & (np.isnan(matrix_phi))
        matrix_phi[idx] = 0
        idx = (y<self.ycen) & (np.isnan(matrix_phi))
        matrix_phi[idx] = np.pi
        idx = (x<self.xcen)
        matrix_phi[idx] = (np.pi-matrix_phi[idx])+np.pi
#        matrix_phi[idx] = temp+n.pi
        self.matrix_phi = matrix_phi
        self.msg("...done")

        self.msg("calculating pol matrix...",cr=0)
        Pout = 1-Pplane
        pol = Pout*(1-(np.sin(matrix_phi)*np.sin(matrix_theta))**2)+\
                Pplane*(1-(np.cos(matrix_phi)*np.sin(matrix_theta))**2)

        self.msg("... done")
        self.pol=pol
        theta_max = np.nanmax(matrix_theta[~self.mask])

        self.msg("calculating digitize")
        if isinstance(phiBins, list):
            if max(phiBins)<(2*np.pi-0.01):
                phiBins.append(2*np.pi)
            if min(phiBins)>0:
                phiBins.append(0)
            phiBins.sort()
            self.nphi = len(phiBins)
            pbm = self.matrix_phi + (phiBins[1]-phiBins[0])/2
            pbm[pbm>=2*np.pi] -= 2*np.pi
            self.phiVec = np.array(phiBins)
        else:
            self.nphi = phiBins
            phiint = 2*np.pi/phiBins
            pbm = self.matrix_phi + phiint/2
            pbm[pbm>=2*np.pi] -= 2*np.pi
            self.phiVec = np.linspace(0,2*np.pi+np.spacing(np.min(pbm)),phiBins+1)

        self.pbm = pbm #added for debugging of epix10k artifacts.
        self.idxphi = np.digitize(pbm.ravel(),self.phiVec)-1
        self.matrix_q = 4*np.pi/lam*np.sin(self.matrix_theta/2)
        q_max = np.nanmax(self.matrix_q[~self.mask])
        q_min = np.nanmin(self.matrix_q[~self.mask])
        qbin = np.array(qbin)
        if qbin.size==1:
            if rank==0:
                print('qmax: ',q_max,' qbin ',qbin)
            #self.qbins = np.arange(0,q_max+qbin,qbin)
            self.qbins = np.arange(q_min-qbin,q_max+qbin,qbin)
        else:
            self.qbins = qbin
        self.q = (self.qbins[0:-1]+self.qbins[1:])/2
        self.theta = 2*np.arcsin(self.q*lam/4/np.pi)
        self.nq = self.q.size
        self.idxq    = np.digitize(self.matrix_q.ravel(),self.qbins)-1
        last_idx = self.idxq.max()
        self.idxq[self.mask.ravel()] = 0; # send the masked ones in the first bin

        # 2D binning!
        self.Cake_idxs = np.ravel_multi_index((self.idxphi,self.idxq),(self.nphi,self.nq))
        self.Cake_idxs[self.mask.ravel()] = 0; # send the masked ones in the first bin
        
        #print("last index",last_idx)
        self.msg("...done")
        # include geometrical corrections
        geom    = (self.d/r) ; # pixels are not perpendicular to scattered beam
        geom *= (self.d/r**2); # scattered radiation is proportional to 1/r^2
        self.msg("calculating normalization...",cr=0)
        self.geom = geom
        self.geom /= self.geom.max()
        self.correction = self.geom*self.pol
        self.Npixel = np.bincount(self.idxq,minlength=self.nq); self.Npixel = self.Npixel[:self.nq]
        self.norm     = self.Npixel
        self.Cake_Npixel = np.bincount(self.Cake_idxs,minlength=self.nq*self.nphi)
        #self.Cake_Npixel = self.Npixel[:self.nq*self.nphi]
        self.Cake_norm=np.reshape(self.Cake_Npixel,(self.nphi,self.nq));#/self.correction1D
        #self.correction1D    =self.correction1D[:self.nq]/self.Npixel
        self.header    = "# Parameters for data reduction\n"
        self.header += "# xcen, ycen = %.2f m %.2f m\n" % (self.xcen,self.ycen)
        self.header += "# sample det distance = %.4f m\n" % (self.d)
        self.header += "# wavelength = %.4f Ang\n" % (lam)
        self.header += "# detector angles x,y = %.3f,%.3f deg\n" % (np.rad2deg(tx),np.rad2deg(ty))
        self.header += "# fraction of inplane pol %.3f\n" % (Pplane)
        if isinstance(qbin,float):
            self.header += "# q binning : %.3f Ang-1\n" % (qbin)
        return 

    def msg(self,s,cr=True):
        if (self.debug):
            if (cr):
                print(s)
            else:
                print(s,)
        sys.stdout.flush()

    def doAzimuthalAveraging(self,img,applyCorrection=True):
        if self.darkImg is not None: img-=self.darkImg
        if self.gainImg is not None: img/=self.gainImg
        if applyCorrection:
            I=np.bincount(self.idxq, weights = img.ravel()/self.correction.ravel(), minlength=self.nq); I=I[:self.nq]
        else:
            I=np.bincount(self.idxq, weights = img.ravel()                        , minlength=self.nq); I=I[:self.nq]
        self.sig = np.sqrt(1./self.ADU_per_photon)*np.sqrt(I)/self.norm
        self.I = I/self.norm
        return self.I


    def doCake(self,img,applyCorrection=True):
        if self.darkImg is not None: img-=self.darkImg
        if self.gainImg is not None: img/=self.gainImg
        if applyCorrection:
            I=np.bincount(self.Cake_idxs, weights = img.ravel()/self.correction.ravel(), minlength=self.nq*self.nphi); I=I[:self.nq*self.nphi]
        else:
            I=np.bincount(self.Cake_idxs, weights = img.ravel()                        , minlength=self.nq*self.nphi); I=I[:self.nq*self.nphi]
        I = np.reshape(I,(self.nphi,self.nq))
        self.sig = 1./np.sqrt(self.ADU_per_photon)*np.sqrt(I)/self.Cake_norm    # ??? where comes this sqrt from? Ah I see...
        self.Icake = I/self.Cake_norm
        return self.Icake

    def saveChiFile(self,fname):
        header = "q(Ang-1) I sig"
        #mc.writev(fname,self.q,n.vstack((self.I,self.sig)),header=header)
        #n.savetxt(fname,n.vstack((self.q,self.I,self.sig)),header=header)
        np.savetxt(fname,np.transpose(np.vstack((self.q,self.I,self.sig))),
            fmt=["%.3f","%.4f","%.4f"])


def test():
    mask=np.ones( (2000,2000) )
    az=azimuthal_averaging(mask,-80,1161,pixelsize=82e-6,d=4.7e-2,tx=0,ty=90-28.,thetabin=1e-1,lam=1,debug=1)
    plt.subplot("121")
    displayimg(np.rad2deg(az.matrix_theta))
    print(az.matrix_theta.min())
    plt.clim(0,180)
    plt.colorbar()
    plt.subplot("122")
    displayimg(np.rad2deg(az.matrix_phi))
    print(az.matrix_phi.min(),az.matrix_phi.max())
    plt.colorbar()
    plt.clim(0,360)
