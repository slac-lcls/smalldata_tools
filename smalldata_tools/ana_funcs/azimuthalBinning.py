import numpy as np
import sys
from scipy import hypot,arcsin,arccos
import time
import h5py
from scipy.interpolate import griddata
import smalldata_tools.utilities as util
from smalldata_tools.DetObjectFunc import DetObjectFunc
from mpi4py import MPI
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpiSize = comm.Get_size()

class azimuthalBinning(DetObjectFunc):
    def __init__(self, **kwargs):
    #new are: ADU/photon, gainImg. darkImg,phiBins
        """ 
        This function azumithally averages images into q & phi bins
        it applies geometrical & (X-ray) polarization corrections
        correctedImage = (Image-darkImg)/gainImg/geom_correction/pol_correction

        Parameters
        ----------
        center (x,y)   = pixel coordinate (1D array each); note: they should be the center of the pixels
        xcen,ycen = center beam position (in micron as derived from psana-detector geometry)
        tx,ty = angle of detector normal with respect to incoming beam (in deg)
                        zeros are for perpendicular configuration
        darkImg    = darkImage to subtract, by default input is pedestal subtracted.
        ADU_per_photon : used to estimate errors
        qbin = rebinning q (def 0.01)
        phiBins = bin in azimuthal angle (def: one bin)
        Pplane = Polarization (1 = horizontal, 0 = vertical)
        dis_to_sam = distance of center of detector to sample (in mm)
        lam = wavelength in Ang
        userMask = userMask as array (detector data shaped)
        thresADU = lower threshold in ADU
        thresADUhigh = high threshold in ADU
        thresRms = lower threshold in RMS
        geomCorr: apply geometry correction (def True)
        polCorr: apply polarization correction (def True)
        """
        # save parameters for later use
        self._name = kwargs.get('name','azav')
        super(azimuthalBinning, self).__init__(**kwargs)
        self._mask = kwargs.pop("userMask",None)
        self.gainImg = kwargs.pop("gainImg",None)
        self.darkImg = kwargs.pop("darkImg",None)
        self._debug = kwargs.pop("debug",False)
        self.ADU_per_photon = kwargs.pop("ADU_per_Photon",1.)
        self.dis_to_sam = kwargs.pop("dis_to_sam",100e-3)
        self.eBeam =  kwargs.pop("eBeam",9.5)
        self.lam = util.E2lam(self.eBeam)*1e10
        self.phiBins = kwargs.pop("phiBins",1) 
        self.Pplane = kwargs.pop("Pplane",0) #0/1
        self.tx = kwargs.pop("tx",0.)
        self.ty = kwargs.pop("ty",0.)
        self.qbin = kwargs.pop("qbin",5e-3)
        self.rbin = kwargs.pop("rbin", None) #in micron
        self.thresRms =  kwargs.pop("thresRms",None)
        self.thresADU =  kwargs.pop("thresADU",None)
        self.thresADUhigh = kwargs.pop("thresADUhigh",None)
        self.x = kwargs.pop("x",None)
        self.y = kwargs.pop("y",None)
        self.geomCorr = kwargs.pop("geomCorr",True)
        self.polCorr = kwargs.pop("polCorr",True)
        self.square = kwargs.pop("square",False)

        center = kwargs.pop("center",None)
        if center is not None:
            self.xcen = center[0]/1e3
            self.ycen = center[1]/1e3 
        else:
            print('no center has been given, will return None')
            return None
            
        if self._mask is not None: self._mask = np.asarray(self._mask,dtype=np.bool)

    def setFromDet(self, det):
        if det.mask is not None and det.cmask is not None:
            if self._mask is not None and self._mask.flatten().shape == det.mask.flatten().shape:
                self._mask = ~(self._mask.flatten().astype(bool)&det.mask.astype(bool).flatten())
            else:
                self._mask = ~(det.cmask.astype(bool)&det.mask.astype(bool))
        self._mask = self._mask.flatten()
        #if self._mask is None and det.mask is not None:
        #    setattr(self, '_mask', det.mask.astype(np.uint8))
        if det.x is not None:
            self.x = det.x.flatten()/1e3
            self.y = det.y.flatten()/1e3
            self.z = det.z.flatten()/1e3
            self.z_off = np.nanmean(self.z)-self.z
            #z_off defined such that z_off >0 is downstream

    def setFromFunc(self, func=None):
        super(azimuthalBinning,self).setFromFunc()
        if func is None:
            self._setup()
            return
        if getattr(func,'_x', None) is not None: self.x = func._x.flatten()/1e3
        if getattr(func,'_y', None) is not None: self.y = func._y.flatten()/1e3
        maskattr='_mask'
        if not hasattr(func, maskattr):
            maskattr='mask'
            if not hasattr(func, maskattr):
                maskattr=None
        if maskattr is not None and getattr(func,maskattr) is not None:
            if isinstance(func, ROIFunc): 
                self._mask = getattr(func, maskattr).astype(bool).flatten()
            else:
                self._mask = (~(getattr(func,maskattr).astype(bool))).flatten()
        #elif func._rms is not None: 
        #    self._mask = np.ones_like(func._rms).flatten()
        self._setup()

    def _setup(self):

        if rank==0:
            if self._mask is not None: 
                print('initialize azimuthal binning, mask %d pixel for azimuthal integration'%self._mask.sum())
            else:
                print('no mask has been passed, will return None')
                return None
            if self.x is None:
                print('no x/y array have been passed, will return None')
                return None
                
        tx = np.deg2rad(self.tx)
        ty = np.deg2rad(self.ty)
        self.xcen = float(self.xcen)
        self.ycen = float(self.ycen)

        # equations based on J Chem Phys 113, 9140 (2000) [logbook D30580, pag 71]
        (A,B,C) = (-np.sin(ty)*np.cos(tx),-np.sin(tx),-np.cos(ty)*np.cos(tx))
        (a,b,c) = (self.xcen+(self.dis_to_sam+self.z_off)*np.tan(ty),float(self.ycen)-(self.dis_to_sam+self.z_off)*np.tan(tx),(self.dis_to_sam+self.z_off))

        x = self.x
        y = self.y
        r = np.sqrt( (x-a)**2+(y-b)**2+c**2)
        self.r = r
        
        self.msg("calculating theta...",cr=0)
        matrix_theta = np.arccos( (A*(x-a)+B*(y-b)-C*c )/r )
        self.matrix_theta = matrix_theta
        self.msg("...done")

        if self._debug:
            print('matrix theta: ',self.matrix_theta.shape)
        self.msg("calculating phi...",cr=0)
        matrix_phi = np.arccos( ((A**2+C**2)*(y-b)-A*B*(x-a)+B*C*c )/ \
                np.sqrt((A**2+C**2)*(r**2-(A*(x-a)+B*(y-b)-C*c)**2)))
        idx = (y>=self.ycen) & (np.isnan(matrix_phi))
        matrix_phi[idx] = 0
        idx = (y<self.ycen) & (np.isnan(matrix_phi))
        matrix_phi[idx] = np.pi
        idx = (x<self.xcen)
        matrix_phi[idx] = (np.pi-matrix_phi[idx])+np.pi
#        matrix_phi[idx] = temp+n.pi
        self.matrix_phi = matrix_phi
        self.msg("...done")

        self.msg("calculating pol matrix...",cr=0)
        Pout = 1-self.Pplane
        pol = Pout*(1-(np.sin(matrix_phi)*np.sin(matrix_theta))**2)+\
                self.Pplane*(1-(np.cos(matrix_phi)*np.sin(matrix_theta))**2)

        self.msg("... done")
        self.pol=pol
        theta_max = np.nanmax(matrix_theta[~self._mask])

        self.msg("calculating digitize")
        if isinstance(self.phiBins, np.ndarray):
            self.phiBins = self.phiBins.tolist()
        if isinstance(self.phiBins, list):
            if max(self.phiBins)<(2*np.pi-0.01):
                #self.phiBins.append(2*np.pi)
                self.phiBins.append(np.array(self.phiBins).max()+0.001)
            if min(self.phiBins)>0:
                #self.phiBins.append(0)
                self.phiBins.append(np.array(self.phiBins).min()-0.001)
            self.phiBins.sort()
            self.nphi = len(self.phiBins)
            pbm = self.matrix_phi + (self.phiBins[1]-self.phiBins[0])/2
            pbm[pbm>=2*np.pi] -= 2*np.pi
            self.phiVec = np.array(self.phiBins)
        else:
            self.nphi = self.phiBins
            #phiint = 2*np.pi/self.phiBins
            phiint = (self.matrix_phi.max()-self.matrix_phi.min())/self.phiBins
            pbm = self.matrix_phi + phiint/2
            pbm[pbm>=2*np.pi] -= 2*np.pi
            #self.phiVec = np.linspace(0,2*np.pi+np.spacing(np.min(pbm)),self.phiBins+1)
            self.phiVec = np.linspace(self.matrix_phi.min(),self.matrix_phi.max()+np.spacing(np.min(pbm)),self.phiBins+1)
            #self.phiVec = np.linspace(0,2*np.pi+np.spacing(np.min(pbm)),self.phiBins+1)

        self.pbm = pbm #added for debugging of epix10k artifacts.
        self.idxphi = np.digitize(pbm.ravel(),self.phiVec)-1
        if self.idxphi.min()<0:
            print('pixels will underflow, will put all pixels beyond range into first bin in phi')
            self.idxphi[self.idxphi<0]=0 #put all 'underflow' bins in first bin.
        if self.idxphi.max()>=self.nphi:
            print('pixels will overflow, will put all pixels beyond range into first bin in phi')
            self.idxphi[self.idxphi==self.nphi]=0 #put all 'overflow' bins in first bin.

        #print('DEBUG phi ',self.phiVec)
        #print('DEBUG phi ',np.unique(self.idxphi))
        
        # include geometrical corrections
        geom  = ((self.dis_to_sam+self.z_off)/r) ; # pixels are not perpendicular to scattered beam
        geom *= ((self.dis_to_sam+self.z_off)/r**2); # scattered radiation is proportional to 1/r^2
        self.msg("calculating normalization...",cr=0)
        self.geom = geom
        self.geom /= self.geom.max()
        if not self.geomCorr:
            self.geom = np.ones_like(self.geom).astype(float)
        if not self.polCorr:
            self.pol = np.ones_like(self.pol).astype(float)
        self.correction = self.geom*self.pol

        #coordinates in q
        self.matrix_q = 4*np.pi/self.lam*np.sin(self.matrix_theta/2)
        q_max = np.nanmax(self.matrix_q[~self._mask])
        q_min = np.nanmin(self.matrix_q[~self._mask])
        qbin = np.array(self.qbin)
        if qbin.size==1:
            if rank==0 and self._debug:
                print('q-bin size has been given: qmax: ',q_max,' qbin ',qbin)
            #self.qbins = np.arange(0,q_max+qbin,qbin)
            self.qbins = np.arange(q_min-qbin,q_max+qbin,qbin)
        else:
            self.qbins = qbin
        self.q = (self.qbins[0:-1]+self.qbins[1:])/2
        self.nq = self.q.size
        self.idxq    = np.digitize(self.matrix_q.ravel(),self.qbins)-1
        self.idxq[self._mask.ravel()] = 0; # send the masked ones in the first bin

        self.theta = 2*np.arcsin(self.q*self.lam/4/np.pi)
        # 2D binning!
        self.Cake_idxs = np.ravel_multi_index((self.idxphi,self.idxq),(self.nphi,self.nq))
        self.Cake_idxs[self._mask.ravel()] = 0; # send the masked ones in the first bin
        self.Cake_Npixel = np.bincount(self.Cake_idxs,minlength=self.nq*self.nphi)
        #self.Cake_Npixel = self.Npixel[:self.nq*self.nphi]
        self.Cake_norm=np.reshape(self.Cake_Npixel,(self.nphi,self.nq));#/self.correction1D
        #self.correction1D    =self.correction1D[:self.nq]/self.Npixel

        #coordinates in r
        if self.rbin is not None:
            x = self.x
            y = self.y
            rl = np.sqrt( (x-self.xcen)**2+(y-self.ycen)**2 )
            self.rlocal = rl
            #r_max = np.nanmax(self.rlocal)
            #r_min = np.nanmin(self.rlocal)
            r_max = np.nanmax(self.rlocal[~self._mask])
            r_min = np.nanmin(self.rlocal[~self._mask])
            rbin = np.array(self.rbin)
            if rbin.size==1:
                if rank==0 and self._debug:
                    print('q-bin size has been given: rmax: ',r_max,' rbin ',rbin)
                #self.rbins = np.arange(0,q_max+rbin,rbin)
                self.rbins = np.arange(r_min-rbin,r_max+rbin,rbin)
            else:
                self.rbins = rbin
            self.rbinsbound = (self.rbins[0:-1]+self.rbins[1:])/2
            #print('self.rbinsbound ',self.rbinsbound)
            self.nr = self.rbinsbound.size
            #print('nr ',self.nr)
            #here, the bin-range is > thew actual range unlike for q where the bined range is < than the actual range...
            if ( (np.nanmax(self.rlocal)-np.nanmin(self.rlocal)) < (self.rbins.max()-self.rbins.min()) ):
                self.nr+=1
            self.idxr = np.digitize(self.rlocal.ravel(),self.rbins)-1
            self.idxr[self._mask.ravel()] = 0; # send the masked ones in the first bin
            #print('2 r', self.idxr.min(), self.idxr.max(), self.idxphi.min(), self.idxphi.max())
            #print('2 ns' , self.nphi, self.nr, self.rbins.size, self.rbinsbound.size)
            self.Cake_idxs = np.ravel_multi_index((self.idxphi,self.idxr),(self.nphi,self.nr))
            #print('3', self.Cake_idxs.shape, self.Cake_idxs.max())
            self.Cake_idxs[self._mask.ravel()] = 0; # send the masked ones in the first bin
            self.Cake_Npixel = np.bincount(self.Cake_idxs,minlength=self.nr*self.nphi)
            #self.Cake_Npixel = self.Npixel[:self.nq*self.nphi]
            self.Cake_norm=np.reshape(self.Cake_Npixel,(self.nphi,self.nr));#/self.correction1D
            #print('nrend ',self.nr, self.Cake_idxs.max())

        #last_idx = self.idxq.max()
        #print("last index",last_idx)
        self.msg("...done")


        self.header    = "# Parameters for data reduction\n"
        self.header += "# xcen, ycen = %.2f m %.2f m\n" % (self.xcen,self.ycen)
        self.header += "# sample det distance = %.4f m\n" % (self.dis_to_sam)
        self.header += "# wavelength = %.4f Ang\n" % (self.lam)
        self.header += "# detector angles x,y = %.3f,%.3f deg\n" % (np.rad2deg(tx),np.rad2deg(ty))
        self.header += "# fraction of inplane pol %.3f\n" % (self.Pplane)
        if isinstance(qbin,float):
            self.header += "# q binning : %.3f Ang-1\n" % (qbin)
        #remove idx & correction values for masked pixels. Also remove maks pixels from image in process fct
        self.Cake_idxs = self.Cake_idxs[self._mask.ravel()==0]
        self.correction = self.correction.flatten()[self._mask.ravel()==0]
        #print('return ', self.Cake_idxs.shape, self.Cake_idxs.max())
        return 

    def msg(self,s,cr=True):
        if (self._debug):
            if (cr):
                print(s)
            else:
                print(s,)
        sys.stdout.flush()

    def doCake(self,img,applyCorrection=True):
        if self.darkImg is not None: img-=self.darkImg
        if self.gainImg is not None: img/=self.gainImg

        img = img.ravel()[self._mask.ravel()==0]
        #print('img:', img.shape, self.Cake_idxs.max())

        nradial=self.nq
        if self.rbin is not None:
            nradial=self.nr

        if applyCorrection:
            #I=np.bincount(self.Cake_idxs, weights = img.ravel()/self.correction.ravel(), minlength=self.nq*self.nphi); I=I[:self.nq*self.nphi]
            I=np.bincount(self.Cake_idxs, weights = img/self.correction.ravel(), minlength=nradial*self.nphi); I=I[:nradial*self.nphi]
        else:
            #I=np.bincount(self.Cake_idxs, weights = img.ravel()                        , minlength=self.nq*self.nphi); I=I[:self.nq*self.nphi]
            I=np.bincount(self.Cake_idxs, weights = img, minlength=self.nq*self.nphi); I=I[:self.nq*self.nphi]
        I = np.reshape(I,(self.nphi,self.nq))
        #don't calculate this - I don't think this is stored.
        #self.sig = 1./np.sqrt(self.ADU_per_photon)*np.sqrt(I)/self.Cake_norm    # ??? where comes this sqrt from? Ah I see...
        self.Icake = I/self.Cake_norm
        return self.Icake

    def process(self, data):
        data=data.copy()
        if self.thresADU is not None:
            data[data<self.thresADU]=0.
        if self.thresADUhigh is not None:
            data[data>self.thresADUhigh]=0.
        if self.thresRms is not None:
            data[data>self.thresRms*self.rms]=0.
        if self.square:
            data=data*data
        return {'azav': self.doCake(data)}
    
        
#make this a real test class w/ assertions.
def test():
    mask=np.ones( (2000,2000) )
    az=azimuthal_averaging(mask,-80,1161,pixelsize=82e-6,d=4.7e-2,tx=0,ty=90-28.,thetabin=1e-1,lam=1,debug=1)
    print(az.matrix_theta.min())
    print(az.matrix_phi.min(),az.matrix_phi.max())
