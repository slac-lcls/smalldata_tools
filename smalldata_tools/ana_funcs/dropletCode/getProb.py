import numpy as np
import scipy.special as sp

def NB_dist(k,M,kavg,I0=1.):
    temp1 = sp.gammaln(k+M)-sp.gammaln(k+1)-sp.gammaln(M)
    temp2 = -k*np.log(1 + M/kavg)
    temp3 = -M*np.log(1 + kavg/M)
    return I0*np.exp(temp1+temp2+temp3)

def chisqs(p,kavg,M,nroi):
    N = np.size(kavg)
    k = np.reshape((np.arange((4))),(4,1))
    k = np.tile(k,(1,N))
    kavg = np.tile(kavg,(4,1))
    return -2*np.nansum((p*nroi*np.log(1/p*NB_dist(k,M,kavg,1.))))

def getContrast(ps, nroi,low, high, Mmax):
    ps = np.transpose(ps)
    kavg = ps[-1]
    kavg_filter = (kavg>=low)&(kavg<=high)
    kavg = kavg[kavg_filter]
    ps = ps[:,kavg_filter]
    
    nn = (Mmax-1)*1000+1
    Ms = np.linspace(1,Mmax,nn)
    chi2 = np.zeros(Ms.size)
    for ii,M0 in enumerate(Ms):
        chi2[ii] = chisqs(p = ps[:4], kavg = ps[-1],M = M0, nroi = nroi)
    pos = np.argmin(chi2)
    M0 = Ms[pos]
    #curvature as error analysis
    dM = Ms[1] - Ms[0]
    try:
        delta_M = np.sqrt(dM**2/(chi2[pos+1]+chi2[pos-1]-2*chi2[pos]))
    except:
        delta_M = 0.
    return M0, delta_M


def getProb(photonlist, i0s, roi, Np = 12):
    nx,ny = roi.shape
    flag = np.zeros_like(i0s,dtype = "int")
    p = np.zeros((i0s.size,Np))
    #note the x and y in python and in c - photonlist has both & is correctly shaped
    #photon_arr = np.append(ones[:,[0,2,1]],photonlist,axis = 0)
    photon_arr = photonlist
    sequence = np.argsort(photon_arr[:,0]) # this would not sort in tiles?
    photon_arr_sorted = photon_arr[sequence,:]
    #aa is the number of frames with droplets in it
    aa,bb = np.unique(photon_arr_sorted[:,0],return_counts = True)
    aa = np.int32(aa)
    flag[aa] = 1
    # frame i starts from ppos[i]
    ppos = np.append(0,np.cumsum(bb))
    nn = len(aa)
    for j in xrange(nn):
        fr = photon_arr_sorted[ppos[j]:ppos[j+1]]
        ave0,xedges, yedges = np.histogram2d(fr[:,1]+0.5,fr[:,2]+0.5,bins = [np.arange(nx+1),np.arange(ny+1)])
#         sumimgs += ave0
        p[aa[j],:] = np.bincount(np.int32(ave0[roi>0].ravel()),minlength = Np)[:Np]
        p[aa[j],-1] = ave0[roi>0].sum()
    return p/float(roi.sum())

def getProb_img(photonlist, mask, Np=12):
    p = np.zeros((1,Np))
    if mask is None or len(photonlist)==0:
        return p
    #print('prob2 - shape', isinstance(mask, np.ndarray))
    if len(mask.shape)>2:
        mask = mask.sum(axis=0)
    nx,ny = mask.shape
    #note the x and y in python and in c 
#     fr = np.append(ones[:,[0,2,1]],photonlist,axis = 0)
#     ave0,xedges, yedges = np.histogram2d(fr[:,1]+0.5,fr[:,2]+0.5,bins = [np.arange(nx+1),np.arange(ny+1)])
    ave0, xedges, yedges = np.histogram2d(photonlist[:,1]+0.5, photonlist[:,2]+0.5, bins=[np.arange(nx+1),np.arange(ny+1)])
    p[0,:] = np.bincount(np.int32(ave0[mask>0].ravel()),minlength = Np)[:Np]
    p[0,-1] = ave0[mask>0].sum()
    return p/float(mask.sum())
    
