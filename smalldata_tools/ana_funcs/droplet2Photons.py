import numpy as np
from smalldata_tools.ana_funcs.dropletCode.convert_img import convert_img
from smalldata_tools.ana_funcs.dropletCode.loopdrops import *
from smalldata_tools.ana_funcs.dropletCode.getProb import *
import scipy.ndimage.measurements as measurements
import skimage.measure as measure
import scipy.ndimage.filters as filters
from scipy import sparse
import time
from smalldata_tools.DetObjectFunc import DetObjectFunc

class droplet2Photons(DetObjectFunc):
    '''
    return_img: whether or not to return the img with how many photons at each coordinate
    threshold: # (noise sigma for) threshold
    mask: pass a mask in here (array-form), is None: use mask stored in DetObject
    aduspphot: 
    offset: 
    
    uses convert_img to make droplets to analyze 
    uses loopdrops to find the photons in the droplets (don't forget to append the ones)
    
    counts the number of photons at each (rounded) coordinate
    returns either photonlist or img depending on return_img
    '''
    def __init__(self, **kwargs):
        self.return_img = kwargs.get('return_img',False)
        self._name =  kwargs.get('name','droplet')
        super(droplet2Photons, self).__init__(**kwargs)
        self.threshold = kwargs.get('threshold', None)
        self.mask = kwargs.get('mask', None)
        self.aduspphot = kwargs.get('aduspphot', 0)
        self.offset = kwargs.get('offset', 0)
        self.photpts = np.arange(1000000)*self.aduspphot-self.aduspphot+self.offset
        # self.Np = kwargs.get('Np', None)
        self.cputime = kwargs.get('cputime', False)
        
    def setFromDet(self, det):
        super(droplet2Photons, self).setFromDet(det)
        if self.mask is None:
            self.mask = det.mask
        else:
            self.mask = np.logical_and(self.mask, det.mask)
        
    def process(self, data):
        time0 = time.time()
        sum_img = None
        img = data
        
        #make droplets
        ones,ts,pts,h,b = convert_img(img,self.threshold,self.photpts,self.mask)
        timed = time.time()
        #find photons
        photonlist = loopdrops(ones,ts,pts,self.aduspphot,self.photpts)
        photonlist = np.append(ones[:,[0,2,1]], photonlist, axis=0) # indexes are inverted for ones because of c vs python indexing
        if sum_img is None:
            sum_img = img.copy()
            hh = h.copy()
        else:
            sum_img += img
            hh += h.copy()
            
        nx, ny = img.shape

        timep = time.time()        
        phot_img, xedges, yedges = np.histogram2d(photonlist[:,1]+0.5, photonlist[:,2]+0.5, bins=[np.arange(nx+1),np.arange(ny+1)])
        
        # look at this
        p = getProb_img(photonlist, self.mask, 12)
        
        # output dictionary
        ret_dict = {'prob': np.squeeze(p)}
        if self.cputime: ret_dict['cputime'] = np.array([time.time()-time0, timed-time0, timep-timed, time.time()-timep]) 
      
        self.dat = np.ma.masked_array(phot_img, mask=(~(self.mask.astype(bool))).astype(np.uint8))
        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]

        return ret_dict
