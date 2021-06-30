import psana
import sys
sys.path.append('XPCS_analysis/dropletCode/')
from convert_img import convert_img
from loopdrops import *
from getProb import *
from utilities import printMsg
import resource
import socket
import numpy as np
import scipy.ndimage.measurements as measurements
import skimage.measure as measure
import scipy.ndimage.filters as filters
from scipy import sparse
from smalldata_tools.DetObject import DetObjectFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc

class droplet2Func(DetObjectFunc):
    def __init__(self, **kwargs):
        self._name = kwargs.get('name', 'droplet')
        super(droplet2Func, self).__init__(**kwargs)
        self.threshold = kwargs.get('threshold', None)
        self.mask = kwargs.get('mask', None)
        self.aduspphot = kwargs.get('aduspphot', 0)
        self.offset = kwargs.get('offset', 0)
        self.photpts = np.arange(1000000)*self.aduspphot-self.aduspphot+self.offset
        # self.Np = kwargs.get('Np', None)
        
    def setFromDet(self, det):
        super(droplet2Func, self).setFromDet(det)
        self.mask = det.mask
        
    def process(self, data):
        sum_img = None
        img = data
            
        #make droplets
        ones,ts,pts,h,b = convert_img(img,self.threshold,self.photpts,self.mask)
        #find photons
        photonlist = loopdrops(ones,ts,pts,self.aduspphot,self.photpts)
        if sum_img is None:
            sum_img = img.copy()
            hh = h.copy()
        else:
            sum_img += img
            hh += h.copy()
            
        # look at this
        # p = getProb_img(ones, photonlist, self.mask, self.Np)
        
        # save per-event data
        
        #array of arrays to dictionary 
        d = {}
        photonList = photonlist.tolist()
            
        for same, key, value in photonList:
            if key in d:
                d[key].append(value)
            else:
                d[key] = list()
                d[key].append(value)
                
        return d