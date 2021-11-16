from pathlib import Path
import h5py as h5
import glob
import os
import numpy as np
import sys
#sys.path.append('/reg/d/psdm/xcs/xcsx41018/results/smalldata_tools/')
#sys.path.append('/reg/d/psdm/xcs/xcsx41018/results/vesp/')

from smalldata_tools.ana_funcs.correlations import correlation as corr
from smalldata_tools.ana_funcs.correlations import utils
from smalldata_tools.DetObject import DetObjectFunc


class Autocorrelation(DetObjectFunc):
    """ 
    """
    def __init__(self, **kwargs):
        """ 
        Args:
            name (str): DetObjectName, default: autocorr
            thresh (list or tuple): low and high pixel intensity tresholds [low, high]
            roi (list or array): [roi0, roi1, roi3, roi4] rectangular ROI coordinates
            mask (str or Path object): path to npy file containing mask
        """
        self._name = kwargs.get('name','autocorr')
        super(Autocorrelation, self).__init__(**kwargs)
        self.thresholds = kwargs.get('thresADU', [-1e6, 1e6])
        self.save_range = kwargs.get('save_range', None)
        self.save_lineout = kwargs.get('save_lineout', False)
        self.correct_illumination = kwargs.get('correct_illumination', False) # not implemented
        self.roi = kwargs.get('roi', None)
        if 'mask' in kwargs:
            self.mask = np.load(kwargs['mask'])
        else:
            self.mask = None
        if self.mask is not None:
            print('A mask is given, ROI set to None')
            self.roi = None
        return
    
            
    def setFromDet(self, det):
        """ """
        super(Autocorrelation, self).setFromDet(det)
        return
    
    
    def process(self, img):
        """
        Perform autocorrelation on masked detector images
        """
#         if self.roi is not None:
#             img = img[self.roi[0][0]:self.roi[0][1], self.roi[1][0]:self.roi[1][1]]
        img[img<self.thresholds[0]] = 0
        img[img>self.thresholds[1]] = 0
        
        if self.mask is not None:
            if self.mask.ndim==3:
                autocorr = [corr.spatial_correlation_fourier(img, mask=mask) for mask in self.mask]
                autocorr = np.asarray(autocorr)
                cr = [img[mask].mean() for mask in self.mask]
                cr = np.asarray(cr)
            else:
                autocorr = corr.spatial_correlation_fourier(img, mask=self.mask)
                cr = img[self.mask].mean()
        else:
            autocorr = corr.spatial_correlation_fourier(img, mask=self.mask)
            cr = img.mean()
        
        if self.save_range is not None:
            cx,cy = utils.get_center(autocorr)
            rr = self.save_range
            autocorr = autocorr[cx-rr[0]:cx+rr[0], cy-rr[1]:cy+rr[1]]
        
        if self.save_lineout:
            cv,ch = utils.get_center(autocorr)
            lineout_h = autocorr[cv,:]
            lineout_v = autocorr[:,ch]
            output = {
                'lineout_h': lineout_h,
                'lineout_v': lineout_v,
                'cr': cr
            }
        else:
            output = {
#                 'img': img,
                'autocorr': autocorr,
                'cr': cr,
            }
        return output
