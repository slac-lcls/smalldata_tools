from pathlib import Path
import h5py as h5
import glob
import os
import numpy as np
import sys
import time

from smalldata_tools.ana_funcs.correlations import correlation as corr
from smalldata_tools.ana_funcs.correlations import utils
from smalldata_tools.common.detector_base import DetObjectFunc


class Autocorrelation(DetObjectFunc):
    """ 
    """
    def __init__(self, **kwargs):
        """ 
        Parameters
        ----------
            name (str):
                DetObjectName, default: autocorr
            thresh (list or tuple):
                low and high pixel intensity tresholds [low, high]
            roi (list or array):
                [roi0, roi1, roi3, roi4] rectangular ROI coordinates
            mask (str or Path object):
                path to an npy file containing mask(s).
                The mask can be used to define non-rectangular region of interest. 
                Several mask can be passed by as the first dimension of the array.
                If several masks are passed, save_range must be set.
            save_lineout (bool):
                Save autocorr image or only vertical and horizontal lineouts.
                Default: False
            save_range (2-tuple):
                Size (square) of the autocorr image to save. Default: (50, 50)
            illumination_correction (dict):
                Correction arrays for each mask.
                Dict items:
                    'correction': path to correction arrays. One array per mask /
                                  ROI.
                    'kernel': kernel size used in the creation of the illumination
                              correction
                Default: None
        """
        self._name = kwargs.get('name','autocorr')
        super(Autocorrelation, self).__init__(**kwargs)
        self.thresholds = kwargs.get('thresADU', [-1e6, 1e6])
        self.save_range = kwargs.get('save_range', (50, 50))
        self.save_lineout = kwargs.get('save_lineout', False)
        self.roi = kwargs.get('roi', None)
        if 'mask' in kwargs:
            self.mask = np.load(kwargs['mask']).astype(bool)
            if self.mask.ndim == 2:
                self.mask = np.asarray([self.mask])
        else:
            self.mask = None
        if self.mask is not None:
            print('A mask is given, ROI set to None')
            self.roi = None

        if 'illumination_correction' in kwargs:
            self.illumination_correction = np.load(
                kwargs['illumination_correction']['correction'],
                allow_pickle=True
            )
            self.illumination_kernel = kwargs['illumination_correction']['kernel']
        else:
            self.illumination_correction = None
        
        # check save_range vs data shape to avoid wrapping issues
        self.check_mask_save_range() 
        return
    
            
    def setFromDet(self, det):
        """ """
        super(Autocorrelation, self).setFromDet(det)
        return


    def check_mask_save_range(self):
        min_size = 1e6
        for ii, mask in enumerate(self.mask):
            _, c_mask = utils.box_to_roi(mask, mask)
            print(f"Cropped mask shape: {c_mask.shape}")
            if np.min(c_mask.shape) < min_size:
                min_size = np.min(c_mask.shape)
        if min_size//2 < np.max(self.save_range):
            print(f"save_range is bigger than the autocorr array. Reducing it to ({min_size//2-1}, {min_size//2-1})")
            self.save_range = (min_size//2-1, min_size//2-1)
        return

    
    def process(self, img):
        """
        Perform autocorrelation on masked detector images
        """
        if self.roi is not None:
            img = img[self.roi[0][0]:self.roi[0][1], self.roi[1][0]:self.roi[1][1]]
        img[img<self.thresholds[0]] = 0
        img[img>self.thresholds[1]] = 0
        
        if self.mask is not None:
            autocorr = []
            cr = []
            for ii, mask in enumerate(self.mask):
                if self.illumination_correction is not None:
                    img_box, mask = utils.box_to_roi_extend(
                        img,
                        mask,
                        extend = 2*self.illumination_kernel
                    )
                    img_corr = img_box / self.illumination_correction[ii]
                    autocorr.append(corr.spatial_correlation_fourier(img_box, mask=mask))
                    cr.append(img_box[mask].mean()) # careful if mask is not bool, this does weird things
                else:
                    autocorr.append(corr.spatial_correlation_fourier(img, mask=mask))
                    cr.append(img[mask].mean()) # careful if mask is not bool, this does weird things
            cr = np.asarray(cr)
        else:
            autocorr = corr.spatial_correlation_fourier(img, mask=self.mask)
            cr = img.mean()

        if self.save_range is not None:
            if not isinstance(autocorr, list):
                autocorr = [autocorr]
            for ii, ac in enumerate(autocorr):
                cx,cy = utils.get_center(ac)
                rr = self.save_range
                autocorr[ii] = ac[cx-rr[0]:cx+rr[0], cy-rr[1]:cy+rr[1]]
        autocorr = np.squeeze(np.asarray(autocorr))
        
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
                'autocorr': autocorr,
                'cr': cr,
            }
        return output
