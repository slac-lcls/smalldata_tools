import numpy as np
import time
import scipy.ndimage as ndi
import skimage.measure as measure
import scipy.ndimage.filters as filters
from scipy import sparse
from smalldata_tools.DetObjectFunc import DetObjectFunc

class dropletFunc(DetObjectFunc):
    """ 
    Parameters
    ----------
    threshold : float (default = 5)
         Treshold for pixel to be part of a droplet in sigma or ADU, depending on the 
         value of the useRms parameters.
    thresholdLow : float (default = same as threshold)
        Lower threshold: this is to make the spectrum sharper, but not find peaks 
        out of pixels with low significance.
    mask: np.ndarray (default = None)
        Pass a mask in here, is None: use mask stored in DetObject
    name: str (default: 'droplet') 
        Name used in hdf5 for data field
    thresADU: float (default = None)
        Threshold on droplets' ADU (sum of all pixels in a droplet) for droplet 
        to be further processed. Rejects droplets that are considered too low, i.e.
        that do not contain enough intensity for a single photon.
    useRms (def True): 
        If True, threshold and thresholdLow are # of rms of data, otherwise ADU are used.
    relabel (def True): 
        After initial droplet finding and allowing pixels above the lower threshold, 
        relabel the image (so that droplets merge). This allows pixel below the first 
        threshold that are neigboring to existing droplet to be accounted for and perhaps 
        round the intensity to a full photon ADU for example.

    By default, only total number of droplets is returned by process(data)
    
    Many more information about the droplets can be saved there.
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name', 'droplet')
        super(dropletFunc, self).__init__(**kwargs)
        self.threshold = kwargs.get('threshold', 5.)
        self.thresholdLow = kwargs.get('thresholdLow', self.threshold)
        self.thresADU = kwargs.get('thresADU', None)
        self.useRms = kwargs.get('useRms', True)
        self.relabel = kwargs.get('relabel', True)
        self._mask = kwargs.get('mask', None)
        self._debug = False
        self.footprint = np.array([
            [0,1,0],
            [1,1,1],
            [0,1,0]
        ])
        self._footprint2d = np.array([
            [0,1,0],
            [1,1,1],
            [0,1,0]
        ])
        self._saveDrops = False
        self._flagMasked = False
        self._needProps = False
        self._nMaxPixels = 15
 

    def setFromDet(self, det):
        super(dropletFunc, self).setFromDet(det)
        if self._mask is None and det.mask is not None:
            setattr(self, '_mask', det.mask.astype(np.uint8))
        setattr(self, '_rms', det.rms)
        self._compData = np.ones_like(self._mask).astype(float)
        if self.useRms:
            if (len(self._rms.shape) > len(det.mask.shape)):
                self._compData *= self._rms[0]/det.gain[0]
            else:
                self._compData *= self._rms/det.gain
        self._is_tiled = False
        if len(det.ped.shape)>2:
            self.footprint = np.array([ 
                [[0,0,0],[0,0,0],[0,0,0]], 
                [[0,1,0],[1,1,1],[0,1,0]],  
                [[0,0,0],[0,0,0],[0,0,0]] 
            ])
            self._is_tiled = True
        return

            
    def addFunc(self, func):
        """ 
        If a function is added, set _saveDrops to True so that the droplet 
        image is passed to the next function for further analysis. Also sets
        other parameters in case the parameters are set in the sub-function.
        """
        super().addFunc(func)
        self._saveDrops = True
        self._needProps = getattr(func, '_needProps', self._needProps)
        return


    def applyThreshold(self, img, donut=False, invert=False, low=False):
        if not donut:
            if low:
                threshold = self.thresholdLow
            else:
                threshold = self.threshold
            if not invert:
                img[img < self._compData*threshold] = 0.0
            else:
                img[img > self._compData*threshold] = 0.0
        else:
            img[img < self._compData*self.thresholdLow] = 0.0
            img[img > self._compData*self.threshold] = 0.0
        return

            
    def neighborImg(self, img):
        return filters.maximum_filter(img, footprint=self._footprint2d)

    
    def prepareImg(self, img, donut=False, invert=False, low=False):
        imgIn = img.copy()
        if self._mask is not None:
            imgIn[self._mask==0] = 0
        self.applyThreshold(imgIn, donut, invert, low) # modifies in place
        return imgIn

    
    def returnEmpty(self):
        ret_dict = {'nDroplets': -1}
        for dropSave in self.dropletSaves:
            dropDict = dropSave.initArray()
            for key in dropDict.keys():
                ret_dict[key] = dropDict[key]
        if self.photonize:
            ret_dict['nPhotons'] = -1
            ret_dict['photons'] = np.zeros([self.maxPhotons,2])     
        return ret_dict


    def process(self, data):
        ret_dict = self.dropletize(data)

        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict[f'{k}_{kk}'] = subfuncResults[k][kk]
        return ret_dict

 
    def dropletize(self, data):
        tstart=time.time()
        if data is None:
            print('img is None!')
            self.ret_dict = self.returnEmpty()
            return
        time_start = time.time()
        img = self.prepareImg(data)
        img_drop = ndi.label(img, structure=self.footprint)
        time_label = time.time()
        
        # get all neighbors
        if (self.threshold != self.thresholdLow):
            if (len(img_drop[0].shape) == 2):
                imgDrop = self.neighborImg(img_drop[0])
            else:
                imgDrop = np.array([self.neighborImg(imgd_tile) for imgd_tile in img_drop[0]])
            img = self.prepareImg(data, low=True)

            if self.relabel:
                imgDrop[img==0]=0
                img_drop_relabel = ndi.label(imgDrop, structure=self.footprint)
                imgDrop = img_drop_relabel[0]
        else:
            imgDrop = img_drop[0]

        drop_ind = np.arange(1, img_drop[1]+1)
        ret_dict = {'nDroplets_all': len(drop_ind)} # number of droplets before ADU cut.
        adu_drop = ndi.sum_labels(img, labels=imgDrop, index=drop_ind)
        tfilled = time.time()

        if self.thresADU is not None:
            # Apply threshold to droplet total ADU (i.e. all pixels in droplet)
            vThres = np.where(adu_drop<self.thresADU)[0]
            vetoed = np.in1d(imgDrop.ravel(), (vThres+1)).reshape(imgDrop.shape)
            imgDrop[vetoed] = 0
            drop_ind_thres = np.delete(drop_ind, vThres)
            ret_dict['nDroplets'] = len(drop_ind_thres)
        else:
            drop_ind_thres = drop_ind
            ret_dict['nDroplets'] = ret_dict['nDroplets_all']

        if not self._saveDrops:
            return ret_dict


        # Get more info on the droplets if requested
        if not self._needProps:
            # Faster option
            drop_adu = np.array(ndi.sum_labels(img, labels=imgDrop, index=drop_ind_thres))
            pos_drop = np.array(ndi.center_of_mass(
                img,
                imgDrop,
                drop_ind_thres)
                )
            npix_drop = (ndi.sum_labels(
                img.astype(bool).astype(int),
                labels=imgDrop,
                index=drop_ind_thres
            )).astype(int)

            dat_dict = {'data': drop_adu}  # adu_drop}
            dat_dict['npix'] = npix_drop
            if drop_adu.shape[0] == 0:
                dat_dict['row'] = np.array([])
                dat_dict['col'] = np.array([])
            else:
                dat_dict['row'] = pos_drop[:, pos_drop.shape[1] - 2]
                dat_dict['col'] = pos_drop[:, pos_drop.shape[1] - 1]
                if self._is_tiled:
                    dat_dict['tile'] = pos_drop[:, 0]
                else:
                    dat_dict['tile'] = np.zeros_like(pos_drop[:, 0])
        else:
            # Use region props - this is not particularly performant on busy data.
            # this should be tested for tiled detectors!
            pos_drop = []
            moments = []
            bbox = []
            adu_drop = []
            npix_drop = []
            images = []
        
            self.regions = measure.regionprops(imgDrop,
                                               intensity_image=img,
                                               cache=True)
            dropSlices = ndi.find_objects(imgDrop)
            for droplet, ds in zip(self.regions, dropSlices):
                pos_drop.append(droplet['weighted_centroid'])
                moments.append(droplet['weighted_moments_central'])
                bbox.append(droplet['bbox'])
                adu_drop.append(droplet['intensity_image'].sum())
                npix_drop.append((droplet['intensity_image'] > 0).sum())
                pixelArray = droplet['intensity_image'].flatten()

                if pixelArray.shape[0] > self._nMaxPixels:
                    images.append(pixelArray[:self._nMaxPixels])
                else:
                    images.append(np.append(
                        pixelArray,
                        np.zeros(self._nMaxPixels-pixelArray.shape[0])
                        ))

            dat_dict = {'data': np.array(adu_drop)}
            dat_dict['npix'] = np.array(npix_drop)
            dat_dict['bbox'] = np.array(bbox)
            dat_dict['moments'] = np.array(moments)
            dat_dict['pixels'] = np.array(images)
            dat_dict['row'] = np.array(pos_drop)[:, 0]
            dat_dict['col'] = np.array(pos_drop)[:, 1]

        if self._flagMasked:
            maxImg = filters.maximum_filter(imgDrop, footprint=self.footprint)
            maskMax = ndi.sum_labels(self._mask, labels=maxImg, index=drop_ind)
            imgDropMin = imgDrop.copy()
            imgDropMin[imgDrop == 0] = (imgDrop.max() + 1)
            minImg = filters.minimum_filter(
                imgDropMin,
                footprint=self.footprint
                )
            minImg[minImg == (imgDrop.max()+1)] = 0
            maskMin = ndi.sum(self._mask, label=maxImg, index=drop_ind)
            maskDrop = maskMax + maskMin
            dat_dict['masked'] = maskDrop

        dat_dict['_mask'] = self._mask
        dat_dict['_imgDrop'] = imgDrop
        dat_dict['_image'] = img
        self.dat = dat_dict

        return ret_dict
