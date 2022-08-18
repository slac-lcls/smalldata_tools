import numpy as np
from scipy import sparse
import time
try:
    from ImgAlgos.PyAlgos import photons
except:
    pass
import scipy.ndimage.filters as filters
import itertools
from smalldata_tools.DetObjectFunc import DetObjectFunc

def fcn(buffer):
    if len(buffer[buffer<buffer[2]])>0:
        return max(buffer[buffer<buffer[2]])
    else:
        return 0.

class photonFunc(DetObjectFunc):
    """
    Wrapper for the psana algorithms described at
    
    Parameters: ADU_per_photon, mask, name, nphotMax, nphotRet, thresADU, retImg, ROI
    ADU_per_photon (def:154): expected value for a single photon in the detector in question
    mask (def:None): pass a mask in here, is None: use mask stored in DetObject
    name (def:'photon'): name used in hdf5 for data
    thresADU (def: 0.9): fraction of ADU_per_photon in the most neighboring pixels needed to count as a photon 
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','photon')
        super(photonFunc, self).__init__(**kwargs)
        self.ADU_per_photon = kwargs.get('ADU_per_photon',154)
        self.thr_fraction = kwargs.get('thr_fraction',0.9)
        self.thresADU = kwargs.get('thresADU',None)
        self._mask = kwargs.get('mask',None)
        print('photon init mask ',self._mask)

    def setFromDet(self, det):
        super(photonFunc, self).setFromDet(det)
        #photon uses the mask multiplicatively!
        if self._mask is None and det.mask is not None:
            mask = (det.mask.astype(bool))
            setattr(self, '_mask', mask.astype(np.uint8))

    def prepImage(self,image):
        """
        convert image from ADU to (fractional) photons
        #apply the ROI if applicable (to both image&mask)
        """
        if self.thresADU is not None:
            image[image<self.thresADU]=0
        image = image/self.ADU_per_photon
        return image

    def process(self,image):
        """
        use PyAlgo method to find photons.

        """
        tstart=time.time()
        fimage = self.prepImage(image)
        locMask = self._mask
        #photons_nda_nomask = photons(fimage, np.ones_like(self._mask),thr_fraction=self.thresADU)
        #note: this works on tiled detectors.
        photons_nda = photons(fimage, locMask,thr_fraction=self.thr_fraction)
        #photons_nda = photons(fimage, locMask)

        #make filling of histogram a function to be run on dict with nphot(ADU?),row,col 
        ret_dict = {'nPhot': photons_nda.sum()}
        #ret_dict['nPhot_mask']=photons_nda[locMask>0].sum() 
        #ret_dict = {'pHist': np.histogram(photons_nda.flatten(), np.arange(0,self.nphotMax))[0]}
        #the masked array will _mask_ the 1 in the mask!
        #given that the mask is used above, it should not do anything...
        self.dat = np.ma.masked_array(photons_nda, mask=(~(locMask.astype(bool))).astype(np.uint8))
        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]

        return ret_dict

#
#
#
import scipy.ndimage.measurements as measurements
class photon2(DetObjectFunc):
    """ deprecated: attemps to code pyalgo algo using pure python"""
    def __init__(self, **kwargs):
        self.ADU_per_photon = kwargs.get('ADU_per_photon',154)
        self.thresRms = kwargs.get('thresRms',3.)
        self.thresADU = kwargs.get('thresADU',0.7)
        self._name = kwargs.get('name','photon2')
        self._mask = kwargs.get('mask',None)
        self.rms = kwargs.get('rms',None)
        self.nphotMax = kwargs.get('nphotMax',25)
        self.nphotRet = kwargs.get('nphotRet',100)
        self.retImg = kwargs.get('retImg',0)

    def setFromDet(self, det):
        super(photon2, self).setFromDet(det)
        if self._mask is None and det.mask is not None:
            setattr(self, 'mask', det.mask.astype(bool))
        if self.rms is None and det.rms is not None:
            setattr(self, 'mask', det.rms)
            
    def prepImage(self,image):
        img = image.copy()
        img[self._mask==0]=0
        if self.rms is not None:
            img[img<self.rms*self.thresRms]=0
        return img/self.ADU_per_photon

    #just change name.
    def process(self,image):
    #    retDict = self.photon(image)
    #    return retDict
    #
    #def photon(self,image):
        """
        use a method has on the PyAlgo photon method w/ flexible threshold
    
        """
        tstart=time.time()
        fphotons = self.prepImage(image)
        fpi=fphotons//1
        fpf=fphotons%1
        # label the fractional image
        img_drop = measurements.label(fpf)
        drop_ind = np.arange(1,img_drop[1]+1)
        adu_drop = measurements.sum(fpf,img_drop[0], drop_ind)
        pos_drop = np.array(measurements.center_of_mass(fpf,img_drop[0], drop_ind)).astype(int)
        # for sum ADU > threshold: 1 @ COM(x,y) as int (sparse array -> add to full photon image)
        vThres = np.where(adu_drop>self.thresADU)[0]
        adu_phot = adu_drop[vThres]
        com_phot = (pos_drop[vThres]+0.5)//1
        # save nphot & sparse full img as for photon
        #deal with 0.9 vs 1.1 etc.
        sImgSparse = sparse.coo_matrix((adu_phot//1+(adu_phot%1>0.5).astype(int), (com_phot[:,0], com_phot[:,1])),shape=fpi.shape).todense()
        photImg = fpi+sImgSparse
        ret_dict = {'pHist': np.histogram(photImg.flatten(), np.arange(0,self.nphotMax))[0]}
        ret_dict['nPhot']=photImg.sum()

        if self.retImg>0:
            #now sparsify the image (data, xIdx, yIdx)
            if self.retImg>1:
                #print('DEBUG: shape ',photImg.shape, fphotons.shape
                ret_dict['img']=photImg.copy()
            else:
                sImage = sparse.coo_matrix(photImg)
                if sImage.data.shape[0] >= self.nphotRet:
                    ret_dict['data']=sImage.data[:self.nphotRet]
                    ret_dict['row']=sImage.row[:self.nphotRet]
                    ret_dict['col']=sImage.col[:self.nphotRet]
                else:
                    ret_dict['data']=(np.append(sImage.data[:self.nphotRet], np.zeros(self.nphotRet-len(sImage.data)))).astype(int)
                    ret_dict['row']=(np.append(sImage.row[:self.nphotRet], np.zeros(self.nphotRet-len(sImage.row)))).astype(int)
                    ret_dict['col']=(np.append(sImage.col[:self.nphotRet], np.zeros(self.nphotRet-len(sImage.col)))).astype(int)
        ret_dict['evtTime']=time.time() - tstart

        # re-label this img, and make hist of relabelled ADU histo as well.
        img_drop_phot = measurements.label(photImg)
        adu_drop = measurements.sum(photImg,img_drop_phot[0], np.arange(1,img_drop_phot[1]+1))
        ret_dict['pHist_RL'] = np.histogram(adu_drop, np.arange(0,self.nphotMax))[0]
        return ret_dict



class photon3(DetObjectFunc):
    """deprecated: recode pyalgo photon algorithms w/ flexible threshold before that was introduced."""
    def __init__(self, **kwargs):
        self.ADU_per_photon = kwargs.get('ADU_per_photon',154)
        self.thresRms = kwargs.get('thresRms',3.)
        self.thresADU = kwargs.get('thresADU',0.9)
        self._name = kwargs.get('name','photon3')
        self._mask = kwargs.get('mask', None)
        self.rms = kwargs.get('rms', None)
        self.nphotMax = kwargs.get('nphotMax', 25)
        self.nphotRet = kwargs.get('nphotRet', 100)
        self.retImg = kwargs.get('retImg', 0)
        self.maxMethod = kwargs.get('maxMethod', 0)

    def setFromDet(self, det):
        super(photon3, self).setFromDet(det)
        if self._mask is None and det.mask is not None:
            setattr(self, 'mask', det.mask.astype(bool))
        if self.rms is None and det.rms is not None:
            setattr(self, 'mask', det.rms)

    def prepImage(self,image):
        img = image.copy()
        img[self._mask==0]=0
        if self.rms is not None:
            img[img<self.rms*self.thresRms]=0
        return img/self.ADU_per_photon

    #just change name.
    def process(self,image):
        retDict = self.photon(image)
        return retDict

    def photon(self,image):
        """
        use a method has on the PyAlgo photon method w/ flexible threshold
    
        """
        tstart=time.time()    
        fphotons = self.prepImage(image)
        fpi=fphotons//1
        #ignore partial photons.
        if self.retImg<0:
            ret_dict={'img':fpi}
            return ret_dict

        fpf=fphotons%1
        #
        fpf_seed=fpf.copy()
        fpf_seed[fpf<0.5]=0
        if self.maxMethod==-1:
            #generic filter is not fast: likely loop in python as pre-optimization is not possible
            fpf_max=filters.generic_filter(fpf, fcn, footprint=np.array([[0,1,0],[1,1,1],[0,1,0]]),mode='constant')
        else:
            fpf_max=filters.maximum_filter(fpf,footprint=np.array([[0,1,0],[1,0,1],[0,1,0]]),mode='constant')
        fpf_max[fpf<0.5]=0
        if self.maxMethod>0:
            if self.maxMethod>=2:
                fpf_max2nd=filters.rank_filter(fpf,rank=-2,footprint=np.array([[0,1,0],[1,0,1],[0,1,0]]),mode='constant')
                fpf_max = np.where(fpf_max>fpf_seed,fpf_max2nd,fpf_max)
            fpf_max[fpf_max>fpf_seed]=0

        fpf_max+=fpf_seed
        fpf_max[fpf_max<self.thresADU]=0
        fpf_max[fpf_max>0]=1
            
        photImg = fpf_max+fpi
        ret_dict = {'pHist': np.histogram(photImg.flatten(), np.arange(0,self.nphotMax))[0]}
        ret_dict['nPhot']=photImg.sum()

        if self.retImg>0:
            #now sparsify the image (data, xIdx, yIdx)
            if self.retImg>1:
                ret_dict['img']=photImg
            else:
                sImage = sparse.coo_matrix(photImg)
                if sImage.data.shape[0] >= self.nphotRet:
                    ret_dict['data']=sImage.data[:self.nphotRet]
                    ret_dict['row']=sImage.row[:self.nphotRet]
                    ret_dict['col']=sImage.col[:self.nphotRet]
                else:
                    ret_dict['data']=(np.append(sImage.data, np.zeros(self.nphotRet-len(sImage.data)))).astype(int)
                    ret_dict['row']=(np.append(sImage.row, np.zeros(self.nphotRet-len(sImage.row)))).astype(int)
                    ret_dict['col']=(np.append(sImage.col, np.zeros(self.nphotRet-len(sImage.col)))).astype(int)
        ret_dict['evtTime']=time.time() - tstart

        # re-label this img, and make hist of relabelled ADU histo as well.
        img_drop_phot = measurements.label(photImg)
        adu_drop = measurements.sum(photImg,img_drop_phot[0], np.arange(1,img_drop_phot[1]+1))
        ret_dict['pHist_RL'] = np.histogram(adu_drop, np.arange(0,self.nphotMax))[0]
        return ret_dict
