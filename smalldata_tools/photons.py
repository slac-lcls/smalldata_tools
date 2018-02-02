import numpy as np
from scipy import sparse
import time
from ImgAlgos.PyAlgos import photons
import scipy.ndimage.filters as filters

def fcn(buffer):
    if len(buffer[buffer<buffer[2]])>0:
        return max(buffer[buffer<buffer[2]])
    else:
        return 0.

class photon:
    def __init__(self, ADU_per_photon=154, mask=None, rms=None, name='photon', nphotMax=25, nphotRet=100, thresADU=0.9, retImg=0):
        self.ADU_per_photon = ADU_per_photon
        self.name = name
        self.mask = mask
        self.rms = rms
        self.nphotMax = nphotMax
        self.nphotRet = nphotRet
        self.retImg = retImg
        self.thresADU = thresADU

    def prepImage(self,image):
        return image/self.ADU_per_photon

    def photon(self,image):
        """
        use PyAlgo method to find photons.
    
        """
        tstart=time.time()
        fimage = self.prepImage(image)
        photons_nda_nomask = photons(fimage, np.ones_like(self.mask),thr_fraction=self.thresADU)
        photons_nda = photons(fimage, self.mask)

        ret_dict = {'pHist': np.histogram(photons_nda.flatten(), np.arange(0,self.nphotMax))[0]}
        ret_dict['nPhot']=photons_nda.sum()
        ret_dict['nPhot_mask']=photons_nda[self.mask>0].sum()
        ret_dict['nPhot_nomask']=photons_nda_nomask.sum()

        if self.retImg>0:
            #next line: work around fact that mask sometimes gets ignored.
            photonsImg = photons_nda.copy()
            photonsImg[self.mask==0]=0
            if self.retImg>1:
                ret_dict['img']=photonsImg
            else:    
                sImage = sparse.coo_matrix(photonsImg)
                #back to normal code
                data = sImage.data
                if data.shape[0] >= self.nphotRet:
                    ret_dict['data']=data[:self.nphotRet]
                    ret_dict['row']=sImage.row[:self.nphotRet]
                    ret_dict['col']=sImage.col[:self.nphotRet]
                else:
                    ret_dict['data']=(np.append(data[:self.nphotRet], np.zeros(self.nphotRet-len(data)))).astype(int)
                    ret_dict['row']=(np.append(sImage.row[:self.nphotRet], np.zeros(self.nphotRet-len(sImage.row)))).astype(int)
                    ret_dict['col']=(np.append(sImage.col[:self.nphotRet], np.zeros(self.nphotRet-len(sImage.col)))).astype(int)
        ret_dict['evtTime']=time.time() - tstart
        #print 'shapes: ',ret_dict['data'].shape,ret_dict['row'].shape,ret_dict['col'].shape
        #print ret_dict['data'][:10]
        #print ret_dict['row'][:10]
        #print ret_dict['col'][:10]
        return ret_dict

#
#
#
import scipy.ndimage.measurements as measurements
class photon2:
    def __init__(self, ADU_per_photon=154, thresADU=0.7, thresRms=3., mask=None, rms=None, name='photon2', nphotMax=25, retImg=0, nphotRet=100):
        self.ADU_per_photon = ADU_per_photon
        self.thresRms=thresRms
        self.thresADU=thresADU
        self.name = name
        self.mask = mask
        self.rms = rms
        self.nphotMax = nphotMax
        self.nphotRet = nphotRet
        self.retImg = retImg

    def prepImage(self,image):
        img = image.copy()
        img[self.mask==0]=0
        img[img<self.rms*self.thresRms]=0
        return img/self.ADU_per_photon

    def photon(self,image):
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
                #print 'DEBUG: shape ',photImg.shape, fphotons.shape
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



class photon3:
    def __init__(self, ADU_per_photon=154, thresADU=0.9, thresRms=3., mask=None, rms=None, name='photon3', nphotMax=25, retImg=0, nphotRet=100, maxMethod=0):
        self.ADU_per_photon = ADU_per_photon
        self.thresRms=thresRms
        self.thresADU=thresADU
        self.name = name
        self.mask = mask
        self.rms = rms
        self.nphotMax = nphotMax
        self.nphotRet = nphotRet
        self.retImg = retImg
        self.maxMethod = maxMethod

    def prepImage(self,image):
        img = image.copy()
        img[self.mask==0]=0
        img[img<self.rms*self.thresRms]=0
        return img/self.ADU_per_photon

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
