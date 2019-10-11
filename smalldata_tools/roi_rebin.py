import os
import numpy as np
from scipy import sparse

import time
from smalldata_tools.utilities import rebin, getBins
from smalldata_tools.DetObject import DetObjectFunc
from smalldata_tools.droplet import dropletFunc
#from droplet import dropletFunc

#
# for now, this works in "raw" coordinates.
# enable ROI in other coordinates: e.g. x, y: save ROI in x/y save across tiles.
# add asImg. Is this an option or a separate function?
#
class ROIFunc(DetObjectFunc):
    """
    apply as ROI to input data
    ROI: boundary for ROI
    writeArea (default False): if true, save pixels in ROI
    calcPars (defalt True): if True, calculate sum, center-of-mass and highest entry in ROI
    userMask (default None): apply mask for calculation or subfunctions
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','ROI')
        super(ROIFunc, self).__init__(**kwargs)
        ROI = kwargs.get('ROI',None)
        if ROI is None:
            if self._name == 'ROI': self._name='full'
            ROI=[0, 1e6]
        if isinstance(ROI, list):
            ROI = np.array(ROI)
        else:
            ROI=ROI
        if ROI.ndim==1 and ROI.shape[0]>2:
            ROI = ROI.reshape(ROI.shape[0]/2,2)
        self.bound = ROI.squeeze()
        #print 'DEBUG: def ROI', self.bound
        self.writeArea = kwargs.get('writeArea',False)
        self._calcPars = kwargs.get('calcPars',True)
        self._mask =  kwargs.get('userMask',None)

    def setFromDet(self, det):
        super(ROIFunc, self).setFromDet(det)
        self._rms = det.rms
        if self._rms is not None:
            if self._rms.ndim==4:
                self._rms = self._rms[0]#.squeeze()
            self._rms = self.applyROI(self._rms)
        if self._mask is not None and det.mask is not None and self._mask.shape == det.mask.shape:
            self._mask = ~(self._mask.astype(bool)&det.mask.astype(bool))
        else:
            try:
                self._mask = ~(det.cmask.astype(bool)&det.mask.astype(bool))
            except:
                try:
                    self._mask = ~(np.ones_like(det.ped).astype(bool))
                except:
                    self._mask = None

        if self._mask is not None:
            try:
                self._mask = self.applyROI(self._mask)
            except:
                pass

        try:
            self._x = self.applyROI(det.x)
            self._y = self.applyROI(det.y)
        except:
            pass

    def addFunc(self, func):
        if isinstance(func, dropletFunc) and isinstance(self, ROIFunc):
            try:
                func.setKeyData('_mask', (~self._mask).astype(np.uint8))
                func.setKeyData('_compData', self._rms)
            except:
                print('failed to set parameters needed to run droplets on ROI')
                return
        self.__dict__[func._name] = func

    def applyROI(self, array):
        #array = np.squeeze(array) #added for jungfrau512k. Look here if other detectors are broken now...
        if array.ndim < self.bound.ndim:
            print('array has fewer dimensions that bound: ',array.ndim,' ',len(self.bound))
            return array
        #ideally would check this from FrameFexConfig and not on every events
        if array.ndim == self.bound.ndim:
            new_bound=[]
            if array.ndim==1:
                new_bound=[max(min(self.bound), 0), min(max(self.bound), array.shape[0])]
            else:
                for dim,arsz in zip(self.bound, array.shape):
                    if max(dim) > arsz:
                        new_bound.append([min(dim), arsz])
                    else:
                        new_bound.append([min(dim), max(dim)])
            self.bound = np.array(new_bound)
        elif self.bound.ndim==1:
            self.bound = np.array([min(self.bound), min(max(self.bound), array.shape[0])]).astype(int)
        #this needs to be more generic....of maybe just ugly spelled out for now.
        if self.bound.shape[0]==2 and len(self.bound.shape)==1:
            subarr = array[self.bound[0]:self.bound[1]]
        elif self.bound.shape[0]==2 and len(self.bound.shape)==2:
            subarr = array[self.bound[0,0]:self.bound[0,1],self.bound[1,0]:self.bound[1,1]]
        elif self.bound.shape[0]==3:
            subarr = array[self.bound[0,0]:self.bound[0,1],self.bound[1,0]:self.bound[1,1],self.bound[2,0]:self.bound[2,1]]
        return subarr.copy()
    #calculates center of mass of first 2 dim of array within ROI using the mask
    def centerOfMass(self, array):        
        array=array.squeeze()
        if array.ndim<2:
            return np.nan, np.nan
        #why do I need to invert this?
        X, Y = np.meshgrid(np.arange(array.shape[0]), np.arange(array.shape[1]))
        try:
            imagesums = np.nansum(np.nansum(array,axis=0),axis=0)
            centroidx = np.nansum(np.nansum(array*X.T,axis=0),axis=0)
            centroidy = np.nansum(np.nansum(array*Y.T,axis=0),axis=0)
            return centroidx/imagesums,centroidy/imagesums 
        except:
            return np.nan,np.nan
    def addNsat(self,highLim=None):
        self.Nsat = highLim
    def process(self, data):
        ret_dict = {}
        ROIdata=self.applyROI(data)
        if self._mask is not None:
            if ROIdata.dtype == np.uint16:
                ROIdata[self._mask]=0
            else:
                ROIdata[self._mask]=np.nan
###
### this is why processFuns works. For droplet & photon, store dict of ADU, x, y (& more?)
### photons: row & col. make droplet use that too.
###
        self.dat = ROIdata
###
        if self.writeArea:
            ret_dict['area'] = ROIdata.squeeze()
        if self._calcPars:
            ret_dict['max'] = np.nanmax(ROIdata.astype(np.float64))
            ret_dict['sum'] = np.nansum(ROIdata.astype(np.float64))
            ret_dict['com'] = self.centerOfMass(ROIdata)
        if 'Nsat' in self.__dict__.keys():
            ret_dict['nsat'] =  (ROIdata >= self.Nsat).astype(int).sum()

        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                if isinstance(subfuncResults[k][kk], list):
                    try:
                        ret_dict['%s_%s'%(k,kk)] = np.array(subfuncResults[k][kk])
                    except:
                        print 'issue with: ',subfuncResults[k][kk], '%s_%s'%(k,kk), len(subfuncResults[k][kk])

                else:
                    ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]
        #print 'ret_dict ',ret_dict.keys()
        return ret_dict

class rebinFunc(DetObjectFunc):
    """
    function to rebin input data to new shape
    shape: desired shape of array
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','rebin')
        super(rebinFunc, self).__init__(**kwargs)
        self.shape =  kwargs.get('shape',None)
        if self.shape is None:
            print('need shape to rebin data to!')
            return None
        self.shape = np.array(self.shape)

    def process(self, data):
        #masked array????
        newArray = rebin(data, self.shape)
        ret_dict = {'data': newArray}
        return ret_dict

class projectionFunc(DetObjectFunc):
    def __init__(self, **kwargs):
        """
        add a projection of the ROI to the user data. Parameters are pjName, axis, singlePhoton, mean, thresADU and thresRms
        pjName (def ''): name of projection data in hdf5 files. Needed to save >1 projection/ROI.
        axis (def -1): axis to project onto. axis=-1 will save a single number
        mean (def False): if False, save sum, save mean otherwise
        thresADU (def 1e-6): pixels with ADU < thresADU will be set to 0 before summing.
        thresRms (def 1e-6): pixels with ADU < thresRms*rms will be set to 0 before summing
        singlePhoton (def False): set all pixels > threshold to 1
        """
        self.axis =  kwargs.get('axis',-1)
        self._name = kwargs.get('name','pj_ax_%d'%abs(self.axis))
        super(projectionFunc, self).__init__(**kwargs)
        self.thresADU =  kwargs.get('thresADU',1e-6)
        self.thresRms =  kwargs.get('thresRms',1e-6)
        self.singlePhoton =  kwargs.get('singlePhoton',False)
        self.mean =  kwargs.get('mean',False)
    def process(self,data):
        if isinstance(data, np.ma.masked_array):
            array = data.data.copy().squeeze()
        else:
            array = data.copy().squeeze()
        array[array<self.thresADU]=0
        if 'rms' in self.__dict__.keys() and self.rms is not None:
            array[array<self.thresRms*self.rms.squeeze()]=0
        if self.singlePhoton:
            array[array>0]=1
        #array.data = array.data.astype(np.float64)
        if self.mean:
            if self.axis<0:
                retDict={'data': np.nanmean(array.squeeze())}
            else:
                retDict={'data': np.nanmean(array.squeeze(),axis=self.axis)}
        if self.axis<0:
            retDict={'data': np.nanmean(array.squeeze())}
        else:
            retDict={'data': np.nanmean(array.squeeze(),axis=self.axis)}
        return retDict

#effectitely a projection onto a non-spatial coordinate.
class spectrumFunc(DetObjectFunc):
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','spectrum')
        super(spectrumFunc, self).__init__(**kwargs)
        #need to get bin range
        #need to get ROI ? ROI w/ ROI func before!
        self.bins =  kwargs.get('bins',None)
        if self.bins is None or isinstance(self.bins, np.ndarray):
            print('need bins as a list to rebin data to!')
            return None
        #self.bins = np.array(getBins(self.bins))
        self.bins = getBins(self.bins)

    def process(self, data):
        if isinstance(data, np.ma.masked_array):
            his=np.histogram(data.compressed(), self.bins)
        elif  isinstance(data, np.ndarray):
            his=np.histogram(data[~np.isnan(data)].flatten(), bins=self.bins)
        elif  isinstance(data, dict) and 'data' in data:
            his=np.histogram(data['data'], self.bins)
        else:
            print('cannot make a spectrum of input data', data)
            
        ret_dict = {'histogram': his[0]}

        #store for further processing
        self.dat = his[0]
        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]

        return ret_dict

#effectitely a projection onto a non-spatial coordinate.
class sparsifyFunc(DetObjectFunc):
    """
    Function to sparisify a passed array (2 or 3-d input)
    nData: if passed, make output array rectangular (for storing in event based smlData)
    if a dictionary w/ data, row, col is passed, only make rectangular
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','sparse')
        super(sparsifyFunc, self).__init__(**kwargs)
        self.nData =  kwargs.get('nData',None)
        self._flagMasked = kwargs.get('flagMasked',False)
        self._needProps = kwargs.get('needProps',False)

    def process(self, data):
        #apply mask - set to zero, so pixels will fall out in sparify step.
        if isinstance(data, np.ma.masked_array):            
            data = data.filled(data, fill_value=0)

        #already have dict w/ data
        if  isinstance(data, dict):
            if 'data' not in data or 'row' not in data or 'col' not in data:
                print('cannot make a make a sparse, rectangular array of', data)
                return
            ret_dict = data

        #sparsify image
        if  isinstance(data, np.ndarray):
            photonsImg = data.copy()
            if len(photonsImg.shape)>2: #tiled detector!
                data=[]
                row=[]
                col=[]
                tile=[]
                for iTile,photonTile in enumerate(photonsImg):
                    sImage = sparse.coo_matrix(photonTile)
                    data = list(itertools.chain.from_iterable([data, sImage.data.tolist()]))
                    row = list(itertools.chain.from_iterable([row, sImage.data.tolist()]))
                    col = list(itertools.chain.from_iterable([col, sImage.data.tolist()]))
                    tile = list(itertools.chain.from_iterable([tile, (np.ones_like(sImage.data)*iTile).tolist()]))
                data = np.array(data)
                row = np.array(row)
                col = np.array(col)
                tile= np.array(tile)
            else:
                sImage = sparse.coo_matrix(photonsImg)
                data = sImage.data
                row = sImage.row
                col = sImage.col
                tile = np.zeros_like(data)
            ret_dict={'data':data}
            ret_dict['row']=row
            ret_dict['col']=col
            ret_dict['tile']=tile
            
        #now fix shape of data in dict.
        if self.nData > 0:
            for key in ret_dict.keys():
                if ret_dict[key].shape[0] >= self.nData:
                    ret_dict[key]=ret_dict[key][:self.nData]
                else:
                    ret_dict[key]=(np.append(ret_dict[key], np.zeros(self.nData-len(ret_dict[key])))).astype(int)

        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]

        return ret_dict


class imageFunc(DetObjectFunc):
    """
    function to cast detector data into a different coordinate system
    examples: 3-d tiled detector to x/y image.
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','image')
        super(imageFunc, self).__init__(**kwargs)
        self.coords = kwargs.get('coords',None)
        self.imgShape = kwargs.get('imgShape',None)
        self.correction = kwargs.get('correction',None)
        self.mask = kwargs.get('mask',None)
        if self.mask is not None: self.mask = np.asarray(self.mask,dtype=np.bool).flatten()

    def setFromDet(self, det):
        super(imageFunc, self).setFromDet(det)
        #set imgShape if obvious (and not previously set)
        if self.imgShape is None:
            #cspad
            if det.x.shape==(32,185,388): self.imgShape=[1689,1689]
            #cs140k
            #elif det.x.shape==(2,185,388): self.imgShape=[391,371] #at least for one geometry
            elif det.x.shape==(2,185,388): self.imgShape=[371,391]
            #epix100a
            elif det.x.shape==(704,768): self.imgShape=[709,773]
            #jungfrau512k
            elif det.x.shape==(1,512,1024): self.imgShape=[514,1030]
            elif det.x.shape==(512,1024): self.imgShape=[514,1030]
            #jungfrau1M
            elif det.x.shape==(2,512,1024): self.imgShape=[1064,1030]
            #epix10a
            elif det.x.shape==(704,768): self.imgShape=[709,773]

        if self.mask is None:
            self.mask = ~(det.cmask.astype(bool)&det.mask.astype(bool)).flatten()

        if isinstance(self.coords, list):
            for coord in self.coords:
                self.__dict__[coord] = getattr(det, coord)
                try:
                    self.__dict__['i%s'%coord] = getattr(det, 'i%s'%coord).flatten()
                except:
                    #
                    # DEBUG ME: looks right, but not sure
                    #
                    intcoord = self.__dict__[coord].copy()
                    intcoord = intcoord - np.nanmin(intcoord)
                    intcoord = (intcoord/np.nanmax(intcoord)*self.imgShape[0]).astype(int)
                    self.__dict__['i%s'%coord] = intcoord

            self._coordTuple = tuple( self.__dict__['i%s'%coord] for coord in self.coords)
            self._n_coordTuple = tuple( np.max(self.__dict__['i%s'%coord]+1) for coord in self.coords)
            try:
                self._multidim_idxs = np.ravel_multi_index(self._coordTuple, self._n_coordTuple)
                #DEBUG ME HERE....
                #self._multidim_idxs[self.mask.ravel()] = 0; # send the masked ones in the first bin
                self._n_multidim_idxs = 1
                for nBin in self._n_coordTuple:
                    self._n_multidim_idxs *= int(nBin)
                self._npix_flat = np.bincount(self._multidim_idxs,minlength=int(self._n_multidim_idxs))
                self.npix = self._npix_flat.reshape(self._n_coordTuple)
                self._npix_div = None
                if self.npix.max()>1: #if we map multiple pixels into one image pixel, we need to normalize
                    self._npix_div = self.npix.copy()
                    self._npix_div[self.npix>1] = 1./self.npix[self.npix>1]
                    #print self._npix_div.max()
                
            except:
                pass

        if self.imgShape is None:
            self.imgShape = det.imgShape

        if len(self.coords)==2:
            self.imgShape = (int(max(np.max(self.__dict__['i%s'%self.coords[0]])+1, \
                                     self.imgShape[0])), \
                             int(max(np.max(self.__dict__['i%s'%self.coords[1]])+1, \
                                     self.imgShape[1])))
            if self.mask is not None:
                self.mask_img = np.array(sparse.coo_matrix((self.mask.flatten(),(self.__dict__['i%s'%self.coords[0]],self.__dict__['i%s'%self.coords[1]])), shape=self.imgShape).todense())
                self.mask_img[self.mask_img!=0]=1
                ones_mask = np.array(sparse.coo_matrix((np.ones_like(self.mask.flatten()),(self.__dict__['i%s'%self.coords[0]],self.__dict__['i%s'%self.coords[1]])), shape=self.imgShape).todense())
                self.mask_ones = ones_mask
                self.mask_img[ones_mask==0]=1
                self.mask_img = self.mask_img.astype(int)

    def process(self, data):
        #already have dict w/ data
        if  isinstance(data, dict):
            if 'data' not in data or 'row' not in data or 'col' not in data:
                print('dictionay, but does not have sparsified data', data)
                return
            #
            # should ideally specify the output shape if possible.
            #
            #img = sparse.coo_matrix((d.flatten(), (ix.flatten(), iy.flatten())), shape=outShape).todense()
            if 'tile' not in data or max(data['tile'])==0:
                data = sparse.coo_matrix(data['data'],(data['row'],data['col'])).todense()
            else:
                data = np.array([ sparse.coo_matrix(data['data'][data['tile']==itile],(data['row'][data['tile']==itile],data['col'][data['tile']==itile])).todense() for itile in range(1+max(data['tile'])) ])

        elif not isinstance(data, np.ndarray):
            print('cannot work with this data, need array or sparsified array')

        if self.correction is not None:
            data /= self.correction

        if self.coords is None:
            return {'img': data}

        #use sparse matrix method to create image in new coordinate space. 
        #also: check if this is special case of multidim.
        if len(self.coords)==2:
            retDict={}
            ##using the sparse array os slower (1ms versus 0.55ms for bincount)
            ##also, mapping of multiple pixels into one final pixel will not be normalizable
            ##as a notel the normalization costs about 0.2ms
            #data2d = sparse.coo_matrix((data.flatten(),(self.__dict__['i%s'%self.coords[0]],self.__dict__['i%s'%self.coords[1]])), shape=self.imgShape).todense()            
            #retDict['img_sparse'] = np.array(data2d)
            I=np.bincount(self._multidim_idxs, weights = data.flatten(), minlength=int(self._n_multidim_idxs))
            I=np.reshape(I[:self._n_multidim_idxs], self._n_coordTuple)
            retDict['img_bc'] = I
            if self._npix_div is not None:
                Inorm=I*self._npix_div
                retDict['img_bc'] = Inorm
            else:
                retDict['img_bc'] = I
            return retDict

        elif len(self.coords)==1: #this might be a special case of the multi dim thing....
            data = np.bincount(data, self.__dict__['i%s'%self.coord[0]])
            return {'img_1d': data/self.npix}            
        else:
            I=np.bincount(self._multidim_idxs, weights = data, minlength=self._n_multidim_idxs)
            I=I[:self._n_multidim_idxs]
        return {'img_n': np.reshape(I, self._n_coordTuple)}

