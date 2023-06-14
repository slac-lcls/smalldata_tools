import os
import numpy as np
import numpy.ma as ma
import itertools
from scipy import sparse

import time
from smalldata_tools.utilities import rebin, getBins
from smalldata_tools.DetObjectFunc import DetObjectFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc

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
        self.writeArea = kwargs.get('writeArea',False)
        self._calcPars = kwargs.get('calcPars',True)
        self.mask =  kwargs.get('userMask',None)
        self.thresADU = kwargs.get('thresADU',None)

    def setFromDet(self, det):
        super(ROIFunc, self).setFromDet(det)
        self._rms = det.rms
        if self._rms is not None:
            if self._rms.ndim==4:
                self._rms = self._rms[0]#.squeeze()
            self._rms = self.applyROI(self._rms)
        if self.mask is not None and det.mask is not None and self.mask.shape == det.mask.shape:
            self.mask = ~(self.mask.astype(bool)&det.mask.astype(bool))
        else:
            try:
                self.mask = ~(det.cmask.astype(bool)&det.mask.astype(bool))
            except:
                if det.ped is None:
                    self.mask = None
                else:
                    try:
                        self.mask = ~(np.ones_like(det.ped).astype(bool))
                    except:
                        self.mask = None
        if self.mask is not None:
            try:
                self.mask = self.applyROI(self.mask)
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
                func.setKeyData('_mask', (~self.mask).astype(np.uint8))
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
            imagesums = np.sum(np.sum(array,axis=0),axis=0)
            centroidx = np.sum(np.sum(array*X.T,axis=0),axis=0)
            centroidy = np.sum(np.sum(array*Y.T,axis=0),axis=0)
            return centroidx/imagesums,centroidy/imagesums 
        except:
            return np.nan,np.nan

    def addNsat(self,highLim=None):
        self.Nsat = highLim

    def process(self, data):
        ret_dict = {}
        ROIdata=self.applyROI(data)
        if self.mask is not None:
            ROIdata = ma.array(ROIdata, mask=self.mask)
        else:
            ROIdata = ma.array(ROIdata)
            #if ROIdata.dtype == np.uint16:
            #    ROIdata[self.mask]=0
            #else:
            #    ROIdata[self.mask]=np.nan
###
### this is why processFuns works. For droplet & photon, store dict of ADU, x, y (& more?)
### photons: row & col. make droplet use that too.
###
        if self.thresADU is not None:
            if isinstance(self.thresADU, list):
                ROIdata[ROIdata<self.thresADU[0]] = 0
                ROIdata[ROIdata>self.thresADU[1]] = 0
            else:
                ROIdata[ROIdata<self.thresADU] = 0
        self.dat = ROIdata
###
        if self.writeArea:
            ret_dict['area'] = ROIdata.data.squeeze()
        if self._calcPars:
            ret_dict['sum'] = ROIdata.filled(fill_value=0).sum()
            ret_dict['mean'] = ROIdata.filled(fill_value=0).mean()
            ret_dict['max'] = ROIdata.filled(fill_value=0).max()
            ret_dict['com'] = self.centerOfMass(ROIdata)
        if 'Nsat' in self.__dict__.keys():
            ret_dict['nsat'] =  (ROIdata.filled(fill_value=0) >= self.Nsat).astype(int).sum()

        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                if isinstance(subfuncResults[k][kk], list):
                    try:
                        ret_dict['%s_%s'%(k,kk)] = np.array(subfuncResults[k][kk])
                    except:
                        print('issue with: ',subfuncResults[k][kk], '%s_%s'%(k,kk), len(subfuncResults[k][kk]))

                else:
                    ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]
        return ret_dict

#DEBUG ME WITH MASKED ARRAY#
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

#TEST ME WITH MASKED ARRAY#
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
        self.thresADU =  kwargs.get('thresADU',None)
        self.thresRms =  kwargs.get('thresRms',None)
        self.singlePhoton =  kwargs.get('singlePhoton',False)
        self.mean =  kwargs.get('mean',False)

    def process(self,data):
        array = np.ma.array(data.copy().squeeze())
        if self.thresADU is not None:
            array.data[array.data<self.thresADU]=0
        if self.thresRms is not None and 'rms' in self.__dict__.keys() and self.rms is not None:
            array.data[array.data<self.thresRms*self.rms.squeeze()]=0
        if self.singlePhoton:
            array.data[array.data>0]=1
        #array.data = array.data.astype(np.float64)
        if self.mean:
            if self.axis<0:
                retDict={'data': np.nanmean(array)}
            else:
                meanRes = np.nanmean(array,axis=self.axis)
                if isinstance(meanRes, np.ma.masked_array):
                    meanRes = meanRes.data
                retDict={'data': meanRes}
        else:
            if self.axis<0:
                retDict={'data': np.nansum(array)}
            else:
                sumRes = np.nansum(array,axis=self.axis)
                if isinstance(sumRes, np.ma.masked_array):
                    sumRes = sumRes.data
                retDict={'data': sumRes}
        return retDict

#effectitely a projection onto a non-spatial coordinate.
#TEST ME WITH MASKED ARRAY - WAS ALEADY WRITTEN WITH MASKED ARRAY IN MIND#
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


#DEBUG ME WITH MASKED ARRAY#
class imageFunc(DetObjectFunc):
    """
    function to cast detector data into a different coordinate system
    examples: 3-d tiled detector to x/y image.
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','image')
        super(imageFunc, self).__init__(**kwargs)
        self._coords = kwargs.get('coords',['x','y'])
        self.imgShape = kwargs.get('imgShape',None)
        self.correction = kwargs.get('correction',None)
        self.mask = kwargs.get('mask',None)
        if self.mask is not None: self.mask = np.asarray(self.mask,dtype=np.bool).flatten()

    def setFromDet(self, det):
        super(imageFunc, self).setFromDet(det)
        self._pedShape = det.ped.shape
        #set imgShape if obvious (and not previously set)
        if self.imgShape is None:
            try:
                self.imgShape = [det.ix.max()-det.ix.min()+1,
                                 det.iy.max()-det.iy.min()+1]
            except:
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

        if isinstance(self._coords, list) or isinstance(self._coords, dict):
            if isinstance(self._coords, list):
                for coord in self._coords:
                    if not hasattr(det, coord): 
                        print('Could not get information for coordinate ',coord,' from detector')
                        continue
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
            else:
                coordNames=[]
                for coord, values in self._coords.items():
                    coordNames.append(coord)
                    self.__dict__[coord] = values
                    if isinstance(values.flatten()[0], (int, np.uint64)):
                        self.__dict__['i%s'%coord] = values
                    else:
                        intcoord = self.__dict__[coord].copy()
                        intcoord = intcoord - np.nanmin(intcoord)
                        intcoord = (intcoord/np.nanmax(intcoord)*self.imgShape[len(coordNames)-1]).astype(int)
                        self.__dict__['i%s'%coord] = intcoord
                self._coords = coordNames

            self._coordTuple = tuple( self.__dict__['i%s'%coord].flatten() for coord in self._coords)
            self._n_coordTuple = tuple( np.max(self.__dict__['i%s'%coord]+1) for coord in self._coords)
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
                
            except:
                pass

        if self.imgShape is None:
            self.imgShape = det.imgShape

        if self._coords is not None and len(self._coords)==2:
            self.imgShape = (int(max(np.max(self.__dict__['i%s'%self._coords[0]])+1, \
                                     self.imgShape[0])), \
                             int(max(np.max(self.__dict__['i%s'%self._coords[1]])+1, \
                                     self.imgShape[1])))
            if self.mask is not None:
                maskArray = np.array(self.mask).flatten()
                self.mask_img = np.array(
                        sparse.coo_matrix((
                            maskArray,
                            (
                                self.__dict__['i%s'%self._coords[0]].flatten(),
                                self.__dict__['i%s'%self._coords[1]].flatten()
                            )
                        ), 
                            shape=self.imgShape).toarray()
                    )
                self.mask_img[self.mask_img!=0]=1
                ones_mask = np.array(
                        sparse.coo_matrix((
                            np.ones_like(maskArray),
                            (
                                self.__dict__['i%s'%self._coords[0]].flatten(),
                                self.__dict__['i%s'%self._coords[1]].flatten()
                            )
                        ), 
                            shape=self.imgShape).toarray()
                    )
                self.mask_ones = ones_mask
                self.mask_img[ones_mask==0]=1
                self.mask_img = self.mask_img.astype(int)

    #THIS NEEDS DEBUGGING.....MASKED ARRAY? ONLY DIRECT DATA
    def process(self, data):
        #already have dict w/ data
        if  isinstance(data, dict):
            if 'data' not in data or 'row' not in data or 'col' not in data:
                print('dictionay, but does not have sparsified data', data)
                return
            #
            # should ideally specify the output shape if possible.
            #
            #img = sparse.coo_matrix((d.flatten(), (ix.flatten(), iy.flatten())), shape=outShape).toarray()
            if 'tile' not in data or max(data['tile'])==0:
                aduData = data['data']
                rowData = data['row'].astype(int)[aduData>0]
                colData = data['col'].astype(int)[aduData>0]
                if self.imgShape is not None:
                    data = sparse.coo_matrix((aduData[aduData>0],(rowData, colData)), 
                                                 shape=(self.imgShape[1], self.imgShape[0])).toarray()
                else:
                    data = sparse.coo_matrix((aduData[aduData>0],(rowData, colData))).toarray()
            else:
                rowData = data['row'].astype(int)[data['data']>0]
                colData = data['col'].astype(int)[data['data']>0]
                tileData = data['tile'].astype(int)[data['data']>0]
                taduData = data['data'].copy()[data['data']>0]
                dataList=[]
                for itile in range(self._pedShape[-3]):
                    if len(taduData[tileData==itile])==0:
                        dataList.append(np.zeros((self._pedShape[-2], self._pedShape[-1])))
                    else:
                        tileD = taduData[tileData==itile]
                        tileR = rowData[tileData==itile]
                        tileC = colData[tileData==itile]
                        dataTile = sparse.coo_matrix((tileD,(tileR, tileC)), 
                                                     shape=(self._pedShape[-2], self._pedShape[-1])).toarray()
                        dataList.append(dataTile)
                #    print('itile data ',taduData[tileData==itile])
                #data = np.array([ sparse.coo_matrix(data['data'][data['tile']==itile],(data['row'][data['tile']==itile],data['col'][data['tile']==itile])).toarray() for itile in range(1+max(data['tile'])) ])
                data = np.array(dataList)

        elif not isinstance(data, np.ndarray):
            print('cannot work with this data, need array or sparsified array')

        if self.correction is not None:
            data /= self.correction

        retDict={}
        if self._coords is None:
            if isinstance(data, np.ma.masked_array):
                if data.dtype==np.uint16:
                    data = data.filled(fill_value=0)
                else:
                    #data from droplet2Func already has fill_value defined.
                    try:
                        data = data.filled(data, fill_value=np.nan)
                    except:
                        data = data.filled(data)

            if isinstance(data, np.matrix):
                data = np.asarray(data)
            self.dat = data
            retDict['img']=data

        #DEBUG ME: masked pixels should be in extra pixel!
        ##data will be used as weights in bincount. Thus masked pixels should be set to 0.
        #if isinstance(data, np.ma.masked_array):            
        #    print 'fill'
        #    data = data.filled(data, fill_value=0)
        #    print 'filled'
        elif len(self._coords)==2:
            ##using the sparse array is slower (1ms versus 0.55ms for bincount)
            ##also, mapping of multiple pixels into one final pixel will not be normalizable
            ##as a note the normalization costs about 0.2ms
            #data2d = sparse.coo_matrix((data.flatten(),(self.__dict__['i%s'%self._coords[0]],self.__dict__['i%s'%self._coords[1]])), shape=self.imgShape).toarray()            
            #retDict['img_sparse'] = np.array(data2d)
            I = np.bincount(self._multidim_idxs, weights=data.flatten(), minlength=int(self._n_multidim_idxs))
            I = np.reshape(I[:self._n_multidim_idxs], self._n_coordTuple)
            if self._npix_div is not None:
                Inorm=I*self._npix_div
                retDict['img'] = Inorm
            else:
                retDict['img'] = I
            #cast to same type that input array was.
            retDict['img'] = retDict['img'].astype(data.dtype)
            self.dat = retDict['img']

        elif len(self._coords)==1: #this might be a special case of the multi dim thing....
            data = np.bincount(data, self.__dict__['i%s'%self._coords[0]])
            retDict['img_1d']=data/self.npix
            self.dat = retDict['img_1d']
        else:
            I = np.bincount(self._multidim_idxs, weights=data, minlength=self._n_multidim_idxs)
            I = I[:self._n_multidim_idxs]
            retDict['img_n'] = np.reshape(I, self._n_coordTuple)
            self.dat = retDict['img_n']

        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                retDict['%s_%s'%(k,kk)] = subfuncResults[k][kk]

        return retDict

