import os
import copy
import numpy as np
from smalldata_tools.DetObject import event
from smalldata_tools.DetObject import DetObjectFunc
from future.utils import iteritems
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
try:
    basestring
except NameError:
    basestring = str

import psana
#from collections import Counter

def DetObject_lcls2(srcName, run, **kwargs):
    print('Getting the detector for: ',srcName)
    det = None
    try:
        det = run.Detector(srcName)
    except:
        return None
    det.alias = srcName
    detector_classes = {
        'opal': OpalObject_lcls2,
        'hsd':  HsdObject,
        'pv':   PVObject_lcls2
    }
    cls = detector_classes[det._dettype]
    return cls(det, run, **kwargs)
    ##should throw an exception here.
    #return None

class DetObjectClass_lcls2(object):
    def __init__(self,det, run, **kwargs):#name=None, common_mode=None, applyMask=0):
        self.det=det
        self._detid=det._detid
        self._name = kwargs.get('name', self.det._det_name)#srcName)

        self.run=run
        self._storeSum = {}
        self.applyMask = kwargs.get('applyMask', 0)

        self.dataAccessTime=0.

    def params_as_dict(self):
        """returns parameters as dictionary to be stored in the hdf5 file (once/file)"""
        parList =  {key:self.__dict__[key] for key in self.__dict__ if (key[0]!='_' and isinstance(getattr(self,key), (basestring, int, float, np.ndarray)) )}
        parList.update({key: np.array(self.__dict__[key]) for key in self.__dict__ if (key[0]!='_' and isinstance(getattr(self,key), tuple) and isinstance(getattr(self,key)[0], (basestring, int, float, np.ndarray)))}) 
        parList.update({key: np.array(self.__dict__[key]) for key in self.__dict__ if (key[0]!='_' and isinstance(getattr(self,key), list) and isinstance(getattr(self,key)[0], (basestring, int, float, np.ndarray)))}) 
        #add parameters of function to dict with composite keyt(sf._name, key)
        subFuncs = [ self.__dict__[key] for key in self.__dict__ if isinstance(self.__dict__[key], DetObjectFunc) ]
        for sf in subFuncs:
            sfPars = sf.params_as_dict()
            parList.update({'%s__%s'%(sf._name,key): value for key,value in iteritems(sfPars) if (key[0]!='_' and isinstance(value, (basestring, int, float, np.ndarray, tuple)))})

        #remKeys = [key for key in self.__dict__ if (key not in parList)]
        #print('DEBUG: keys which are not parameters:',remKeys)
        #for k in remKeys:
        #    if k[0]!='_':
        #        print k, self.__dict__[k]
        return parList

    #here we should get masks from the calibstore once we learn how to.
    def _getMasks(self):               
        self.mask = None
        self.cmask = None

    def _applyMask(self):
        try:
            if self.applyMask==1:
                self.evt.dat[self.mask==0]=0
            if self.applyMask==2:
                self.evt.dat[self.cmask==0]=0
        except:
            print('Could not apply mask to data for detector ',self.name)

    def storeSum(self, sumAlgo=None):
        if sumAlgo is not None:
            self._storeSum[sumAlgo]=None
        else:
            return self._storeSum

    def setPed(self, ped):
        self.ped = ped

    def setMask(self, mask):
        self.mask = mask

    def setcMask(self, mask):
        self.cmask = np.amin(np.array([self.mask,mask]),axis=0)

    def setGain(self, gain):
        """
        set a local gain. 
        This file is supposed to be applied on top of whatever corrections DetObject will apply, given the common mode
        """ 
        self.local_gain = gain

    def getData(self, evt):
        try:
            getattr(self, 'evt')
        except:
            self.evt = event()
        self.evt.dat = None

    def addFunc(self, func):
        func.setFromDet(self) #pass parameters from det (rms, geometry, .....)
        try:
            func.setFromFunc() #pass parameters from itself to children (rms, bounds, .....)            
        except:
            print('Failed to pass parameters to children of ',func._name)
        self.__dict__[func._name] = func

    def processFuncs(self):
        if self.evt.dat is None:
            print('This event has no data to be processed for %s'%self._name)
            return 
        for func in [self.__dict__[k] for k in  self.__dict__ if isinstance(self.__dict__[k], DetObjectFunc)]:
            try:
                retData=func.process(self.evt.dat)
                self.evt.__dict__['_write_%s'%func._name] = retData
            except:
                print('Could not run function %s on data of detector %s of shape'%(func._name, self._name), self.evt.dat.shape)

    def processSums(self):
        for key in self._storeSum.keys():
            asImg=False
            thres=-1.e9
            for skey in key.split('_'):
                if skey.find('img')>=0:
                    asImg=True
                else:
                    if skey.find('thresADU')>=0:
                        thres=float(skey.replace('thresADU',''))
            
            if self.evt.dat is None:
                return
            dat_to_be_summed = self.evt.dat
            if thres>1e-9:
                dat_to_be_summed[self.evt.dat<thres]=0.

            if key.find('nhits')>=0:
                dat_to_be_summed[dat_to_be_summed>0]=1
            
            if key.find('square')>=0:
                dat_to_be_summed = np.square(dat_to_be_summed)

            if asImg:
                try:
                    dat_to_be_summed = self.det.image(self.run,dat_to_be_summed)
                except:
                    pass

            if self._storeSum[key] is None:
                if dat_to_be_summed is not None:
                    self._storeSum[key] = dat_to_be_summed.copy()
            else:
                try:
                    self._storeSum[key] += dat_to_be_summed
                except:
                    print('could not add ',dat_to_be_summed)
                    print('could not to ',self._storeSum[key])
            #print('%s'%key, self._storeSum[key] )

class CameraObject_lcls2(DetObjectClass_lcls2): 
    def __init__(self, det,run,**kwargs):
        super(CameraObject_lcls2, self).__init__(det,run, **kwargs)
        self._common_mode_list = [0,-1, 30] #none, raw, calib
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list and type(self) is CameraObject:
            print('Common mode %d is not an option for a CameraObject, please choose from: '%self.common_mode, self._common_mode_list)
        self.pixelsize=[25e-6]
        self.isGainswitching=False

        self.ped = None
        self.rms = None
        self.mask = None
        #self.det.calibconst['pop_rbfs'][1] return meta data for the calib data.
        #self.rms = self.det.rms(run)
        #self.gain_mask = self.det.gain_mask(run)
        #self.gain = self.det.gain(run)
        #self.common_mode_pars=self.det.common_mode(run)
        self.local_gain = None
        self._getImgShape() #sets self.imgShape
        #self._getMasks() #sets mask, cmask, statusMask
        self._gainSwitching = False
        self.x = None
        self.y = None

    def getData(self, evt):
        super(CameraObject_lcls2, self).getData(evt)

    #this is not important until we get to tiled detectors w/ geometry.
    def _getImgShape(self):
        self.imgShape = None


class OpalObject_lcls2(CameraObject_lcls2): 
    def __init__(self, det, run,**kwargs):
        super(OpalObject_lcls2, self).__init__(det, run, **kwargs)

    def getData(self, evt):
        super(OpalObject_lcls2, self).getData(evt)
        if self.common_mode<0:
            self.evt.dat = self.det.raw.raw(evt)
        elif self.common_mode==0: #we need to figure out how to do this. Don't implement for, return raw
            self.evt.dat = self.det.raw.raw(evt)
            #if not self._gainSwitching:
            #    if self.ped is not None:
            #        self.evt.dat = self.det.raw_data(evt)-self.ped
            #    else:
            #        self.evt.dat = self.det.raw_data(evt)
            #    self._applyMask()
            #    if (self.gain is not None) and self.gain.std() != 0 and self.gain.mean() != 1.:
            #        if self.gain.shape == self.evt.dat.shape:
            #            self.evt.dat*=self.gain   
        elif self.common_mode%100==30:
            self.evt.dat = self.det.raw.calib(evt)

        #override gain if desired
        if self.local_gain is not None and self.common_mode in [0,30] and self._gainSwitching is False and self.local_gain.shape == self.evt.dat.shape:
            self.evt.dat*=self.local_gain   #apply own gain

class PVObject_lcls2(CameraObject_lcls2): 
    def __init__(self, det, run,**kwargs):
        super(PVObject_lcls2, self).__init__(det, run, **kwargs)

    def getData(self, evt):
        super(PVObject_lcls2, self).getData(evt)
        if self.common_mode<0:
            self.evt.dat = self.det.raw.value(evt)
        elif self.common_mode==0: #we need to figure out how to do this. Don't implement for, return raw
            self.evt.dat = self.det.raw.value(evt)
        elif self.common_mode%100==30:
            self.evt.dat = self.det.raw.value(evt)


class WaveformObject_lcls2(DetObjectClass_lcls2): 
    def __init__(self, det,run,**kwargs):
        super(WaveformObject_lcls2, self).__init__(det,run, **kwargs)
        self.common_mode = kwargs.get('common_mode', -1)
        self.rms = None
        self.mask = None
        self.wfx = None
        #if self.det.dettype == 16: #acqiris
        #  cfg = env.configStore().get(psana.Acqiris.ConfigV1, self._src)
        #  self.interval = [cfg.horiz().sampInterval()]
        #  self.delayTime = [cfg.horiz().delayTime()]
        #  self.nSamples =  [cfg.horiz().nbrSamples()]
        #  self.fullScale =[]
        #  self.offset = []

        #  for c in cfg.vert():
        #    self.fullScale.append(c.fullScale())
        #    self.offset.append(c.offset())
        self.gain = None
    def getData(self, evt):
        super(WaveformObject_lcls2, self).getData(evt)
        #self.evt.dat = self.det.raw.waveform(evt)
        #if self.evt.dat is None:
        #    self.wfx = self.det.wftime(evt)
        
class HsdObject(WaveformObject_lcls2): 
    def __init__(self, det,run,**kwargs):
        super(HsdObject, self).__init__(det,run, **kwargs)
        self.cidx = [k for k in self.det.raw._seg_chans()]

    def getData(self, evt):
        super(HsdObject, self).getData(evt)
        datadict = self.det.raw.waveforms(evt)
        #if self.channels == [-1]:
        #    self.channels = [k for k in datadict.keys()]
        if self.wfx is None:
            self.wfxlen = np.array([datadict[k]['times'].shape[0] for k in self.cidx])
            self.wfx = np.zeros(self.wfxlen.sum())
            startidx=0
            for ik,k in enumerate(self.cidx):
                self.wfx[startidx:(startidx+self.wfxlen[ik])] = datadict[k]['times']
                startidx+=self.wfxlen[ik]

        if self.evt.dat is None:
            startidx=0
            self.evt.dat = np.zeros(self.wfxlen.sum())
            for ik,k in enumerate(self.cidx):
                self.evt.dat[startidx:(startidx+self.wfxlen[ik])] = datadict[k][0]
                startidx+=self.wfxlen[ik]
            try:
                self.evt.dat = np.array(self.evt.dat)
            except:
                print('HsdObject: could not cast waveform times to array ',self.evt.dat)

