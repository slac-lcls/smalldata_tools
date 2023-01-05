import os
import copy
import numpy as np
from psana.pscalib.calib.MDBWebUtils import calib_constants
from smalldata_tools.DetObjectFunc import DetObjectFunc
from smalldata_tools.DetObjectFunc import event
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
        'epix10ka': Epix10kObject_lcls2,
        'opal': OpalObject_lcls2,
        'hsd':  HsdObject,
        'wave8':  Wave8Object,
        'pv':   PVObject_lcls2
    }
    #print('dettype: ', det._dettype)
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
        #we cannot save arraus of chars/strings. There should be a proper test instead of this explicit rejection of coords.
        #I can't see how to do that in a single line, but that is not a great reason...
        for sf in subFuncs:
            sfPars = sf.params_as_dict()
            parList.update({'%s__%s'%(sf._name,key): value for key,value in iteritems(sfPars) if (key[0]!='_' and isinstance(value, (basestring, int, float, np.ndarray, tuple)) and key.find('coords')<0)})

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

        #try calibconst...
        #detrawid = det.raw._uniqueid
        #self.peds = calib_constants(detrawid, exp=run.expt, ctype='pedestals', run=run.runnum)[0]
        #self.rms = calib_constants(detrawid, exp=run.expt, ctype='pixel_rms', run=run.runnum)[0]

        try:
            self.ped = det.raw._pedestals()
        except:
            self.ped = None
        try:
            self.rms = det.raw._rms()
        except:
            self.rms = None
        try:
            self.gain = det.raw._gain()
        except:
            self.gain = None
        try:
            self.mask = det.raw._mask(calib=False, status=True, edges=True)
            self.cmask = det.raw._mask(calib=True, status=True, edges=True)
        except:
            self.mask = None
            self.cmask = None
        #self.det.calibconst['pop_rbfs'][1] return meta data for the calib data.
        #self.rms = self.det.rms(run)
        #self.gain_mask = self.det.gain_mask(run)
        #self.gain = self.det.gain(run)
        #self.common_mode_pars=self.det.common_mode(run)
        self.local_gain = None
        self._getImgShape() #sets self.imgShape
        #self._getMasks() #sets mask, cmask, statusMask
        self._gainSwitching = False
        try:
            self.x, self.y, self.z = det.raw._pixel_coords(do_tilt=True, cframe=0)
            self.x = self.x.squeeze()
            self.y = self.y.squeeze()
            self.z = self.z.squeeze()
        except:
            self.x, self.y, self.z = None, None, None

    def getData(self, evt):
        super(CameraObject_lcls2, self).getData(evt)

    #this is not important until we get to tiled detectors w/ geometry.
    def _getImgShape(self):
        self.imgShape = None

class OpalObject_lcls2(CameraObject_lcls2): 
    def __init__(self, det, run, **kwargs):
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

class TiledCameraObject_lcls2(CameraObject_lcls2): 
    def __init__(self, det, run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(TiledCameraObject_lcls2, self).__init__(det,run, **kwargs)
        try:
            self.ix, self.iy = det.raw._pixel_coord_indexes()
        except:
            if rank==0:
                print('failed to get geometry info, likely because we do not have a geometry file')
            self.ix=self.x #need to change this so ix & iy are integers!
            self.iy=self.y
        self.ix = self.ix.squeeze()
        self.iy = self.iy.squeeze()
        self._needsGeo=True #FIX ME: not sure it should be here.            
    def getData(self, evt):
        super(TiledCameraObject_lcls2, self).getData(evt)

class Epix10kObject_lcls2(TiledCameraObject_lcls2): 
    def __init__(self, det, run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(Epix10kObject_lcls2, self).__init__(det,run, **kwargs)
        self._common_mode_list = [80, 0, -1, -2, 30, 84, 85] # calib-noCM, ped sub, raw, raw_gain, calib, calib-CMrow, calib-CMrowcol
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print('Common mode %d is not an option for as Epix detector, please choose from: '%self.common_mode, self._common_mode_list)
        self.pixelsize=[100e-6]
        self.isGainswitching=True

        if self.rms is None or self.rms.shape!=self.ped.shape:
            self.rms=np.ones_like(self.ped)
        try:
            #ok, this call does NOT work (yet) for LCLS2. Needs an event...
            #self.imgShape=det.raw.image(run.runnum, self.ped[0]).shape
            self.imgShape=(self.ix.max(), self.iy.max())
        except:
            if len(self.ped[0].squeeze().shape)==2:
                self.imgShape=self.ped[0].squeeze().shape
            else:
                self.imgShape=None
        self._gainSwitching = True                

    def getData(self, evt):
        super(Epix10kObject_lcls2, self).getData(evt)
        mbits=0 #do not apply mask (would set pixels to zero)
        #mbits=1 #set bad pixels to 0
        if self.common_mode<0:
          if self.common_mode==-2:
            self.evt.dat = self.evt.dat
          else:
            self.evt.dat = self.evt.dat&0x3fff
          #self.evt.gainbit = (self.evt.dat&0xc000>0)
        elif self.common_mode==0:
            ##########
            ### FIX ME epix10ka
            #will need to read gain bit from data and use right pedestal.
            #will hopefully get calib function for this.
            ##########
            if len(self.ped.shape)>3:
                self.evt.dat = (self.det.raw.raw(evt)&0x3fff)-self.ped[0]
            else:
                self.evt.dat = (self.det.raw.raw(evt)&0x3fff)-self.ped
        elif self.common_mode==30:
            self.evt.dat = self.det.raw.calib(evt)
        elif self.common_mode==80:
            self.evt.dat = self.det.raw.calib(evt, cmpars=(7,0,100))
        elif self.common_mode==84:
            self.evt.dat = self.det.raw.calib(evt, cmpars=(7,2,100,10))
        elif self.common_mode==85:
            self.evt.dat = self.det.raw.calib(evt, cmpars=(7,3,100,10))


        #override gain if desired -- this looks like CsPad.
        if self.local_gain is not None and self.local_gain.shape == self.evt.dat.shape and self.common_mode in [1,5,55,10]:
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
        self.ped = None
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

class Wave8Object(WaveformObject_lcls2): 
    def __init__(self, det,run,**kwargs):
        super(Wave8Object, self).__init__(det,run, **kwargs)
        self._chan_names = [name for name in dir(det.raw) if name[0]!='_']
        self.chan_names = ' '.join(self._chan_names)

    def getData(self, evt):
        super(Wave8Object, self).getData(evt)
        vetolist=['config']
        self.evt.dat = [ getattr(self.det.raw, name)(evt) for name in self._chan_names if name not in vetolist]
        try:
            self.evt.dat = np.array(self.evt.dat)
        except:
            print('Wave8: could not cast waveform times to array ',self.evt.dat)

