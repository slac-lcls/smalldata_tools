import os
import copy
import numpy as np
import tables

from smalldata_tools.utilities import cm_epix
from smalldata_tools.utilities import cm_uxi
from smalldata_tools.read_uxi import getDarks
from smalldata_tools.DetObjectFunc import DetObjectFunc
from smalldata_tools.DetObjectFunc import event
from future.utils import iteritems
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
try:
    basestring
except NameError:
    basestring = str
from read_uxi import read_uxi, get_uxi_timestamps, getDarks
#for common mode.
#from read_uxi import getUnbonded, getZeroPeakFit 

import Detector as DetTypes
import psana
#from collections import Counter

def DetObject(srcName, env, run, **kwargs):
    print('Getting the detector for: ',srcName)
    det = None
    try:
        det = psana.Detector(srcName)
    except:
        try:
            det = run.Detector(srcName)
        except:
            pass
    if det is None:
        if srcName.find('uxi')>=0:
            return UxiObject(run, **kwargs)
        else:
            return None
    det.alias = srcName

    detector_classes = {
        1: CsPad2MObject,
        2: CsPadObject,
        5: PulnixObject,
        6: OpalObject,
        7: OpalObject,
        8: OpalObject,
        9: OpalObject,
       10: OrcaObject, #fee_spec
       13: EpixObject,
       15: AndOrObject,
       16: AcqirisObject,
       19: RayonixObject,
       23: FliObject,
       26: JungfrauObject,
       27: ZylaObject,
       28: ControlsCameraObject,
       29: Epix10kObject, #this is a small epix10k
       30: IcarusObject,
       32: Epix10k2MObject, #this is a 2M
       33: Epix10k2MObject, #this is a quad-lookds like a 2M
       34: CameraObject,
       36: iStarObject,
       37: AlviumObject,
       98: OceanOpticsObject,
       99: ImpObject,
    }

    detector_lookup={DetTypes.GenericWFDetector.GenericWFDetector: GenericWFObject}

    if hasattr(det,'dettype'):
        cls = detector_classes[det.dettype]
    else:
        cls = detector_lookup[type(det)]

    return cls(det, env, run, **kwargs)
    ##should throw an exception here.
    #return None

class DetObjectClass(object):
    # These constants are essentially the inverse of the detector_classes dict above.
    # Where are these really defined?  PSCalib!  In src/GlobalUtils.py.
    CsPad2M = 1           # CSPAD
    CsPad = 2             # CSPAD2X2
    Pulnix = 5
    Opal = 6              # OPAL1000
    Opal2k = 7
    Opal4k = 8
    Opal8k = 9
    Orca = 10
    Epix = 13
    AndOr = 15
    Acqiris = 16
    Rayonix = 19
    Fli = 23
    Jungfrau = 26
    Zyla = 27
    ControlsCamera = 28
    Epix10k = 29
    Icarus = 30
    Epix10k2M = 32
    Epix10k2M_quad = 33
    Camera = 34           # STREAK?
    iStar = 36
    Alvium = 37
    OceanOptics = 98      # ?
    Imp = 99              # ?
    
    def __init__(self,det,env,run, **kwargs):#name=None, common_mode=None, applyMask=0):
        self.det=det
        self._src=det.source
        self._name = kwargs.get('name', self.det.alias)#srcName)

        self.run=run
        self._storeSum = {}
        self.applyMask = kwargs.get('applyMask', 0)
        self.maskCentral = kwargs.get('maskCentral', True)

        self.dataAccessTime=0.
        self.rawShape = None

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

    def _getImgShape(self):
        if self.ped is not None:
            self.imgShape = self.ped.shape
            try:
                self.imgShape = self.det.image(self.run, self.ped).shape
            except:
                if len(self.ped.shape)>2:
                    try:
                        self.imgShape = self.det.image(ped[0], self.run).shape
                    except:         
                        print('multi dim pedestal & image function does do nothing: multi gain detector.....')
                        self.imgShape=self.ped.shape[1:]
        else:
            self.imgShape = None

    #
    # In what follows, we should call this to set the raw shape whenever we set self.evt.dat.
    # It would be nice to do this once at the beginning of time, but probably better not to
    # search for an event with data.
    #
    def _getRawShape(self):
        if self.rawShape is None and self.evt is not None and self.evt.dat is not None:
            self.rawShape = self.evt.dat.shape

    def _get_coords_from_ped(self):
        self.x = np.arange(0,self.ped.shape[-2]*self.pixelsize[0]*1e6, self.pixelsize[0]*1e6)
        if len(self.pixelsize)==1:
            self.y = np.arange(0,self.ped.shape[-1])*self.pixelsize[0]*1e6
        else:
            self.y = np.arange(0,self.ped.shape[-1])*self.pixelsize[1]*1e6
        self.x, self.y = np.meshgrid(self.x, self.y, indexing='ij')
        self.ix = self.x.copy()
        self.ix = self.ix - np.min(self.ix)
        self.ix = (self.ix/np.max(self.ix)*self.imgShape[0]).astype(int)
        self.iy = self.y.copy()
        self.iy = self.iy - np.min(self.iy)
        self.iy = (self.iy/np.max(self.iy)*self.imgShape[1]).astype(int)

    def _getMasks(self):               
        try:
            self.statusMask = self.det.mask(self.run, status=True)
            self.mask = self.det.mask(self.run, unbond=True, unbondnbrs=True, status=True,  edges=True, central=self.maskCentral)

            if rank==0 and self.mask is not None:
                print('masking %d pixel (status & edge,..) of %d'%(np.ones_like(self.mask).sum()-self.mask.sum(), np.ones_like(self.mask).sum()))
            self.cmask = self.det.mask(self.run, unbond=True, unbondnbrs=True, status=True,  edges=True, central=self.maskCentral,calib=True)
            if self.cmask is not None and self.cmask.sum()!=self.mask.sum() and rank==0:
                print('found user mask, masking %d pixel'%(np.ones_like(self.mask).sum()-self.cmask.sum()))
        except:
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


    #make this a private function?
    def getMask(self, ROI):
        ROI = np.array(ROI)
        #print 'DEBUG getMask: ',ROI.shape
        if ROI.shape != (2,2):
            return np.ones_like(self.ped) 

        mask_roi=np.zeros(self.imgShape)#<-- no, like image. Need image size.
        #print 'DEBUG getMask: img shape ',self.imgShape, self.ped.shape
        mask_roi[ROI[0,0]:ROI[0,1],ROI[1,0]:ROI[1,1]]=1
        if self._needsGeo:
            mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(self.ix,self.iy)] )
        else:
            mask_r_nda = mask_roi
        #print 'mask from rectangle (shape):',mask_r_nda.shape
        return mask_r_nda

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
            print('This event has no data to be processed')
            return 
        for func in [self.__dict__[k] for k in  self.__dict__ if isinstance(self.__dict__[k], DetObjectFunc)]:
            if func._proc == False: # so that we don't process all events in cube
                continue
            try:
                retData=func.process(self.evt.dat)
                self.evt.__dict__['_write_%s'%func._name] = retData
            except Exception as E:
                print('Could not run function %s on data of detector %s of shape'%(func._name, self._name), self.evt.dat.shape)
                print('Error message: {}'.format(E))

    def processSums(self, **kwargs):
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
            dat_to_be_summed = self.evt.dat.copy()
            if thres>-1e9:
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

            #loop over kwargs & look for them in the key
            for k,v in kwargs.items():
                # the key indicates a given selection is to be respected
                if key.find(k)>=0:
                    # is the selection is not true, set array to zeros
                    if not v: dat_to_be_summed = np.zeros_like(dat_to_be_summed)
            
            if self._storeSum[key] is None:
                self._storeSum[key] = dat_to_be_summed.copy()
            else:
                try:
                    self._storeSum[key] += dat_to_be_summed
                except:
                    print('could not add ',dat_to_be_summed)
                    print('could not to ',self._storeSum[key])

class WaveformObject(DetObjectClass): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(WaveformObject, self).__init__(det,env,run, **kwargs)
        self.common_mode = kwargs.get('common_mode', -1)
        self.ped = None
        self.rms = None
        self.mask = None
        self.wfx = None
        self.gain = None

    def getData(self, evt):
        super(WaveformObject, self).getData(evt)
        if self.wfx is None:
            self.wfx = self.det.wftime(evt)

class GenericWFObject(WaveformObject): 
    def __init__(self, det,env,run,**kwargs):
        super(GenericWFObject, self).__init__(det,env,run, **kwargs)

    def getData(self, evt):
        super(GenericWFObject, self).getData(evt)
        self.evt.dat = self.det.raw(evt)
        if isinstance(self.evt.dat, list):
            self.evt.dat = np.array(self.evt.dat)
        self._getRawShape()

class AcqirisObject(WaveformObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(AcqirisObject, self).__init__(det,env,run, **kwargs)
        cfg = env.configStore().get(psana.Acqiris.ConfigV1, self._src)
        self.interval = [cfg.horiz().sampInterval()]
        self.delayTime = [cfg.horiz().delayTime()]
        self.nSamples =  [cfg.horiz().nbrSamples()]
        self.fullScale =[]
        self.offset = []
        
        for c in cfg.vert():
            self.fullScale.append(c.fullScale())
            self.offset.append(c.offset())

    def getData(self, evt):
        super(AcqirisObject, self).getData(evt)
        self.evt.dat = self.det.waveform(evt)
        self._getRawShape()

class OceanOpticsObject(WaveformObject): 
    def __init__(self, det,env,run,**kwargs):
        super(OceanOpticsObject, self).__init__(det,env,run, **kwargs)

    def getData(self, evt):
        super(OceanOpticsObject, self).getData(evt)
        self.evt.dat = self.det.intensity(evt)
        if self.wfx is None:
            self.wfx = self.det.wavelength(evt)
        self._getRawShape()

class ImpObject(WaveformObject): 
    def __init__(self, det,env,run,**kwargs):
        super(ImpObject, self).__init__(det,env,run, **kwargs)

    def getData(self, evt):
        super(ImpObject, self).getData(evt)
        self.evt.dat = self.det.waveform(evt)
        self._getRawShape()
        
class CameraObject(DetObjectClass): 
    def __init__(self, det,env,run,**kwargs):
        super(CameraObject, self).__init__(det,env,run, **kwargs)
        self._common_mode_list = [0,-1, 30] #none, raw, calib
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list and type(self) is CameraObject:
            print('Common mode %d is not an option for a CameraObject, please choose from: '%self.common_mode, self._common_mode_list)
        self.pixelsize=[25e-6]
        self.isGainswitching=False

        self.rms = self.det.rms(run)
        self.ped = self.det.pedestals(run)
        self.gain_mask = self.det.gain_mask(run)
        self.gain = self.det.gain(run)
        
        self.common_mode_pars=self.det.common_mode(run)
        self.local_gain = None
        self._getImgShape() #sets self.imgShape
        self._getMasks() #sets mask, cmask, statusMask
        self._gainSwitching = False
        try:
            self.x = self.det.coords_x(run)
            self.ix=self.x
            self.y = self.det.coords_y(run)
            self.iy=self.y
            self.z = self.det.coords_z(run)
        except:
            self.x = None
            self.y = None

    def getData(self, evt):
        super(CameraObject, self).getData(evt)
        if self.common_mode<0:
            self.evt.dat = self.det.raw_data(evt)
        elif self.common_mode==0:
            if not self._gainSwitching:
                if self.ped is not None:
                    self.evt.dat = self.det.raw_data(evt)-self.ped
                else:
                    self.evt.dat = self.det.raw_data(evt)
                self._applyMask()
                if (self.gain is not None) and self.gain.std() != 0 and self.gain.mean() != 1.:
                    if self.gain.shape == self.evt.dat.shape:
                        self.evt.dat*=self.gain   
        elif self.common_mode%100==30:
            self.evt.dat = self.det.calib(evt)
        #override gain if desired
        if self.local_gain is not None and self.common_mode in [0,30] and self._gainSwitching is False and self.local_gain.shape == self.evt.dat.shape:
            self.evt.dat*=self.local_gain   #apply own gain
        self._getRawShape()

class OpalObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(OpalObject, self).__init__(det,env,run, **kwargs)
        if env.configStore().get(psana.TimeTool.ConfigV2, self._src) is not None and env.configStore().get(psana.TimeTool.ConfigV2, self._src).write_image()==0:
            self.ped=np.array([[-1,-1],[-1,-1]])
        else:
            fexCfg = env.configStore().get(psana.Camera.FrameFexConfigV1, self._src)
            if fexCfg.forwarding() == fexCfg.Forwarding.values[2]: #make sure we are only doing ROI
                self.ped = np.zeros([fexCfg.roiEnd().column()-fexCfg.roiBegin().column(), fexCfg.roiEnd().row()-fexCfg.roiBegin().row()]).T
            if self.ped is None or self.ped.shape==(0,0): #this is the case for e.g. the xtcav recorder but can also return with the DAQ. Assume Opal1k for now.
                #if srcName=='xtcav':
                #  self.ped = np.zeros([1024,1024])
                if det.dettype == DetObjectClass.Opal:
                    self.ped = np.zeros([1024,1024])
                else:
                    self.ped = None
        if det.dettype == DetObjectClass.Opal:
            self.imgShape = self.ped.shape
            if self.x is None:
                self._get_coords_from_ped()
            if self.mask is None or self.mask.shape!=self.imgShape:
                self.mask = np.ones(self.imgShape)
                self.cmask = np.ones(self.imgShape)

class ZylaObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(ZylaObject, self).__init__(det,env,run, **kwargs)
        zylaCfg = env.configStore().get(psana.Zyla.ConfigV1, self._src)
        #needs to be this way around to match shape of data.....
        self.imgShape = (zylaCfg.numPixelsY(), zylaCfg.numPixelsX())
        if self.ped is None:
            self.ped = np.zeros([zylaCfg.numPixelsY(), zylaCfg.numPixelsX()])
        if self.imgShape != self.ped.shape:            
            print('Detector %s does not have a pedestal taken for run %s, this might lead to issues later on!'%(self._name,run))
            if self.common_mode != -1:
                print('We will use the raw data for Detector %s in run %s'%(self._name,run))
            self.common_mode = -1
        self.imgShape = self.ped.shape
        if self.x is None:
            self._get_coords_from_ped()
        if self.mask is None or self.mask.shape!=self.imgShape:
            self.mask = np.ones(self.imgShape)
            self.cmask = np.ones(self.imgShape)

class iStarObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(iStarObject, self).__init__(det,env,run, **kwargs)
        iStarCfg = env.configStore().get(psana.iStar.ConfigV1, self._src)
        #needs to be this way around to match shape of data.....
        self.imgShape = (iStarCfg.numPixelsY(), iStarCfg.numPixelsX())
        if self.ped is None:
            self.ped = np.zeros([iStarCfg.numPixelsY(), iStarCfg.numPixelsX()])
        if self.imgShape != self.ped.shape:            
            print('Detector %s does not have a pedestal taken for run %s, this might lead to issues later on!'%(self._name,run))
            if self.common_mode != -1:
                print('We will use the raw data for Detector %s in run %s'%(self._name,run))
            self.common_mode = -1
        self.imgShape = self.ped.shape
        if self.x is None:
            self._get_coords_from_ped()
        if self.mask is None or self.mask.shape!=self.imgShape:
            self.mask = np.ones(self.imgShape)
            self.cmask = np.ones(self.imgShape)

class AlviumObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(AlviumObject, self).__init__(det,env,run, **kwargs)
        alviumCfg = env.configStore().get(psana.Vimba.AlviumConfigV1, self._src)
        #needs to be this way around to match shape of data.....
        self.imgShape = (alviumCfg.numPixelsY(), alviumCfg.numPixelsX())
        if self.ped is None:
            self.ped = np.zeros([alviumCfg.numPixelsY(), alviumCfg.numPixelsX()])
        if self.imgShape != self.ped.shape:            
            print('Detector %s does not have a pedestal taken for run %s, this might lead to issues later on!'%(self._name,run))
            if self.common_mode != -1:
                print('We will use the raw data for Detector %s in run %s'%(self._name,run))
            self.common_mode = -1
        self.imgShape = self.ped.shape
        if self.x is None:
            self._get_coords_from_ped()
        if self.mask is None or self.mask.shape!=self.imgShape:
            self.mask = np.ones(self.imgShape)
            self.cmask = np.ones(self.imgShape)

class ControlsCameraObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(ControlsCameraObject, self).__init__(det,env,run, **kwargs)
        camrecCfg = env.configStore().get(psana.Camera.ControlsCameraConfigV1, self._src)
        self.imgShape = (camrecCfg.height(), camrecCfg.width())
        if self.ped is None:
            self.ped = np.zeros([camrecCfg.height(), camrecCfg.width()])
        if self.imgShape != self.ped.shape:            
            self.common_mode = -1
        pedImg = self.det.image(run, self.ped)
        self.imgShape = self.ped.shape
        if self.x is None:
            self._get_coords_from_ped()
        if self.mask is None:
            self.mask = np.ones(self.imgShape)
            self.cmask = np.ones(self.imgShape)

class PulnixObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(PulnixObject, self).__init__(det,env,run, **kwargs)
        yag2Cfg = env.configStore().get(psana.Pulnix.TM6740ConfigV2,self._src)
        self.ped = np.zeros([yag2Cfg.Row_Pixels, yag2Cfg.Column_Pixels])
        self.imgShape = self.ped.shape
        if self.x is None:
            self._get_coords_from_ped()
        if self.mask is None or self.mask.shape!=self.imgShape:
            self.mask = np.ones(self.imgShape)
            self.cmask = np.ones(self.imgShape)

class OrcaObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(OrcaObject, self).__init__(det,env,run, **kwargs)
        oCfg = env.configStore().get(psana.Orca.ConfigV1,self._src)
        try:
            ofCfg = env.configStore().get(psana.Camera.FrameFexConfigV1,self._src)
            self.ped = np.zeros([ofCfg.roiEnd().row()-ofCfg.roiBegin().row(),
                                 ofCfg.roiEnd().column()-ofCfg.roiBegin().column()])
        except:
            ofCfg = None
            self.ped = np.zeros([oCfg.Row_Pixels, oCfg.Column_Pixels])
        self.imgShape = self.ped.shape
        if self.x is None:
            self._get_coords_from_ped()
        if self.mask is None or self.mask.shape!=self.imgShape:
            self.mask = np.ones(self.imgShape)
            self.cmask = np.ones(self.imgShape)

class AndOrObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(AndOrObject, self).__init__(det,env,run, **kwargs)
        andorCfg = env.configStore().get(psana.Andor.ConfigV2,self._src)
        self.imgShape = (andorCfg.numPixelsY(), andorCfg.numPixelsX())
        if self.ped is None:
            self.ped = np.zeros([andorCfg.numPixelsY(), andorCfg.numPixelsX()])
        if self.imgShape != self.ped.shape:            
            print('Detector %s does not have a pedestal taken for run %s, this might lead to issues later on!'%(self._name,run))
            if self.common_mode != -1:
                print('We will use the raw data for Detector %s in run %s'%(self._name,run))
            self.common_mode = -1
        self.imgShape = self.ped.shape
        if self.x is None:
            self._get_coords_from_ped()
        if self.mask is None or self.mask.shape!=self.imgShape:
            self.mask = np.ones(self.imgShape)
            self.cmask = np.ones(self.imgShape)

class FliObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        super(FliObject, self).__init__(det,env,run, **kwargs)
        fliCfg = env.configStore().get(psana.Fli.ConfigV1,self._src)
        self.imgShape = (fliCfg.numPixelsY(), fliCfg.numPixelsX())
        if self.ped is None:
            self.ped = np.zeros([fliCfg.numPixelsY(), fliCfg.numPixelsX()])
        if self.imgShape != self.ped.shape:            
            print('Detector %s does not have a pedestal taken for run %s, this might lead to issues later on!'%(self._name,run))
            if self.common_mode != -1:
                print('We will use the raw data for Detector %s in run %s'%(self._name,run))
            self.common_mode = -1
        self.imgShape = self.ped.shape
        if self.x is None:
            self._get_coords_from_ped()
        if self.mask is None or self.mask.shape!=self.imgShape:
            self.mask = np.ones(self.imgShape)
            self.cmask = np.ones(self.imgShape)

class TiledCameraObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(TiledCameraObject, self).__init__(det,env,run, **kwargs)
        try:
            ix, iy = self.det.indexes_xy(run)
            self.ix=np.array(ix)
            self.iy=np.array(iy)
        except:
            if rank==0:
                print('failed to get geometry info, likely because we do not have a geometry file')
            self.ix=self.x #need to change this so ix & iy are integers!
            self.iy=self.y
        self._needsGeo=True #FIX ME: not sure it should be here.            

        self.areas = None
        if not self.maskCentral:
            try:
                self.areas = self.det.areas(run)
            except:
                pass

    def getData(self, evt):
        super(TiledCameraObject, self).getData(evt)

class IcarusObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(IcarusObject, self).__init__(det,env,run, **kwargs)
        self._common_mode_list = [0,80,81,98,99,-1] #none, unb-calc, unb-applied, raw
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print('Common mode %d is not an option for an Icarus detector, please choose from: '%self.common_mode, self._common_mode_list)
        if self.x is None and self.ped is not None:
            self.x = np.arange(0,self.ped.shape[-2]*self.pixelsize[0], self.pixelsize[0])*1e6
            self.y = np.arange(0,self.ped.shape[-1]*self.pixelsize[0], self.pixelsize[0])*1e6
            self.y, self.x = np.meshgrid(self.y, self.x)
            self.x=np.array([self.x for i in range(self.ped.shape[0])])
            self.y=np.array([self.y for i in range(self.ped.shape[0])])
        if self.mask is None and self.ped is not None:
            self.mask = np.ones(self.ped.shape)
            self.cmask = np.ones(self.ped.shape)
        #common mode stuff
        self.cm_maskedROI =  kwargs.get('cm_maskedROI', None)
        if self.cm_maskedROI is not None and isinstance(self.cm_maskedROI, list):
            self.cm_maskedROI = np.array(self.cm_maskedROI).flatten()
        self.cm_photonThres =  kwargs.get('cm_photonThres', 30)
        self.cm_maskNeighbors =  kwargs.get('cm_maskNeighbors', 0)
        self.cm_maxCorr =  kwargs.get('cm_maxCorr', self.cm_photonThres)
        self.cm_minFrac =  kwargs.get('cm_minFrac', 0.05)
        self.cm_mask_unb = np.zeros([1024,512])
        self.cm_mask_unb[:,255:257]=1
        self.cm_mask_unb=self.cm_mask_unb.astype(bool)

    def getData(self, evt):
        super(IcarusObject, self).getData(evt)
        if self.common_mode >= 80 and self.common_mode < 100: #this here is some icarus stuff.
            self.evt.dat = self.det.calib(evt)

            if self.cm_maskedROI is not None:
                cm_mask_mean = [ np.nanmean(frame[self.cm_maskedROI[0]:self.cm_maskedROI[1],self.cm_maskedROI[2]:self.cm_maskedROI[3]]) for frame in self.evt.dat ]
                cm_mask_med = [ np.nanmedian(frame[self.cm_maskedROI[0]:self.cm_maskedROI[1],self.cm_maskedROI[2]:self.cm_maskedROI[3]]) for frame in self.evt.dat ]
                cm_mask_std = [ np.nanstd(frame[self.cm_maskedROI[0]:self.cm_maskedROI[1],self.cm_maskedROI[2]:self.cm_maskedROI[3]]) for frame in self.evt.dat ]
                self.evt.__dict__['env_cm_masked_med'] = np.array(cm_mask_med)
                self.evt.__dict__['env_cm_masked_mean'] = np.array(cm_mask_mean)
                self.evt.__dict__['env_cm_masked_std'] = np.array(cm_mask_std)
                dataFrame = np.array([df-corr for df, corr in zip(self.evt.dat,cm_mask_med)])

            #only subtract the masked ROI
            if self.common_mode==80:
                self.evt.dat=dataFrame
                return

            if self.common_mode==81:
                corrImg, cmDict = cm_uxi(dataFrame, self.cm_photonThres, self.cm_maxCorr, self.cm_minFrac, self.cm_maskNeighbors)
                for k, value in iteritems(cmDict):
                    self.evt.__dict__['env_%s'%k] = value
                self.evt.dat = corrImg

            cmVals=[]
            for tile in self.evt.dat:
                cmVals.append(np.nanmedian(tile[self.cm_mask_unb>0]))
                if self.common_mode==99:
                    tile-=cmVals[-1] #apply the common mode.
            self.evt.__dict__['write_cmUnb'] = cmVals
            self._applyMask()
            #gain is ignored here for now
            self._getRawShape()

class JungfrauObject(TiledCameraObject):
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(JungfrauObject, self).__init__(det,env,run, **kwargs)
        self._common_mode_list = [0, 7,71,72,-1, 30] #none, epix-style corr on row*col, row, col, raw, calib
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print('Common mode %d is not an option for Jungfrau, please choose from: '%self.common_mode, self._common_mode_list)
        self.pixelsize=[75e-6]
        self.isGainswitching=True
        self.offset=self.det.offset(run)
        self.pixel_status=self.det.status(run)
        try:
            self.imgShape = self.det.image(run, self.ped[0]).shape
        except:
            pass
        self._gainSwitching = True
        #self._common_mode_list.append()

    def getData(self, evt):
        super(JungfrauObject, self).getData(evt)
        mbits=0 #do not apply mask (would set pixels to zero)
        #mbits=1 #set bad pixels to 0
        if self.common_mode==0:
            self.evt.dat = self.det.calib(evt, cmpars=(7,0,10), mbits=mbits)
        elif self.common_mode%100==71:
            self.evt.dat = self.det.calib(evt, cmpars=(7,1,10,10), mbits=mbits) #correction in rows
        elif self.common_mode%100==72:
            self.evt.dat = self.det.calib(evt, cmpars=(7,2,10,10), mbits=mbits) #correction in columns
        elif self.common_mode%100==7:
            self.evt.dat = self.det.calib(evt, cmpars=(7,3,10,10), mbits=mbits) #correction in rows&columns

        #override gain if desired
        if self.local_gain is not None and self.local_gain.shape == self.evt.dat.shape and self.common_mode in [7,71,72,0]:
            self.evt.dat*=self.local_gain   #apply own gain

        #correct for area of pixels.
        if self.areas is not None:
            self.evt.dat/=self.areas
        self._getRawShape()

class CsPadObject(TiledCameraObject):  
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(CsPadObject, self).__init__(det,env,run, **kwargs)
        self._common_mode_list = [1,5,55,10,0,-1, 30] #zero-peak, unbonded, unbonded (high no-correcion threshold), mixed, none, raw data, calib
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print('Common mode %d is not an option for a CsPad detector, please choose from: '%self.common_mode, self._common_mode_list)
        self.pixelsize=[110e-6]
        if self.ped is not None:
            try:
                 self.imgShape = self.det.image(run, self.ped).shape
            except:
                 self.imgShape = self.ped.shape
        if self.gain_mask.std()!=0:
          self.gain_mask = self.gain_mask.astype(float)
          self.gain_mask*=6.85
          self.gain_mask[self.gain_mask==0]=1.
          self.gain*=self.gain_mask

    def getData(self, evt):
        super(CsPadObject, self).getData(evt)
        mbits=0 #do not apply mask (would set pixels to zero)
        #mbits=1 #set bad pixels to 0
        if self.common_mode%100==1:
            self.evt.dat = self.det.calib(evt, cmpars=(1,25,40,100,0), mbits=mbits)
            try:
                cmValues = self.det.common_mode_correction(self.run, self.det.calib(evt, cmpars=(0), mbits=mbits), cmpars=(1,25,40,100,0), mbits=mbits)
                self.evt.__dict__['env_cmValues'] = np.array([tile[0][0] for tile in cmValues])
            except:
                pass

            needGain=False
        elif self.common_mode%100==5:
            self.evt.dat = self.det.calib(evt, cmpars=(5,100), mbits=mbits)
            needGain=False
        elif self.common_mode%100==55:
            self.evt.dat = self.det.calib(evt, cmpars=(5,5000), mbits=mbits)
            needGain=False
        elif self.common_mode%100==10:
            self.evt.dat = self.det.calib(evt, cmpars=(1,25,40,100,1), mbits=mbits)
            needGain=False

        #override gain if desired
        if self.local_gain is not None and self.local_gain.shape == self.evt.dat.shape and self.common_mode in [1,5,55,10]:
            self.evt.dat*=self.local_gain   #apply own gain

        #correct for area of pixels.
        if self.areas is not None:
            self.evt.dat/=self.areas
        self._getRawShape()

class CsPad2MObject(CsPadObject):  
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(CsPad2MObject, self).__init__(det,env,run, **kwargs)
        print('cspad alias',det.alias)
        self._common_mode_list = [0, 1,5,10,-1,30] #none, zero-peak, unbonded, mixed, raw data, calib
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print('Common mode %d is not an option for a CsPad2M detector, please choose from: '%self.common_mode, self._common_mode_list)
        if self.ped is not None:
            try:
                 self.imgShape = self.det.image(run, self.ped).shape
            except:
                 self.imgShape = self.ped.shape

    def getData(self, evt):
        super(CsPad2MObject, self).getData(evt)

class EpixObject(TiledCameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(EpixObject, self).__init__(det,env,run, **kwargs)
        self._common_mode_list = [6, 36, 4, 34, 45, 46, 47, 0,-1, 30] # Jacob (norm), Jacob, def, def(ami-like), mine, mine (norm), mine (norm-bank), none, raw, calib
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print('Common mode %d is not an option for as Epix detector, please choose from: '%self.common_mode, self._common_mode_list)
        self.pixelsize=[50e-6]

        epixCfg = env.configStore().get(psana.Epix.Config100aV2, self.det.source)
        self.carrierId0 = epixCfg.carrierId0()
        self.carrierId1 = epixCfg.carrierId1()
        self.digitalCardId0 = epixCfg.digitalCardId0()
        self.digitalCardId1 = epixCfg.digitalCardId1()
        self.analogCardId0 = epixCfg.analogCardId0()
        self.analogCardId1 = epixCfg.analogCardId1()

        if self.ped is not None:
            try:
                 self.imgShape = self.det.image(run, self.ped).shape
            except:
                 self.imgShape = self.ped.shape

        if self.common_mode==47:
            for i in range(0,16):
                bmask = np.zeros_like(self.rms)
                bmask[(i%2)*352:(i%2+1)*352,768/8*(i/2):768/8*(i/2+1)]=1
                self.bankMasks.append(bmask.astype(bool))
                
    def getData(self, evt):
        super(EpixObject, self).getData(evt)
        mbits=0 #do not apply mask (would set pixels to zero)
        #mbits=1 #set bad pixels to 0
        needGain=True
        if self.common_mode%100==6:
            self.evt.dat = self.det.calib(evt, cmpars=[6], mbits=mbits, rms = self.rms, normAll=True)
            needGain=False
        elif self.common_mode%100==36:
            self.evt.dat = self.det.calib(evt, cmpars=[6], mbits=mbits, rms = self.rms)
            needGain=False
        elif self.common_mode%100==34:
            self.evt.dat = self.det.calib(evt, cmpars=(4,6,100,100), mbits=mbits)
            needGain=False
        elif self.common_mode%100==4:
            self.evt.dat = self.det.calib(evt, cmpars=(4,6,30,10), mbits=mbits)
            needGain=False
        elif self.common_mode%100==45:
            self.evt.dat = self.det.raw_data(evt)-self.ped
            self.evt.dat = cm_epix(self.evt.dat,self.rms,mask=self.statusMask)
        elif self.common_mode%100==46:
            self.evt.dat = self.det.raw_data(evt)-self.ped
            self.evt.dat = cm_epix(self.evt.dat,self.rms,normAll=True,mask=self.statusMask)
        elif self.common_mode%100==47:
            self.evt.dat = self.det.raw_data(evt)-self.ped
            for ib,bMask in enumerate(self.bankMasks):
                #print ib, np.median(self.evt.dat[bMask]), bMask.sum()
                self.evt.dat[bMask]-=np.median(self.evt.dat[bMask])
            self.evt.dat = cm_epix(self.evt.dat,self.rms,mask=self.statusMask)

        #override gain if desired
        if self.local_gain is not None and self.local_gain.shape == self.evt.dat.shape and self.common_mode in [6,36,34,3,4,45,46,47]:
            self.evt.dat*=self.local_gain   #apply own gain
        elif self.local_gain is None and self.gain is not None and self.gain.shape == self.evt.dat.shape and self.common_mode in [45,46,47]:
            self.evt.dat*=self.gain   #apply gain after own common mode correction

        #correct for area of pixels.
        if self.areas is not None:
            self.evt.dat/=self.areas

        #store environmental row 
        try:
            envRow = evt.get(psana.Epix.ElementV3,psana.Source(self.det.alias)).environmentalRows()[1]
            self.evt.__dict__['env_temp0']  = envRow[0]*0.01
            self.evt.__dict__['env_temp1']  = envRow[1]*0.01
            self.evt.__dict__['env_humidity']  = envRow[2]*0.01

            self.evt.__dict__['env_AnalogI']  = envRow[3]*0.001
            self.evt.__dict__['env_DigitalI']  = envRow[4]*0.001
            self.evt.__dict__['env_GuardI']  = envRow[5]*0.000001
            self.evt.__dict__['env_BiasI']  = envRow[6]*0.000001
            self.evt.__dict__['env_AnalogV']  = envRow[7]*0.001
            self.evt.__dict__['env_DigitalV']  = envRow[8]*0.001
        except:
            pass
        self._getRawShape()

#
# as a
#
class Epix10kObject(TiledCameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(Epix10kObject, self).__init__(det,env,run, **kwargs)
        self._common_mode_list = [80, 84, 85, 0, -1, -2, 30] # official, ped sub, raw, raw_gain, calib
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print('Common mode %d is not an option for as Epix detector, please choose from: '%self.common_mode, self._common_mode_list)
        self.pixelsize=[100e-6]
        self.isGainswitching=True

        epixCfg = env.configStore().get(psana.Epix.Config10kaV1, det.source)              
        if epixCfg is None:
            epixCfg = env.configStore().get(psana.Epix.Config10kaV2, det.source)    
        self.pixelGain=[]
        self.nomGain=[1.,3.,100.,100.,3.,100.] #H,M,L,(HL auto), (ML auto - 2 values)
        #asicList=[0,3,1,2]
        asicList=[0,1,3,2]
        self.carrierId0 = epixCfg.carrierId0()
        self.carrierId1 = epixCfg.carrierId1()
        self.pixelConfig = epixCfg.asicPixelConfigArray()
        trbits=[]
        for ia in range(epixCfg.asics_shape()[0]):
            trbits.append(epixCfg.asics(ia).trbit())
        self.trbit = np.array(trbits)

        cfgShape=epixCfg.asicPixelConfigArray().shape
        cfgReshape=epixCfg.asicPixelConfigArray().reshape(int(cfgShape[0]/2), cfgShape[1]*2,order='F')
        pixelGain=np.ones_like(cfgReshape).astype(float)
        for ia in asicList:
            if epixCfg.asics(ia).trbit()==1:
                continue
            asicGainConfig=cfgReshape[:,ia*int(cfgShape[1]/2):(ia+1)*int(cfgShape[1]/2)]
            pixelGain[:,ia*int(cfgShape[1]/2):(ia+1)*int(cfgShape[1]/2)]=((asicGainConfig&0x4)/4).astype(float)*self.nomGain[1] + ((np.ones_like(asicGainConfig)-(asicGainConfig&0x4)/4)).astype(float)*self.nomGain[2]
            #pixelGain[:,ia*cfgShape[1]/2:(ia+1)*cfgShape[1]/2]=((asicGainConfig&0x4)/4).astype(float)*100./3. + ((np.ones_like(asicGainConfig)-(asicGainConfig&0x4)/4)).astype(float)*100.
        self.pixelGain = pixelGain.reshape(cfgShape,order='F')

        self.gainSetting = 0
        if self.pixelGain.mean()==0:
            if self.trbit.mean()==1:
                self.gainSetting = 1 #gain switch HL
            else:
                self.gainSetting = 2 #gain switch ML

        if self.rms is None or self.rms.shape!=self.ped.shape:
            self.rms=np.ones_like(self.ped)
        try:
            self.imgShape=self.det.image(run, self.ped[0]).shape
        except:
            if len(self.ped[0].squeeze().shape)==2:
                self.imgShape=self.ped[0].squeeze().shape
            else:
                self.imgShape=None
        self._gainSwitching = True                

    def getData(self, evt):
        super(Epix10kObject, self).getData(evt)
        mbits=0 #do not apply mask (would set pixels to zero)
        #mbits=1 #set bad pixels to 0
        if self.common_mode<0:
          if self.common_mode==-2:
            self.evt.dat = self.evt.dat
          else:
            self.evt.dat = self.evt.dat&0x3fff
          #self.evt.gainbit = (self.evt.dat&0xc000>0)
          try:
            epixCalRows = evt.get(psana.Epix.ArrayV1, psana.Source(self.det.alias)).calibrationRows()
            self.evt.__dict__['env_calibRows']  = epixCalRows
          except:
            epixCalRows = None
        elif self.common_mode==0:
            ##########
            ### FIX ME epix10ka
            #will need to read gain bit from data and use right pedestal.
            #will hopefully get calib function for this.
            ##########
            if len(self.ped.shape)>3:
                self.evt.dat = (self.det.raw_data(evt)&0x3fff)-self.ped[0]
            else:
                self.evt.dat = (self.det.raw_data(evt)&0x3fff)-self.ped
        elif self.common_mode%100==80: #no common mode
            #self.evt.dat = self.det.calib(evt, mbits=mbits)
            evt_dat = self.det.calib(evt, cmpars=(7,0,100))
            self.evt.dat = np.empty(evt_dat.shape, dtype=np.float32) #can't reember why I'm doing this. I think for the data type.
            self.evt.dat[:,:,:] = evt_dat[:,:,:]
        elif self.common_mode%100==84: #correct in row
            evt_dat = self.det.calib(evt, cmpars=(7,2,10,10))
            self.evt.dat = np.empty(evt_dat.shape, dtype=np.float32) 
            self.evt.dat[:,:,:] = evt_dat[:,:,:]
        elif self.common_mode%100==85: #correct in row&cols
            evt_dat = self.det.calib(evt, cmpars=(7,3,10,10))
            self.evt.dat = np.empty(evt_dat.shape, dtype=np.float32) 
            self.evt.dat[:,:,:] = evt_dat[:,:,:]
            
        #override gain if desired -- this looks like CsPad.
        if self.local_gain is not None and self.local_gain.shape == self.evt.dat.shape and self.common_mode in [1,5,55,10]:
            self.evt.dat*=self.local_gain   #apply own gain

        #correct for area of pixels.
        if self.areas is not None:
            self.evt.dat/=self.areas

        #store environmental row 
        try:
            envRows = evt.get(psana.Epix.ArrayV1,psana.Source(self.det.alias)).environmentalRows()
            #envRows = evt.get(psana.Epix.ArrayV1,psana.Source(self.det.alias)).environmentalRows().astype(np.uint16)
            #print 'envRow: iunt16', envRows.shape
            temp1    = []
            temp1H    = []
            temp2    = []
            temp3    = []
            hum      = []
            humH      = []
            anacur   = []
            digcur   = []
            anav     = []
            digv     = []
            anatemp  = []
            digtemp  = []
            for moduleRow in envRows:
                tempVal= np.array([ divmod(temp, 1<<16) for temp in moduleRow[0]]).flatten()
                humH.append(tempVal[16]/65535.0 * 100)
                temp1H.append(getThermistorTemp(tempVal[26]))
                moduleRow = moduleRow.astype(np.uint16)
                #print moduleRow.shape, moduleRow[1][26], moduleRow[1][27]
                temp1.append(getThermistorTemp(moduleRow[1][26]))
                temp2.append(getThermistorTemp(moduleRow[1][27]))
                temp3.append(moduleRow[1][17]/65535.0  * 175 - 45)
                hum.append(moduleRow[1][16]/65535.0 * 100)
                anacur.append(moduleRow[1][31]*1024.0/4095/0.2)
                digcur.append(moduleRow[1][28]*1024.0/4095/0.2)
                anav.append(moduleRow[1][32]*1024.0/4095 * 100)
                digv.append(moduleRow[1][29]*1024.0/4095 * 100)
                anatemp.append(moduleRow[1][33]*2.048/4095*(130/(0.882-1.951)) + (0.882/0.0082+100))
                digtemp.append(moduleRow[1][30]*2.048/4095*(130/(0.882-1.951)) + (0.882/0.0082+100))

            #print 'temp1',temp1
            #print 'temp1H',temp1H
            #print 'hum',hum
            #print 'humH',humH
            self.evt.__dict__['env_temp1']  = np.array(temp1) 
            self.evt.__dict__['env_temp2']  = np.array(temp2) 
            self.evt.__dict__['env_temp3']  = np.array(temp3) 
            self.evt.__dict__['env_humidity']  = np.array(hum) 

            self.evt.__dict__['env_AnalogI']  = np.array(anacur) 
            self.evt.__dict__['env_DigitalI']  = np.array(digcur) 
            self.evt.__dict__['env_AnalogV']  = np.array(anav) 
            self.evt.__dict__['env_DigitalV']  = np.array(digv) 
            self.evt.__dict__['env_AnalogTemp']  =  np.array(anatemp) 
            self.evt.__dict__['env_DigitalTemp']  =  np.array(digtemp) 
        except:
            pass
        self._getRawShape()

class Epix10k2MObject(TiledCameraObject): 
    def __init__(self, det, env, run, **kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(Epix10k2MObject, self).__init__(det,env,run, **kwargs)
        self._common_mode_list = [80, 84, 81, 82, 85, 180, 181, 0, -1, -2, 30] # official, sn kludge, ped sub, raw, calib
        self.common_mode = kwargs.get('common_mode', self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print('Common mode %d is not an option for as Epix detector, please choose from: '%self.common_mode, self._common_mode_list)
        self.pixelsize=[100e-6]
        self.isGainswitching=True

        configs = ['Config10ka2MV1', 'Config10ka2MV2', 'Config10kaQuadV1', 'Config10kaQuadV2']
        for config in configs:
            epixCfg = env.configStore().get(getattr(psana.Epix, config), det.source)
            if epixCfg is not None:
                print('Loaded config {}'.format(config))
                break

        self.carrierId0 = []
        self.carrierId1 = []
        self.pixelConfig = []
        self.trbit = []
        self.pixelGain=[]
        self.nomGain=[1.,3.,100.,100.,3.,100.] #H,M,L,(HL auto), (ML auto - 2 values)
        #asicList=[0,3,1,2]
        asicList=[0,1,3,2]
        for i in range(epixCfg.elemCfg_shape()[0]):
            elemCfg=epixCfg.elemCfg(i)
            self.carrierId0.append(elemCfg.carrierId0())
            self.carrierId1.append(elemCfg.carrierId1())
            self.pixelConfig.append(elemCfg.asicPixelConfigArray())
            trbits=[]
            for ia in range(elemCfg.asics_shape()[0]):
                trbits.append(elemCfg.asics(ia).trbit())
            self.trbit.append(trbits)
            cfgShape=elemCfg.asicPixelConfigArray().shape
            cfgReshape=elemCfg.asicPixelConfigArray().reshape(int(cfgShape[0]/2), cfgShape[1]*2,order='F')
            pixelGain=np.ones_like(cfgReshape).astype(float)
            for ia in asicList:
                if elemCfg.asics(ia).trbit()==1:
                    continue
                asicGainConfig=cfgReshape[:,ia*int(cfgShape[1]/2):(ia+1)*int(cfgShape[1]/2)]
                pixelGain[:,ia*int(cfgShape[1]/2):(ia+1)*int(cfgShape[1]/2)]=((asicGainConfig&0x4)/4).astype(float)*self.nomGain[1] + ((np.ones_like(asicGainConfig)-(asicGainConfig&0x4)/4)).astype(float)*self.nomGain[2]
                #pixelGain[:,ia*cfgShape[1]/2:(ia+1)*cfgShape[1]/2]=((asicGainConfig&0x4)/4).astype(float)*100./3. + ((np.ones_like(asicGainConfig)-(asicGainConfig&0x4)/4)).astype(float)*100.
            pixelGain=pixelGain.reshape(cfgShape,order='F')
            self.pixelGain.append(pixelGain)
        self.pixelGain=np.array(self.pixelGain)
        self.trbit = np.array(self.trbit)
        self.gainSetting = 0
        if self.pixelGain.mean()==0:
            if self.trbit.mean()==1:
                self.gainSetting = 1 #gain switch HL
            else:
                self.gainSetting = 2 #gain switch ML
        self.rms = self.rms.squeeze()
        self.gain = self.gain.squeeze()
        self.ped = self.ped.squeeze()
        if self.rms is None or self.rms.shape!=self.ped.shape:
            self.rms=np.ones_like(self.ped)
        self.imgShape=self.det.image(run, self.ped[0]).shape
        self._gainSwitching = True                

        ##stuff for ghost correction
        if self.common_mode > 100:
            #read p0&p1 from h5 files.
            calibDir = '/reg/d/psdm/%s/%s/calib/Epix10ka2M::CalibV1/XcsEndstation.0:Epix10ka2M.0/ghost'%(env.experiment()[:3], env.experiment())
            ghostFile = tables.open_file('%s/GhostFits_%s_Run%03d.h5'%(calibDir,env.experiment(),run))
            self.p0 = ghostFile.root.epix10ka2m_OLS_p0.read()
            self.p1 = ghostFile.root.epix10ka2m_OLS_p1.read()
            #self.ghostMask = (abs(self.p1-0.0073)<0.002)
            self.dsIdxRun = psana.DataSource('exp=%s:run=%d:idx'%(env.experiment(),run)).runs().next()
            self.dsIdxSecs=[]
            self.dsIdxFiducials=[]
            self.dsIdxTimes=self.dsIdxRun.times()
            for evtt in self.dsIdxTimes:
              self.dsIdxFiducials.append(evtt.fiducial())
              self.dsIdxSecs.append(evtt.time()>>32)
            self.dsIdxSecs=np.array(self.dsIdxSecs)
            self.dsIdxFiducials=np.array(self.dsIdxFiducials)
            self.preEvtDet=DetObject(self.det.alias,env,run, common_mode=self.common_mode%100)
            ghostmask = (~(np.isnan(self.p0)) & ~(np.isnan(self.p1)))
            self.mask = (ghostmask&self.mask)

    def getData(self, evt):
        super(Epix10k2MObject, self).getData(evt)
        mbits=0 #do not apply mask (would set pixels to zero)
        #mbits=1 #set bad pixels to 0
        if self.common_mode<0:
          if self.common_mode==-2:
            self.evt.dat = self.evt.dat
          else:
            self.evt.dat = self.evt.dat&0x3fff
          #self.evt.gainbit = (self.evt.dat&0xc000>0)
          try:
            epixCalRows = evt.get(psana.Epix.ArrayV1, psana.Source(self.det.alias)).calibrationRows()
            self.evt.__dict__['env_calibRows']  = epixCalRows
          except:
            epixCalRows = None
        elif self.common_mode==0:
            ##########
            ### FIX ME epix10ka
            #will need to read gain bit from data and use right pedestal.
            #will hopefully get calib function for this.
            ##########
            if len(self.ped.shape)>3:
                self.evt.dat = (self.det.raw_data(evt)&0x3fff)-self.ped[0]
            else:
                self.evt.dat = (self.det.raw_data(evt)&0x3fff)-self.ped
        elif self.common_mode%100==30: #default calib
            evt_dat = self.det.calib(evt)
            self.evt.dat = np.empty(evt_dat.shape, dtype=np.float32) #can't reember why I'm doing this. I think for the data type.
            self.evt.dat[:,:,:] = evt_dat[:,:,:]
        elif self.common_mode%100==80: #placeholder for specific treatment of epix10k
            evt_dat = self.det.calib(evt, cmpars=(7,0,100))
            self.evt.dat = np.empty(evt_dat.shape, dtype=np.float32) #can't reember why I'm doing this. I think for the data type.
            self.evt.dat[:,:,:] = evt_dat[:,:,:]
        elif self.common_mode%100==84: #rows 
            evt_dat = self.det.calib(evt, cmpars=(7,2,10,10))
            self.evt.dat = np.empty(evt_dat.shape, dtype=np.float32) 
            self.evt.dat[:,:,:] = evt_dat[:,:,:]
        elif self.common_mode%100==85: #rows&cols
            evt_dat = self.det.calib(evt, cmpars=(7,3,10,10))
            self.evt.dat = np.empty(evt_dat.shape, dtype=np.float32) 
            self.evt.dat[:,:,:] = evt_dat[:,:,:]
        elif self.common_mode%100==81: #my horrible kludge that might work for fixed gains.
            if len(self.ped.shape)>3:
                print('not a supported mode, return None')
                self.evt.dat = None
            else:
                raw_dat = self.det.raw_data(evt)
                gain_dat = raw_dat&0xc000
                dat_dat =  (raw_dat&0x3fff)-self.ped
                if self.gainSetting == 0:
                    self.evt.dat = dat_dat*self.pixelGain 
                elif self.gainSetting == 1:
                    self.evt.dat = dat_dat[gain_dat==0]+dat_dat[gain_dat>1]*self.nomGain[3]
                else:
                    self.evt.dat = dat_dat[gain_dat==0]*self.nomGain[4]+dat_dat[gain_dat>1]*self.nomGain[5]
                    #self.evt.dat = dat_dat[gain_dat==0]*100./3.+dat_dat[gain_dat>1]*100.
        elif self.common_mode%100==82: #use new common mode.
            raw_dat = self.det.raw_data(evt)
            gain_dat = raw_dat&0xc000
            #get the calibration rows.
            try:
                epixCalRows = evt.get(psana.Epix.ArrayV1, psana.Source(self.det.alias)).calibrationRows()
                self.evt.__dict__['env_calibRows']  = epixCalRows
            except:
                epixCalRows = None
                
            if len(self.ped.shape)==3:
                print('own pedestal...')
                dat_dat =  (raw_dat&0x3fff)-self.ped
            else:
                ped = np.zeros_like(self.ped[0])
                #make pedestal from normal calibration 7-array based on gain mode.
                if self.gainSetting == 1:
                    ped=self.ped[3][gain_dat==0]+self.ped[4][gain_dat>1]
                elif self.gainSetting == 2:
                    ped=self.ped[5][gain_dat==0]+self.ped[6][gain_dat>1]
                else:
                    if np.abs(self.pixelGain.mean()-self.nomGain[0])<0.01:
                        ped=self.ped[0]
                    elif np.abs(self.pixelGain.mean()-self.nomGain[1])<0.01:
                        ped=self.ped[1]
                    elif np.abs(self.pixelGain.mean()-self.nomGain[2])<0.01:
                        ped=self.ped[2]
                dat_dat =  (raw_dat&0x3fff)-ped

            if self.gainSetting == 0:
                evt_dat = dat_dat*self.pixelGain 
            elif self.gainSetting == 1:
                evt_dat = dat_dat[gain_dat==0]+dat_dat[gain_dat>1]*self.nomGain[3]
            else:
                evt_dat = dat_dat[gain_dat==0]*self.nomGain[4]+dat_dat[gain_dat>1]*self.nomGain[5]
                #evt_dat = dat_dat[gain_dat==0]*100./3.+dat_dat[gain_dat>1]*100.
            
            self.evt.dat = np.empty(evt_dat.shape, dtype=np.float32) #can't reember why I'm doing this. I think for the data type.
            self.evt.dat[:,:,:] = evt_dat[:,:,:]
            
        if self.common_mode>100: #correct for the ghost image after official calibration
            #current event info
            thisFid = evt.get(psana.EventId).fiducials()
            thisSec = evt.get(psana.EventId).time()[0]
            #find this in indedx dataset.
            sameFid = np.argwhere(self.dsIdxFiducials==thisFid)
            #if we have a choice, check the second field.
            sameIdx=-1
            if sameFid.shape[0]==1:
                sameIdx=sameFid[0][0]
            else:
                for thisFid in sameFid:
                    if thisSec==self.dsIdxSecs[thisFid[0]]:
                        sameIdx=thisFid[0]
            if sameIdx<1:
                print('did not find event or this event is first event...',sameIdx)
                return None
            #now get the event
            try:
                preEvt = self.dsIdxRun.event(self.dsIdxTimes[sameIdx-1])
            except:
                print('DEBUG: failed to get previous event')
                return None
            #now get data for this event
            orgData = self.evt.dat.copy()
            self.preEvtDet.getData(preEvt)
            #preData = self.preEvtDet.evt.dat.copy()
            corrData = self.p0+self.p1*self.preEvtDet.evt.dat
            self.evt.dat-=corrData

        #override gain if desired -- this looks like CsPad.
        if self.local_gain is not None and self.local_gain.shape == self.evt.dat.shape and self.common_mode in [1,5,55,10]:
            self.evt.dat*=self.local_gain   #apply own gain

        #correct for area of pixels.
        if self.areas is not None:
            self.evt.dat/=self.areas

        #store environmental row 
        try:
            envRows = evt.get(psana.Epix.ArrayV1,psana.Source(self.det.alias)).environmentalRows()
            #envRows = evt.get(psana.Epix.ArrayV1,psana.Source(self.det.alias)).environmentalRows().astype(np.uint16)
            #print 'envRow: iunt16', envRows.shape
            temp1    = []
            temp1H    = []
            temp2    = []
            temp3    = []
            hum      = []
            humH      = []
            anacur   = []
            digcur   = []
            anav     = []
            digv     = []
            anatemp  = []
            digtemp  = []
            for moduleRow in envRows:
                tempVal= np.array([ divmod(temp, 1<<16) for temp in moduleRow[0]]).flatten()
                humH.append(tempVal[16]/65535.0 * 100)
                temp1H.append(getThermistorTemp(tempVal[26]))
                moduleRow = moduleRow.astype(np.uint16)
                #print moduleRow.shape, moduleRow[1][26], moduleRow[1][27]
                temp1.append(getThermistorTemp(moduleRow[1][26]))
                temp2.append(getThermistorTemp(moduleRow[1][27]))
                temp3.append(moduleRow[1][17]/65535.0  * 175 - 45)
                hum.append(moduleRow[1][16]/65535.0 * 100)
                anacur.append(moduleRow[1][31]*1024.0/4095/0.2)
                digcur.append(moduleRow[1][28]*1024.0/4095/0.2)
                anav.append(moduleRow[1][32]*1024.0/4095 * 100)
                digv.append(moduleRow[1][29]*1024.0/4095 * 100)
                anatemp.append(moduleRow[1][33]*2.048/4095*(130/(0.882-1.951)) + (0.882/0.0082+100))
                digtemp.append(moduleRow[1][30]*2.048/4095*(130/(0.882-1.951)) + (0.882/0.0082+100))

            #print 'temp1',temp1
            #print 'temp1H',temp1H
            #print 'hum',hum
            #print 'humH',humH
            self.evt.__dict__['env_temp1']  = np.array(temp1) 
            self.evt.__dict__['env_temp2']  = np.array(temp2) 
            self.evt.__dict__['env_temp3']  = np.array(temp3) 
            self.evt.__dict__['env_humidity']  = np.array(hum) 

            self.evt.__dict__['env_AnalogI']  = np.array(anacur) 
            self.evt.__dict__['env_DigitalI']  = np.array(digcur) 
            self.evt.__dict__['env_AnalogV']  = np.array(anav) 
            self.evt.__dict__['env_DigitalV']  = np.array(digv) 
            self.evt.__dict__['env_AnalogTemp']  =  np.array(anatemp) 
            self.evt.__dict__['env_DigitalTemp']  =  np.array(digtemp) 
        except:
            pass
        self._getRawShape()

#needs to be retrofitted to work with both rayonix cameras.
class RayonixObject(CameraObject): 
    def __init__(self, det,env,run,**kwargs):
        #super().__init__(det,env,run, **kwargs)
        super(RayonixObject, self).__init__(det,env,run, **kwargs)
        self.common_mode = kwargs.get('common_mode', -1)
        if self.common_mode is None:
            self.common_mode = -1
        #this is NOT correct for the new, bigger rayonix. Fix me.
        rcfg = env.configStore().get(psana.Rayonix.ConfigV2,psana.Source('%s'%det.name))
        binning=rcfg.binning_s()
        try:
          self.pixelsize=[rcfg.pixelWidth(),rcfg.pixelHeight()]
        except:
          self.pixelsize=[170e-3/3840*binning] #needs to change for bigger rayonix.
        try:
          self.ped = np.zeros([rcfg.width(), rcfg.height()])
          self.imgShape = (rcfg.width(), rcfg.height())
        except: 
          npix = int(170e-3/self.pixelsize[0])
          self.ped = np.zeros([npix, npix])
          self.imgShape = (npix, npix)
        if self.x is None:
            self.x = np.arange(0,self.ped.shape[0]*self.pixelsize[0], self.pixelsize[0])*1e6
            self.y = np.arange(0,self.ped.shape[0]*self.pixelsize[0], self.pixelsize[0])*1e6
            self.x, self.y = np.meshgrid(self.x, self.y)
        if self.mask is None:
            self.mask = np.ones(self.imgShape)
            self.cmask = np.ones(self.imgShape)


class UxiObject(DetObjectClass):
    def __init__(self, run, **kwargs):
        self.det = None
        self._name = kwargs.get('name', 'uxi')
        self._uxiDict =  kwargs.get('uxiDict', None)
        uxiConfigDict =  kwargs.get('uxiConfigDict', None)
        self.common_mode =  kwargs.get('common_mode', 0)
        if self.common_mode is None:
            self.common_mode = 0
        for k, value in iteritems(uxiConfigDict):
          setattr(self, k, value)
        #now get the pedestals from run unless we stuff this into the configDict before.
        iDark, darkA, darkB =  getDarks(int(run))
        setattr(self, 'pedestal_run', iDark)
        setattr(self, 'ped', np.array([darkA, darkB]))
        setattr(self, 'rms', np.ones_like(getattr(self, 'ped')))
        setattr(self, 'mask', None)#  np.ones_like(getattr(self, 'ped')))
        setattr(self, 'cmask', None)# getattr(self, 'mask'))
        #geometry information. all frames have same x/y
        self.pixelsize=[25e-6]
        self.x = np.arange(0,self.ped.shape[-2]*self.pixelsize[0], self.pixelsize[0])*1e6
        self.y = np.arange(0,self.ped.shape[-1]*self.pixelsize[0], self.pixelsize[0])*1e6
        self.y, self.x = np.meshgrid(self.y, self.x)
        self.x=np.array([self.x for i in range(self.ped.shape[0])])
        self.y=np.array([self.y for i in range(self.ped.shape[0])])
        #common mode stuff
        self.cm_maskedROI =  kwargs.get('cm_maskedROI', None)
        if self.cm_maskedROI is not None and isinstance(self.cm_maskedROI, list):
            self.cm_maskedROI = np.array(self.cm_maskedROI)
            if len(self.cm_maskedROI.shape)==1:
              self.cm_maskedROI=np.array([self.cm_maskedROI.tolist()]*2)
        self.cm_photonThres =  kwargs.get('cm_photonThres', 50)
        self.cm_maskNeighbors =  kwargs.get('cm_maskNeighbors', 0)
        self.cm_maxCorr =  kwargs.get('cm_maxCorr', self.cm_photonThres)
        self.cm_minFrac =  kwargs.get('cm_minFrac', 0.05)

    def getData(self, evt):
        try:
            getattr(self, 'evt')
        except:
            self.evt = event()
        self.evt.dat = None

        #get the timestamp. 
        evttime = evt.get(psana.EventId).time()
        evtfid = evt.get(psana.EventId).fiducials()
        #get the right frame from self._uxiDict

        #find uxi pictures taken the same second & get nsec & fiducials
        uxiEvent=False; uxiIdx=-1
        if (len(np.argwhere(self._uxiDict['lcls_ts_secs']==evttime[0]))>0):
            if (len(np.argwhere(self._uxiDict['lcls_ts_secs']==evttime[0]))>1):
                print('more than one UXI match???')
            secIdx=np.argwhere(self._uxiDict['lcls_ts_secs']==evttime[0])[0][0]
            if len(np.argwhere(self._uxiDict['lcls_ts_necs']==evttime[1]))>0 and len(np.argwhere(self._uxiDict['lcls_ts_fids']==evtfid))>0:
                #evtNr_withUxi=eventNr
                #evtFid_withUxi=evtfid
                if evttime[1]==self._uxiDict['lcls_ts_necs'][secIdx] and evtfid==self._uxiDict['lcls_ts_fids'][secIdx]:
                    uxiEvent=True
                    uxiIdx=secIdx
        #print uxiEvent, uxiIdx
        if not uxiEvent:
            return None

        #now add uxi data.
        evtUxiDict={'uxi':{}}
        dataFrame=[]
        dataFrame.append(self._uxiDict['frameA'][uxiIdx])
        dataFrame.append(self._uxiDict['frameB'][uxiIdx])
        for key in self._uxiDict.keys():
            if key in ['frameA','frameB']: continue
            if key.find('lcls_ts')>=0: continue
            #print 'append key to event: ',key
            if isinstance(self._uxiDict[key][uxiIdx],basestring):
                if self._uxiDict[key][uxiIdx].find('.')>=0:
                    self.evt.__dict__['env_%s'%key]=float(self._uxiDict[key][uxiIdx])
                else:
                    self.evt.__dict__['env_%s'%key]=int(self._uxiDict[key][uxiIdx])
            else:
                self.evt.__dict__['env_%s'%key]=self._uxiDict[key][uxiIdx]

        if self.common_mode<0:
            #return raw data 
            self.evt.dat = np.array(dataFrame)
            return

        #subtract the pedestal.
        try:
            dataFrame-=self.ped
        except:
            print('Uxi pedestal subtraction failed: shapes:',dataFrame.shape, self.ped.shape)

        if self.common_mode==0:
            #return pedestal subtraction.
            self.evt.dat = dataFrame
            return

        if self.cm_maskedROI is not None:
            #print('subtract mean from masked area is passed., use ROI: ',self.cm_maskedROI)
            cm_mask_mean = [ np.nanmean(frame[maskedROI[0]:maskedROI[1],maskedROI[2]:maskedROI[3]]) for frame,maskedROI in zip(dataFrame, self.cm_maskedROI )]
            cm_mask_med = [ np.nanmedian(frame[maskedROI[0]:maskedROI[1],maskedROI[2]:maskedROI[3]]) for frame,maskedROI in zip(dataFrame, self.cm_maskedROI )]
            cm_mask_std = [ np.nanstd(frame[maskedROI[0]:maskedROI[1],maskedROI[2]:maskedROI[3]]) for frame,maskedROI in zip(dataFrame, self.cm_maskedROI )]
            self.evt.__dict__['env_cm_masked_med'] = np.array(cm_mask_med)
            self.evt.__dict__['env_cm_masked_mean'] = np.array(cm_mask_mean)
            self.evt.__dict__['env_cm_masked_std'] = np.array(cm_mask_std)
            dataFrame = np.array([df-corr for df, corr in zip(dataFrame,cm_mask_med)])

        #only subtract the masked ROI
        if self.common_mode==80:
            self.evt.dat=dataFrame
            return

        corrImg, cmDict = cm_uxi(dataFrame, self.cm_photonThres, self.cm_maxCorr, self.cm_minFrac, self.cm_maskNeighbors)
        for k, value in iteritems(cmDict):
          self.evt.__dict__['env_%s'%k] = value
        self.evt.dat = corrImg

        #corrImg, cmValues, cmValues_used = cm_uxi(dataFrame, self.cm_photonThres, self.cm_maxCorr, self.cm_minFrac, self.cm_maskNeighbors)
        #self.evt.__dict__['env_cm_RowMed'] = cmValues
        #self.evt.__dict__['env_cm_RowMed_used'] = cmValues_used
        #self.evt.dat = corrImg
        self._getRawShape()
        return

##---------------------------------------------------------------------

def _test_cs140k():
    """
    test of the DetObject via a factory class for the cs140k
    """
    import time
    import psana
    import numpy as np

    run=211; detname='cs140_0'
    ds = psana.DataSource('exp=xpptut15:run=%d'%run)
    evt = ds.events().next()

    from smalldata_tools.DetObject import DetObject
    det = DetObject(detname,ds.env(), run)

    try:
        assert(det.ped.dtype == np.float32)
        assert(det.ped.shape == (2, 185,388))
    except:
        print('test of pedestal for cs140 in xpptut15 failed')

    try:
        assert(det.mask.shape == (2, 185,388))
        assert(det.mask.sum() == 140120)
    except:
        print('test of mask for cs140 in xpptut15 failed')

    try:
        assert(det.det.image(run, det.ped).shape == (392,398))
    except:
        print('test of geometry for cs140 in xpptut15 failed')

    try:
        assert(det.common_mode == 1)
    except:
        print('test of 140k: default common_mode is not 1')

    try:
        dat = det.getData(evt)
    except:
        print('test of 140k: could not get raw data')


def _test_opal():
    """
    test of the DetObject via a factory class for the cs140k
    """
    import time
    import psana
    import numpy as np

    run=211; detname='opal_1'
    ds = psana.DataSource('exp=xpptut15:run=%d'%run)
    evt = ds.events().next()

    from smalldata_tools.DetObject import DetObject
    det = DetObject(detname,ds.env(), run)

    try:
        assert(det.ped.shape == (451, 1022))
    except:
        print('test of pedestal for opal in xpptut15 failed')

    try:
        assert(det.mask.shape == (451, 20122))
        assert(det.mask.sum() == 460992)
    except:
        print('test of mask for opal in xpptut15 failed')

    try:
        assert(det.det.image(run, det.ped).shape == (451, 1022))
    except:
        print('test of geometry for opal in xpptut15 failed')

    try:
        assert(det.common_mode == 0)
    except:
        print('test of opal: default common_mode is not 1')

    try:
        det.getData(evt)
        assert(det.evt.dat.shape == det.ped.shape)
    except:
        print('test of opal: could not get raw data')
