import os
import copy
import numpy as np
import h5py
import fnmatch

from matplotlib import pyplot as plt
from utilities import E2lam
from utilities import cm_epix
from utilities import rebin
from utilities import dropObject
import azimuthalBinning as ab
import droplet as droplet
import acf
import fitCenter
import photons

from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()

import psana
from collections import Counter
#import time

class ROIObject(dropObject):
  def __init__(self, ROI_limit, name='noname', writeArea=False, rms=None, mask=None):
    self.name = name
    if isinstance(ROI_limit, list):
      ROI = np.array(ROI_limit)
    else:
      ROI=ROI_limit
    if ROI.ndim==1 and ROI.shape[0]>2:
      ROI = ROI.reshape(ROI.shape[0]/2,2)        
    self.addField('bound', ROI.astype(int))
    self.writeArea = writeArea
    self.rebin = False
    if rms is not None:
      self.rms = self.applyROI(rms)
  def applyROI(self, array):
    if array.ndim < self.bound.ndim:
      print 'array has fewer dimensions that bound: ',array.ndim,' ',len(self.bound)
      return array
    #ideally would check this from FrameFexConfig and not on every events
    self.bound = np.array(self.bound)
    if array.ndim == self.bound.ndim:
      new_bound=[]
      for dim,arsz in zip(self.bound, array.shape):
        if max(dim) > arsz:
          new_bound.append([min(dim), arsz])
        else:
          new_bound.append([min(dim), max(dim)])
      self.bound = np.array(new_bound)
    elif self.bound.ndim==1:
      self.bound = np.array([min(self.bound), min(max(self.bound), array.shape[0])])
    #this needs to be more generic....of maybe just ugly spelled out for now.
    if self.bound.ndim==1:
      subarr = array[self.bound[0]:self.bound[1]]
    elif self.bound.ndim==2:
      subarr = array[self.bound[0,0]:self.bound[0,1],self.bound[1,0]:self.bound[1,1]]
    elif self.bound.ndim==3:
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
  def addRebin(self,shape=[],asImage=False):
    self.rebin = True
    self.__dict__['rebin'] = shape
  def Rebin(self, array):
    return rebin(array, self.__dict__['rebin'])
  def addProj(self, pjName = '', axis=0, singlePhoton=False, mean=False, cutADU=-1e6, cutRms=-1e6):
    self.__dict__['proj'+pjName]=[axis, singlePhoton, mean, cutADU, cutRms, 'proj'+pjName]
  def projection(self,arrayIn,params=[]):
    if len(params)<4:
      print 'projection does not have enough parameters'
      return np.nan
    singlePhoton=params[1]
    mean=params[2]
    cutADU=params[3]
    cutRms=params[4]
    array = arrayIn.copy()
    array.data[array.data<cutADU]=0
    if 'rms' in self.__dict__.keys() and self.rms is not None:
      array.data[array.data<cutRms*self.rms]=0
    if singlePhoton:
      array.data[array.data>0]=1
    #array.data = array.data.astype(np.float64)
    if mean:
      if params[0]<0:
        return array.squeeze().mean()
      else:
        return array.squeeze().mean(axis=params[0])
    if params[0]<0:
      return array.squeeze().sum()
    else:
      return array.squeeze().sum(axis=params[0])
  def getProjs(self):
    retkey = [ self[key] for key in self.__dict__.keys() if key[:4]=='proj' ]
    #print retkey
    return retkey
    
#idea is to get rid of "common" lines of all detectors as they come automatically.
class DetObject(dropObject):
    def __init__(self,srcName,env,run,name=None, common_mode=1, applyMask=0):
      self._srcName=srcName
      if name is not None:
        self._name = name
      else:
        self._name = srcName
      self.det=psana.Detector(srcName)
      self.run=run
      self.common_mode = common_mode
      self.storeEnv = False
      #det.dettype
      #1->CsPad, 2->Cs2x2, 13->epix100a
      #6->opal
      #16 -> aqiris, ?->oceanOptics
      self.applyMask = applyMask
      #default to CsPad
      self.pixelsize=[110e-6]
      if srcName.find('ayon')>=0:
        binning=env.configStore().get(psana.Rayonix.ConfigV2,psana.Source('rayonix')).binning_s()
        self.pixelsize=[170e-3/3840*binning]
      if srcName.find('epix')>=0:
        self.pixelsize=[75e-6]
        epixCfg = env.configStore().get(psana.Epix.Config100aV2, self.det.source)
        self.carrierId0 = epixCfg.carrierId0()
        self.carrierId1 = epixCfg.carrierId1()
        self.digitalCardId0 = epixCfg.digitalCardId0()
        self.digitalCardId1 = epixCfg.digitalCardId1()
        self.analogCardId0 = epixCfg.analogCardId0()
        self.analogCardId1 = epixCfg.analogCardId1()
      if self.det.name.__str__().find('OceanOptics')<0 and self.det.dettype != 16:
        self.rms = self.det.rms(run)
        self.ped = self.det.pedestals(run)
        self.gain = self.det.gain(run)
        self.imgShape = None
        if self.ped is not None:
          self.imgShape = self.ped.shape
        elif self.det.dettype == 19:
          if self.common_mode > 0:
            self.common_mode = -1
          npix = int(170e-3/self.pixelsize[0])
          self.ped = np.zeros([npix, npix])
          self.imgShape = [npix, npix]
        elif self.det.dettype == 6:
          if env.configStore().get(psana.TimeTool.ConfigV2, psana.Source(srcName)) is not None and env.configStore().get(psana.TimeTool.ConfigV2, psana.Source(srcName)).write_image()==0:
            self.ped=np.array([[-1,-1],[-1,-1]])
          else:
            fexCfg = env.configStore().get(psana.Camera.FrameFexConfigV1, psana.Source(srcName))
            self.ped = np.zeros([fexCfg.roiEnd().column()-fexCfg.roiBegin().column(), fexCfg.roiEnd().row()-fexCfg.roiBegin().row()])
          self.imgShape = self.ped.shape
        try:
          pedImg = self.det.image(run, self.ped)
          if pedImg is not None:
            self.imgShape = pedImg.shape
        except:
          pass
        try:
          self.mask = self.det.mask(run, unbond=True, unbondnbrs=True, status=True,  edges=True, central=True).squeeze()
          if rank==0:
            print 'masking %d pixel (status & edge,..) of %d'%(np.ones_like(self.mask).sum()-self.mask.sum(), np.ones_like(self.mask).sum())
          self.cmask = self.det.mask(run, unbond=True, unbondnbrs=True, status=True,  edges=True, central=True,calib=True).squeeze()
          if self.cmask.sum()!=self.mask.sum() and rank==0:
            print 'found user mask, masking %d pixel'%(np.ones_like(self.mask).sum()-self.cmask.sum())
        except:
          self.mask = None
          self.cmask = None
        if self.det.dettype == 6:
          self.common_mode = -1
        #geometry
        try:
          self.x = self.det.coords_x(run).squeeze()
          self.y = self.det.coords_y(run).squeeze()
        except:
          if self.det.dettype == 19:
            self.x = np.arange(0,self.ped.shape[0]*self.pixelsize[0], self.pixelsize[0])*1e6
            self.y = np.arange(0,self.ped.shape[0]*self.pixelsize[0], self.pixelsize[0])*1e6
            self.x, self.y = np.meshgrid(self.x, self.y)
          if self.det.dettype == 6:
            self.x = np.arange(0,self.ped.shape[0]*self.pixelsize[0], self.pixelsize[0])*1e6
            self.y = np.arange(0,self.ped.shape[1]*self.pixelsize[0], self.pixelsize[0])*1e6
            self.x, self.y = np.meshgrid(self.x, self.y)
        if self.det.dettype == 1 or self.det.dettype == 2 or self.det.dettype == 13: 
            try:
                iX,iY = self.det.indexes_xy(run)
                self.iX=np.array(iX)
                self.iY=np.array(iY)
            except:
                if rank==0:
                  print 'failed to get geometry info, likely because we do not have a geometry file'
                self.iX=self.x
                self.iY=self.y
        else:          
          self.iX=self.x
          self.iY=self.y

      else:
        self.rms = None
        self.mask = None
        self.wfx = None
        if self.det.dettype == 16: #acqiris
          cfg = env.configStore().get(psana.Acqiris.ConfigV1, psana.Source(srcName))
          self.interval = [cfg.horiz().sampInterval()]
          self.delayTime = [cfg.horiz().delayTime()]
          self.nSamples =  [cfg.horiz().nbrSamples()]
          self.fullScale =[]
          self.offset = []
          for c in cfg.vert():
            self.fullScale.append(c.fullScale())
            self.offset.append(c.offset())
        self.gain = None
      #common mode timing
      self.dataAccessTime=0.

      self.bankMasks=[]
      if self.common_mode==47:
        for i in range(0,16):
          bmask = np.zeros_like(self.rms)
          bmask[(i%2)*352:(i%2+1)*352,768/8*(i/2):768/8*(i/2+1)]=1
          self.bankMasks.append(bmask.astype(bool))
          #print i,'|', (i%2)*352,(i%2+1)*352,768/8*(i/2),768/8*(i/2+1),'|',bmask.sum(), self.bankMasks[i].sum(),'---',np.ones_like(self.rms)[352:704,0:96].sum()
      self.needsGeo=False
      if 'ped' in self.__dict__.keys() and self.ped.shape != self.imgShape:
        self.needsGeo=True

    def storeEnvironment(self):
      self.storeEnv = True
    def setMask(self, mask):
      self.mask = mask
    def setcMask(self, mask):
      self.cmask = np.amin(np.array([self.mask,mask]),axis=0)
    def setGain(self, gain):
      self.gain = gain
    def getAzAvs(self):
      return [ self[key] for key in self.__dict__.keys() if isinstance(self[key], ab.azimuthalBinning) ]
    def getAzAvKeys(self):
      return [ key for key in self.__dict__.keys() if isinstance(self[key], ab.azimuthalBinning) ]
    def saveFull(self):
      self.__dict__['full']=ROIObject([0,1e6], name='full', writeArea=True, rms=self.rms)
    def addROI(self, ROIname, ROI_limit, writeArea=False):
      self.__dict__[ROIname]=ROIObject(ROI_limit, name=ROIname, writeArea=writeArea, rms=self.rms)
    def getROIs(self):
      return [ self[key] for key in self.__dict__.keys() if isinstance(self[key], ROIObject) ]
    def processROIs(self):
      if self.evt.dat is None:
        print 'no data for', self._name,' , let mpiDataSource take care of this'
        return
      for ROI in self.getROIs():
        #print 'ROI: ',ROI.name
        if self.mask is not None:
          ROI.area=ROI.applyROI(np.ma.masked_array(self.evt.dat, ~(self.mask.astype(bool))))
        else:
          ROI.area=ROI.applyROI(np.ma.masked_array(self.evt.dat, ~(np.ones_like(self.evt.dat).astype(bool))))
        if ROI.writeArea:
          self.evt.__dict__['write_'+ROI.name] = ROI.area.squeeze()
        self.evt.__dict__['write_'+ROI.name+'_max'] = ROI.area.astype(np.float64).max()
        self.evt.__dict__['write_'+ROI.name+'_sum'] = ROI.area.astype(np.float64).sum()
        self.evt.__dict__['write_'+ROI.name+'_com'] = ROI.centerOfMass(ROI.area)
        for pj in ROI.getProjs():
          self.evt.__dict__['write_'+ROI.name+'_'+pj[-1]] = ROI.projection(ROI.area, pj[:-1])
          self.evt.__dict__['write_'+ROI.name+'_'+pj[-1]+'_max'] = self.evt.__dict__['write_'+ROI.name+'_'+pj[-1]].max()
          self.evt.__dict__['write_'+ROI.name+'_'+pj[-1]+'_sum'] = self.evt.__dict__['write_'+ROI.name+'_'+pj[-1]].sum()
        if ROI.rebin:
          self.evt.__dict__['write_'+ROI.name+'_rebin'] = ROI.Rebin(ROI.area)

    def getDroplets(self):
      return [ self[key] for key in self.__dict__.keys() if isinstance(self[key], droplet.droplet)  or isinstance(self[key], droplet.droplet) ]

    def processDroplets(self):
      for Drops in self.getDroplets():
        return_vals = Drops.dropletize(self.evt.dat)
        if return_vals is None:
          return_vals = Drops.ret_dict
        for key in return_vals:
          #print 'key in return val: write_%s'%(key)
          #try:
          #  print 'has shape: ',return_vals[key].shape, return_vals[key].dtype
          #except:
          #  pass
          if key.find('ragged')>=0:
            self.evt.__dict__['ragged_write_%s'%(key.replace('ragged_',''))] = return_vals[key]
          else:
            self.evt.__dict__['write_%s'%(key)] = return_vals[key]
        del Drops.ret_dict

    def getAcfs(self):
      return [ self[key] for key in self.__dict__.keys() if isinstance(self[key], acf.acf) ]

    def processAcf(self):
      for acf in self.getAcfs():
        if len(self.evt.dat.shape)==2:
          return_vals = acf.speckle_profile_image(self.evt.dat, acf.resolution)
        else:
          return_vals = acf.speckle_profile_image(self.det.image(int(run),self.evt.dat), acf.resolution)
        for key in return_vals.keys():
          self.evt.__dict__['write_%s_%s'%(acf.name,key)] = return_vals[key]

    def getFitCenters(self):
      return [ self[key] for key in self.__dict__.keys() if isinstance(self[key], fitCenter.fitCenter) ]

    def processFitCenters(self):
      for fitCenter in self.getFitCenters():
        return_vals = fitCenter.getCenter(self.evt.dat)
        for key in return_vals:
          self.evt.__dict__['write_%s_%s'%(fitCenter.name, key)] = return_vals[key]

    def getPhotons(self):
      return [ self[key] for key in self.__dict__.keys() if isinstance(self[key], photons.photon) ]

    def processPhotons(self):
      for photon in self.getPhotons():
        return_vals = photon.photon(self.evt.dat)
        for key in return_vals.keys():
          self.evt.__dict__['write_%s_%s'%(photon.name, key)] = return_vals[key]

    def getPhotons2(self):
      return [ self[key] for key in self.__dict__.keys() if isinstance(self[key], photons.photon2) ]

    def processPhotons2(self):
      for photon in self.getPhotons2():
        return_vals = photon.photon(self.evt.dat)
        for key in return_vals.keys():
          self.evt.__dict__['write_%s_%s'%(photon.name, key)] = return_vals[key]

    def getPhotons3(self):
      return [ self[key] for key in self.__dict__.keys() if isinstance(self[key], photons.photon3) ]

    def processPhotons3(self):
      for photon in self.getPhotons3():
        return_vals = photon.photon(self.evt.dat)
        for key in return_vals.keys():
          self.evt.__dict__['write_%s_%s'%(photon.name, key)] = return_vals[key]


    def processEnvironment(self):
      if self.storeEnv:
        if self.evt.envRow is not None:
          self.evt.__dict__['env_temp0']  = self.evt.envRow[0]*0.01
          self.evt.__dict__['env_temp1']  = self.evt.envRow[1]*0.01
          self.evt.__dict__['env_humidity']  = self.evt.envRow[2]*0.01

          self.evt.__dict__['env_AnalogI']  = self.evt.envRow[3]*0.001
          self.evt.__dict__['env_DigitalI']  = self.evt.envRow[4]*0.001
          self.evt.__dict__['env_GuardI']  = self.evt.envRow[5]*0.000001
          self.evt.__dict__['env_BiasI']  = self.evt.envRow[6]*0.000001
          self.evt.__dict__['env_AnalogV']  = self.evt.envRow[7]*0.001
          self.evt.__dict__['env_DigitalV']  = self.evt.envRow[8]*0.001

    def processDetector(self):
      self.processROIs()
      self.processDroplets()
      self.processAcf()
      self.processFitCenters()
      self.processPhotons()
      self.processPhotons2()
      self.processPhotons3()
      self.processEnvironment()
      # calculate azimuthal average if requested
      for thisAzavName,thisAzav in zip(self.getAzAvKeys(),self.getAzAvs()):
        if self.evt.dat is not None:
          self.evt.__dict__['write_%s'%thisAzavName] = thisAzav.doCake(self.evt.dat)
        else:
          #now let mpiDataSource take care of this.
          continue
          #nanShp = [ shpi-1 for shpi in self.azav.qbins.shape ]
          #nanAr = np.empty(nanShp)
          #nanAr[:]=np.nan
          #self.evt.write_azav = nanAr
      #add stuff for binned image.
      if 'binnedImgShp' in self.__dict__.keys():
        try:
          if len(self.evt.dat.shape)>2:
            img = self.det.image(self.evt.dat)
          else:
            img = self.evt.dat
          if self.orgImgShp is None:
            self.orgImgShp = img.shape
          self.evt.write_binnedImg = rebin(img, self.binnedImgShp)
        except:
          #now let mpiDataSource take care of this.
          pass
          #if self.orgImgShp is not None:
          #  if len(self.orgImgShp)==1:
          #    self.evt.write_binnedImg = np.empty([self.orgImgShp,self.orgImgShp])
          #  else:
          #    self.evt.write_binnedImg = np.empty(self.orgImgShp)
          #  self.evt.write_binnedImg[:] = np.nan
          #else:
          #  print 'the first event must have no image, come up with better code'

    def getMask(self, ROI):
      ROI = np.array(ROI)
      print 'DEBUG getMask: ',ROI.shape
      if ROI.shape != (2,2):
        return np.ones_like(self.ped) 

      mask_roi=np.zeros(self.imgShape)#<-- no, like image. Need image size.
      print 'DEBUG getMask: img shape ',self.imgShape, self.ped.shape
      mask_roi[ROI[0,0]:ROI[0,1],ROI[1,0]:ROI[1,1]]=1
      if self.needsGeo:
        mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(self.iX,self.iY)] )
      else:
        mask_r_nda = mask_roi
      print 'mask from rectangle (shape):',mask_r_nda.shape
      return mask_r_nda

    def addBinnedImg(self,shape):
      self.binnedImgShp = shape
      self.orgImgShp = None

    def addAzAv(self,phiBins=1,qBin=1e-2,center=None,dis_to_sam=None, eBeam=None, azavName='azav'):
      azavMask = ~(self.cmask.astype(bool)&self.mask.astype(bool))
      if rank==0:
        print 'mask %d pixel for azimuthal integration'%azavMask.sum()
      if dis_to_sam is None:
        dis_to_sam = self.__dict__[azavName+'_dis_to_sam']
      if center is None:
        center = self.__dict__[azavName+'_center']
      if eBeam is None:
        eBeam = self.__dict__[azavName+'_eBeam']
      self.__dict__[azavName] = ab.azimuthalBinning(x=self.x.flatten()/1e3,y=self.y.flatten()/1e3,xcen=center[0]/1e3,ycen=center[1]/1e3,d=dis_to_sam,mask=azavMask.flatten(),lam=E2lam(eBeam)*1e10,Pplane=1,phiBins=phiBins,qbin=qBin)
      self.__dict__[azavName+'_q'] = self.__dict__[azavName].q
      self.__dict__[azavName+'_correction'] = self.__dict__[azavName].correction
      self.__dict__[azavName+'_norm'] = self.__dict__[azavName].norm
      self.__dict__[azavName+'_normPhi'] = self.__dict__[azavName].Cake_norm
      self.__dict__[azavName+'_phi'] = self.__dict__[azavName].phiVec

    def addDroplet(self,threshold=5.0, thresholdLow=3., thresADU=71., name='droplet', useRms=True, ROI=[], relabel=True, boolMask=False):
      if len(ROI)>0:
        print 'ROI DEBUG: ',self.cmask.shape,self.getMask(ROI).shape
        mask = ( self.mask.astype(bool) & self.cmask.astype(bool) & self.getMask(ROI).astype(bool) )
      else:
        mask = ( self.cmask.astype(bool) & self.mask.astype(bool) )
      self.__dict__[name] = droplet.droplet(threshold=threshold, thresholdLow=thresholdLow, thresADU=thresADU, rms=self.rms, mask=mask, name=name,useRms=useRms,relabel=relabel)
      self.__dict__[name+'_threshold'] = threshold
      self.__dict__[name+'_thresholdLow'] = thresholdLow
      self.__dict__[name+'_thresADU'] = thresADU
      
      for aduHist in self.__dict__[name].aduHists:
        self.__dict__[name+'_'+aduHist.name] = aduHist.bins

    def addAcf(self,resolution=0.1):
      if rank==0:
        print 'defined autocorrelation function'
      self.acf = acf.acf(resolution=0.1)
      self.acf_resolution = resolution

    def addFitCenter(self,maskName='', threshold=90, name=''):
      if rank==0:
        print 'define event-by-event center fitting'
      self.fitCenter = fitCenter.fitCenter(maskName=maskName, threshold=threshold, name=name, xArr=self.x, yArr =self.y )
      self.fitCenter_threshold = threshold
      self.fitCenter_maskName = maskName

    def addPhotons(self, ADU_per_photon=154, mask=None, rms=None, name='photon', nphotMax=25, retImg=False, nphotRet=100):
      self.__dict__[name] = photons.photon(ADU_per_photon=ADU_per_photon, mask=self.cmask, rms=self.rms, nphotMax=nphotMax, retImg=retImg, nphotRet=nphotRet,name=name)

    def addPhotons2(self, ADU_per_photon=154, thresADU=0.7, thresRms=3., mask=None, rms=None, name='photon', nphotMax=25, retImg=False, nphotRet=100): 
      self.__dict__[name] = photons.photon2(ADU_per_photon=ADU_per_photon, thresADU=thresADU, thresRms=thresRms, mask=self.cmask, rms=self.rms, nphotMax=nphotMax, retImg=retImg, nphotRet=nphotRet,name=name)

    def addPhotons3(self, ADU_per_photon=154, thresADU=0.9, thresRms=3., mask=None, rms=None, name='photon', nphotMax=25, retImg=False, nphotRet=100, maxMethod=0):
      self.__dict__[name] = photons.photon3(ADU_per_photon=ADU_per_photon, thresADU=thresADU, thresRms=thresRms, mask=self.cmask, rms=self.rms, nphotMax=nphotMax, retImg=retImg, nphotRet=nphotRet,name=name,maxMethod=maxMethod)

    def getData(self,evt):
      time_start = MPI.Wtime()
      if self.det.name.__str__().find('OceanOptics')>=0:
        self.evt.dat = self.det.intensity(evt)
        if self.wfx is None:
            self.wfx = self.det.wavelength(evt)
        return
      if self.det.dettype == 16:
        self.evt.dat = self.det.waveform(evt)
        if self.wfx is None:
            self.wfx = self.det.wftime(evt)
        return
      #print 'still trying to get data for ',self._name
      if self.common_mode<0:
          self.evt.dat = self.det.raw_data(evt)
      elif self.common_mode==0:
        self.evt.dat = self.det.raw_data(evt)-self.ped
        #apply mask if requested
        if self.applyMask==1:
          self.evt.dat[self.mask==0]=0
        if self.applyMask==2:
          self.evt.dat[self.cmask==0]=0
      elif self.common_mode%100==1:
        self.evt.dat = self.det.calib(evt, cmpars=(1,25,40,100))
      elif self.common_mode%100==5:
        self.evt.dat = self.det.calib(evt, cmpars=(5,100))
      elif self.common_mode%100==4:
        self.evt.dat = self.det.calib(evt, cmpars=(4,6,30,30))
      elif self.common_mode%100==45:
        self.evt.dat = self.det.raw_data(evt)-self.ped
        self.evt.dat[self.mask==0]=0
        self.evt.dat = cm_epix(self.evt.dat,self.rms)
      elif self.common_mode%100==46:
        self.evt.dat = self.det.raw_data(evt)-self.ped
        self.evt.dat[self.mask==0]=0
        self.evt.dat = cm_epix(self.evt.dat,self.rms,normAll=True)
      elif self.common_mode%100==47:
        self.evt.dat = self.det.raw_data(evt)-self.ped
        self.evt.dat[self.mask==0]=0
        for ib,bMask in enumerate(self.bankMasks):
          #print ib, np.median(self.evt.dat[bMask]), bMask.sum()
          self.evt.dat[bMask]-=np.median(self.evt.dat[bMask])
        self.evt.dat = cm_epix(self.evt.dat,self.rms)
      elif self.common_mode%100==55:
        self.evt.dat = self.det.calib(evt, cmpars=(5,5000))
      elif self.common_mode%100==10:
        #data = self.det.raw_data(evt)-self.det.pedestals(evt)        
        data = self.det.raw_data(evt)-self.ped
        if self.applyMask==1:
          data[self.mask==0]=0
        if self.applyMask==2:
          data[self.cmask==0]=0
        data_def = self.det.calib(evt, cmpars=(1,25,40,200))
        data_unb = self.det.calib(evt, cmpars=(5,100))
        data_diff = data_unb-data_def
        data_diff[(data_def-data)!=0]=0
        self.evt.dat = data_def + data_diff
        #tileAvs = [ tile[tile!=0].flatten().mean() for tile in data]
      else:
        self.evt.dat = self.det.calib(evt)
      if self.gain is not None and self.gain.shape == self.evt.dat.shape and self.common_mode>=0 and self.common_mode<100:
        self.evt.dat*=self.gain
      self.dataAccessTime+=MPI.Wtime()-time_start
      #store environmental row if desired.
      if self.storeEnv:
        try:
          e = evt.get(psana.Epix.ElementV3,psana.Source(self._srcName))
          self.evt.envRow = e.environmentalRows()[1]
        except:
          self.evt.envRow = None
        
