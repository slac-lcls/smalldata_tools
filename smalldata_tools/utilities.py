import os
import json
import copy
import numpy as np
import time
import h5py
import tables
from scipy import optimize
from scipy import ndimage
from matplotlib import pyplot as plt
import resource

import psana
import RegDB.experiment_info


def MAD(a, c=0.6745, axis=None):
  """
  Median Absolute Deviation along given axis of an array:
  
  median(abs(a - median(a))) / c
  
  c = 0.6745 is the constant to convert from MAD to std; it is used by
  default  
  """
  
  a = np.ma.masked_where(a!=a, a)
  if a.ndim == 1:
    d = np.ma.median(a)
    m = np.ma.median(ma.fabs(a - d) / c)
  else:
    d = np.ma.median(a, axis=axis)
    # I don't want the array to change so I have to copy it?
    if axis > 0:
      aswp = np.ma.swapaxes(a,0,axis)
    else:
      aswp = a
      m = np.ma.median(ma.fabs(aswp - d) / c, axis=0)
      
  return m
    
def nanmedian(arr, **kwargs):
      """
      Returns median ignoring NAN
      """
      return np.ma.median( np.ma.masked_where(arr!=arr, arr), **kwargs )
    
def rebinFactor(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def rebin(a, shape):
    if isinstance(shape, float) or isinstance(shape, int):
      shape = [shape, shape]
    if (a.shape[0]%shape[0]) == 0 and (a.shape[1]%shape[1]) == 0:
      rebinFactor(a, shape)
    else:
      factor = [ float(int(a.shape[0]/shape[0])+1), float(int(a.shape[1]/shape[1])+1)]
      bigImg = ndimage.zoom(a, [shape[0]*factor[0]/a.shape[0],shape[1]*factor[1]/a.shape[1]])
      img = rebinFactor(bigImg, shape)
    return img

def reduceVar(vals, sigROI,threshold=-1e25):
    if threshold!=-1e25:
      vals = vals[vals<threshold]=0
    print 'array shape: ',len(vals.shape)
    if len(vals.shape)>1 and sigROI!=[]:
      if len(vals.shape)==2:
        if not isinstance(sigROI, list):
          return vals[:,sigROI]
        elif len(sigROI)>1:
          return vals[:,sigROI[0]:sigROI[1]]
        else:
          return vals[:,sigROI[0]]
      elif len(vals.shape)==3:
        if not isinstance(sigROI, list):
          return vals[:,sigROI,:]
        elif len(sigROI)==1:
          return vals[:,sigROI[0],:]
        elif len(sigROI)==2:
          return vals[:,sigROI[0]:sigROI[1],:]
        else:
          return vals[:,sigROI[0]:sigROI[1],sigROI[2]:sigROI[3]]
      elif len(vals.shape)==4:
        if not isinstance(sigROI, list):
          return vals[:,sigROI,:,:]
        elif len(sigROI)==1:
          return vals[:,sigROI[0],:,:]
        elif len(sigROI)==2:
          return vals[:,sigROI[0]:sigROI[1],:,:]
        elif len(sigROI)==4:
          return vals[:,sigROI[0]:sigROI[1],sigROI[2]:sigROI[3],:]
        else:
          return vals[:,sigROI[0]:sigROI[1],sigROI[2]:sigROI[3],sigROI[4]:sigROI[5]]

    print 'this dimension is not yet implemented:',vals.shape,' ROI: ',sigROI
    return vals
    
#from ixppy
def E2lam(E,o=0):
    """ Computes photon wavelength in m
        E is photon energy in eV or keV
        set o to 0 if working at sub-100 eV energies
    """
    if o:
      E=E
    else:
      E=eV(E)
    lam=(12398.4/E)*1e-10
    return lam
  #lam = 12.39842 /E
  #return lam

def eV(E):
  if E < 100:
    E=E*1000.0;
    return E*1.0

def checkDet(env, detname):
  for key in env.configStore().keys():
    if key.alias()==detname:
      return True
  for key in env.configStore().keys():
    #print key.src().__repr__(),detname,key.src().__repr__()==detname
    if key.src().__repr__()==detname:
      return True
  return False
        
def printMsg(eventNr, run, rank):
  printFreq = 10
  #if eventNr > 10000:
  #  printFreq = 10000
  if eventNr > 1000:
    printFreq = 1000
  elif eventNr > 120:
    printFreq = 100
    
  if eventNr%printFreq == 0:
    if rank == 0:
      usage = resource.getrusage(resource.RUSAGE_SELF)
      print "*** In Event: run", run, ",event# =", eventNr,' memory used: ',usage[2]*resource.getpagesize()/1000000.,' at ',time.strftime('%X')

def getExpName(env):  
  if (env.jobName()).find('shmem')>=0:
    return RegDB.experiment_info.active_experiment('XPP')[1]
  else:
    return env.experiment()

##########################################################################################
###  helper classes & functions
##########################################################################################


def fitCircle(x,y):
  def calc_R(x,y, xc, yc):
    return np.sqrt((x-xc)**2 + (y-yc)**2)

  def f(c,x,y):
    Ri = calc_R(x,y,*c)
    return Ri - Ri.mean()

  x_m = np.mean(x)
  y_m = np.mean(y)
  center_estimate = x_m, y_m
  center, ier = optimize.leastsq(f, center_estimate, args=(x,y))
  xc, yx = center
  Ri     = calc_R(x, y, *center)
  R      = Ri.mean()
  residu = np.sum((Ri - R)**2)
  return xc, yx, R, residu
    
# --------------------------------------------------------------    

def hasKey(inkey, inh5=None, printThis=False):
    hasKey = False
    if inh5 is None:
      print 'no input file given'
      return hasKey
    if not isinstance(inh5, tables.file.File):
      print 'input file is not a tables h5file'
      return hasKey

    try:
      if inkey[0]==['/']:
        inkey=inkey[1:]
      if inkey.find('/')>=0:
        inh5.get_node('/'+inkey.split('/')[0],inkey.split('/')[1])
      else:
        inh5.get_node('/'+inkey)
      hasKey=True
    except:
      pass
    return hasKey

def getVar(fh5, plotvar):
  if not hasKey(plotvar,fh5):
    print 'signal variable %s not in littleData file'%(plotvar)
    return

  if plotvar[0]!=['/']:
    plotvar='/'+plotvar
  try:
    if len(plotvar.split('/'))>2:
      vals = fh5.get_node('/'.join(plotvar.split('/')[:-1]),plotvar.split('/')[-1]).read()
    else:
      vals = fh5.get_node(plotvar).read()
    return vals.squeeze()
  except:
    print 'failed to get data for ',plotvar
    return

def getTTstr(fh5):
  """
  function to determine the string for the timetool variables in the desired run
  necessary as naming scheme evolved over time
  """
  ttCorr = None
  ttBaseStr = 'tt/'
  if hasKey('tt/ttCorr',fh5):
    ttCorr = 'tt/ttCorr'
  elif hasKey('ttCorr/tt',fh5):
    ttCorr = 'ttCorr/tt'
  if not hasKey(ttBaseStr+'AMPL',fh5):
    if hasKey('tt/XPP_TIMETOOL_AMPL',fh5):
      ttBaseStr = 'tt/XPP_TIMETOOL_'
    elif hasKey('tt/TIMETOOL_AMPL',fh5):
      ttBaseStr = 'tt/TIMETOOL_'
    elif hasKey('tt/TTSPEC_AMPL',fh5):
      ttBaseStr = 'tt/TTSPEC_'
  return ttCorr, ttBaseStr

def getDelay(fh5, use_ttCorr=True, addEnc=False):
  """
  function to get the xray-laser delay from the data
  usage:
  getDelay(fh5): get the delay from lxt and/or encoder stage, add the timetool correction
  getDelay(fh5, use_ttCorr=False): get the delay from lxt and/or encoder stage, NO timetool correction
  getDelay(fh5, addEnc=True): get the delay from lxt, add encoder stage and timetool correction
  fh5 is the pointer to the h5py or tables objext
  """
  ttCorrStr, ttBaseStr = getTTstr(fh5)
  if ttCorrStr is not None:
    ttCorr=getVar(fh5,ttCorrStr)
    if (np.nanstd(ttCorr)==0):
      ttCorr=getVar(fh5,ttBaseStr+'FLTPOS_PS')
  nomDelay=np.zeros_like(ttCorr)

  isDaqDelayScan=False
  for node in fh5.get_node('/scan')._f_list_nodes():
    if node.name.find('var')<0 and node.name.find('none')<0 and node.name.find('lxt')>=0:
      isDaqDelayScan=True
      #print 'DEBUG: found that we have a delay scan'
      nomDelay=node.read()*1e12

  if not isDaqDelayScan:
    if hasKey('enc/lasDelay',fh5):
      encVal = getVar(fh5,'enc/lasDelay')
      #print 'DEBUG: encoder info',encVal.std()
      if encVal.std()>1e-9:
        nomDelay=encVal
        addEnc=False
      elif encVal.std()>1e-15:
        nomDelay=encVal*1e12
        addEnc=False
    elif hasKey('enc/ch0',fh5):
      encVal = getVar(fh5,'enc/ch0')
      if encVal.std()>1e-15 and encVal.std()<1e-9:
        nomDelay=encVal*1e12
        #now look at the EPICS PV if everything else has failed.
    else:
      epics_delay = getVar(fh5,'epics/lxt_ttc')
      if epics_delay.std()!=0:
        nomDelay = epics_delay

  if addEnc and not hasKey('enc/lasDelay',fh5):
    print 'required to add encoder value, did not find encoder!'
  if addEnc and hasKey('enc/lasDelay',fh5):            
    if getVar(fh5,'enc/lasDelay').std()>1e-6:
      nomDelay+=getVar(fh5,'enc/lasDelay').value

  if use_ttCorr:
    #print 'DEBUG adding ttcorr,nomdelay mean,std: ',ttCorr.mean(),nomDelay.mean(),ttCorr.std(),nomDelay.std()
    return (ttCorr+nomDelay)
  else:
    return nomDelay


def addToHdf5(fh5, key, npAr):
    arShape=()
    for i in range(0,len(npAr.shape)):
        arShape+=(npAr.shape[i],)
    dset = fh5.create_dataset(key, arShape)
    dset[...] = npAr.astype(float)


def getBins(bindef=[], fh5=None, debug=False):
  #have full list of bin boundaries, just return it
  if len(bindef)>3:
    return bindef
  
  if len(bindef)==2:
    if debug:
      print 'only two bin boundaries, so this is effectively a cut...cube will have a single image'
    Bins = np.array([min(bindef[0],bindef[1]),max(bindef[0],bindef[1])])
    return Bins

  if len(bindef)==3:
    if type(bindef[2]) is int:
        Bins=np.linspace(min(bindef[0],bindef[1]),max(bindef[0],bindef[1]),bindef[2]+1,endpoint=True)
    else:
        Bins=np.arange(min(bindef[0],bindef[1]),max(bindef[0],bindef[1]),bindef[2])
    if Bins[-1]<bindef[1]:
      Bins = np.append(Bins,max(bindef[0],bindef[1]))
    return Bins

  #have no input at all, assume we have unique values in scan. If not, return empty list 
  scanVarName = ''
  if len(bindef)<=1:
    for key in fh5['scan'].keys():
      if key.find('var')<0 and key.find('none')<0:
        scanVarName = key

  if len(bindef)==0:
    if scanVarName=='':
      print 'this run is no scan, will need bins as input, quit now'
      return []
    print 'no bins as input, we will use the scan variable %s '%scanVarName

    Bins = np.unique(fh5['scan/var0'])
    if scanVarName.find('lxt')>=0:
      Bins*=1e12
    if debug:
      print 'Bins: ',Bins
    return Bins

  #give a single number (binwidth or numBin)
  if len(bindef)==1:
    #this is a "normal" scan, use scanVar
    if scanVarName!='':
      Bins = np.unique(fh5['scan/var0'])
      if scanVarName.find('lxt')>=0:
        Bins*=1e12
      valBound = [ min(Bins), max(Bins)]
      if type(bindef[0]) is int:
        Bins=np.linspace(valBound[0],valBound[1],bindef[0],endpoint=True)
      else:
        Bins=np.arange(valBound[0],valBound[1],bindef[0])
        Bins=np.append(Bins,valBound[1])
      return Bins
      
    else:
      if hasKey('enc/ch0', fh5) or hasKey('enc/lasDelay',fh5):
        if hasKey('enc/ch0', fh5):
          vals = fh5['enc/ch0'].value.squeeze()*1e12
        else:
          vals = fh5['enc/lasDelay'].value.squeeze()
        minEnc = (int(min(vals)*10.))
        if minEnc<0:
            minEnc+=-1
        minEnc /= 10.
        maxEnc = (int(max(vals)*10.)+1)/10.
        print minEnc,maxEnc,bindef[0]
        if minEnc!=maxEnc and abs(minEnc)<101 and abs(maxEnc<1001):
          if type(bindef[0]) is int:
            Bins=np.linspace(minEnc, maxEnc, bindef[0],endpoint=True)
          else:
            Bins=np.arange(minEnc, maxEnc, bindef[0])
            if Bins[-1]< maxEnc:
              Bins=np.append(Bins,maxEnc)
          if debug:
            print 'Bins....',Bins
          return Bins
        else:
          print 'you passed only one number and the this does not look like a new delay scan or a normal scan'
          return []

    print 'why am I here? I should have hit one of the if/else staments....'
    print 'you passed: ',bindef
###
# utility functions for droplet stuff
###
def gaussian(x, amp, cen, wid):
                return amp * exp(-(x-cen)**2 /(2*wid**2))
def lorentzian(x,p0,p1):
		return (p0**2)/(p0**2 + (x-p1)**2)

def neighborImg(img):
    img_up = np.roll(img,1,axis=0); img_up[0,:]=0
    img_down = np.roll(img,-1,axis=0); img_down[-1,:]=0
    img_left = np.roll(img,1,axis=1); img_left[:,0]=0
    img_right = np.roll(img,-1,axis=1); img_right[:,-1]=0
    return np.amax(np.array([img_up, img_down, img_left, img_right]),axis=0)

def cm_epix(img,rms,maxCorr=30, histoRange=30, colrow=3, minFrac=0.25, normAll=False):
    #make a mask: all pixels > 10 rms & neighbors & pixels out of historange
    imgThres = img.copy()
    imgThres[img>=rms*10]=1
    imgThres[img<rms*10]=0
    imgThres+=neighborImg(imgThres)
    imgThres+=imgThres+(abs(img)>histoRange)

    maskedImg = np.ma.masked_array(img, imgThres)
    if normAll:
      #this should be done in the 16 blocks.
      maskedImg -= maskedImg.mean()
    if colrow%2==1:
        rs = maskedImg.reshape(704/2,768*2,order='F')
        rscount = np.ma.count_masked(rs,axis=0)
        rsmed = np.ma.median(rs,axis=0)
        rsmed[abs(rsmed)>maxCorr]=0
        rsmed[rscount>((1.-minFrac)*352)]=0
        imgCorr = np.ma.masked_array((rs.data-rsmed[None,:]).data.reshape(704,768,order='F'),imgThres)
    else:
        imgCorr = maskedImg.copy()

    if colrow>=2:
        rs = imgCorr.reshape(704*8,96)
        rscount = np.ma.count_masked(rs,axis=1)
        rsmed = np.ma.median(rs,axis=1)
        rsmed[abs(rsmed)>maxCorr]=0
        rsmed[rscount>((1.-minFrac)*96)]=0
        imgCorr = (np.ma.masked_array(rs.data-rsmed[:,None], imgThres)).reshape(704,768)
                                                
    return imgCorr.data  

###
# utility functions for plotting data as 2-d histogram
###
def hist2d(ar1, ar2,limits=[1,99.5],numBins=[100,100],histLims=[np.nan,np.nan, np.nan, np.nan],weights=None, doPlot=True):
        pmin0 = np.nanmin(ar1); pmin1 = np.nanmin(ar2)
        pmax0 = np.nanmax(ar1); pmax1 = np.nanmax(ar2)
        if not np.isnan(np.percentile(ar1,limits[0])):
            pmin0 = np.percentile(ar1,limits[0])
        if not np.isnan(np.percentile(ar2,limits[0])):
            pmin1 = np.percentile(ar2,limits[0])
        if limits[1]<100:
            if not np.isnan(np.percentile(ar1,limits[1])):
                pmax0 = np.percentile(ar1,limits[1])
            if not np.isnan(np.percentile(ar2,limits[1])):
                pmax1 = np.percentile(ar2,limits[1])
        if histLims[0] is not np.nan:
            pmin0 = histLims[0]
            pmax0 = histLims[1]
            pmin1 = histLims[2]
            pmax1 = histLims[3]
        v0 = ar1
        v1 = ar2
        binEdges0 = np.linspace(pmin0, pmax0, numBins[0])
        binEdges1 = np.linspace(pmin1, pmax1, numBins[1])
        ind0 = np.digitize(v0, binEdges0)
        ind1 = np.digitize(v1, binEdges1)
        ind2d = np.ravel_multi_index((ind0, ind1),(binEdges0.shape[0]+1, binEdges1.shape[0]+1)) 
        if weights is None:
                iSig = np.bincount(ind2d, minlength=(binEdges0.shape[0]+1)*(binEdges1.shape[0]+1)).reshape(binEdges0.shape[0]+1, binEdges1.shape[0]+1) 
        else:
                iSig = np.bincount(ind2d, weights=weights, minlength=(binEdges0.shape[0]+1)*(binEdges1.shape[0]+1)).reshape(binEdges0.shape[0]+1, binEdges1.shape[0]+1)    
        if doPlot:
          plt.imshow(iSig,aspect='auto', interpolation='none',origin='lower',extent=[binEdges1[1],binEdges1[-1],binEdges0[1],binEdges0[-1]],clim=[np.percentile(iSig,limits[0]),np.percentile(iSig,limits[1])])
          plt.colorbar()
        return iSig

###
def dictToHdf5(filename, indict):
  f = h5py.File(filename,'w')
  for key in indict.keys():
    npAr = np.array(indict[key])
    dset = f.create_dataset(key, npAr.shape, dtype='f')
    dset[...] = npAr
  f.close()


#
# should work on removing it
#
class dropObject(object):
  def __init__(self,name='noname',parent=None):
    self._name = name
    self._parent = parent

  def add(self,name,data):
    if (name not in self.__dict__):
      self._add(name,data)
    else:
      self.__dict__[name].append(data)
  def _add(self,name,data):
    self.__dict__[name]=[data]
  def addField(self,name,data):
    if (name not in self.__dict__):
      if isinstance(data, list) or  isinstance(data, np.ndarray):
        self.__dict__[name]=data
      else:
        self.__dict__[name]=[data]
    else:
      print 'field ',name,' already in dropObject: ',self._name
  def __repr__(self):
    return "dropObject with fields: "+str(self.__dict__.keys())
  def __getitem__(self,x):
    return self.__dict__[x]
  def __setitem__(self,name,var,setParent=True):
    self._add(name,var)
    if setParent:
      try:
        self[name]._parent = self
      except:
	pass
  def _get_keys(self):
    return [tk for tk in self.__dict__.keys() if not tk[0]=='_']


def shapeFromKey_h5(fh5, thiskey):
    if len(thiskey.split('/')) > 2:
        return fh5.get_node('/'.join(thiskey.split('/')[:-1]),thiskey.split('/')[-1]).shape
    else:
        return fh5.get_node(thiskey).shape

