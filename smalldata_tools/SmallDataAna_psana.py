import numpy as np
try:
    from pylab import ginput
except:
    print('could not import ginput')
    pass
try:
    from matplotlib import pyplot as plt
except:
    print('could not import pyplot')
    pass
from matplotlib import gridspec
from matplotlib import path
import itertools
import os
import socket
import holoviews as hv
import bokeh.plotting as bp
from matplotlib import gridspec
import psana
import logging
import time
import requests
try:
    basestring
except NameError:
    basestring = str
try:
    raw_input
except NameError:
    raw_input = input


import smalldata_tools.SmallDataAna as sda

from smalldata_tools.DetObject import DetObject
from smalldata_tools.SmallDataUtils import getUserData
from smalldata_tools.utilities import printR
from smalldata_tools.utilities import addToHdf5
from smalldata_tools.utilities import rename_reduceRandomVar
from smalldata_tools.utilities_plotting import plotImageBokeh
from smalldata_tools.utilities_plotting import hv_image
from smalldata_tools.utilities_plotting import hv_image_ctl
from smalldata_tools.utilities_plotting import hv_3dimage
from smalldata_tools.utilities_FitCenter import FindFitCenter
from smalldata_tools.utilities_FitCenter import fitCircle
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, sparsifyFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning

import smalldata_tools.cube.cube_mpi_fun as mpi_fun

from mpi4py import MPI
import h5py
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SmallDataAna_psana(object):
    def __init__(self, expname='', run=-1, dirname='', filename='', plotWith='matplotlib'):
        self.run=int(run)
        self.expname=expname
        self.hutch=expname[0:3]
        self.plotWith=plotWith
        self.sda=None

        ws_url = "https://pswww.slac.stanford.edu/ws/lgbk"
        if self.hutch == 'cxi':
            print('Will assume the first CXI station, if this is wrong, make sure to add  -e <expname> on commandline')
        resp = requests.get(ws_url + "/lgbk/ws/activeexperiment_for_instrument_station", {"instrument_name": self.hutch.upper(), 
                                                                                          "station": 0})
        try:
            currExpname = resp.json().get("value", {}).get("name")
            print('Current experiment for %s is %s'%(self.hutch, currExpname))
            rundoc = requests.get(ws_url + "/lgbk/" + expname  + "/ws/current_run").json()["value"]
            if not rundoc:
                logger.error("Invalid response from server")
            lastRun = int(rundoc['num'])
            if self.run > lastRun:
                printR(rank, 'experiment %s does only have %d runs, requested %d'%(expname, lastRun, self.run))
                return None
        except:
            lastRun = -1
    
        isLive = False
        if self.run == lastRun:
            end_time = rundoc.get('end_time', None)
            if end_time is None:
                isLive = True
        self.dsname='exp=%s:run=%i:smd'%(expname,self.run)
        hostname = socket.gethostname()
        if hostname.find('drp-srcf')>=0:
            xtcdirname='/cds/data/drpsrcf/%s/%s/xtc'%(self.hutch.lower(),expname)
            self.dsname+=':dir=%s'%xtcdirname
        self.dsnameIdx=self.dsname.replace('smd','idx').replace(':live','')
        if isLive:
            self.dsname=self.dsname+':live:stream=0-79'
            self.dsnameIdx=None
        comm.bcast(self.dsnameIdx, root=0)

        printR(rank, 'make SmallDataAna_psana from dsname: %s'%self.dsname)
        try:
            self.ds = psana.DataSource(str(self.dsname))
        except:
            printR(rank, 'Failed to set up small data psana dataset!')
            self.ds = None
        if self.dsnameIdx is None:
            printR(rank, 'Failed to set up index based psana dataset, likely because this run is still live')
        else:
            try:
                self.dsIdx = psana.DataSource(str(self.dsnameIdx))
            except:
                self.dsIdx = None
                printR(rank, 'Failed to set up index based psana dataset from dsname: %s'%self.dsnameIdx)
            if self.dsIdx is not None:
                try:
                    self.dsIdxRun = self.dsIdx.runs().next()
                except:
                    try:
                        self.dsIdxRun = next(self.dsIdx.runs())
                    except:
                        self.dsIdxRun = None
                        printR(rank, 'Failed to get run from psana dataset from dsname: %s'%self.dsnameIdx)
        printR(rank, 'try to make SmallDataAna using dirname %s, for exp %s and run %s'%(dirname,expname,run))
        try:
            printR(rank, 'setting up SmallData ana from anaps ')
            self.sda = sda.SmallDataAna(expname,self.run,dirname=dirname, filename=filename,plotWith=plotWith)
        except:
            printR(rank, 'failed, set anaps.sda to None')
            self.sda = None
        self.calibhisto={}
        self.jobsIds = []
        self.commonModeStrings=['raw','pedSub','unb','hist','histAmiLike','median','medianNoNorm','medianSN','median45','cm47','cm71','cm72','cm10','cm145','cm146','cm147','cm110','calib','cm80','cm81']

    def commonModeStr(self, common_mode=0):
        if common_mode<0:
            return 'raw_'
        elif common_mode==5:
            return 'unb_'
        elif common_mode==4 or common_mode==1:
            return 'hist_'
        elif common_mode==6:
            return 'median_'
        elif common_mode==30:
            return 'calib_'
        elif common_mode==34:
            return 'histAmiLike_'
        elif common_mode==36:
            return 'medianNoNorm_'
        elif common_mode==45:
            return 'median45_'
        elif common_mode==46:
            return 'medianSN_'
        elif common_mode==47:
            return 'cm47_'
        elif common_mode==71:
            return 'cm71_'
        elif common_mode==72:
            return 'cm72_'
        elif common_mode==80:
            return 'cm80_'
        elif common_mode==81:
            return 'cm81_'
        elif common_mode==10:
            return 'cm10_'
        elif common_mode==105:
            return 'cm105_'
        elif common_mode==145:
            return 'cm145_'
        elif common_mode==146:
            return 'cm146_'
        elif common_mode==147:
            return 'cm147_'
        elif common_mode==110:
            return 'cm110_'
        elif common_mode==0:
            return 'pedSub_'
        else:
            return ''

    def resetDs(self, idx=True):
        if idx:
            del self.dsIdx
            self.dsIdx = psana.DataSource(self.sda.dsnameIdx)
            self.dsIdxRun = self.dsIdx.runs().next()
        else:
            del self.ds
            self.ds = psana.DataSource(self.sda.dsname)

    def Keys(self,printthis=False):
        try:
            keys=[]
            for k in self.dsIdxRun.env().configStore().keys():
                if k.alias() is not None and k.alias()!='':
                    newKey=True
                    for oldKey in keys:
                        if oldKey.alias() == k.alias():
                            newKey=False
                    if newKey:
                        keys.append(k)
                else:
                    newKey=True
                    for oldKey in keys:
                        if oldKey.src().__repr__() == k.src().__repr__():
                            newKey=False
                    if newKey:
                        keys.append(k)
            #keys=evt.keys()
            if printthis: 
                print(keys)
                return
            else:
                return keys
        except:
            printR(rank, 'failed to use last event in idx mode, try smlData')
            try:
                keys=self.ds.events().next().keys()
                if printthis: 
                    print(keys)
                    return
                else:
                    return keys
            except:
                printR(rank, 'WARNING: smd data also did not work, give up. No keys returned! ')

    def CfgKeys(self,idx=True,printthis=False):
        if idx:
            keys=self.dsIdx.env().configStore().keys()
        else:
            keys=self.ds.env().configStore().keys()
        if printthis: 
            print(keys)
        else:
            return keys

    def EpicsAliases(self,idx=True,printthis=False):
        if idx:
            keys=self.dsIdx.env().epicsStore().aliases()
        else:
            keys=self.ds.env().epicsStore().aliases()
        if printthis: 
            print(keys)
        else:
            return keys
        
    def plotVar(self, plotvar, numBins=[100], useFilter=False, limits=[1,99],fig=None,asHist=False):
        self.sda.plotVar(plotvar=plotvar, numBins=numBins, useFilter=useFilter, limits=limits,fig=fig,asHist=asHist)

    def plotScan(self, ttCorr=False, sig='diodeU/channels', sigROI=[], i0='ipm3/sum', numBins=100):
        self.sda.plotScan(ttCorr=ttCorr, sig=sig, sigROI=sigROI, i0=i0, numBins=numBins)

#
# these functions need psana as well, make separate class that imports SmallDataAna?
#
    def addDetInfo(self, detname='None', common_mode=None):
        if detname=='None':
            aliases = self._getDetName()
            if len(aliases)==1:
                detname = aliases[0]
            else:
                detsString='detectors in event: \n'
                for alias in aliases:
                    detsString+=alias+', '
                print(detsString)
                detname = raw_input("Select detector to get detector info for?:\n")
        printR(rank, 'try to make psana Detector with: %s'%detname)

        detnameDict=None
        if isinstance(detname, dict):
            detnameDict = detname
            if 'common_mode' in detnameDict.keys() and common_mode is not None:
                common_mode=detnameDict['common_mode']
            detname = detnameDict['source']

        #only do this if information is not in object yet
        if detname in self.__dict__.keys() and self.__dict__[detname].common_mode==common_mode and detnameDict is not None:
            return detname

        if detname in self.__dict__.keys() and common_mode is None:
            return detname

        if detname in self.__dict__.keys():
            printR(rank, 'redefine detector object with different common mode: %d instead of %d'%( common_mode,self.__dict__[detname].common_mode))
        det = DetObject(detname , self.dsIdx.env(), self.run, name=detname, common_mode=common_mode)
        self.__dict__[detname]=det
        if (detname+'_pedestals') in self.__dict__.keys():
            return detname
        self.__dict__[detname+'_pedestals'] = det.ped
        self.__dict__[detname+'_rms'] = det.rms
        self.__dict__[detname+'_x']=det.x
        self.__dict__[detname+'_y']=det.y
        self.__dict__[detname+'_ix']=det.ix
        self.__dict__[detname+'_iy']=det.iy

        if detnameDict is not None:
            for key in detnameDict.keys():
                if key=='full' or key=='Full':
                    self.__dict__[detname].addFunc(ROIFunc(writeArea=True))
                if key.find('ROI')==0:
                    self.__dict__[detname].addFunc(ROIFunc(name=key,ROI=detnameDict[key],writeArea=True))
                if key.find('photon')==0:
                    thres=0.9
                    if key.find('_')>0:
                        try:
                            thres=float(key.split('_')[1].replace('p','.'))
                        except:
                            pass
                    self.__dict__[detname].addFunc(photon(ADU_per_photon=detnameDict[key], thresADU=thres, retImg=2))

        return detname

    def _getDetName(self, detname=None):
        detNameList = ['cs','Cs','epix','Epix','opal','Opal','zyla','Zyla','jungfrau','Jungfrau','gige','Camera','icarus','ayonix']
        #look for detector
        aliases=[]
        for key in self.Keys():
            if key.alias()!='':
                kname = key.alias()
            else:
                kname = key.src().__repr__()
            for detName in detNameList:
                if kname.find(detName)>=0 and kname.find('Epics')<0:
                    aliases.append(kname)
        if len(aliases)==1:
            return aliases[0]
        elif detname is not None:
            alias_match = [alias for alias in aliases if alias.find(detname)>=0]
            return alias_match
        else:
            return aliases

    def getImages(self, detname='None', numEvts=-1, useFilter=None, nSkip=0, common_mode=0):
        if not isinstance(detname, basestring):
            print('please give parameter name unless specifying arguments in right order. detname is first')
            return
        #look for detector
        if detname=='None':
            detname = self._getDetName()
            if isinstance(detname, list):
                print('detectors in event: \n')
                for alias in detname:
                    print(alias)
                detname = raw_input("Select detector to get detector info for?:\n")

        if not detname in self.__dict__.keys() or self.__dict__[detname].common_mode!=common_mode:
            self.addDetInfo(detname=detname, common_mode=common_mode)
        det=self.__dict__[detname]
        rms = self.__dict__[detname+'_rms']
        pedestals = self.__dict__[detname+'_pedestals']
        ix = self.__dict__[detname+'_ix']
        iy = self.__dict__[detname+'_iy']

        if detname.find('opal')>=0:
            common_mode = -1
            if pedestals[0,0]==-1:
                print('time tool opal, image was not saved. Should ideally exclude from detector list')
                return
        else:
            print('done setting up the geometry')

        #now get the non-image data
        imgAr = []
        run = self.dsIdxRun
        times=[]
        if self.sda is None or 'fh5' not in self.sda.__dict__.keys():
            useFilter = None
        if useFilter is not None:
            evttsSel = self.sda.getSelIdx(useFilter)
            print('using smd based selection, have %s events for selection %s'%(len(evttsSel),useFilter))
            if numEvts==-1:
                numEvts = len(evttsSel)-nSkip
                if numEvts<0:
                    print('have no events, quit')
                    return
            for evtts in evttsSel[nSkip:min(nSkip+numEvts, len(evttsSel))]:
                times.append(psana.EventTime(int(evtts[1]),evtts[0]))
        else:
            times = run.times()[nSkip:]
        print('requested ',numEvts,' used ',min(len(times),numEvts), ' now actually get events')
        if (min(len(times),numEvts) < numEvts*0.5):
            if raw_input('too few events, quit?') in ['y','Y','yes','Yes']:
                return
            else:
                numEvts = len(times)

        for tm in times:
            #print 'numEvts ',numEvts
            if numEvts<=0:
                break
            try:
                evt=run.event(tm)
            except:
                print('Could not get this event, skip ')
                continue
            if evt is None:
                print('Returned event is None, skip')
                continue
            aliases = [ k.alias() for k in evt.keys() ]
            if not detname in aliases:
                continue

            det.getData(evt)
            imgAr.append(det.evt.dat.copy())
            numEvts-=1

        return np.array(imgAr)

    def _getEventTimestamps(self, numEvts=100, useFilter=None, nSkip=0, uniform=False, printFid=False):
        #now get the non-image data
        run = self.dsIdxRun
        times=[]
        if self.sda is None or 'fh5' not in self.sda.__dict__.keys():
            useFilter = None
        if useFilter is not None:
            evttsSel = self.sda.getSelIdx(useFilter)
            print('using smd based selection, have %s events for selection %s'%(len(evttsSel),useFilter))
            if numEvts==-1:
                numEvts = len(evttsSel)-nSkip
                if numEvts<0:
                    print('have no events, quit')
                    return
            for evtts in evttsSel[nSkip:min(nSkip+numEvts, len(evttsSel))]:
                times.append(psana.EventTime(int(evtts[1]),evtts[0]))
        else:
            times = run.times()[nSkip:]
        print('requested %d, used %d now actually get events'%(numEvts,min(len(times),numEvts)))
        if (min(len(times),numEvts) < numEvts*0.5):
            if raw_input('too few events, quit?') in ['y','Y','yes','Yes']:
                return
            else:
                numEvts = len(times)

        if len(times)>numEvts: 
            if uniform:
                tr = np.arange(len(times))
                evtFraction = float(numEvts)/len(times)
                times_sel = np.array(times)[(tr*evtFraction)%1<(evtFraction)]
                times = times_sel
            else:
                times = times[:numEvts]

        if printFid:
            for tm in times:
                print(tm)
                try:
                    evt=run.event(tm)
                except:
                    print('Could not get this event, skip ')
                    continue
                if evt is None:
                    print('Returned event is None, skip')
                    continue
                print((evt.get(psana.EventId)).fiducials())
        return times
    
    
    def pixelHistogram(self, detname=None, numEvts=100, nBins=180, useFilter=None, nSkip=0,minIpm=-1., common_mode=None, std=False, median=False, printFid=False,useMask=True, uniform=False, returnEnv=False):
        """ Return a pixel histogram
        
        Arguments:
            detname: detector name
            numEvts: number of events to be used
            nSkip: number of events to be skipped
            common_mode: calibration applied, including optional correction for common mode noise. Default choice will be used if none are passed if not supplied
            useFilter: use a filter defined using anaps
            uniform: pick events respecting nSkip&useFilter across whole run (default: False, pick first numEvts events)
            printFid: print the picked fiducials (debug)
        """
        
        if detname=='None':
            detname = self._getDetName()
            if isinstance(detname, list):
                detsString='detectors in event: \n'
                for alias in detname:
                    detsString+=alias+', '
                print(detsString)
                detname = raw_input("Select detector to get detector info for?:\n")
                
        if not detname in self.__dict__.keys() or self.__dict__[detname].common_mode!=common_mode:
            self.addDetInfo(detname=detname, common_mode=common_mode)
        det=self.__dict__[detname]
        
        imgAr = []
        run = self.dsIdxRun
        times = self._getEventTimestamps(numEvts=numEvts, useFilter=useFilter, nSkip=nSkip, uniform=uniform, printFid=printFid)
        envDict={'timestamp':[]}
        
        data = []
        for tm in times[:10]:
            if numEvts<=0:
                break
            try:
                evt=run.event(tm)
            except:
                print('Could not get this event, skip ')
                continue
            if evt is None:
                print('Returned event is None, skip')
                continue

            det.getData(evt)
            data.append(det.evt.dat.copy())
        data = np.asarray(data).ravel()
        
        bmin, bmax = np.percentile(data,0.1), 2.5*np.percentile(data,99.5)
        bin_centers = np.linspace(bmin,bmax,nBins)
        hist, bin_edges = np.histogram(data, bins=bin_centers)
        
        for tm in times[10:]:
            #print('numEvts ',numEvts)
            #this loop here is used for the minIpm option when no smallData file is available.
            if numEvts<=0:
                break
            try:
                evt=run.event(tm)
            except:
                print('Could not get this event, skip ')
                continue
            if evt is None:
                print('Returned event is None, skip')
                continue

            det.getData(evt)
            data = det.evt.dat.copy()
            hist+=np.histogram(data, bins=bin_centers)[0]
        return hist, bin_centers, bin_edges

    
    def AvImage(self, detname='None', numEvts=100, thresADU=0., thresRms=0., useFilter=None, nSkip=0,minIpm=-1., common_mode=None, std=False, median=False, printFid=False,useMask=True, uniform=False, returnEnv=False):
        """
        make an average (summed) image for a given detector
        if a detector name is not passed, you will be presented with a choice.
        
        Arguments:
            numEvts: number of events to be used
            nSkip: number of events to be skipped
            common_mode: calibration applied, including optional correction for common mode noise. Default choice will be used if none are passed if not supplied
            thresADU: threshold  each images at thresADU (in units after calibration)
            thresRms: threshhold in multtiples of pixels noise. Currently does not work for gain switching detectors.
            median: get the median, rather than mean image
            std: get the standard deviation as image, rather than mean image
            useMask: apply the user mask to the average image.
            useFilter: use a filter definied using anaps
            uniform: pick events respecting nSkip&useFilter across whole run (default: False, pick first numEvts events)
            printFid: print the picked fiducials (debug)
            minIpm: require a minimal IPM value (if hdf5 file needed to use a filter does not exist). Currently IPM setup only (not wave8)!
            returnEnv: return dictionary of timestamps & enviroment variables if detector provides them
        """
        if not isinstance(detname, basestring):
            print('please give parameter name unless specifying arguments in right order. detname is first')
            return
        #look for detector
        if detname=='None':
            detname = self._getDetName()
            if isinstance(detname, list):
                detsString='detectors in event: \n'
                for alias in detname:
                    detsString+=alias+', '
                print(detsString)
                detname = raw_input("Select detector to get detector info for?:\n")

        if not detname in self.__dict__.keys() or self.__dict__[detname].common_mode!=common_mode:
            self.addDetInfo(detname=detname, common_mode=common_mode)
        det=self.__dict__[detname]
        rms = self.__dict__[detname+'_rms']
        pedestals = self.__dict__[detname+'_pedestals']
        ix = self.__dict__[detname+'_ix']
        iy = self.__dict__[detname+'_iy']

        if detname.find('opal')>=0:
            common_mode = -1
            if pedestals[0,0]==-1:
                print('time tool opal, image was not saved. Should ideally exclude from detector list')
                return
        else:
            print('done setting up the geometry')

        #now get the non-image data
        imgAr = []
        run = self.dsIdxRun
        times = self._getEventTimestamps(numEvts=numEvts, useFilter=useFilter, nSkip=nSkip, uniform=uniform, printFid=printFid)
        envDict={'timestamp':[]}
        for tm in times:
            #print('numEvts ',numEvts)
            #this loop here is used for the minIpm option when no smallData file is available.
            if numEvts<=0:
                break
            try:
                evt=run.event(tm)
            except:
                print('Could not get this event, skip ')
                continue
            if evt is None:
                print('Returned event is None, skip')
                continue
            if minIpm!=-1 and ( (self.hutch=='xpp' and evt.get(psana.Lusi.IpmFexV1, psana.Source('BldInfo(XppSb2_Ipm)')).sum() < minIpm) or (self.hutch=='xcs' and evt.get(psana.Lusi.IpmFexV1, psana.Source('BldInfo(XCS-IPM-05)')).sum() < minIpm)):
                continue
            #aliases = [ k.alias() for k in evt.keys() ]
            #if not detname in aliases:
            #    continue

            det.getData(evt)
            data = det.evt.dat.copy()
            if thresADU != 0:
                data[data<abs(thresADU)]=0
                if thresADU < 0:
                    data[data>=abs(thresADU)]=1
            if thresRms != 0:
                data[data<abs(thresRms)*rms]=0
                if thresRms < 0:
                    data[data>=abs(thresRms)*rms]=1
            imgAr.append(data)                      
            numEvts-=1
            ts = evt.get(psana.EventId).time()
            timesec = float(ts[0])+float(ts[1]/1e6)*1e-3
            envDict['timestamp'].append(timesec)
            for dval in det.evt.__dict__:
                if dval.find('env')>=0: 
                    if dval.replace('env_','') not in envDict:
                        envDict[dval.replace('env_','')]=[]
                    if isinstance(det.evt.__dict__[dval], list) or isinstance(det.evt.__dict__[dval], np.ndarray):
                        envDict[dval.replace('env_','')].append(det.evt.__dict__[dval][0])
                    else:
                        envDict[dval.replace('env_','')].append(det.evt.__dict__[dval])

        #make array
        data='AvImg_';
        if useFilter is not None:
            data+='Filter'+useFilter+'_'
        if thresADU!=0:
            data+='thresADU%d_'%int(thresADU)
        if thresRms!=0:
            data+='thresRms%d_'%int(thresRms*10.)

        if common_mode is None:
            common_mode = det.common_mode
        print('use common mode: ',common_mode)
        data+=self.commonModeStr(common_mode)
        data+=detname

        img = np.array(imgAr)
        if thresADU >= 0 and thresRms >=0:
            imgA = img.mean(axis=0)#.squeeze()
        else:
            imgA = img.sum(axis=0)#.squeeze()
        if det.mask is not None and useMask:
            try:
                imgA[det.mask==0]=0
            except:
                np.squeeze(imgA)[np.squeeze(det.mask)==0]=0
        self.__dict__[data]=imgA
        if std:
            imgS = img.std(axis=0)#.squeeze()
            self.__dict__[data.replace('AvImg_','AvImg_std_')]=imgS
        if median:
            imgM = np.median(img,axis=0)#.squeeze()
            self.__dict__[data.replace('AvImg_','AvImg_median_')]=imgM
            
        if returnEnv: return envDict

    def getAvImage(self,detname=None, imgName=None):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('AvImg')!=0:
                continue
            if key.find('_mask_')>=0:
                continue
            if key.find('_azint_')>=0:
                continue
            if imgName is not None and key.find(imgName)<0:
                continue
            if detname is not None and key.find(detname)<0:
                continue
            avImages.append(key)
        if len(avImages)==0:
            print('please create the AvImage first!')
            return
        elif len(avImages)>1:
            print('we have the following options: ',avImages)
            avImage=raw_input('type the name of the AvImage to use:')
        else:
            avImage=avImages[0]
        detname = self._getDetName_from_AvImage(avImage)
        img = self.__dict__[avImage]
        return detname, img, avImage

    def getAzInt(self,detname=None):
        azInts=[]
        for key in self.__dict__.keys():
            if key.find('_mask_')>=0:
                continue
            print('key ',key)
            if key.find('azint')>=0:
                if detname is not None and key.find(detname)>=0:
                    azInts.append(key)
                elif detname is None:
                    azInts.append(key)
        if len(azInts)==0:
            print('please create an azimuthal integral first!')
            return
        elif len(azInts)>1:
            print('we have the following options: ',azInts)
            avInt=raw_input('type the name of the AvImage to use:')
        else:
            avInt=azInts[0]
        detname = self._getDetName_from_AvImage(avInt.replace('_azint_',''))
        values = self.__dict__[avInt]
        #azavName = 'azav'
        #bins = self.__dict__[detname].__dict__[azavName+'_q']
        return detname, values, avInt

    def _getDetName_from_AvImage(self,avimage):
        detname=''
        avimage=avimage.replace('std_','')
        avimage=avimage.replace('AvImg_','')
        for thisCmString in self.commonModeStrings:
            avimage=avimage.replace(thisCmString+'_','')
        dns = avimage.split('_')
        for ddns in dns:
            if ddns.find('thres')<0 and ddns.find('Filter')<0:
                detname+=ddns;detname+='_'
        if detname[-1]=='_':
            detname = detname[:-1]
        return detname

    def plotAvImage(self,detname=None,imgName=None, use_mask=False, ROI=[], limits=[5,99.5], returnIt=False, plotWith=None,debugPlot=-1, inImage={}):
        if inImage!={}:
            img=inImage['image']
            if 'mask' in inImage:
                mask=inImage['mask']
            else:
                mask=np.ones_like(img)
            detname=inImage['detname']
            avImage='input Image'
            if 'name' in inImage:
                avImage=inImage['name']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)
        else:
            detname, img, avImage = self.getAvImage(detname=detname, imgName=imgName)

        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
            img = (img*mask)

        plotMax = np.nanpercentile(img, limits[1])
        plotMin = np.nanpercentile(img, limits[0])
        print('plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        if len(img.shape)>2:
            image = self.__dict__[detname].det.image(self.run, img)
            if image is None:
                if self.__dict__[detname].det.dettype != 30:
                    image = img.squeeze()
                else:
                    imgNew = img[0]
                    for iFrame, frame in enumerate(img):
                        if iFrame>0:
                            imgNew = np.append(imgNew, frame, axis=1)
                    image = imgNew
        else:
            image = img

        #have issues with Jupyter notebook io size.
        if plotWith is None:
            plotWith=self.plotWith
        if plotWith=='bokeh_notebook':
            plot_title="%s in Run %d"%(avImage, self.run)
            if debugPlot>0:
                img = image[:debugPlot, :debugPlot]
            else:
                img = image
            layout, p, im = plotImageBokeh(img, plot_title=plot_title, plotMaxP=np.nanpercentile(img,99),plotMinP=np.nanpercentile(img,1), plotWidth=700, plotHeight=700)
            bp.output_notebook()
            bp.show(layout)
        else:
            fig=plt.figure(figsize=(10,6))
            if ROI!=[]:
                gs=gridspec.GridSpec(1,2,width_ratios=[2,1])        
                im1 = plt.subplot(gs[1]).imshow(img[ROI[0][0],ROI[1][0]:ROI[1][1],ROI[2][0]:ROI[2][1]],clim=[plotMin,plotMax],interpolation='None')
                cbar1 = plt.colorbar(im1)
            else:
                gs=gridspec.GridSpec(1,2,width_ratios=[99,1])        
            im0 = plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')
            cbar1 = plt.colorbar(im0)
            plt.title(avImage)
            plt.show()

        if returnIt:
            return image

    def saveAvImage(self,detname=None,dirname=''):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('_mask_')>=0:
                continue
            if key.find('_azint_')>=0:
                continue
            if detname is not None and key.find(detname)<0:
                continue
            if key.find('AvImg')==0:
                if detname is not None and key.find(detname)>=0:
                    avImages.append(key)
                elif detname is None:
                    avImages.append(key)

        detname, img, avImage = self.getAvImage(detname=None)

        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True
        data_to_save={}
        data_to_save['ped'] = self.__dict__[detname].ped
        data_to_save['rms'] = self.__dict__[detname].rms
        data_to_save['mask'] = self.__dict__[detname].mask
        data_to_save['cmask'] = self.__dict__[detname].cmask
        if self.__dict__[detname].x is not None:
            data_to_save['x'] = self.__dict__[detname].x
            data_to_save['y'] = self.__dict__[detname].y
        for rawImg in avImages:
            data_to_save[rawImg] = self.__dict__[rawImg]
            if needsGeo:
                data_to_save[rawImg+'_image'] = self.__dict__[detname].det.image(self.run,self.__dict__[rawImg])
        if dirname=='':
            dirname=self.sda.dirname
        if len(avImages)<=0:
            return
        fname = '%s%s_Run%03d_%s.h5'%(dirname,self.expname,self.run,avImages[0])
        print('now save information to file: ',fname)
        imgFile = h5py.File(fname, "w")
        for key in data_to_save:
            addToHdf5(imgFile, key, data_to_save[key])
        imgFile.close()
            
#bokeh, use BoxSelectTool to get selected region
#https://stackoverflow.com/questions/34164587/get-selected-data-contained-within-box-select-tool-in-bokeh
    def SelectRegion(self,detname=None, limits=[5,99.5]):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('AvImg')==0:
                if detname is not None and key.find(detname)>=0:
                    avImages.append(key)
                elif detname is None:
                    avImages.append(key)
        if len(avImages)==0:
            print('please create the AvImage first!')
            return
        elif len(avImages)>1:
            print('we have the following options: ',avImages)
            avImage=raw_input('type the name of the AvImage to use:')
        else:
            avImage=avImages[0]
        detname = self._getDetName_from_AvImage(avImage)
        img = self.__dict__[avImage]

        plotMax = np.nanpercentile(img, limits[1])
        plotMin = np.nanpercentile(img, limits[0])
        print('plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        fig=plt.figure(figsize=(10,6))
        gs=gridspec.GridSpec(1,2,width_ratios=[2,1])

        needsGeo=False
        if self.__dict__[detname].det.dettype==26:
            if self.__dict__[detname].ped[0].shape != self.__dict__[detname].imgShape:
                needsGeo=True
        elif self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if needsGeo:
            image = self.__dict__[detname].det.image(self.run, img)
        else:
            image = img.squeeze()

        plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')
        plt.pause(0.0001)

        ix = self.__dict__[detname+'_ix']
        iy = self.__dict__[detname+'_iy']

        happy = False
        while not happy:
            p =np.array(ginput(2))
            p_axis1=[int(p[:,1].min()),int(p[:,1].max())+1]
            p_axis0=[int(p[:,0].min()),int(p[:,0].max())+1]
            print('points:',p)
            plt.subplot(gs[1]).imshow(image[p_axis1[0]:p_axis1[1],p_axis0[0]:p_axis0[1]],clim=[plotMin,plotMax],interpolation='None')
            plt.pause(0.0001)
            if needsGeo:
                mask_roi=np.zeros_like(image)
                mask_roi[p_axis1[0]:p_axis1[1],p_axis0[0]:p_axis0[1]]=1
                mask_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(ix,iy)] )
            if raw_input("Happy with this selection:\n") in ["y","Y"]:
                happy = True
                    
        if not needsGeo:
            print('ROI: [[%i,%i], [%i,%i]]'%(p[:,1].min(),p[:,1].max(),p[:,0].min(),p[:,0].max()))
            return

        if len(mask_nda.shape)>2:
            for itile,tile in enumerate(mask_nda):
                if tile.sum()>0:
                    ax0 = np.arange(0,tile.sum(axis=0).shape[0])[tile.sum(axis=0)>0]
                    ax1 = np.arange(0,tile.sum(axis=1).shape[0])[tile.sum(axis=1)>0]
                    print('ROI: [[%i,%i], [%i,%i], [%i,%i]]'%(itile,itile+1,ax1.min(),ax1.max(),ax0.min(),ax0.max()))
                    fig=plt.figure(figsize=(6,6))
                    plt.imshow(img[itile,ax1.min():ax1.max(),ax0.min():ax0.max()],interpolation='none')
                    plt.pause(0.0001)
                    plt.show()
        else:
            tile=mask_nda
            if tile.sum()>0:
                ax0 = np.arange(0,tile.sum(axis=0).shape[0])[tile.sum(axis=0)>0]
                ax1 = np.arange(0,tile.sum(axis=1).shape[0])[tile.sum(axis=1)>0]
                print('ROI: [[%i,%i], [%i,%i]]'%(ax1.min(),ax1.max(),ax0.min(),ax0.max()))
            

    def FitCircleAuto(self, detname=None, plotRes=True, forceMask=False, inImage={},inParams={}):
        if inImage!={}:
            img=inImage['image']
            if 'mask' in inImage:
                mask=inImage['mask']
            else:
                mask=np.ones_like(img)
            detname=inImage['detname']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)

        else:
            detname, img, avImage = self.getAvImage(detname=None)
            try:
               mask = self.__dict__[detname].det.cmask.astype(bool)
            except:
               mask = self.__dict__[detname].det.mask(self.run, calib=True, status=True).astype(bool)
            nPixRaw=1.
            for idim in mask.shape:
                nPixRaw=nPixRaw*idim
            print('check calib mask: ',mask.sum(),nPixRaw,(mask.sum()/nPixRaw)>0.5,(mask.sum().astype(float)/nPixRaw))
            if (mask.sum()/nPixRaw)<0.5 and not forceMask:
                mask=~mask
            try:
                maskgeo = self.__dict__[detname].det.mask_geo(self.run).astype(bool)
                mask = mask*maskgeo
            except:
                pass
        #apply mask to image.
        img = (img*mask)

        needsGeo=False
        if self.__dict__[detname].ped is not None and self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if needsGeo:
            image = self.__dict__[detname].det.image(self.run, img)
            mask =  self.__dict__[detname].det.image(self.run, mask)
            #limits =  [[self.__dict__[detname+'_x'].min(),self.__dict__[detname+'_x'].max()],
            #           [self.__dict__[detname+'_y'].min(),self.__dict__[detname+'_y'].max()]]
        else:
            image = img
        limits=[[0,image.shape[0]],[0,image.shape[1]]]

        #plot image in grayscale
        combRes, ringInfo, arSparse = FindFitCenter(image, mask, inParams=inParams)

        if not plotRes:
            return combRes
        #plot image in grayscale
        plt.figure(figsize=[12,12])
        plt.imshow(image, interpolation='none', cmap='gray',clim=[0,np.nanpercentile(img.flatten(),99.5)])
        #plot sparse images in blue.
        plt.plot(arSparse.col, arSparse.row,markersize=5,color='#ff9900',marker='.',linestyle='None')
        if combRes==-1:
            return -1
        greens=['#666600','#669900','#66cc00','#66cc99','#6699cc','#6633ff']
        print('center ',combRes['xCen'],combRes['yCen'])
        for ir,thisRingInfo,rFit in itertools.izip(itertools.count(),ringInfo,combRes['R']):
            #plot ransac selected data in green <-- consider different greens for first circles.
            plt.plot(arSparse.col[thisRingInfo['pointsInCircle']], arSparse.row[thisRingInfo['pointsInCircle']],marker='.',color=greens[ir],linestyle='None', markersize=4)
            circle1 = plt.Circle((combRes['xCen'],combRes['yCen']),rFit,color=greens[ir],fill=False,linestyle='dashed',linewidth=1)
            plt.gca().add_artist(circle1)
        #circle1 = plt.Circle((combRes['xCen'],combRes['yCen']),ringInfo[0]['rInput'],color='r',fill=False,linestyle='dashed')
        #plt.gca().add_artist(circle1)
        plt.plot([combRes['xCen'],combRes['xCen']],[combRes['yCen']-20,combRes['yCen']+20],color='m')
        plt.plot([combRes['xCen']-20,combRes['xCen']+20],[combRes['yCen'],combRes['yCen']],color='m')
        plt.xlim(limits[0])
        plt.ylim(limits[1])
        plt.show()
        if needsGeo:
            geo = self.__dict__[detname].det.geometry(self.run)
            d=160000.
            ix0,iy0 = geo.point_coord_indexes(p_um=(0,0))
            ixd,iyd = geo.point_coord_indexes(p_um=(d,d))
            combRes['yCen'] = d*(combRes['yCen']-ix0)/(ixd-ix0)
            combRes['xCen'] = d*(combRes['xCen']-iy0)/(iyd-iy0)
            print('center Final mid',combRes['xCen'],combRes['yCen'])
            helpVar =  combRes['xCen']
            combRes['xCen'] = combRes['yCen']
            combRes['yCen'] = helpVar
            for r in combRes['R']:
                print('aradii: ',r,d*r/(ixd-ix0))
        else:
            for r in combRes['R']:
                print('aradii: ',r)
        print('center Final ',combRes['xCen'],combRes['yCen'])
        return combRes
 
    def FitCircleMouse(self, detname=None, use_mask=False, limits=[5,99.5], use_mask_local=False):
        self.FitCircle(detname=detname, use_mouse=True, use_mask=use_mask, use_mask_local=use_mask_local, limits=[5,99.5])

    def FitCircleThreshold(self, detname=None, use_mask=False, use_mask_local=False, limits=[5,99.5],plotIt=False, thresIn=None):
        self.FitCircle(detname=detname, use_mouse=False, use_mask=use_mask, use_mask_local=use_mask_local, limits=[5,99.5], thresIn=thresIn, plotIt=plotIt)

    def FitCircle(self, **kwargs):
        """
        function to fit (single) circles in images
        options are to fit point selected by mouse or by threshold
        
        Parameters
        ----------
        detname: name of detector we have created an average image for
        use_mouse: select points to be fitted by mouse (if false, use threshold)
        use_mask: use the mask of detector as saved in calibration data, default is False
        use_mask_local: use a locally defined mask, default is False
        limits: range of z-axisof plot in percent of image data, default is [5, 99.5]
        plotIt: show plots of results & for selection of parameters. Default is True. Turn off for running this in batch
        thresIn: supply a minimum threshold for points selected for fit in percentile of image values
        singleTile: select a single tile of image (e.g. useful for Icarus detector), default -1 (use whole detector)
        """
        detname = kwargs.pop("detname",None)
        use_mouse = kwargs.pop("use_mouse",None)
        use_mask = kwargs.pop("use_mask",False)
        use_mask_local = kwargs.pop("use_mask_local",False)
        plotIt = kwargs.pop("plotIt",True)
        thresIn = kwargs.pop("thresIn",None)
        limits = kwargs.pop("limits",[5, 99.5])
        singleTile = kwargs.pop("singleTile",-1)

        try: 
            detname, img, avImage = self.getAvImage(detname=detname)
        except:
            detname, img, avImage = self.getAvImage(detname=None)

        if use_mouse is None or (not isinstance(use_mouse, bool) and not isinstance(use_mouse, basestring)):
            #set pplotting to yes
            plotIt=True
            if raw_input("Select Circle Points by Mouse?:\n") in ["y","Y"]:
                use_mouse = True
            else:
                use_mouse = False

        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
            img = (img*mask)

        if use_mask_local and self.__dict__.has_key('_mask_'+avImage):
                mask = self.__dict__['_mask_'+avImage]
                img = (img*mask)
        elif use_mask_local and not use_mouse:
            maskName = raw_input('no local mask defined for %s, please enter file name '%avImage)
                
            if maskName!='' and os.path.isfile(maskName):
                locmask=np.loadtxt(maskName)
                if self.__dict__[detname].x is not None and locmask.shape != self.__dict__[detname].x.shape:
                    if locmask.shape[1] == 2:
                        locmask = locmask.transpose(1,0)
                    locmask = locmask.reshape(self.__dict__[detname].x.shape)
                self.__dict__['_mask_'+avImage] = locmask
                img = (img*locmask)

        plotMax = np.nanpercentile(img[img!=0], 99.5)
        plotMin = np.nanpercentile(img[img!=0], 5)
        printR(rank, 'plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        needsGeo=False
        if self.__dict__[detname].ped is not None and self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        x = self.__dict__[detname+'_x']
        y = self.__dict__[detname+'_y']
        if x is None:
            x = np.arange(0, image.shape[1])
            y = np.arange(0, image.shape[0])
            x,y = np.meshgrid(x,y)
        extent=[x.min(), x.max(), y.min(), y.max()]

        if needsGeo:
            if self.__dict__[detname].det.dettype==30:
                if singleTile<0:
                    print('we need to select a tile for the icarus')
                    return
                elif singleTile>=img.shape[0]:
                    print('requested tile %d, detector %s only has %d tiles'%(singleTile, detname, image.shape[0]))
                    return
                    
                image = img[singleTile]
                x = x[singleTile]
                y = y[singleTile]
                extent=[x.min(), x.max(), y.min(), y.max()]
                print('icarus image: ',img.shape, image.shape, extent)
            else:
                image = self.__dict__[detname].det.image(self.run, img)
        else:
            image = img

        happy = False
        if use_mouse:
            while not happy:
                fig=plt.figure(figsize=(10,10))
                if needsGeo:
                    plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
                else:
                    plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None')
                points=ginput(n=0)
                parr=np.array(points)
                if parr.shape[0]==1:
                    print('[x,y]: [%f, %f]'%(parr[0][0],parr[0][1]))
                    happy=True
                    break
                res = fitCircle(parr[:,0],parr[:,1])
                #draw the circle. now need to redraw the whole image thanks to the conda matplotlib
                fig=plt.figure(figsize=(10,10))
                if needsGeo:
                    plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
                else:
                    plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None')
                circle = plt.Circle((res['xCen'],res['yCen']),res['R'],color='b',fill=False)
                plt.gca().add_artist(circle)
                plt.plot([res['xCen'],res['xCen']],[y.min(),y.max()],'r')
                plt.plot([x.min(),x.max()],[res['yCen'],res['yCen']],'r')
                print('[x,y] (micron): [%f, %f] R (in mm): %f '%(res['xCen'],res['yCen'],res['R']/1000.))
                if raw_input("Happy with this selection:\n") in ["y","Y"]:
                    happy = True

        else:
            while not happy:
                if thresIn is not None:
                    thres = thresIn
                    thresP = np.nanpercentile(img[img!=0], thres)
                    happy = True
                else:                    
                    if not plotIt:
                        print('this is not going to work, you either need to specify a threshold or require plots')
                        return
                    fig=plt.figure(figsize=(10,10))
                    if needsGeo:
                        plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
                    else:
                        plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None')
                    thres = float(raw_input("min percentile % of selected points:\n"))

                    thresP = np.nanpercentile(img[img!=0], thres)
                    print('thresP',thresP)
                    imageThres=image.copy()
                    imageThres[image>thresP]=1
                    imageThres[image<thresP]=0

                    fig=plt.figure(figsize=(5,5))
                    plt.imshow(imageThres,clim=[-0.1,1.1])
                    if raw_input("Happy with this threshold (y/n):\n") in ["y","Y"]:
                        happy=True

            if singleTile<0:
                res = fitCircle(x.flatten()[img.flatten()>thresP],y.flatten()[img.flatten()>thresP])
            else:
                if len(x.shape)==len(img.shape) and len(x.shape)>2:
                    res = fitCircle(x[singleTile].flatten()[img[singleTile].flatten()>thresP],y[singleTile].flatten()[img[singleTile].flatten()>thresP])
                else:
                    res = fitCircle(x.flatten()[img[singleTile].flatten()>thresP],y.flatten()[img[singleTile].flatten()>thresP])
            circleM = plt.Circle((res['xCen'],res['yCen']),res['R'],color='b',fill=False)
            print('[x,y] (micron): [%f, %f] R (in mm): %f '%(res['xCen'],res['yCen'],res['R']/1000.))
            if plotIt:
                fig=plt.figure(figsize=(10,10))
                if needsGeo:
                    plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
                else:
                    plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None')
                plt.gca().add_artist(circleM)
                plt.plot([res['xCen'],res['xCen']],[y.min(),y.max()],'r')
                plt.plot([x.min(),x.max()],[res['yCen'],res['yCen']],'r')
                plt.show()

    def MakeMask(self, detname=None, limits=[5,99.5], singleTile=-1):
        detname, img, avImage = self.getAvImage(detname=None)

        plotMax = np.nanpercentile(img, limits[1])
        plotMin = np.nanpercentile(img, limits[0])
        print('plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if self.__dict__[detname].det.dettype==30:
            needsGeo=False
            if singleTile<0 or singleTile>= img.shape[0]:
                image = img.sum(axis=0)
            else:
                image = img[singleTile]
        elif needsGeo:
            image = self.__dict__[detname].det.image(self.run, img)
            if image is None and self.__dict__[detname].ped.shape[1]==1:
                image = img.squeeze()
        else:
            image = img

        det = self.__dict__[detname].det
        x = self.__dict__[detname+'_x']
        y = self.__dict__[detname+'_y']
        if x is None:
            xVec = np.arange(0, image.shape[1])
            yVec = np.arange(0, image.shape[0])
            x, y = np.meshgrid(xVec, yVec)
            self.__dict__[detname+'_x'] = x
            self.__dict__[detname+'_y'] = y
        ix = self.__dict__[detname+'_ix']
        iy = self.__dict__[detname+'_iy']
        extent=[x.min(), x.max(), y.min(), y.max()]
        #print('DEBUG: extent(x,y min,max)',extent)
        if self.__dict__[detname].det.dettype==30:
            if singleTile<0 or singleTile>= img.shape[0]:
                x=x[0]
                y=y[0]
            else:
                x=x[singleTile]
                y=y[singleTile]
        
        mask=[]
        mask_r_nda=None
        select=True
        while select:
            fig=plt.figure(figsize=(12,10))
            gs=gridspec.GridSpec(1,2,width_ratios=[2,1])
            #print("rectangle(r-click, R-enter), circle(c), polygon(p), dark(d), noise(n) or edgepixels(e)?:")
            #needs to be pixel coordinates for rectable selection to work.
            plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')

            #edgepixels(e)?:\n")
            #polygon(p), 
            #rectangle(r-click,
            #rectangle(R-enter), 

            #circle(c), 
            #dark(d), 
            #noise(n) or 

            shape = raw_input("rectangle(r-click, R-enter), circle(c), polygon(p), dark(d), noise(n) or edgepixels(e)?:\n")
            #shape = raw_input()
            #this definitely works for the rayonix...
            if shape=='r':
                print('select two corners: ')
                p =np.array(ginput(2))
                mask_roi=np.zeros_like(image)
                p=p.astype(int)
                mask_roi[p[:,1].min():p[:,1].max(),p[:,0].min():p[:,0].max()]=1
                if needsGeo:
                    mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(ix,iy)] )
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    mask_r_nda = mask_roi
                    plt.subplot(gs[1]).imshow(mask_r_nda)
                    if self.__dict__[detname].det.dettype==30:
                        maskTuple=[mask_r_nda for itile in range(self.__dict__[detname].ped.shape[0])]
                        mask_r_nda = np.array(maskTuple)
                print('mask from rectangle (shape):',mask_r_nda.shape)
            elif shape=='R':
                print('coordinates to select: ')
                htot = raw_input("horizontal,vertical h1 h2 v1 v2 ?\n")
                h = htot.split(' ');h1=float(h[0]);h2=float(h[1]);v1=float(h[2]);v2=float(h[3]);
                mask_roi=np.zeros_like(image)
                mask_roi[int(h1):int(h2),int(v1):int(v2)]=1
                if needsGeo:
                    mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(ix,iy)] )
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    mask_r_nda = mask_roi
                    plt.subplot(gs[1]).imshow(mask_r_nda)
                print('mask from rectangle (shape):',mask_r_nda.shape)
            elif shape=='c':
                fig=plt.figure(figsize=(12,10))
                gs=gridspec.GridSpec(1,2,width_ratios=[2,1])
                if det.dettype==30:
                    plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None',extent=(y.min(),y.max(),x.min(),x.max()))
                else:
                    plt.subplot(gs[0]).imshow(np.rot90(image),clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
                if raw_input("Select center by mouse?\n") in ["y","Y"]:
                    c=ginput(1)
                    cx=c[0][0];cy=c[0][1]
                    print('Corrdinates of selected center: ',cx,' ',cy)
                else:
                    ctot = raw_input("center (x y)?\n")
                    c = ctot.split(' ');cx=float(c[0]);cy=float(c[1]);
                if raw_input("Select outer radius by mouse?\n") in ["y","Y"]: 
                    r=ginput(1)
                    rox=r[0][0];roy=r[0][1]
                    print('outer radius point: ',r[0])
                    ro=np.sqrt((rox-cx)**2+(roy-cy)**2)
                    if raw_input("Select inner radius by mouse (for donut-shaped mask)?\n") in ["y","Y"]:
                        r=ginput(1)
                        rix=r[0][0];riy=r[0][1]
                        print('inner radius point: ',r[0])
                        ri=np.sqrt((rix-cx)**2+(riy-cy)**2)
                    else:
                        ri=0
                    print('radii: ',ro,' ',ri)
                else:
                    rtot = raw_input("radii (r_outer r_inner)?\n")
                    r = rtot.split(' ');ro=float(r[0]);ri=max(0.,float(r[1]));        
                mask_router_nda = np.array( [(ix-cx)**2+(iy-cy)**2<ro**2 for ix, iy in zip(x,y)] )
                mask_rinner_nda = np.array( [(ix-cx)**2+(iy-cy)**2<ri**2 for ix, iy in zip(x,y)] )
                mask_r_nda = mask_router_nda&~mask_rinner_nda
                print('mask from circle (shape):',mask_r_nda.shape)
                if raw_input("Invert [y/n]? (n/no inversion: masked pixels will get rejected)?\n") in ["y","Y"]:
                    mask_r_nda = ~mask_r_nda

                if needsGeo:
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    plt.subplot(gs[1]).imshow(mask_r_nda)
            elif shape=='p':
                fig=plt.figure(figsize=(12,10))
                gs=gridspec.GridSpec(1,2,width_ratios=[2,1])
                if needsGeo:
                    plt.subplot(gs[0]).imshow(np.rot90(image),clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
                elif det.dettype==13:
                    plt.subplot(gs[0]).imshow(np.rot90(image),clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
                elif det.dettype==19:
                    #plt.subplot(gs[0]).imshow(np.rot90(image),clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
                    plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
                elif det.dettype==30:
                    plt.subplot(gs[0]).imshow(np.rot90(image),clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
                else:
                    plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
                nPoints = int(raw_input("Number of Points (-1 until middle mouse click)?\n"))
                p=np.array(ginput(nPoints))
                print(p)
                mpath=path.Path(p)
                all_p = np.array([ (ix,iy) for ix,iy in zip(x.flatten(),y.flatten()) ] )
                mask_r_nda = np.array([mpath.contains_points(all_p)]).reshape(x.shape)
                if needsGeo:
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    plt.subplot(gs[1]).imshow(mask_r_nda)
                if self.__dict__[detname].det.dettype==30:
                    maskTuple=[mask_r_nda for itile in range(self.__dict__[detname].ped.shape[0])]
                    mask_r_nda = np.array(maskTuple)
                print('mask from polygon (shape):',mask_r_nda.shape)
            elif shape=='d' or shape=='n':
                figDark=plt.figure(figsize=(12,10))
                gsPed=gridspec.GridSpec(1,2,width_ratios=[1,1])
                if shape=='d':
                    pedResult = det.pedestals(self.run)
                    if det.dettype in [26,29,32,33]:
                        pedResult = pedResult[0]
                    hstPed = np.histogram(pedResult.flatten(), np.arange(0, pedResult.max()*1.05))
                else:
                    pedResult = det.rms(self.run)
                    if det.dettype in [26,29,32,33]:
                        pedResult = pedResult[0]
                    hstPed = np.histogram(pedResult.flatten(), np.arange(0, pedResult.max()*1.05,0.05))
                if needsGeo:
                    pedResultImg = det.image(self.run,pedResult)
                else:
                    pedResultImg = pedResult.copy()
                plt.subplot(gsPed[0]).imshow(pedResultImg,clim=[np.nanpercentile(pedResult,1),np.nanpercentile(pedResult,99)])
                plt.subplot(gsPed[1]).plot(hstPed[1][:-1],np.log(hstPed[0]),'o')
                ctot=raw_input("Enter allowed pedestal range (min max)")
                c = ctot.split(' ');pedMin=float(c[0]);pedMax=float(c[1]);
                mask_r_nda=np.zeros_like(pedResult)
                mask_r_nda[pedResult<pedMin]=1
                mask_r_nda[pedResult>pedMax]=1
                mask_r_nda = (mask_r_nda.astype(bool)).astype(int)
                if needsGeo:
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    plt.subplot(gs[1]).imshow(mask_r_nda)
            if shape=='e':
                ctot=raw_input("Enter number of edge rows that should be masked:")
                try:
                    nEdge = int(ctot)
                except:
                    print('please enter an integer')
                    continue
                mask_r_nda = np.zeros_like(det.coords_x(self.run))
                if needsGeo:
                    for tile in mask_r_nda:
                        tile[0:nEdge,:]=1
                        tile[tile.shape[0]-nEdge:,:]=1
                        tile[:,0:nEdge]=1
                        tile[:,tile.shape[1]-nEdge:]=1
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda.astype(bool)))
                else:
                    tile=mask_r_nda
                    tile[0:nEdge,:]=1
                    tile[nEdge:-1,:]=1
                    tile[:,0:nEdge]=1
                    tile[:,nEdge:-1]=1
                    plt.subplot(gs[1]).imshow(mask_r_nda.astype(bool))

            if shape=='cen' or shape=='center':
                ctot=raw_input("Enter number of center rows that should be masked:")
                try:
                    nEdge = int(ctot)
                except:
                    print('please enter an integer')
                    continue
                mask_r_nda = np.zeros_like(det.coords_x(self.run))
                if needsGeo:
                    for tile in mask_r_nda:
                        tile[int(tile.shape[0]/2-nEdge):int(tile.shape[0]/2+nEdge),:]=1
                        tile[:,int(tile.shape[1]/2-nEdge):int(tile.shape[1]/2+nEdge)]=1
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda.astype(bool)))
                else:
                    tile=mask_r_nda
                    tile[int(tile.shape[0]/2-nEdge):int(tile.shape[0]/2+nEdge),:]=1
                    tile[:,int(tile.shape[1]/2-nEdge):int(tile.shape[1]/2+nEdge)]=1
                    plt.subplot(gs[1]).imshow(mask_r_nda.astype(bool))

            if mask_r_nda is not None:
                print('created a mask....',len(mask))
                mask.append(mask_r_nda.astype(bool).copy())
                countMask=1
                totmask_nm1 = mask[0]
                for thismask in mask[1:-1]:
                    totmask_nm1 = np.logical_or(totmask_nm1,thismask)
                    countMask+=1
                if len(mask)>1:
                    totmask = np.logical_or(totmask_nm1,mask[-1])
                    countMask+=1
                else:
                    totmask = totmask_nm1
            #print('DEBUG: ',mask_r_nda.shape, totmask_nm1.shape, x.shape)
            print('masked in this step: ',np.ones_like(self.__dict__[detname].x)[mask_r_nda.astype(bool)].sum())
            print('masked up to this step: ',np.ones_like(self.__dict__[detname].x)[totmask_nm1].sum())
            print('masked tot: ',np.ones_like(self.__dict__[detname].x)[totmask].sum())

            if len(mask)>1:
                fig=plt.figure(figsize=(15,9))
                gs2=gridspec.GridSpec(1,2,width_ratios=[1,1])
                plt.show()
                image_mask = img.copy(); image_mask[totmask]=0;
                #image_mask_nm1 = img.copy(); image_mask_nm1[totmask_nm1]=0;
                #print('DEBUG: ',(img-image_mask_nm1).sum(), (img-image_mask).sum())
                plt.subplot(gs2[0]).imshow(image,clim=[plotMin,plotMax])
                if needsGeo:
                    plt.subplot(gs2[1]).imshow(det.image(self.run,image_mask),clim=[plotMin,plotMax])
                else:
                    plt.subplot(gs2[1]).imshow(image_mask,clim=[plotMin,plotMax])
            else:
                fig=plt.figure(figsize=(12,10))
                plt.show()
                image_mask = img.copy(); image_mask[totmask]=0;
                if needsGeo:
                    plt.imshow(det.image(self.run,image_mask),clim=[plotMin,plotMax])
                else:
                    if (self.__dict__[detname].det.dettype==30):
                        if singleTile<0 or singleTile>= img.shape[0]:
                            image = image_mask.sum(axis=0)
                        else:
                            image = [singleTile]
                        plt.imshow(image,clim=[plotMin,plotMax])
                    else:
                        plt.imshow(image_mask,clim=[plotMin,plotMax])
                

            if raw_input("Add this mask?\n") in ["n","N"]:
                mask = mask[:-1]

            if len(mask)>0:
                #remake the image with mask up to here.
                totmask = mask[0]
                for thismask in mask[1:]:
                    totmask = np.logical_or(totmask,thismask)
                image_mask = img.copy(); image_mask[totmask]=0;
                if needsGeo:
                    image = self.__dict__[detname].det.image(self.run, image_mask)
                else:
                    image = image_mask

            if raw_input("Done?\n") in ["y","Y"]:
                select = False

        #end of mask creating loop
        if len(mask)==0:
            return
        totmask = mask[0]
        for thismask in mask[1:]:
            totmask = np.logical_or(totmask,thismask)

        if raw_input("Invert [y/n]? (n/no inversion: masked pixels will get rejected)?\n") in ["y","Y"]:
            totmask = (totmask.astype(bool)).astype(int)
        else:
            totmask = (~(totmask.astype(bool))).astype(int)

        self.__dict__['_mask_'+avImage]=totmask

        if  det.dettype == 2:
            mask=totmask.reshape(2,185*388).transpose(1,0)
        elif det.dettype == 1:
            mask=totmask.reshape(32*185,388)
        elif det.dettype == 13:
            mask=totmask.reshape(704,768)
        else:
            mask=totmask
        #2x2 save as 71780 lines, 2 entries
        #cspad save as 5920 lines, 388 entries
        if raw_input("Save to calibdir?\n") in ["y","Y"]:
            srcStr=det.source.__str__().replace('Source("DetInfo(','').replace(')")','')
            if det.dettype==2:
                dirname='/reg/d/psdm/%s/%s/calib/CsPad2x2::CalibV1/%s/pixel_mask/'%(self.expname[:3],self.expname,srcStr)
            elif det.dettype==1:
                dirname='/reg/d/psdm/%s/%s/calib/CsPad::CalibV1/%s/pixel_mask/'%(self.expname[:3],self.expname,srcStr)        
            elif det.dettype==13:
                dirname='/reg/d/psdm/%s/%s/calib/Epix100a::CalibV1/%s/pixel_mask/'%(self.expname[:3],self.expname,srcStr)
            elif det.dettype==19:
                dirname='/reg/d/psdm/%s/%s/calib/Camera::CalibV1/%s/pixel_mask/'%(self.expname[:3],self.expname,srcStr)        
            elif det.dettype==32:
                dirname='/reg/d/psdm/%s/%s/calib/Epix10ka2M::CalibV1/%s/pixel_mask/'%(self.expname[:3],self.expname,srcStr)        
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            fname='%s-end.data'%self.run
            print('save mask in %s as %s '%(dirname,fname))
            #np.savetxt(dirname+fname,mask)
            det.save_txtnda(dirname+fname,mask.astype(float), fmt='%d',addmetad=True)
        elif raw_input("Save to local?\n") in ["y","Y"]:
            if len(mask.shape)<=2:
                np.savetxt('Mask_%s_%s_Run%03d.data'%(avImage,self.expname,int(self.run)),mask)
            elif len(mask.shape)==3:
                np.savetxt('Mask_%s_%s_Run%03d.data'%(avImage,self.expname,int(self.run)),mask.reshape(mask.shape[0]*mask.shape[1], mask.shape[2]))
                
        return mask
         
    def makePedestal(self, detname, useFilter=None, numEvts=-1, pedRange=[10,10000], rmsRange=[2.,7.], i0Check='ipm', dirname='./'):
        if i0Check=='' or i0Check is None:
            i0List=[]
        elif i0Check=='ipm':
            if self.expname[:3]=='xpp':
                i0List = ['ipm3/sum']
            elif self.expname[:3]=='xcs':
                i0List = ['ipm5/sum']
            else:
                i0List=[]
        else:
            i0List = ['gas_detector/f_22_ENRC']
    
        minFrac = 1e6
        if useFilter is not None:
            for i0 in i0List:
                try:
                    i0Median = np.nanmedian(self.sda.getVar(i0))
                    i0Off = np.nanmedian(self.sda.getVar(i0, useFilter=useFilter))
                except:
                    print('if not smallData file is available, try pass i0Check=\'\'')
                if minFrac > i0Off/i0Median:
                    minFrac=i0Off/i0Median
                print('median value for ',i0,' is ',i0Median,' and for the off events ',i0Off,' offRatio: ',i0Off/i0Median)
            if minFrac > 0.05 and i0List!=[]:
                print('This selection seems to lets too many events with beam through, will quit')
                return

        self.AvImage(detname,numEvts=numEvts,useFilter=useFilter, common_mode=-1, useMask=False, median=True)
        self.AvImage(detname,numEvts=numEvts,useFilter=useFilter, common_mode=-1, useMask=False, std=True)
        fname='%s-end.data'%self.sda.run
        if useFilter is not None:
            pedImg = self.__dict__['AvImg_median_Filter%s_raw_%s'%(useFilter,detname)]
            rmsImg = self.__dict__['AvImg_std_Filter%s_raw_%s'%(useFilter,detname)]
        else:
            pedImg = self.__dict__['AvImg_median_raw_%s'%(detname)]
            rmsImg = self.__dict__['AvImg_std_raw_%s'%(detname)]
        pedStat = np.logical_and(pedImg > min(pedRange), pedImg < max(pedRange))
        rmsStat = np.logical_and(rmsImg > min(rmsRange), rmsImg < max(rmsRange))
        status = (~(np.logical_and(rmsStat, pedStat))).astype(int)

        det = self.__dict__[detname].det     
        #if raw_input("Save to calibdir?\n") in ["y","Y"]:
        if dirname == 'calib':
            #detname, img, avImage = self.getAvImage(detname=None)
            srcStr=det.source.__str__().replace('Source("DetInfo(','').replace(')")','')
            if det.dettype==2:
                dirname='/reg/d/psdm/%s/%s/calib/CsPad2x2::CalibV1/%s/'%(self.expname[:3],self.expname,srcStr)
            elif det.dettype==1:
                dirname='/reg/d/psdm/%s/%s/calib/CsPad::CalibV1/%s/'%(self.expname[:3],self.expname,srcStr)        
            elif det.dettype==13:
                dirname='/reg/d/psdm/%s/%s/calib/Epix100a::CalibV1/%s/'%(self.expname[:3],self.expname,srcStr)
            elif det.dettype==19:
                dirname='/reg/d/psdm/%s/%s/calib/Camera::CalibV1/%s/'%(self.expname[:3],self.expname,srcStr)        
            elif det.dettype==29:
                dirname='/reg/d/psdm/%s/%s/calib/Epix10ka::CalibV1/%s/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            elif det.dettype==32:
                dirname='/reg/d/psdm/%s/%s/calib/Epix10ka2M::CalibV1/%s/'%(self.sda.expname[:3],self.sda.expname,srcStr)        

        if not os.path.exists(dirname+'pedestals'):
            os.makedirs(dirname+'pedestals')
        if not os.path.exists(dirname+'pixel_rms'):
            os.makedirs(dirname+'pixel_rms')
        if not os.path.exists(dirname+'pixel_status'):
            os.makedirs(dirname+'pixel_status')
        print('save pedestal file in %s as %s '%(dirname+'pedestals/',fname))
        if useFilter is not None:
            det.save_txtnda(dirname+'pedestals/'+fname,self.__dict__['AvImg_median_Filter%s_raw_%s'%(useFilter,detname)], fmt='%.1f',addmetad=True)
        else:
            det.save_txtnda(dirname+'pedestals/'+fname,self.__dict__['AvImg_median_raw_%s'%(detname)], fmt='%.1f',addmetad=True)
        print('save noise file in %s as %s '%(dirname+'pixel_rms/',fname))
        print('save status file in %s as %s '%(dirname+'pixel_status/',fname))
        if useFilter is not None:
            det.save_txtnda(dirname+'pixel_rms/'+fname,self.__dict__['AvImg_std_Filter%s_raw_%s'%(useFilter,detname)], fmt='%.3f',addmetad=True)
        else:
            det.save_txtnda(dirname+'pixel_rms/'+fname,self.__dict__['AvImg_std_raw_%s'%(detname)], fmt='%.3f',addmetad=True)
        det.save_txtnda(dirname+'pixel_status/'+fname,status, fmt='%d',addmetad=True)

    def addAzInt(self, detname=None, phiBins=1, qBin=0.01, eBeam=9.5, center=None, dis_to_sam=None, name='azav', Pplane=1,userMask=None,tx=0,ty=0, geomCorr=True, polCorr=True, inImage={}):
        if inImage!={}:
            img=inImage['image']
            if 'mask' in inImage:
                mask=inImage['mask']
            else:
                mask=np.ones_like(img)
            detname=inImage['detname']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)
        else:
            detname, img, avImage = self.getAvImage(detname=detname)

        if dis_to_sam==None:
            dis_to_sam=float(raw_input('please enter the detector distance'))
        if center==None or len(center)!=2:
            centerString=raw_input('please enter the coordinates of the beam center as c1,c2 or [c1,c2]:')
            center=[float(centerString.replace('[','').replace(']','').split(',')[0]),
                    float(centerString.replace('[','').replace(']','').split(',')[1])]

        azav = azimuthalBinning(center=center, dis_to_sam=dis_to_sam,  phiBins=phiBins, eBeam=eBeam, Pplane=Pplane, userMask=userMask, qbin=qBin, tx=tx, ty=ty, geomCorr=geomCorr, polCorr=polCorr, name=name)
        getattr(self,detname).addFunc(azav)

    def getAzAvs(self,detname=None):
        if detname is None:
            detname, img, avImage = self.getAvImage(detname=None)   
            if detname is None:
                return
        azintArray = [ getattr(getattr(self,detname),key) for key in getattr(self, detname).__dict__.keys() if isinstance(getattr(getattr(self, detname),key), azimuthalBinning) ]
        azintNames = [ key for key in getattr(self,detname).__dict__.keys() if isinstance(getattr(getattr(self, detname),key), azimuthalBinning) ]
        return azintNames, azintArray

    def AzInt(self, detname=None, use_mask=False, use_mask_local=False, plotIt=False, azintName=None, inImage={}, imgName=None):
        avImage=None
        if inImage!={}:
            img=inImage['image']
            if 'mask' in inImage:
                mask=inImage['mask']
            else:
                mask=np.ones_like(img)
            avImage='input Image'
            if 'name' in inImage:
                avImage=inImage['name']
            detname=inImage['detname']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)
        else:
            detname, img, avImage = self.getAvImage(detname=detname, imgName=None)
        mask=np.ones_like(img)
        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
        if use_mask_local:
            if avImage is None:
                detname_dump, img_dump, avImage = self.getAvImage(detname=detname)
            if self.__dict__.has_key('_mask_'+avImage):
                mask = self.__dict__['_mask_'+avImage]
            else:
                print('no local mask defined for ',avImage)

        img = (img*mask)
        azIntNames,azIntegrations = self.getAzAvs(detname)
        if len(azIntegrations)>1:
            if azintName is None:
                print('we have the following options: ',azIntNames)
                azintName=raw_input('type the name of the Azimuthal integral to use:')
            azIntegs = [ iiaz for iName, iiaz in zip(azIntNames,azIntegrations) if iName==azintName]
            azInteg = azIntegs[0]
        else:
            azInteg = azIntegrations[0]
            azintName = azIntNames[0]
        azintValues = azInteg.doCake(img).squeeze()
        self.__dict__['_azint_'+azintName] = azintValues
        if plotIt:
            fig=plt.figure(figsize=(8,5))
            if len(azintValues.shape)==1:
                qVals = getattr(getattr(getattr(self,detname), azintName),'q')
                plt.plot(qVals,azintValues,'o')
            elif len(azintValues.shape)==2:
                plt.imshow(azintValues,aspect='auto',interpolation='none')
                plt.colorbar()
        else:
            return azintValues

    def plotAzInt(self, detname=None, azintName=None):
        detname, img, avImage = self.getAvImage(detname=detname)
        azIntNames,azIntegrations = self.getAzAvs(detname)
        azint=''
        if len(azIntegrations)>1:
            if azintName is None:
                print('we have the following options: ',azIntNames)
                azintName=raw_input('type the name of the Azimuthal integral to use:')
            azIntegs = [ iiaz for iName, iiaz in zip(azIntNames,azIntegrations) if iName==azintName]
            azInteg = azIntegs[0]
            azint = azintName
        else:
            azInteg = azIntegrations[0]
            azint = azIntNames[0]
        if azint=='':
            print('did not find azimuthal integral asked for')
            return
        else:
            if ('_azint_'+azint) not in self.__dict__.keys():
                print('did not find azint ',azint,', all keys are: ',self.__dict__.keys())
        try:
            azintValues = self.__dict__['_azint_'+azint]
            fig=plt.figure(figsize=(8,5))
            if len(azintValues.shape)==1:
                qVals = getattr(getattr(getattr(self,detname), azint),'q')
                print(azintValues.shape, qVals.shape)
                plt.plot(qVals,azintValues,'o')
            elif len(azintValues.shape)==2:
                plt.imshow(azintValues,aspect='auto',interpolation='none')
                plt.colorbar()
        except:
            pass

    def AzInt_centerVar(self, detname=None, use_mask=False, center=None, data=None, varCenter=110., zoom=None, qBin=0.001, phiBins=13, dis_to_sam=1000., inImage={}, imgName=None, plotGrid=0):
        if inImage!={}:
            img=inImage['image']
            if 'mask' in inImage:
                mask=inImage['mask']
            else:
                mask=np.ones_like(img)
            detname=inImage['detname']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)
        else:
            detname, img, avImage = self.getAvImage(detname=detname, imgName=imgName)
        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
            img = (img*mask)

        if center==None or len(center)!=2:
            centerString=raw_input('please enter the coordinates of the beam center as c1,c2 or [c1,c2]:')
            center=[float(centerString.replace('[','').replace(']','').split(',')[0]),
                    float(centerString.replace('[','').replace(']','').split(',')[1])]
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, center=[center[0],center[1]], dis_to_sam=dis_to_sam, name='c00', inImage=inImage)
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, center=[center[0]-varCenter,center[1]], dis_to_sam=dis_to_sam, name='cm10', inImage=inImage)
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, center=[center[0]+varCenter,center[1]], dis_to_sam=dis_to_sam, name='cp10', inImage=inImage)
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, center=[center[0],center[1]-varCenter], dis_to_sam=dis_to_sam, name='cm01', inImage=inImage)
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, center=[center[0],center[1]+varCenter], dis_to_sam=dis_to_sam, name='cp01', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='c00', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='cm10', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='cp10', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='cm01', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='cp01', inImage=inImage)

        azintValues_c00 = self.__dict__['_azint_c00']
        azintValues_cm10 = self.__dict__['_azint_cm10']
        azintValues_cp10 = self.__dict__['_azint_cp10']
        azintValues_cm01 = self.__dict__['_azint_cm01']
        azintValues_cp01 = self.__dict__['_azint_cp01']
        
        try:
            fig=plt.figure(figsize=(10,6))
            from matplotlib import gridspec
            p33_11 = plt.subplot2grid((3,3),(1,1))
            p33_10 = plt.subplot2grid((3,3),(1,0))
            p33_12 = plt.subplot2grid((3,3),(1,2))
            p33_01 = plt.subplot2grid((3,3),(0,1))
            p33_21 = plt.subplot2grid((3,3),(2,1))
        except:
            print('failed with getting plot ready')
            return

        while 1:
            ymin=0;ymax=azintValues_c00.shape[1]-1
            maxVal = np.nanpercentile(azintValues_c00, 99.99)
            if zoom is not None:
                if len(zoom)==1:
                    maxVal = np.nanpercentile(azintValues_c00, zoom[0])
                    zoom = None
                elif len(zoom)==3:
                    maxVal = np.nanpercentile(azintValues_c00, zoom[2])
                else:
                    if isinstance(zoom, list) and isinstance(zoom[0], int) and zoom[0]>=0 and zoom[1]>=0 and zoom[0]!=zoom[1]:
                        ymin=zoom[0]
                        ymax=zoom[1]
                    else:
                        yString=raw_input('please enter x-boundaries of the zoomed figure as c1,c2 or [c1,c2]:')
                        ymin=int(yString.replace('[','').replace(']','').split(',')[0])
                        ymax=int(yString.replace('[','').replace(']','').split(',')[1])
            print('we will plot from bin %d to %d '%(ymin, ymax))

            try:
                p33_11.imshow(azintValues_c00[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                p33_10.imshow(azintValues_cm10[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                p33_12.imshow(azintValues_cp10[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                p33_01.imshow(azintValues_cm01[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                p33_21.imshow(azintValues_cp01[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                maxVal_c00 =  np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_c00[:,ymin:ymax])])
                maxVal_cm10 = np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_cm10[:,ymin:ymax])])
                maxVal_cm01 = np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_cm01[:,ymin:ymax])])
                maxVal_cp10 = np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_cp10[:,ymin:ymax])])
                maxVal_cp01 = np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_cp01[:,ymin:ymax])])
                try:
                    p33_11.plot(maxVal_c00[:,1], maxVal_c00[:,0],'r+')
                    p33_10.plot(maxVal_cm10[:,1], maxVal_cm10[:,0],'r+')
                    p33_12.plot(maxVal_cp10[:,1], maxVal_cp10[:,0],'r+')
                    p33_01.plot(maxVal_cm01[:,1], maxVal_cm01[:,0],'r+')
                    p33_21.plot(maxVal_cp01[:,1], maxVal_cp01[:,0],'r+')
                    if plotGrid==0:
                        cMaxVal_c00 = np.nanargmax(azintValues_c00[:,ymin:ymax].sum(axis=0))
                        p33_11.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                        p33_10.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                        p33_12.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                        p33_01.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                        p33_21.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                    else:
                        for il,lineVal in enumerate(np.linspace(ymin,ymax,int(plotGrid))):
                            lw=1
                            if plotGrid>15 and il%5>0:
                                lw=0.5
                            p33_11.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                            p33_10.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                            p33_12.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                            p33_01.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                            p33_21.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                except:
                    print('failed at plotting')

                if ymin==0 and ymax==azintValues_c00.shape[1]-1:
                    break
                if raw_input("done? (y/n):\n") in ["y","Y"]:
                    break
            except:
                print('failed in try-except.')
                break

    def SelectRegionDroplet(self, detname=None, limits=[5,99.5]):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('AvImg')==0:
                if detname is not None and key.find(detname)>=0:
                    avImages.append(key)
                elif detname is None:
                    avImages.append(key)
        if len(avImages)==0:
            print('please create the AvImage first!')
            return
        elif len(avImages)>1:
            print('we have the following options: ',avImages)
            avImage=raw_input('type the name of the AvImage to use:')
        else:
            avImage=avImages[0]
        detname = self._getDetName_from_AvImage(avImage)
        img = self.__dict__[avImage]
        print(img.shape)

        plotMax = np.nanpercentile(img, limits[1])
        plotMin = np.nanpercentile(img, limits[0])
        print('plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if needsGeo:
            image = self.__dict__[detname].det.image(self.run, img)
        else:
            image = img

        det = self.__dict__[detname].det
        x = self.__dict__[detname+'_x']
        y = self.__dict__[detname+'_y']
        ix = self.__dict__[detname+'_ix']
        iy = self.__dict__[detname+'_iy']
        extent=[x.min(), x.max(), y.min(), y.max()]

        fig=plt.figure(figsize=(10,6))
        from matplotlib import gridspec
        gs=gridspec.GridSpec(1,2,width_ratios=[1,1])
        
        mask=None
        mask_r_nda=None

        plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')

        print('select two corners: ')
        p =np.array(ginput(2))
        mask_roi=np.zeros_like(image)
        mask_roi[p[:,1].min():p[:,1].max(),p[:,0].min():p[:,0].max()]=1
        if needsGeo:
            mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(ix,iy)] )
        else:
            mask_r_nda = mask_roi

        if mask_r_nda is not None:
            #print('created a mask....')
            if mask is None:
                mask = mask_r_nda.astype(bool).copy()
            else:
                mask = np.logical_or(mask,mask_r_nda)
        print('masked: ',np.ones_like(x)[mask.astype(bool)].sum())

        image_mask = img.copy(); image_mask[~mask]=0;
        if needsGeo:
            plt.subplot(gs[1]).imshow(det.image(self.run,image_mask),clim=[plotMin,plotMax])
        else:
            plt.subplot(gs[1]).imshow(image_mask,clim=[plotMin,plotMax])

        if  det.dettype == 2:
            mask=mask.reshape(2,185*388).transpose(1,0)
        elif det.dettype == 1:
            mask=mask.reshape(32*185,388)
        elif det.dettype == 13:
            mask=mask.reshape(704,768)

        #2x2 save as 71780 lines, 2 entries
        #cspad save as 5920 lines, 388 entries
        if raw_input("Save to calibdir?\n") in ["y","Y"]:
            mask = (~(mask.astype(bool))).astype(int)
            srcStr=det.source.__str__().replace('Source("DetInfo(','').replace(')")','')
            if det.dettype==2:
                dirname='/reg/d/psdm/%s/%s/calib/CsPad2x2::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)
            elif det.dettype==2:
                dirname='/reg/d/psdm/%s/%s/calib/CsPad::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            elif det.dettype==13:
                dirname='/reg/d/psdm/%s/%s/calib/Epix100a::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            elif det.dettype==29:
                dirname='/reg/d/psdm/%s/%s/calib/Epix10ka::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            elif det.dettype==32:
                dirname='/reg/d/psdm/%s/%s/calib/Epix10ka2M::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            fname='%s-end.data'%self.run
            det.save_txtnda(dirname+fname,mask, fmt='%d',addmetad=True)
        elif raw_input("Save to local?\n") in ["y","Y"]:
            mask = (~(mask.astype(bool))).astype(int)
            np.savetxt('%s_mask_run%s.data'%(self.sda.expname,self.run),mask)
        return mask
         
    def plotCalib(self, detname='None',common_mode=0, plotWith=None):
        if not detname in self.__dict__.keys() or self.__dict__[detname].common_mode!=common_mode:
            detname = self.addDetInfo(detname=detname, common_mode=common_mode)
            if detname == 'None':
                print('need detector name as input! ')
                return
        det=self.__dict__[detname].det
        rms = self.__dict__[detname+'_rms']
        pedestals = self.__dict__[detname+'_pedestals']
        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        np.seterr(divide='ignore', invalid='ignore')
        hstPed=[]
        hstRms=[]
        if det.dettype==26 or det.dettype==30:
            for ped,frms in zip(pedestals, rms):
                hstPed.append(np.histogram(ped.flatten(), np.arange(0, pedestals.max()*1.05)))
                hstRms.append(np.histogram(frms.flatten(), np.arange(0, rms.max()*1.05,0.05)))
        else:
            hstPed.append(np.histogram(pedestals.flatten(), np.arange(0, pedestals.max()*1.05)))
            hstRms.append(np.histogram(rms.flatten(), np.arange(0, rms.max()*1.05,0.05)))

        if needsGeo:
            if det.dettype==26:
                pedImg = det.image(self.run,pedestals[0])
                rmsImg = det.image(self.run,rms[0])
                if pedImg is None and pedestals.shape[1]==1:
                    pedImg = pedestals[0].squeeze()
                    rmsImg = rms[0].squeeze()
            elif det.dettype==30:
                pedImg = pedestals[0]; rmsImg = rms[0]
                for iFrame, pedf, rmsf in zip(itertools.count(), pedestals, rms):
                    if iFrame>0:
                        pedImg = np.append(pedImg, pedf, axis=1)
                        rmsImg = np.append(rmsImg, rmsf, axis=1)
            else:
                pedImg = det.image(self.run,pedestals)
                rmsImg = det.image(self.run,rms)
        else:
            pedImg = pedestals
            rmsImg = rms

        if plotWith is None:
            plotWith=self.plotWith
        if plotWith.find('bokeh')>=0:
            #return hstPed, hstRms, pedImg, rmsImg
            xmin=self.__dict__[detname+'_x'].min()
            xmax=self.__dict__[detname+'_x'].max()
            ymin=self.__dict__[detname+'_y'].min()
            ymax=self.__dict__[detname+'_y'].max()
            bounds = (xmin, ymin, xmax, ymax)
            plotwidth=400
            zPed=(np.nanpercentile(pedImg,5), np.nanpercentile(pedImg,99.5))
            zRms=(np.nanpercentile(rmsImg,5), np.nanpercentile(rmsImg,99.5))
            hstPedO = hv.Overlay([hv.Scatter((0.5*(hst[1][:-1]+hst[1][1:]), np.log(hst[0])),'pedestal','log(nEntries)') for hst in hstPed]).options(width=plotwidth)
            hstRmsO = hv.Overlay([hv.Scatter((0.5*(hst[1][:-1]+hst[1][1:]), np.log(hst[0])),'rms','log(nEntries)') for hst in hstRms]).options(width=plotwidth)
            
            #need to invert the y-axis.
            imgPed = hv.Image(pedImg.T, bounds=bounds, vdims=[hv.Dimension('pedData', range=zPed)]).options( cmap='viridis',colorbar=True, invert_yaxis=True, width=plotwidth)
            imgRms = hv.Image(rmsImg.T, bounds=bounds, vdims=[hv.Dimension('rmsData', range=zRms)]).options( cmap='viridis',colorbar=True, invert_yaxis=True, width=plotwidth) 
            layout = (imgPed + hstPedO + imgRms + hstRmsO).cols(2)
            return layout

        #now for the matplotlib version.
        figDark=plt.figure(figsize=(11,6))
        gsPed=gridspec.GridSpec(2,2,width_ratios=[1,1])
        for i in range(len(hstPed)):
            plt.subplot(gsPed[1]).plot(hstPed[i][1][:-1],np.log(hstPed[i][0]),'o')
            plt.subplot(gsPed[3]).plot(hstRms[i][1][:-1],np.log(hstRms[i][0]),'o')
        if det.dettype==26:
            im0 = plt.subplot(gsPed[0]).imshow(pedImg,clim=[np.nanpercentile(pedestals[0],1),np.nanpercentile(pedestals[0],99)])
            cbar0 = plt.colorbar(im0)
            im2 = plt.subplot(gsPed[2]).imshow(rmsImg,clim=[np.nanpercentile(rms[0],1),np.nanpercentile(rms[0],99)])
        else:
            im0 = plt.subplot(gsPed[0]).imshow(pedImg,clim=[np.nanpercentile(pedestals,1),np.nanpercentile(pedestals,99)])
            cbar0 = plt.colorbar(im0)
            im2 = plt.subplot(gsPed[2]).imshow(rmsImg,clim=[np.nanpercentile(rms,1),np.nanpercentile(rms,99)])
        cbar2 = plt.colorbar(im2)

    def _fillCalibHisto(self, detname='None', printVal=[-1], pedBinWidth=10, rmsBinWidth=0.1):
        if not detname in self.__dict__.keys():
            detname = self.addDetInfo(detname=detname)
            if detname == 'None':
                print('need detector name as input! ')
                return
        det=self.__dict__[detname].det
        rms = self.__dict__[detname+'_rms']
        pedestals = self.__dict__[detname+'_pedestals']
        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True
        #now look at directory and get filenames of darks. Hmm. extra function?
        detNameStr = det.name.__str__()
        if detNameStr.find('Epix')>=0 or detNameStr.find('epix')>=0:
            if detNameStr.find('10k')>=0:
                detTypeStr='Epix10ka2M::CalibV1'
            else:
                detTypeStr='Epix100a::CalibV1'
        elif detNameStr.find('ungfrau')>=0:
            detTypeStr='Jungfrau::CalibV1'
        elif  detNameStr.find('2x2')>=0:
            detTypeStr='CsPad2x2::CalibV1'
        else:
            detTypeStr='CsPad::CalibV1'
        fnames = os.listdir('/reg/d/psdm/%s/%s/calib/%s/%s/pedestals/'%(self.hutch, self.expname, detTypeStr, detNameStr))
        pedNames=[]
        pedRuns=[]
        for fname in fnames:
            if fname[-5:]=='.data':
                pedNames.append(fname)
                pedRuns.append(int(fname.split('-')[0]))

        allPeds=[]
        allRms=[]
        allStatus=[]
        allPedsImg=[]
        allRmsImg=[]
        pedRuns.sort()
        for ipedRun,pedRun in enumerate(pedRuns):
            allPeds.append(det.pedestals(pedRun))
            allRms.append(det.rms(pedRun))
            allStatus.append((det.mask(pedRun, status=True)==0).sum())
            if needsGeo:
                try:
                    if detTypeStr=='Jungfrau::CalibV1' or detTypeStr=='Epix10ka2M::CalibV1':
                        allPedsImg.append(det.image(self.run, allPeds[-1][0]).tolist())
                        allRmsImg.append(det.image(self.run, allRms[-1][0]).tolist())
                    else:
                        allPedsImg.append(det.image(self.run, allPeds[-1]).tolist())
                        allRmsImg.append(det.image(self.run, allRms[-1]).tolist())
                except:
                    #get the first tile/gain/timepoint/...
                    allPedsImg.append(allPeds[-1][0].tolist())
                    allRmsImg.append(allRms[-1][0].tolist())
            else:
                allPedsImg.append(allPeds[-1].tolist())
                allRmsImg.append(allRms[-1].tolist())
            if len(printVal)<2:
                print('getting pedestal from run ',pedRun)
            elif len(printVal)>=4:
                print('run %d, pixel cold/hot, low noise, high noise: %d / %d / %d / %d pixels'%(pedRun,(allPeds[-1]<printVal[0]).sum(), (allPeds[-1]>printVal[1]).sum(),(allRms[-1]<printVal[2]).sum(), (allRms[-1]>printVal[3]).sum()))
            elif len(printVal)>=2:
                print('run %d, pixel cold/hot: %d / %d pixels'%(pedRun,(allPeds[-1]<printVal[0]).sum(), (allPeds[-1]>printVal[1]).sum()))

        allPeds = np.array(allPeds)
        allRms = np.array(allRms)
        ped2d=[] 
        ped2dmax=[] 
        ped2dfwhm=[] 
        rms2d=[] 

        pedHisMin = np.nanmin(allPeds)
        pedHisMax = np.nanmax(allPeds)
        pedHisLow = np.nanpercentile(allPeds,0.5)
        pedHisHigh = np.nanpercentile(allPeds,99.5)
        rmsHisMin = np.nanmin(allRms)
        rmsHisMax = np.nanmax(allRms)
        rmsHisLow = np.nanpercentile(allRms,0.5)
        rmsHisHigh = np.nanpercentile(allRms,99.5)
        #print 'maxmin peds ',allPeds[0].shape, pedHisMax  , pedHisMin , pedHisLow, pedHisHigh
        for thisPed,thisRms in zip(allPeds,allRms):
            ped2d.append(np.histogram(thisPed.flatten(), np.arange(pedHisLow*0.8, pedHisHigh*1.2,pedBinWidth))[0])
            ped2dmax.append(np.argmax(ped2d[-1])*pedBinWidth+pedHisLow*0.8)
            maxwhere = np.argwhere(ped2d[-1]>(np.nanmax(ped2d[-1])*0.5))
            try:
                if isinstance((maxwhere[-1]-maxwhere[0]), np.ndarray):
                    ped2dfwhm.append((maxwhere[-1]-maxwhere[0])[0]*pedBinWidth)
                else:
                    ped2dfwhm.append(maxwhere[-1]-maxwhere[0]*pedBinWidth)
            except:
                ped2dfwhm.append(0)
            rms2d.append(np.histogram(thisRms.flatten(), np.arange(rmsHisLow*0.8, rmsHisHigh*1.2, rmsBinWidth))[0])

        ped2d = np.array(ped2d)
        rms2d = np.array(rms2d)

        #save array to SmallDataAna_psana object to possibly be used again.
        #separate plotting functions later.
        self.calibhisto['status'] = np.array(allStatus)
        self.calibhisto['rms3d'] = np.array(allRmsImg)
        self.calibhisto['ped3d'] = np.array(allPedsImg)
        self.calibhisto['runs'] = np.array(pedRuns)
        self.calibhisto['ped2d'] = np.array(ped2d)
        self.calibhisto['ped2dmax'] = np.array(ped2dmax)
        self.calibhisto['ped2dfwhm'] = np.array(ped2dfwhm)
        self.calibhisto['rms2d'] = np.array(rms2d)
        self.calibhisto['pedBinWidth'] = pedBinWidth
        self.calibhisto['pedBinLow'] = pedHisLow*0.8
        self.calibhisto['rmsBinWidth'] = rmsBinWidth
        self.calibhisto['rmsBinLow'] = rmsHisLow*0.8

    def plotCalibHisto(self, detname='None',  printVal=[-1], pedBinWidth=10, rmsBinWidth=0.1, plotWith=None, showPlot=None):
        if 'rms3d' not in self.calibhisto.keys():
            self._fillCalibHisto(detname=detname, printVal=printVal, pedBinWidth=pedBinWidth, rmsBinWidth=rmsBinWidth)        
        if plotWith is None:
            plotWith=self.plotWith

        if plotWith=='bokeh_notebook':
            plotwidth=400
            boundsPed = (0, np.array(self.calibhisto['runs']).min(),  allPeds[icurrRun].max(), np.array(self.calibhisto['runs']).max())
            pedZ = (np.nanpercentile(self.calibhisto['ped2d'],1),np.nanpercentile(self.calibhisto['ped2d'],99.5))
            boundsRms = (0, np.array(self.calibhisto['runs']).min(),  allRms[icurrRun].max(), np.array(self.calibhisto['runs']).max())
            rmsZ = (np.nanpercentile(self.calibhisto['rms2d'],1),np.nanpercentile(self.calibhisto['rms2d'],99.5))
            runDim=hv.Dimension('run')

            #replace norm plots by zoomed plots.
            if showPlot is None:
                return (hv.Image(self.calibhisto['ped2d'], kdims=[hv.Dimension('pedestal'),runDim], bounds=boundsPed, vdims=[hv.Dimension('pedVals', range=pedZ)]).options( cmap='viridis',colorbar=True, width=plotwidth) + \
                        hv.Image(self.calibhisto['ped2d'], kdims=[hv.Dimension('pedestal'),runDim], bounds=boundsPed, vdims=[hv.Dimension('pedVals', range=pedZ)]).options( cmap='viridis',colorbar=True, width=plotwidth) + \
                        #hv.Image(self.calibhisto['ped2dNorm'], kdims=[hv.Dimension('pedestal/ref'),runDim], bounds=boundsPed, vdims=[hv.Dimension('pedValsNorm', range=pedZNorm)]).options( cmap='viridis',colorbar=True, width=plotwidth) + \
                        hv.Image(self.calibhisto['rms2d'], kdims=[hv.Dimension('rmsestal'),runDim], bounds=boundsRms, vdims=[hv.Dimension('rmsVals', range=rmsZ)]).options( cmap='viridis',colorbar=True, width=plotwidth) + \
                        #hv.Image(self.calibhisto['rms2d'], kdims=[hv.Dimension('rmsestal'),runDim], bounds=boundsRms, vdims=[hv.Dimension('rmsVals', range=rmsZ)]).options( cmap='viridis',colorbar=True, width=plotwidth) + \
                        hv.Image(self.calibhisto['rms2dNorm'], kdims=[hv.Dimension('rmsestal/ref'),runDim], bounds=boundsRms, vdims=[hv.Dimension('rmsValsNorm', range=rmsZNorm)]).options( cmap='viridis',colorbar=True, width=plotwidth) ).cols(2)
            elif showPlot == 'ped3d':
                return hv_3dimage(self.calibhisto['ped3d'])
            elif showPlot == 'rms3d':
                return hv_3dimage(self.calibhisto['rms3d'])
            else:
                print('this is not a defined option')
                return

        #fallback to matplotlib
        figDark=plt.figure(figsize=(11,10))
        gsPed=gridspec.GridSpec(2,2)

        plt.subplot(gsPed[0]).plot(self.calibhisto['runs'],self.calibhisto['ped2dmax'],'bo')
        plt.subplot(gsPed[0]).plot(self.calibhisto['runs'],self.calibhisto['ped2dmax']+self.calibhisto['ped2dfwhm'],'ro')
        plt.subplot(gsPed[0]).plot(self.calibhisto['runs'],self.calibhisto['ped2dmax']-self.calibhisto['ped2dfwhm'],'ro')

        plt.subplot(gsPed[1]).plot(self.calibhisto['runs'],self.calibhisto['status'],'o')

        ch_ped2d = self.calibhisto['ped2d']
        im0 = plt.subplot(gsPed[2]).imshow(ch_ped2d.T,vmin=np.nanpercentile(ch_ped2d,1),vmax=np.nanpercentile(ch_ped2d,99),interpolation='none',aspect='auto',extent=[0, ch_ped2d.shape[0], self.calibhisto['pedBinLow'],self.calibhisto['pedBinLow']+ch_ped2d.shape[1]*self.calibhisto['pedBinWidth']], origin='lower')
        cbar0 = plt.colorbar(im0)

        ch_rms2d = self.calibhisto['rms2d']
        im2 = plt.subplot(gsPed[3]).imshow(ch_rms2d.T,vmin=np.nanpercentile(ch_rms2d,1),vmax=np.nanpercentile(ch_rms2d,99),interpolation='none',aspect='auto',extent=[0, ch_rms2d.shape[0], self.calibhisto['rmsBinLow'],self.calibhisto['rmsBinLow']+ch_rms2d.shape[1]*self.calibhisto['rmsBinWidth']], origin='lower')
        cbar2 = plt.colorbar(im2)


    def compareCommonMode(self, detname='None',common_modes=[], numEvts=100, thresADU=0., thresRms=0., plotWith=None):
        if detname is 'None':
            detname = self.addDetInfo(detname=detname)
            if detname == 'None':
                print('need detector name as input! ')
                return

        if len(common_modes)==0:
            if detname.find('cs')>=0:
                common_modes=[1,5,0]
            elif detname.find('epix')>=0:
                common_modes=[46,45,4,0]
            else:
                common_modes=[0,-1]

        baseName='AvImg_'
        if thresADU!=0:
            baseName+='thresADU%d_'%int(thresADU)
        if thresRms!=0:
            baseName+='thresRms%d_'%int(thresRms*10.)
        imgNames=[]
        imgs=[]

        for cm in common_modes:
            self.AvImage(detname, numEvts=numEvts, thresADU=thresADU, thresRms=thresRms, common_mode=cm)
            imgNames.append('%s%s%s'%(baseName,self.commonModeStr(cm),detname))
            #add needsGeo clause here.
            if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
                imgs.append(self.__dict__[detname].det.image(self.run, self.__dict__[imgNames[-1]]))
            else:
                imgs.append(self.__dict__[imgNames[-1]])
            self.AvImage(detname, numEvts=numEvts, thresADU=thresADU, thresRms=thresRms, common_mode=cm, std=True)
            imgNames.append('%sstd_%s%s'%(baseName,self.commonModeStr(cm),detname))
            if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
                imgs.append(self.__dict__[detname].det.image(self.run, self.__dict__[imgNames[-1]]))
            else:
                imgs.append(self.__dict__[imgNames[-1]])

        if plotWith is None:
            plotWith=self.plotWith

        if plotWith=='bokeh_notebook':
            return hv_3dimage(imgs)

        fig=plt.figure(figsize=(15,10))
        gsCM=gridspec.GridSpec(len(common_modes),2)
        for icm,cm in enumerate(common_modes):
            if cm==0:
                lims=[np.nanpercentile(imgs[icm*2],1),np.nanpercentile(imgs[icm*2],99)]
            else:
                lims=[np.nanpercentile(imgs[0],1),np.nanpercentile(imgs[0],99)]
            limsStd=[np.nanpercentile(imgs[1],1),np.nanpercentile(imgs[1],99)]
            imC = plt.subplot(gsCM[icm*2]).imshow(imgs[icm*2],clim=lims,interpolation='none',aspect='auto')
            plt.colorbar(imC)
            imCS = plt.subplot(gsCM[icm*2+1]).imshow(imgs[icm*2+1],clim=limsStd,interpolation='none',aspect='auto')
            plt.colorbar(imCS)

    def makeCubeData(self, cubeName, dirname='', nEvtsPerBin=-1, offEventsCube=-1, storeMeanStd=False, onoff=2):
        if self.sda is None:
            return
        if dirname=='':
            dirname=self.sda.dirname

        myCube, cubeName_onoff = self.sda.prepCubeData(cubeName)
        
        print('Variables to be read from xtc: ',myCube.targetVarsXtc)

        if isinstance(nEvtsPerBin, str): nEvtsPerBin=int(nEvtsPerBin)

        #smallData only version: only run on single core
        if len(myCube.targetVarsXtc)<=0:
            if rank==0:
                outFileName = dirname+'/Cube_'+self.sda.fname.split('/')[-1].replace('.h5','_%s.h5'%cubeName)
                if (onoff==0):
                    outFileName=outFileName.replace('.h5','_off.h5')
                elif (onoff==1):
                    outFileName=outFileName.replace('.h5','_on.h5')
                fout = h5py.File(outFileName, "w")
                #ADD CNF STUFF HERE!
                printR(rank, 'no big data, bin the data now....be patient')
                cubeData = self.sda.makeCubeData(cubeName,onoff=onoff)
                printR(rank, 'now write outputfile (only small data) to : %s'%outFileName)
                for key in cubeData.variables:
                    addToHdf5(fout, key, cubeData[key].values)

                for cfgVar in myCube.targetVarsCfg:
                    addToHdf5(fout, cfgVar.replace('/','_'), self.sda.getVar(cfgVar))
                    print('add cfgVar to hdf5', cfgVar.replace('/','_'))
                fout.close()
            return cubeData

        printR(rank,'Now make big cube')
        #only run on rank=0 & broadcast.
        outFileName = dirname+'/Cube_'+self.sda.fname.split('/')[-1].replace('.h5','_%s.h5.inprogress'%cubeName)
        if (onoff==0):
            outFileName=outFileName.replace('.h5','_off.h5')
        elif (onoff==1):
            outFileName=outFileName.replace('.h5','_on.h5')
        printR(rank, 'now write outputfile to : %s'%outFileName)
        
        #configuration for cube making
#         cube, cubeName_onoff = self.sda.prepCubeData(cubeName)
        #compare the number of bins to mpi jobs.
#         nBins=myCube.binBounds.shape[0]-1
#         for key in cube.addBinVars:
#             nBins*=cube.addBinVars[key].shape[0]-1
#         if nBins<size:
#             if 'random/random' not in self.sda.Keys('random'):
#                 myrandom=np.random.rand(self.sda.xrData.time.shape[0])
#                 self.sda.addVar('random/random',myrandom)
#             if size/nBins > 1:
#                 cube.add_BinVar({'random/random':[0.,1.,int(size/nBins)]})

        sel = self.sda.Sels[myCube.useFilter]
        selString=''
        for icut,cut in enumerate(sel.cuts):
            selString+=('Cut %i: %f < %s < %f\n'%(icut, cut[1], cut[0],cut[2]))
            
        fout = h5py.File(outFileName, "w")
        dsetcnf = fout.create_dataset('cubeSelection', [1.], dtype='f')
        dsetcnf.attrs['cubeSelection'] = selString

        #back to original code.
        if offEventsCube>0:
            self.sda.getOffVar('fiducials','xon',nNbr=offEventsCube, mean=False)
            self.sda.getOffVar('event_time','xon',nNbr=offEventsCube, mean=False)
            addVars = [ 'offNbrs_event_time_xon_nNbr%02d'%offEventsCube,'offNbrs_fiducials_xon_nNbr%02d'%offEventsCube]
            self.sda.cubes[cubeName].addIdxVar(addVars)
        cubeData, eventIdxDict = self.sda.makeCubeData(cubeName, returnIdx=True, onoff=onoff)
        
        # add small data to hdf5
        for key in cubeData.variables:
            addToHdf5(fout, key, cubeData[key].values)

        # Cube big data
        t0 = time.time()
        bins_info = []
        for f,t in zip(eventIdxDict['fiducial'], eventIdxDict['evttime']):
            if nEvtsPerBin>0:
                bins_info.append([np.asarray(f[:nEvtsPerBin]), np.asarray(t[:nEvtsPerBin])])
            else:
                bins_info.append([np.asarray(f), np.asarray(t)])
        nbins = len(bins_info)
        print('****** Total number of bins: {}'.format(nbins))
        #print('****** BINS: {}'.format(bin))
        print('****** Make big data placeholder dataset and save det config.')
        dets = []
        for detname, detDict in zip(self.detNames, myCube.targetVarsXtc):
            det = self.__dict__[detname]
            if detDict['image']==1:
                det_shape = det.imgShape
            else:
                det_shape = det.mask.shape
            try:
                self.make_det_data_dset(fout, detname, det_shape, nbins)
            except Exception as e:
                logger.warning('Could not make dataset for detector {}. Exit. {}'.format(detname, e))
                comm.Abort()
                
            # save detector config
            logger.info(f'Save config for det {detname}.')
            if det.rms is not None:
                grp = fout.create_group(f'{detname}_cfg')
                if detDict['image']==0:
                    addToHdf5(grp, 'ped', det.ped)
                    addToHdf5(grp, 'rms', det.rms)
                    if det.gain is not None:
                        addToHdf5(grp, 'gain', det.gain)
                    addToHdf5(grp, 'mask', det.mask)
                    addToHdf5(grp, 'calib_mask', det.cmask)
                    if det.x is not None:
                        addToHdf5(grp, 'x', det.x)
                        addToHdf5(grp, 'y', det.y)
                    if det.ix is not None:
                        addToHdf5(grp, 'ix', det.ix)
                        addToHdf5(grp, 'iy', det.iy)
                else:
                    if det.det.dettype==26:
                        addToHdf5(grp, 'ped', det.det.image(self.run,det.ped[0]))
                        addToHdf5(grp, 'rms', det.det.image(self.run,det.rms[0]))
                        addToHdf5(grp, 'gain', det.det.image(self.run,det.gain[0]))
                    else:
                        addToHdf5(grp, 'ped', det.det.image(self.run,det.ped))
                        addToHdf5(grp, 'rms', det.det.image(self.run,det.rms))
                        addToHdf5(grp, 'gain', det.det.image(self.run,det.gain))
                    addToHdf5(grp, 'mask', det.det.image(self.run,det.mask))
                    addToHdf5(grp, 'calib_mask', det.det.image(self.run,det.cmask))
                    if det.x is not None:
                        addToHdf5(grp, 'x', det.x)
                        addToHdf5(grp, 'y', det.y)
                    if det.ix is not None:
                        addToHdf5(grp, 'ix', det.ix)
                        addToHdf5(grp, 'iy', det.iy)
        
        print('Start binning area detectors')                
        save_fct = lambda data=None, bin_idx=None: self.save_bin_to_h5(fout=fout, data=data, bin_idx=bin_idx)
        sum_data = mpi_fun.bin_distribution(bins_info, func=save_fct)
        t3 = time.time()

        print("***** ALL BINS DONE AFTER {:0.2f} min. *****".format((t3-t0)/60))

        print(f'Renaming file from {outFileName} to {outFileName.replace(".inprogress","")}, remove random variable if applicable')
        rename_reduceRandomVar(outFileName)

        bins, nEntries = cubeData.binVar_bins.values, cubeData.nEntries.values
        return cubeName, bins, nEntries

    def _broadcast_xtc_dets(self, cubeName):
        """ Sends the xtc det info to worker so that they can instantiate the DetObjects. Most
        of it is legacy from the old cube.
        """
        myCube, cubeName_onoff = self.sda.prepCubeData(cubeName)
        detInData=[]
        for k in self.Keys():
            if k.alias()!='':
                detInData.append(k.alias())
        self.detNames=[]
        targetVarsXtc=[]
        for idet,det in enumerate(myCube.targetVarsXtc):
            if isinstance(det, dict):
                dName = det['source']
            else:
                dName = det
            if dName in detInData:
                self.detNames.append(dName)
                targetVarsXtc.append(det)
            else:
                printR(rank, 'Detector with alias %s not in data '%det)
        myCube.targetVarsXtc = [ {'source':det, 'full':1} if isinstance(det, str) else det for det in targetVarsXtc]
        for det in myCube.targetVarsXtc:
            self.addDetInfo(det)
        to_worker = myCube.targetVarsXtc
        comm.bcast(to_worker, root=0)
        return
        
    @staticmethod
    def make_det_data_dset(fout, detname, det_shape, nbins):
        # data dset
        dset_name = '{}_data'.format(detname)
        shape = tuple(np.r_[nbins,det_shape])
        dset = fout.create_dataset(dset_name, shape)
        # n_in_bin dset
        dset_name = '{}_nEntries'.format(detname)
        dset = fout.create_dataset(dset_name, (nbins,))
        return
    
    @staticmethod
    def save_bin_to_h5(fout=None, bin_idx=None, data=None):
        """ data[0]: summed_data, data[1]: n_in_bin
        """
        for detname in data[0].keys():
            # data
            dset_name = '{}_data'.format(detname)
            dset = fout[dset_name]
            dset[bin_idx] = data[0][detname]
            # n_in_bin
            dset_name = '{}_nEntries'.format(detname)
            dset = fout[dset_name]
            dset[bin_idx] = data[1][detname]
        return
    
    @staticmethod
    def print_me(data=None, bin_idx=None):
        print(data)
        return
