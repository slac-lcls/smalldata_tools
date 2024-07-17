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
import itertools
import os
import socket
import holoviews as hv
import bokeh.plotting as bp
import psana
import logging
import time
import requests
from pathlib import Path
try:
    basestring
except NameError:
    basestring = str
try:
    raw_input
except NameError:
    raw_input = input


import smalldata_tools.SmallDataAna as sda

from smalldata_tools.BaseSmallDataAna_psana import BaseSmallDataAna_psana
from smalldata_tools.DetObject import DetObject, DetObjectClass
from smalldata_tools.SmallDataUtils import getUserData
from smalldata_tools.utilities import printR
from smalldata_tools.utilities import addToHdf5
from smalldata_tools.utilities import rename_reduceRandomVar
from smalldata_tools.utilities_plotting import hv_image
from smalldata_tools.utilities_plotting import hv_image_ctl
from smalldata_tools.utilities_plotting import hv_3dimage
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, spectrumFunc, projectionFunc
from smalldata_tools.ana_funcs.sparsifyFunc import sparsifyFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning

try:
    import smalldata_tools.cube.cube_mpi_fun as mpi_fun
except:
    print("Can't import smalldata_tools.cube.cube_mpi_fun." \
          + " If you are not on a LCLS-II experiment, consider debugging.")

from mpi4py import MPI
import h5py
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SmallDataAna_psana(BaseSmallDataAna_psana):
    def __init__(self, expname='', run=-1, dirname='', filename='', plotWith='matplotlib'):
        super().__init__(expname=expname, run=run, dirname=dirname, filename=filename, plotWith=plotWith)
    
        self.dsname='exp=%s:run=%i:smd'%(expname,self.run)
        hostname = socket.gethostname()
        if hostname.find('drp-srcf')>=0:
            xtcdirname = '/cds/data/drpsrcf/%s/%s/xtc'%(self.hutch.lower(),expname)
            self.dsname += ':dir=%s'%xtcdirname
        self.dsnameIdx = self.dsname.replace('smd','idx').replace(':live','')
        if self.isLive:
            self.dsname=self.dsname+':live:stream=0-79'
            self.dsnameIdx = None
        comm.bcast(self.dsnameIdx, root=0)

        printR(rank, '\nMake SmallDataAna_psana from dsname: %s'%self.dsname)
        try:
            self.ds = psana.DataSource(str(self.dsname))
        except:
            printR(rank, 'Failed to set up psana datasource!')
            self.ds = None
        if self.dsnameIdx is None:
            printR(rank, 'Failed to set up index-based psana datasource, likely because this run is still live')
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
        return

    def resetDs(self, idx=True):
        if idx:
            del self.dsIdx
            self.dsIdx = psana.DataSource(self.sda.dsnameIdx)
            self.dsIdxRun = self.dsIdx.runs().next()
        else:
            del self.ds
            self.ds = psana.DataSource(self.sda.dsname)

    def Keys(self, printthis=False):
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

    def CfgKeys(self, idx=True, printthis=False):
        if idx:
            keys=self.dsIdx.env().configStore().keys()
        else:
            keys=self.ds.env().configStore().keys()
        if printthis: 
            print(keys)
        else:
            return keys

    def EpicsAliases(self, idx=True, printthis=False):
        if idx:
            keys=self.dsIdx.env().epicsStore().aliases()
        else:
            keys=self.ds.env().epicsStore().aliases()
        if printthis: 
            print(keys)
        else:
            return keys


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
            printR(rank, 'redefine detector object with different common mode:'\
                   '%d instead of %d'%( common_mode,self.__dict__[detname].common_mode))
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
    
    
    def pixelHistogram(self, 
                       detname=None, 
                       numEvts=100, 
                       nBins=180, 
                       useFilter=None, 
                       nSkip=0,maxHis=None,
                       common_mode=None, 
                       std=False,
                       printFid=False,
                       useMask=True, 
                       uniform=False):
        """ Return a pixel histogram
        
        Arguments:
            detname: detector name
            numEvts: number of events to be used
            nSkip: number of events to be skipped
            common_mode: calibration applied, including optional correction for common mode noise. 
                         Default choice will be used if none are passed if not supplied
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
            if det.evt.dat is not None:
                data.append(det.evt.dat.copy())
        data = np.asarray(data).ravel()
        
        # Find histogram boundaries
        # Assumes that most of the signal is around 0,
        # Roughly take 50 times 0-photon peak
        # Can be weird in case of particularly high signal...
        bmin = -2*np.abs(np.percentile(data, 10))
        bmax = 2*np.abs(np.percentile(data, 90))
        bmax = 50*np.std(data[np.logical_and(data>bmin, data<bmax)])
        bmin = -2*np.abs(np.percentile(data, 0.1))
        if maxHis is not None: bmax = np.max(maxHis, bmax)
        bin_centers = np.linspace(bmin,bmax,nBins)
        hist, bin_edges = np.histogram(data, bins=bin_centers)
        
        for tm in times[10:]:
            #print('numEvts ',numEvts)
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

    
    def AvImage(self, 
                detname='None', 
                numEvts=100, 
                thresADU=0., 
                thresRms=0., 
                useFilter=None, 
                nSkip=0,minIpm=-1., 
                common_mode=None, 
                std=False, 
                median=False, 
                printFid=False,
                useMask=True, 
                uniform=False, 
                returnEnv=False):
        """
        make an average (summed) image for a given detector
        if a detector name is not passed, you will be presented with a choice.
        
        Arguments:
            numEvts: number of events to be used
            nSkip: number of events to be skipped
            common_mode: calibration applied, including optional correction for common mode noise. 
                         Default choice will be used if none are passed if not supplied
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
        det = self.__dict__[detname]
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
            if minIpm!=-1 and ( 
                (self.hutch=='xpp' and evt.get(psana.Lusi.IpmFexV1, psana.Source('BldInfo(XppSb2_Ipm)')).sum() < minIpm)\
                or (self.hutch=='xcs' and evt.get(psana.Lusi.IpmFexV1, psana.Source('BldInfo(XCS-IPM-05)')).sum() < minIpm)
            ):
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

         
    def makePedestal(self, 
                     detname, 
                     useFilter=None, 
                     numEvts=-1, 
                     pedRange=[10,10000], 
                     rmsRange=[2.,7.], 
                     i0Check='ipm', 
                     dirname='./'):
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
            calibDir=det.env.calibDir()
            if det.dettype==DetObjectClass.CsPad:
                dirname='/%s/CsPad2x2::CalibV1/%s/'%(calibDir,srcStr)
            elif det.dettype==DetObjectClass.CsPad2M:
                dirname='%s/CsPad::CalibV1/%s/'%(calibDir,srcStr)        
            elif det.dettype==DetObjectClass.Epix:
                dirname='%s/Epix100a::CalibV1/%s/'%(calibDir,srcStr)
            elif det.dettype==DetObjectClass.Rayonix:
                dirname='%s/Camera::CalibV1/%s/'%(calibDir,srcStr)        
            elif det.dettype==DetObjectClass.Epix10k:
                dirname='%s/Epix10ka::CalibV1/%s/'%(calibDir,srcStr)
            elif det.dettype==DetObjectClass.Epix10k2M:
                dirname='%s/Epix10ka2M::CalibV1/%s/'%(calibDir,srcStr)

        if not os.path.exists(dirname+'pedestals'):
            os.makedirs(dirname+'pedestals')
        if not os.path.exists(dirname+'pixel_rms'):
            os.makedirs(dirname+'pixel_rms')
        if not os.path.exists(dirname+'pixel_status'):
            os.makedirs(dirname+'pixel_status')
        print('save pedestal file in %s as %s '%(dirname+'pedestals/',fname))
        if useFilter is not None:
            det.save_txtnda(dirname+'pedestals/'+fname, 
                            self.__dict__['AvImg_median_Filter%s_raw_%s'%(useFilter,detname)], 
                            fmt='%.1f',addmetad=True)
        else:
            det.save_txtnda(dirname+'pedestals/'+fname,
                            self.__dict__['AvImg_median_raw_%s'%(detname)], 
                            fmt='%.1f',addmetad=True)
        print('save noise file in %s as %s '%(dirname+'pixel_rms/',fname))
        print('save status file in %s as %s '%(dirname+'pixel_status/',fname))
        if useFilter is not None:
            det.save_txtnda(dirname+'pixel_rms/'+fname,
                            self.__dict__['AvImg_std_Filter%s_raw_%s'%(useFilter,detname)], 
                            fmt='%.3f',addmetad=True)
        else:
            det.save_txtnda(dirname+'pixel_rms/'+fname, 
                            self.__dict__['AvImg_std_raw_%s'%(detname)], 
                            fmt='%.3f',addmetad=True)
        det.save_txtnda(dirname+'pixel_status/'+fname, status, fmt='%d', addmetad=True)
    

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

        if  det.dettype == DetObjectClass.CsPad:
            mask=mask.reshape(2,185*388).transpose(1,0)
        elif det.dettype == DetObjectClass.CsPad2M:
            mask=mask.reshape(32*185,388)
        elif det.dettype == DetObjectClass.Epix:
            mask=mask.reshape(704,768)

        #2x2 save as 71780 lines, 2 entries
        #cspad save as 5920 lines, 388 entries
        if raw_input("Save to calibdir?\n") in ["y","Y"]:
            mask = (~(mask.astype(bool))).astype(int)
            srcStr=det.source.__str__().replace('Source("DetInfo(','').replace(')")','')
            calibDir=det.env.calibDir()
            if det.dettype==DetObjectClass.CsPad:
                dirname='/%s/CsPad2x2::CalibV1/%s/pixel_mask/'%(calibDir,srcStr)
            elif det.dettype==DetObjectClass.CsPad2M:
                dirname='%s/CsPad::CalibV1/%s/pixel_mask/'%(calibDir,srcStr)        
            elif det.dettype==DetObjectClass.Epix:
                dirname='%s/Epix100a::CalibV1/%s/pixel_mask/'%(calibDir,srcStr)
            elif det.dettype==DetObjectClass.Rayonix:
                dirname='%s/Camera::CalibV1/%s/pixel_mask/'%(calibDir,srcStr)        
            elif det.dettype==DetObjectClass.Epix10k:
                dirname='%s/Epix10ka::CalibV1/%s/pixel_mask/'%(calibDir,srcStr)
            elif det.dettype==DetObjectClass.Epix10k2M:
                dirname='%s/Epix10ka2M::CalibV1/%s/pixel_mask/'%(calibDir,srcStr)

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
        if det.dettype==DetObjectClass.Jungfrau or det.dettype==DetObjectClass.Icarus:
            for ped,frms in zip(pedestals, rms):
                hstPed.append(np.histogram(ped.flatten(), np.arange(0, pedestals.max()*1.05)))
                hstRms.append(np.histogram(frms.flatten(), np.arange(0, rms.max()*1.05,0.05)))
        else:
            hstPed.append(np.histogram(pedestals.flatten(), np.arange(0, pedestals.max()*1.05)))
            hstRms.append(np.histogram(rms.flatten(), np.arange(0, rms.max()*1.05,0.05)))

        if needsGeo:
            if det.dettype==DetObjectClass.Jungfrau:
                pedImg = det.image(self.run,pedestals[0])
                rmsImg = det.image(self.run,rms[0])
                if pedImg is None and pedestals.shape[1]==1:
                    pedImg = pedestals[0].squeeze()
                    rmsImg = rms[0].squeeze()
            elif det.dettype==DetObjectClass.Icarus:
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
        if det.dettype==DetObjectClass.Jungfrau:
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
        calibDir=det.env.calibDir()
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
        fnames = os.listdir('/%s/%s/%s/pedestals/'%(calibDir, detTypeStr, detNameStr))
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
                print(
                    'run %d, pixel cold/hot, low noise, high noise: %d / %d / %d / %d pixels'\
                      %(pedRun, (allPeds[-1]<printVal[0]).sum(), (allPeds[-1]>printVal[1]).sum(),
                        (allRms[-1]<printVal[2]).sum(), (allRms[-1]>printVal[3]).sum())
                )
            elif len(printVal)>=2:
                print('run %d, pixel cold/hot: %d / %d pixels'%(
                    pedRun,(allPeds[-1]<printVal[0]).sum(), (allPeds[-1]>printVal[1]).sum()
                ))

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
            boundsPed = (
                0, np.array(self.calibhisto['runs']).min(),  allPeds[icurrRun].max(), np.array(self.calibhisto['runs']).max()
            )
            pedZ = (
                np.nanpercentile(self.calibhisto['ped2d'],1),np.nanpercentile(self.calibhisto['ped2d'],99.5)
            )
            boundsRms = (
                0, np.array(self.calibhisto['runs']).min(),  allRms[icurrRun].max(), np.array(self.calibhisto['runs']).max()
            )
            rmsZ = (
                np.nanpercentile(self.calibhisto['rms2d'],1),np.nanpercentile(self.calibhisto['rms2d'],99.5)
            )
            runDim=hv.Dimension('run')

            #replace norm plots by zoomed plots.
            if showPlot is None:
                return (
                    hv.Image(self.calibhisto['ped2d'], 
                             kdims=[hv.Dimension('pedestal'),runDim], 
                             bounds=boundsPed, 
                             vdims=[hv.Dimension('pedVals', range=pedZ)]) \
                        .options(cmap='viridis',colorbar=True, width=plotwidth) + \
                    hv.Image(self.calibhisto['ped2d'], 
                             kdims=[hv.Dimension('pedestal'),runDim], 
                             bounds=boundsPed, 
                             vdims=[hv.Dimension('pedVals', range=pedZ)]) \
                        .options(cmap='viridis',colorbar=True, width=plotwidth) + \
                    hv.Image(self.calibhisto['rms2d'], 
                             kdims=[hv.Dimension('rms'),runDim], 
                             bounds=boundsRms, 
                             vdims=[hv.Dimension('rmsVals', range=rmsZ)]) \
                        .options(cmap='viridis',colorbar=True, width=plotwidth) + \
                    hv.Image(self.calibhisto['rms2dNorm'], 
                             kdims=[hv.Dimension('rms/ref'),runDim], 
                             bounds=boundsRms, 
                             vdims=[hv.Dimension('rmsValsNorm', range=rmsZNorm)]) \
                        .options(cmap='viridis',colorbar=True, width=plotwidth) 
                ).cols(2)
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
        im0 = plt.subplot(gsPed[2]).imshow(
            ch_ped2d.T,
            vmin=np.nanpercentile(ch_ped2d,1), vmax=np.nanpercentile(ch_ped2d,99),
            interpolation='none',aspect='auto',
            extent=[0, ch_ped2d.shape[0], self.calibhisto['pedBinLow'],
                    self.calibhisto['pedBinLow']+ch_ped2d.shape[1]*self.calibhisto['pedBinWidth']], 
            origin='lower')
        cbar0 = plt.colorbar(im0)

        ch_rms2d = self.calibhisto['rms2d']
        im2 = plt.subplot(gsPed[3]).imshow(
            ch_rms2d.T,
            vmin=np.nanpercentile(ch_rms2d,1), vmax=np.nanpercentile(ch_rms2d,99),
            interpolation='none',
            aspect='auto',
            extent=[0, ch_rms2d.shape[0], self.calibhisto['rmsBinLow'],
                    self.calibhisto['rmsBinLow']+ch_rms2d.shape[1]*self.calibhisto['rmsBinWidth']],
            origin='lower')
        cbar2 = plt.colorbar(im2)


    def compareCommonMode(self, detname='None',common_modes=[], numEvts=100, thresADU=0., thresRms=0., plotWith=None):
        if detname == 'None':
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
            dirname = self.sda.dirname+'/cube'

        myCube, cubeName_onoff = self.sda.prepCubeData(cubeName)
        
        print('\nVariables to be read from xtc: ',myCube.targetVarsXtc)

        if isinstance(nEvtsPerBin, str): nEvtsPerBin=int(nEvtsPerBin)

        #smallData only version: only run on single core
        if len(myCube.targetVarsXtc)<=0:
            if rank==0:
                outFileName = '/Cube_'+self.sda.fname.split('/')[-1].replace('.h5',f'_{cubeName}.h5')
                outFileName = dirname + outFileName
                if (onoff==0):
                    outFileName=outFileName.replace('.h5','_off.h5')
                elif (onoff==1):
                    outFileName=outFileName.replace('.h5','_on.h5')
                fout = h5py.File(outFileName, "w")
                #ADD CNF STUFF HERE!
                printR(rank, 'No big data, bin the data now....be patient')
                cubeData = self.sda.makeCubeData(cubeName,onoff=onoff)
                printR(rank, f'Now write outputfile (only small data) to : {outFileName}')
                for key in cubeData.variables:
                    addToHdf5(fout, key, cubeData[key].values)

                for cfgVar in myCube.targetVarsCfg:
                    addToHdf5(fout, cfgVar.replace('/','_'), self.sda.getVar(cfgVar))
                    print('add cfgVar to hdf5', cfgVar.replace('/','_'))
                fout.close()
                # stop worker waiting to process area dets
                for worker_id in range(size-1):
                    comm.send('done', dest=worker_id+1)
                bins, nEntries = cubeData.binVar_bins.values, cubeData.nEntries.values
            return cubeName, bins, nEntries

        printR(rank,'Now make big cube')
        #only run on rank=0 & broadcast.
        outFileName = '/Cube_'+self.sda.fname.split('/')[-1].replace('.h5',f'_{cubeName}.h5.inprogress')
        outFileName = dirname + outFileName
        if (onoff==0):
            outFileName = outFileName.replace('.h5','_off.h5')
        elif (onoff==1):
            outFileName = outFileName.replace('.h5','_on.h5')
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
        print('\n\n****** Total number of bins: {}'.format(nbins))
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
                    if det.det.dettype==DetObjectClass.Jungfrau:
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
        # save_fct = lambda data=None, bin_idx=None: self.save_bin_to_h5(fout=fout, data=data, bin_idx=bin_idx)
        # sum_data = mpi_fun.bin_distribution(bins_info, func=save_fct)
        bin_distrib = mpi_fun.Bin_distribution(bins_info, fout)
        bin_distrib.distribute()
        t3 = time.time()

        print("\n***** ALL BINS DONE AFTER {:0.2f} min. *****\n".format((t3-t0)/60))

        print(f'Renaming file from {outFileName} to {outFileName.replace(".inprogress","")}, remove random variable if applicable\n\n')
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
            elif k.src().__repr__()!='':
                detInData.append(k.src().__repr__())
        self.detNames=[]
        targetVarsXtc=[]
        for idet,det in enumerate(myCube.targetVarsXtc):
            if isinstance(det, dict):
                dName = det['source']
            else:
                dName = det
            for d in detInData:
                if dName in d:
                    self.detNames.append(dName)
                    targetVarsXtc.append(det)
                    break
            if det not in targetVarsXtc:
                print(f'Could not find detector {dName} in data.')
        myCube.targetVarsXtc = [ {'source':det, 'full':1} if isinstance(det, str) else det for det in targetVarsXtc]
        for det in myCube.targetVarsXtc:
            self.addDetInfo(det)
        to_worker = myCube.targetVarsXtc
        comm.bcast(to_worker, root=0)
        return
        
    @staticmethod
    def make_det_data_dset(fout, detname, det_shape, nbins):
        """ Is that really necessary? Yes, because of empty bins! """
        # data dset
        dset_name = '{}_data'.format(detname)
        shape = tuple(np.r_[nbins,det_shape])
        dset = fout.create_dataset(dset_name, shape, dtype=float)
        # n_in_bin dset
        dset_name = '{}_nEntries'.format(detname)
        dset = fout.create_dataset(dset_name, (nbins,), dtype=float)
        return
    
    @staticmethod
    def print_me(data=None, bin_idx=None):
        print(data)
        return
