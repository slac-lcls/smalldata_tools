from os import path
import numpy as np
try:
    from pylab import ginput
except:
    print 'could not import ginput'
    pass
try:
    from matplotlib import pyplot as plt
except:
    print 'could not import pyplot'
    pass
from matplotlib import gridspec
from matplotlib import path
import itertools
import os
import psana
import SmallDataAna as sda
from DetObject import DetObject
from SmallDataUtils import getUserData
from utilities import printR
from utilities import addToHdf5
from utilities import dropObject
from utilities import rename_reduceRandomVar
from utilities_plotting import plotImageBokeh
import bokeh.plotting as bp
import azimuthalBinning as ab
from matplotlib import gridspec
import RegDB.experiment_info
from utilities_FitCenter import FindFitCenter
from utilities_FitCenter import fitCircle
from mpi4py import MPI
import h5py
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

class SmallDataAna_psana(object):
    def __init__(self, expname='', run=-1,dirname='', filename=''):
        self.run=run
        self.expname=expname
        self.hutch=expname[0:3]
        currExpname = RegDB.experiment_info.active_experiment(self.hutch.upper())[1]
        if expname==currExpname:
            lastRun = RegDB.experiment_info.experiment_runs(self.hutch.upper())[-1]['num']
            if self.run > lastRun:
                printR(rank, 'experiment %s does only have %d runs, requested %d'%(expname, lastRun, run))
                return None
            
            isLive = (RegDB.experiment_info.experiment_runs(self.hutch.upper())[-1]['end_time_unix'] is None)
            self.dsname='exp=%s:run=%i:smd'%(expname,run)
            xtcdirname = '/reg/d/psdm/%s/%s/xtc'%(self.hutch.lower(),expname)
            idxname=xtcdirname+'/index/e*-r%04d-*'%int(run)
            ffbxtcdirname = '/reg/d/ffb/%s/%s/xtc'%(self.hutch.lower(),expname)
            ffbxtcname=ffbxtcdirname+'/e*-r%04d-*'%int(run)
            ffbidxname=ffbxtcdirname+'/index/e*-r%04d-*'%int(run)
            import glob
            present_ffbXtc=glob.glob('%s'%ffbxtcname)
            present_ffbIdx=glob.glob('%s'%ffbidxname)
            present_Idx=glob.glob('%s'%idxname)
            if len(present_ffbXtc)>0:
                self.dsname='exp=%s:run=%i:smd:dir=/reg/d/ffb/%s/%s/xtc'%(expname,run,self.hutch.lower(),expname)
                if isLive:
                    self.dsname=self.dsname+':live'
            if len(present_ffbIdx)>0:
                self.dsnameIdx='exp=%s:run=%i:idx:dir=/reg/d/ffb/%s/%s/xtc'%(expname,run,self.hutch.lower(),expname)
            elif len(present_Idx)>0:
                self.dsnameIdx='exp=%s:run=%i:idx'%(expname,run)
            else:
                self.dsnameIdx=None

        else:
            xtcdirname = '/reg/d/psdm/%s/%s/xtc/'%(self.hutch, expname)
            haveXtc=False
            dirpath, dirnames, filenames = next(os.walk(xtcdirname))
            for fname in filenames:
                if fname.find('r%04d'%run)>=0 and fname[-3:]=='xtc':
                    haveXtc=True
            if not haveXtc:
                printR(rank, 'Could not find xtc files for SmallDataAna_psana for exp %s and run %03d, return None'%(expname, run))
                return None

            self.dsname='exp=%s:run=%i:smd'%(expname,run)
            self.dsnameIdx='exp=%s:run=%i:idx'%(expname,run)
            xtcdirname = '/reg/d/psdm/%s/%s/xtc'%(self.hutch.lower(),expname)
            idxname=xtcdirname+'/index/e*-r%04d-*'%int(run)
            import glob
            present_Idx=glob.glob('%s'%idxname)
            if len(present_Idx)>0:
                self.dsnameIdx='exp=%s:run=%i:idx'%(expname,run)
            else:
                self.dsnameIdx=None

        printR(rank, 'make SmallDataAna_psana from dsname: %s'%self.dsname)
        try:
            self.ds = psana.DataSource(self.dsname)
        except:
            printR(rank, 'Failed to set up small data psana dataset!')
            self.ds = None
        if self.dsnameIdx is None:
            printR(rank, 'Failed to set up index based psana dataset, likely because no idx files have been produced/moved yet')
        else:
            try:
                self.dsIdx = psana.DataSource(self.dsnameIdx)
                self.dsIdxRun = self.dsIdx.runs().next()
            except:
                printR(rank, 'Failed to set up index based psana dataset')
                self.dsIdx = None
                self.dsIdxRun = None
        printR(rank, 'try to make SmallDataAna using dirname %s, for exp %s and run %s'%(dirname,expname,run))
        try:
            printR(rank, 'setting up SmallData ana from anaps ')
            self.sda = sda.SmallDataAna(expname,run,dirname=dirname, filename=filename)
        except:
            printR(rank, 'failed, set anaps.lda to None')
            self.sda = None
        self.jobsIds = []

    def commonModeStr(self, common_mode=0):
        if common_mode<0:
            return 'raw_'
        elif common_mode==5:
            return 'unb_'
        elif common_mode==6:
            return 'cm6_'
        elif common_mode==45:
            return 'median_'
        elif common_mode==46:
            return 'cm46_'
        elif common_mode==47:
            return 'cm47_'
        elif common_mode==71:
            return 'cm71_'
        elif common_mode==72:
            return 'cm72_'
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
                print keys
                return
            else:
                return keys
        except:
            printR(rank, 'failed to use last event in idx mode, try smlData')
            try:
                keys=self.ds.events().next().keys()
                if printthis: 
                    print keys
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
            print keys
        else:
            return keys

    def EpicsAliases(self,idx=True,printthis=False):
        if idx:
            keys=self.dsIdx.env().epicsStore().aliases()
        else:
            keys=self.ds.env().epicsStore().aliases()
        if printthis: 
            print keys
        else:
            return keys
        
    def plotVar(self, plotvar, numBins=[100], useFilter=False, limits=[1,99],fig=None,asHist=False):
        self.sda.plotVar(plotvar=plotvar, numBins=numBins, useFilter=useFilter, limits=limits,fig=fig,asHist=asHist)

    def plotScan(self, ttCorr=False, sig='diodeU/channels', sigROI=[], i0='ipm3/sum', numBins=100):
        self.sda.plotScan(ttCorr=ttCorr, sig=sig, sigROI=sigROI, i0=i0, numBins=numBins)

#
# these functions need psana as well, make separate class that imports SmallDataAna?
#
    def addDetInfo(self, detname='None', common_mode=0):
        if detname=='None':
            aliases = self._getDetName()
            if len(aliases)==1:
                detname = aliases[0]
            else:
                print 'detectors in event: \n',
                for alias in aliases:
                    print alias
                detname = raw_input("Select detector to get detector info for?:\n")
        printR(rank, 'try to make psana Detector with: %s'%detname)

        detnameDict=None
        if isinstance(detname, dict):
            detnameDict = detname
            if 'common_mode' in detnameDict.keys():
                common_mode=detnameDict['common_mode']
            detname = detnameDict['source']

        #only do this if information is not in object yet
        if detname in self.__dict__.keys() and self.__dict__[detname].common_mode==common_mode and detnameDict is not None:
            return detname

        if detname in self.__dict__.keys():
            printR(rank, 'redefine detector object with different common mode: %d instead of %d'%( common_mode,self.__dict__[detname].common_mode))
        det = DetObject(detname , self.dsIdx.env(), self.run, name=detname,common_mode=common_mode)
        self.__dict__[detname]=det
        if (detname+'_pedestals') in self.__dict__.keys():
            return detname
        self.__dict__[detname+'_pedestals'] = det.ped
        self.__dict__[detname+'_rms'] = det.rms
        self.__dict__[detname+'_x']=det.x
        self.__dict__[detname+'_y']=det.y
        self.__dict__[detname+'_iX']=det.iX
        self.__dict__[detname+'_iY']=det.iY

        if detnameDict is not None:
            for key in detnameDict.keys():
                if key=='full' or key=='Full':
                    self.__dict__[detname].saveFull()
                if key.find('ROI')==0:
                    self.__dict__[detname].addROI(key,detnameDict[key],writeArea=True)
                if key.find('photon')==0:
                    thres=0.9
                    if key.find('_')>0:
                        try:
                            thres=float(key.split('_')[1].replace('p','.'))
                        except:
                            pass
                    self.__dict__[detname].addPhotons(ADU_per_photon=detnameDict[key], thresADU=thres, retImg=2)

        return detname

    def _getDetName(self, detname=None):
        detNameList = ['cs','Cs','epix','Epix','opal','Opal','zyla','Zyla','jungfrau','Jungfrau','gige','Camera']
        #look for detector
        aliases=[]
        for key in self.Keys():
            for detName in detNameList:
                if key.alias().find(detName)>=0:
                    aliases.append(key.alias())
        if len(aliases)<1:
            for key in self.Keys():
                for detName in detNameList:
                    if key.src().__repr__().find(detName)>=0:
                        aliases.append(key.src().__repr__())
        if len(aliases)==1:
            return aliases[0]
        elif detname is not None:
            alias_match = [alias for alias in aliases if alias.find(detname)>=0]
            return alias_match
        else:
            return aliases

    def AvImage(self, detname='None', numEvts=100, thresADU=0., thresRms=0., useFilter='', nSkip=0,minIpm=-1., common_mode=0, std=False, median=False, printFid=False,useMask=True):
        if not isinstance(detname, basestring):
            print 'please give parameter name unless specifying arguments in right order. detname is first'
            return
        #look for detector
        if detname=='None':
            detname = self._getDetName()
            if isinstance(detname, list):
                print 'detectors in event: \n',
                for alias in detname:
                    print alias
                detname = raw_input("Select detector to get detector info for?:\n")

        if not detname in self.__dict__.keys() or self.__dict__[detname].common_mode!=common_mode:
            self.addDetInfo(detname=detname, common_mode=common_mode)
        det=self.__dict__[detname]
        rms = self.__dict__[detname+'_rms']
        pedestals = self.__dict__[detname+'_pedestals']
        iX = self.__dict__[detname+'_iX']
        iY = self.__dict__[detname+'_iY']

        if detname.find('opal')>=0:
            common_mode = -1
            if pedestals[0,0]==-1:
                print 'time tool opal, image was not saved. Should ideally exclude from detector list'
                return
        else:
            print 'done setting up the geometry'

        #now get the non-image data
        imgAr = []
        run = self.dsIdxRun
        times=[]
        if self.sda is None or 'fh5' not in self.sda.__dict__.keys():
            useFilter=''
        if useFilter!='':
            evttsSel = self.sda.getSelIdx(useFilter)
            print 'using ldat base selection, have %s events for selection %s'%(len(evttsSel),useFilter)
            if numEvts==-1:
                numEvts = len(evttsSel)-nSkip
                if numEvts<0:
                    print 'have no events, quit'
                    return
            for evtts in evttsSel[nSkip:min(nSkip+numEvts, len(evttsSel))]:
                times.append(psana.EventTime(evtts[1],evtts[0]))
        else:
            times = run.times()[nSkip:]
        print 'requested ',numEvts,' used ',min(len(times),numEvts), ' now actually get events'
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
                print 'Could not get this event, skip '
                continue
            if evt is None:
                print 'Returned event is None, skip'
                continue
            if minIpm!=-1 and ( (self.hutch=='xpp' and evt.get(psana.Lusi.IpmFexV1, psana.Source('BldInfo(XppSb2_Ipm)')).sum() < minIpm) or (self.hutch=='xcs' and evt.get(psana.Lusi.IpmFexV1, psana.Source('BldInfo(XCS-IPM-05)')).sum() < minIpm)):
                continue
            if printFid:
                print (evt.get(psana.EventId)).fiducials()
            aliases = [ k.alias() for k in evt.keys() ]
            if not detname in aliases:
                continue

            det.evt = dropObject()
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

        #make array
        data='AvImg_';
        if useFilter:
            data+='Filter'+useFilter+'_'
        if thresADU!=0:
            data+='thresADU%d_'%int(thresADU)
        if thresRms!=0:
            data+='thresRms%d_'%int(thresRms*10.)

        print 'use common mode: ',common_mode
        data+=self.commonModeStr(common_mode)
        data+=detname

        img = np.array(imgAr)
        if thresADU >= 0 and thresRms >=0:
            imgA = img.mean(axis=0)#.squeeze()
        else:
            imgA = img.sum(axis=0)#.squeeze()
        if det.mask is not None and useMask:
            imgA[det.mask==0]=0
                
        self.__dict__[data]=imgA
        if std:
            imgS = img.std(axis=0)#.squeeze()
            self.__dict__[data.replace('AvImg_','AvImg_std_')]=imgS
        if median:
            imgM = np.median(img,axis=0)#.squeeze()
            self.__dict__[data.replace('AvImg_','AvImg_median_')]=imgM

    def getAvImage(self,detname=None, imgName=None):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('_mask_')>=0:
                continue
            if key.find('_azint_')>=0:
                continue
            if imgName is not None and key.find(imgName)<0:
                continue
            if key.find('AvImg')==0:
                if detname is not None and key.find(detname)>=0:
                    avImages.append(key)
                elif detname is None:
                    avImages.append(key)
        if len(avImages)==0:
            print 'please create the AvImage first!'
            return
        elif len(avImages)>1:
            print 'we have the following options: ',avImages
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
            print 'key ',key
            if key.find('azint')>=0:
                if detname is not None and key.find(detname)>=0:
                    azInts.append(key)
                elif detname is None:
                    azInts.append(key)
        if len(azInts)==0:
            print 'please create an azimuthal integral first!'
            return
        elif len(azInts)>1:
            print 'we have the following options: ',azInts
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
        dns = avimage.replace('AvImg_','').replace('std_','').replace('median_','').replace('pedSub_','').replace('cm46_','').replace('cm72_','').replace('cm71_','').replace('raw_','').replace('unb_','').split('_')
        for ddns in dns:
            if ddns.find('thres')<0 and ddns.find('Filter')<0:
                detname+=ddns;detname+='_'
        if detname[-1]=='_':
            detname = detname[:-1]
        return detname

    def plotAvImage(self,detname=None, use_mask=False, ROI=[], limits=[5,99.5], returnIt=False, plotWith='',debugPlot=-1):
        detname, img, avImage = self.getAvImage(detname=None)

        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
            img = (img*mask)

        plotMax = np.percentile(img, limits[1])
        plotMin = np.percentile(img, limits[0])
        print 'plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax)

        if len(img.shape)>2:
            image = self.__dict__[detname].det.image(self.run, img)
            if image is None:
                image = img.squeeze()
        else:
            image = img

        #have issues with Jupyter notebook io size.
        if plotWith=='bokeh_notebook':
            plot_title="%s in Run %d"%(avImage, self.run)
            if debugPlot>0:
                img = image[:debugPlot, :debugPlot]
            else:
                img = image
            layout, p, im = plotImageBokeh(img, plot_title=plot_title, plotMaxP=np.percentile(img,99),plotMinP=np.percentile(img,1), plotWidth=700, plotHeight=700)
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
        print 'now save information to file: ',fname
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
            print 'please create the AvImage first!'
            return
        elif len(avImages)>1:
            print 'we have the following options: ',avImages
            avImage=raw_input('type the name of the AvImage to use:')
        else:
            avImage=avImages[0]
        detname = self._getDetName_from_AvImage(avImage)
        img = self.__dict__[avImage]

        plotMax = np.percentile(img, limits[1])
        plotMin = np.percentile(img, limits[0])
        print 'plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax)

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

        iX = self.__dict__[detname+'_iX']
        iY = self.__dict__[detname+'_iY']

        happy = False
        while not happy:
            p =np.array(ginput(2))
            p_axis1=[int(p[:,1].min()),int(p[:,1].max())+1]
            p_axis0=[int(p[:,0].min()),int(p[:,0].max())+1]
            print 'points:',p
            plt.subplot(gs[1]).imshow(image[p_axis1[0]:p_axis1[1],p_axis0[0]:p_axis0[1]],clim=[plotMin,plotMax],interpolation='None')
            plt.pause(0.0001)
            if needsGeo:
                mask_roi=np.zeros_like(image)
                mask_roi[p_axis1[0]:p_axis1[1],p_axis0[0]:p_axis0[1]]=1
                mask_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(iX,iY)] )
            if raw_input("Happy with this selection:\n") in ["y","Y"]:
                happy = True
                    
        if not needsGeo:
            print 'ROI: [[%i,%i], [%i,%i]]'%(p[:,1].min(),p[:,1].max(),p[:,0].min(),p[:,0].max())
            return

        if len(mask_nda.shape)>2:
            for itile,tile in enumerate(mask_nda):
                if tile.sum()>0:
                    ax0 = np.arange(0,tile.sum(axis=0).shape[0])[tile.sum(axis=0)>0]
                    ax1 = np.arange(0,tile.sum(axis=1).shape[0])[tile.sum(axis=1)>0]
                    print 'ROI: [[%i,%i], [%i,%i], [%i,%i]]'%(itile,itile+1,ax1.min(),ax1.max(),ax0.min(),ax0.max())
                    fig=plt.figure(figsize=(6,6))
                    plt.imshow(img[itile,ax1.min():ax1.max(),ax0.min():ax0.max()],interpolation='none')
                    plt.pause(0.0001)
                    plt.show()
        else:
            tile=mask_nda
            if tile.sum()>0:
                ax0 = np.arange(0,tile.sum(axis=0).shape[0])[tile.sum(axis=0)>0]
                ax1 = np.arange(0,tile.sum(axis=1).shape[0])[tile.sum(axis=1)>0]
                print 'ROI: [[%i,%i], [%i,%i]]'%(ax1.min(),ax1.max(),ax0.min(),ax0.max())
            

    def FitCircleAuto(self, detname=None, plotRes=True, forceMask=False, inParams={}):
        detname, img, avImage = self.getAvImage(detname=None)
        try:
            mask = self.__dict__[detname].det.cmask.astype(bool)
        except:
            mask = self.__dict__[detname].det.mask(self.run, calib=True, status=True).astype(bool)
        nPixRaw=1.
        for idim in mask.shape:
            nPixRaw=nPixRaw*idim
        print 'check calib mask: ',mask.sum(),nPixRaw,(mask.sum()/nPixRaw)>0.5,(mask.sum().astype(float)/nPixRaw)
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
        plt.imshow(image, interpolation='none', cmap='gray',clim=[0,np.percentile(img.flatten(),99.5)])
        #plot sparse images in blue.
        plt.plot(arSparse.col, arSparse.row,markersize=5,color='#ff9900',marker='.',linestyle='None')
        if combRes==-1:
            return -1
        greens=['#666600','#669900','#66cc00','#66cc99','#6699cc','#6633ff']
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
        print 'center ',combRes['xCen'],combRes['yCen']
        if needsGeo:
            geo = self.__dict__[detname].det.geometry(self.run)
            d=160000.
            ix0,iy0 = geo.point_coord_indexes(p_um=(0,0))
            ixd,iyd = geo.point_coord_indexes(p_um=(d,d))
            combRes['yCen'] = d*(combRes['yCen']-ix0)/(ixd-ix0)
            combRes['xCen'] = d*(combRes['xCen']-iy0)/(iyd-iy0)
            print 'center Final mid',combRes['xCen'],combRes['yCen']
            helpVar =  combRes['xCen']
            combRes['xCen'] = combRes['yCen']
            combRes['yCen'] = helpVar
            for r in combRes['R']:
                print 'aradii: ',r,d*r/(ixd-ix0)
        else:
            for r in combRes['R']:
                print 'aradii: ',r
        print 'center Final ',combRes['xCen'],combRes['yCen']
        return combRes
 
    def FitCircleMouse(self, detname=None, use_mask=False, limits=[5,99.5]):
        self.FitCircle(detname=None, use_mouse=True, use_mask=False, use_mask_local=False, limits=[5,99.5])

    def FitCircleThreshold(self, detname=None, use_mask=False, limits=[5,99.5],plotIt=False, thresIn=None):
        self.FitCircle(detname=None, use_mouse=False, use_mask=False, use_mask_local=True, limits=[5,99.5], thresIn=thresIn, plotIt=plotIt)

    def FitCircle(self, detname=None, use_mouse=None, use_mask=False, use_mask_local=False, limits=[5,99.5], plotIt=True, thresIn=None):
        detname, img, avImage = self.getAvImage(detname=None)

        if use_mouse is None or (not isinstance(use_mouse, bool) and not isinstance(use_mouse, basestring)):
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

        plotMax = np.percentile(img[img!=0], 99.5)
        plotMin = np.percentile(img[img!=0], 5)
        printR(rank, 'plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        needsGeo=False
        if self.__dict__[detname].ped is not None and self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if needsGeo:
            image = self.__dict__[detname].det.image(self.run, img)
        else:
            image = img

        x = self.__dict__[detname+'_x']
        y = self.__dict__[detname+'_y']
        if x is None:
            x = np.arange(0, image.shape[1])
            y = np.arange(0, image.shape[0])
            x,y = np.meshgrid(x,y)
        extent=[x.min(), x.max(), y.min(), y.max()]

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
                print 'x,y: ',res['xCen'],res['yCen'],' R ',res['R']
                if raw_input("Happy with this selection:\n") in ["y","Y"]:
                    happy = True

        else:
            while not happy:
                if plotIt:
                    fig=plt.figure(figsize=(10,10))
                    plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None')
                if thresIn is None:
                    thres = float(raw_input("min percentile % of selected points:\n"))
                else:
                    thres = thresIn
                    happy = True
                thresP = np.percentile(img[img!=0], thres)
                print 'thresP',thresP
                imageThres=image.copy()
                imageThres[image>thresP]=1
                imageThres[image<thresP]=0
                if plotIt:
                    fig=plt.figure(figsize=(5,5))
                    plt.imshow(imageThres,clim=[-0.1,1.1])
                    if raw_input("Happy with this threshold (y/n):\n") in ["y","Y"]:
                        happy=True

            res = fitCircle(x.flatten()[img.flatten()>thresP],y.flatten()[img.flatten()>thresP])
            circleM = plt.Circle((res['xCen'],res['yCen']),res['R'],color='b',fill=False)
            print 'x,y: ',res['xCen'],res['yCen'],' R ',res['R']
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

    def MakeMask(self, detname=None, limits=[5,99.5]):
        detname, img, avImage = self.getAvImage(detname=None)

        plotMax = np.percentile(img, limits[1])
        plotMin = np.percentile(img, limits[0])
        print 'plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax)

        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if needsGeo:
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
        iX = self.__dict__[detname+'_iX']
        iY = self.__dict__[detname+'_iY']
        extent=[x.min(), x.max(), y.min(), y.max()]
        #print 'DEBUG: extent(x,y min,max)',extent
        
        mask=[]
        mask_r_nda=None
        select=True
        while select:
            fig=plt.figure(figsize=(12,10))
            gs=gridspec.GridSpec(1,2,width_ratios=[2,1])
            #print "rectangle(r-click, R-enter), circle(c), polygon(p), dark(d), noise(n) or edgepixels(e)?:"
            plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')

            shape = raw_input("rectangle(r-click, R-enter), circle(c), polygon(p), dark(d), noise(n) or edgepixels(e)?:\n")
            #shape = raw_input()
            #this definitely works for the rayonix...
            if shape=='r':
                print 'select two corners: '
                p =np.array(ginput(2))
                mask_roi=np.zeros_like(image)
                mask_roi[p[:,1].min():p[:,1].max(),p[:,0].min():p[:,0].max()]=1
                if needsGeo:
                    mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(iX,iY)] )
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    mask_r_nda = mask_roi
                    plt.subplot(gs[1]).imshow(mask_r_nda)
                print 'mask from rectangle (shape):',mask_r_nda.shape
            elif shape=='R':
                print 'coordinates to select: '
                htot = raw_input("horizontal,vertical h1 h2 v1 v2 ?\n")
                h = htot.split(' ');h1=float(h[0]);h2=float(h[1]);v1=float(h[2]);v2=float(h[3]);
                mask_roi=np.zeros_like(image)
                mask_roi[h1:h2,v1:v2]=1
                if needsGeo:
                    mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(iX,iY)] )
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    mask_r_nda = mask_roi
                    plt.subplot(gs[1]).imshow(mask_r_nda)
                print 'mask from rectangle (shape):',mask_r_nda.shape
            elif shape=='c':
                fig=plt.figure(figsize=(12,10))
                gs=gridspec.GridSpec(1,2,width_ratios=[2,1])
                plt.subplot(gs[0]).imshow(np.rot90(image),clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
                if raw_input("Select center by mouse?\n") in ["y","Y"]:
                    c=ginput(1)
                    cx=c[0][0];cy=c[0][1]
                    print 'Corrdinates of selected center: ',cx,' ',cy
                else:
                    ctot = raw_input("center (x y)?\n")
                    c = ctot.split(' ');cx=float(c[0]);cy=float(c[1]);
                if raw_input("Select outer radius by mouse?\n") in ["y","Y"]: 
                    r=ginput(1)
                    rox=r[0][0];roy=r[0][1]
                    ro=np.sqrt((rox-cx)**2+(roy-cy)**2)
                    if raw_input("Select inner radius by mouse (for donut-shaped mask)?\n") in ["y","Y"]:
                        r=ginput(1)
                        rix=r[0][0];riy=r[0][1]
                        ri=np.sqrt((rix-cx)**2+(riy-cy)**2)
                    else:
                        ri=0
                    print 'radii: ',ro,' ',ri
                else:
                    rtot = raw_input("radii (r_outer r_inner)?\n")
                    r = rtot.split(' ');ro=float(r[0]);ri=max(0.,float(r[1]));        
                mask_router_nda = np.array( [(ix-cx)**2+(iy-cy)**2<ro**2 for ix, iy in zip(x,y)] )
                mask_rinner_nda = np.array( [(ix-cx)**2+(iy-cy)**2<ri**2 for ix, iy in zip(x,y)] )
                mask_r_nda = mask_router_nda&~mask_rinner_nda
                print 'mask from circle (shape):',mask_r_nda.shape
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
                else:
                    plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
                nPoints = int(raw_input("Number of Points (-1 until middle mouse click)?\n"))
                p=np.array(ginput(nPoints))
                print p
                mpath=path.Path(p)
                all_p = np.array([ (ix,iy) for ix,iy in zip(x.flatten(),y.flatten()) ] )
                mask_r_nda = np.array([mpath.contains_points(all_p)]).reshape(x.shape)
                if needsGeo:
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    plt.subplot(gs[1]).imshow(mask_r_nda)
                print 'mask from polygon (shape):',mask_r_nda.shape
                print 'not implemented yet....'
            elif shape=='d' or shape=='n':
                figDark=plt.figure(figsize=(12,10))
                gsPed=gridspec.GridSpec(1,2,width_ratios=[1,1])
                if shape=='d':
                    pedResult = det.pedestals(self.run)
                    if det.dettype==26:
                        pedResult = pedResult[0]
                    hstPed = np.histogram(pedResult.flatten(), np.arange(0, pedResult.max()*1.05))
                else:
                    pedResult = det.rms(self.run)
                    if det.dettype==26:
                        pedResult = pedResult[0]
                    hstPed = np.histogram(pedResult.flatten(), np.arange(0, pedResult.max()*1.05,0.05))
                if needsGeo:
                    pedResultImg = det.image(self.run,pedResult)
                else:
                    pedResultImg = pedResult.copy()
                plt.subplot(gsPed[0]).imshow(pedResultImg,clim=[np.percentile(pedResult,1),np.percentile(pedResult,99)])
                plt.subplot(gsPed[1]).plot(hstPed[1][:-1],np.log(hstPed[0]),'o')
                ctot=raw_input("Enter allowed pedestal range (min max)")
                c = ctot.split(' ');pedMin=float(c[0]);pedMax=float(c[1]);
                mask_r_nda=np.zeros_like(pedResult)
                print mask_r_nda.sum()
                mask_r_nda[pedResult<pedMin]=1
                mask_r_nda[pedResult>pedMax]=1
                print mask_r_nda.sum()
                mask_r_nda = (mask_r_nda.astype(bool)).astype(int)
                print mask_r_nda.sum()
                if needsGeo:
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda.astype(bool)))
                else:
                    plt.subplot(gs[1]).imshow(mask_r_nda.astype(bool))
            if shape=='e':
                ctot=raw_input("Enter number of edge rows that should be masked:")
                try:
                    nEdge = int(ctot)
                except:
                    print 'please enter an integer'
                    continue
                mask_r_nda = np.zeros_like(det.rms(self.run))
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

            if mask_r_nda is not None:
                print 'created a mask....',len(mask)
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
            #print 'DEBUG: ',mask_r_nda.shape, totmask_nm1.shape, x.shape
            print 'masked in this step: ',np.ones_like(x)[mask_r_nda.astype(bool)].sum()
            print 'masked up to this step: ',np.ones_like(x)[totmask_nm1].sum()
            print 'masked tot: ',np.ones_like(x)[totmask].sum()

            if len(mask)>1:
                fig=plt.figure(figsize=(15,9))
                gs2=gridspec.GridSpec(1,2,width_ratios=[1,1])
                plt.show()
                image_mask = img.copy(); image_mask[totmask]=0;
                #image_mask_nm1 = img.copy(); image_mask_nm1[totmask_nm1]=0;
                #print 'DEBUG: ',(img-image_mask_nm1).sum(), (img-image_mask).sum()
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
        print 'edited code....'
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
                dirname='/reg/d/psdm/%s/%s/calib/CsPad2x2::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)
            elif det.dettype==1:
                dirname='/reg/d/psdm/%s/%s/calib/CsPad::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            elif det.dettype==13:
                dirname='/reg/d/psdm/%s/%s/calib/Epix100a::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)
            elif det.dettype==19:
                dirname='/reg/d/psdm/%s/%s/calib/Camera::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            fname='%s-end.data'%self.run
            print 'save mask in %s as %s '%(dirname,fname)
            np.savetxt(dirname+fname,mask)
        elif raw_input("Save to local?\n") in ["y","Y"]:
            np.savetxt('Mask_%s_%s_Run%03d.data'%(avImage,self.sda.expname,int(self.run)),mask)
        return mask
         
    def makePedestal(self, detname, filterName='', numEvts=-1, pedRange=[10,10000], rmsRange=[2.,7.], i0Check='ipm', dirname='./'):
        if i0Check=='':
            i0List=[]
        elif i0Check=='ipm':
            if self.sda.expname[:3]=='xpp':
                i0List = ['ipm3/sum']
            elif self.sda.expname[:3]=='xcs':
                i0List = ['ipm5/sum']
            else:
                i0List=[]
        else:
            i0List = ['gas_detector/f_22_ENRC']
    
        minFrac = 1e6
        for i0 in i0List:
            try:
                i0Median = np.nanmedian(self.sda.getVar(i0))
                i0Off = np.nanmedian(self.sda.getVar(i0, filterName))
            except:
                print 'if not smallData file is available, try pass i0Check=\'\''
            if minFrac > i0Off/i0Median:
                minFrac=i0Off/i0Median
            print 'median value for ',i0,' is ',i0Median,' and for the off events ',i0Off,' offRatio: ',i0Off/i0Median
        if minFrac > 0.05 and i0List!=[]:
            print 'This selection seems to lets too many events with beam through, will quit'
            return

        self.AvImage(detname,numEvts=numEvts,useFilter=filterName, common_mode=-1, useMask=False, median=True)
        self.AvImage(detname,numEvts=numEvts,useFilter=filterName, common_mode=-1, useMask=False, std=True)
        fname='%s-end.data'%self.sda.run
        pedImg = self.__dict__['AvImg_median_Filter%s_raw_%s'%(filterName,detname)]
        pedStat = np.logical_and(pedImg > min(pedRange), pedImg < max(pedRange))
        rmsImg = self.__dict__['AvImg_std_Filter%s_raw_%s'%(filterName,detname)]
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

        if not os.path.exists(dirname+'pedestals'):
            os.makedirs(dirname+'pedestals')
        if not os.path.exists(dirname+'pixel_rms'):
            os.makedirs(dirname+'pixel_rms')
        if not os.path.exists(dirname+'pixel_status'):
            os.makedirs(dirname+'pixel_status')
        print 'save pedestal file in %s as %s '%(dirname+'pedestals/',fname)
        det.save_txtnda(dirname+'pedestals/'+fname,self.__dict__['AvImg_median_Filter%s_raw_%s'%(filterName,detname)], fmt='%.1f',addmetad=True)
        print 'save noise file in %s as %s '%(dirname+'pixel_rms/',fname)
        det.save_txtnda(dirname+'pixel_rms/'+fname,self.__dict__['AvImg_std_Filter%s_raw_%s'%(filterName,detname)], fmt='%.1f',addmetad=True)
        print 'save status file in %s as %s '%(dirname+'pixel_status/',fname)
        det.save_txtnda(dirname+'pixel_status/'+fname,status, fmt='%d',addmetad=True)

    def addAzInt(self, detname=None, phiBins=1, qBin=0.01, eBeam=9.5, center=None, dis_to_sam=None, name='azav', Pplane=1,userMask=None,tx=None,ty=None):
        detname, img, avImage = self.getAvImage(detname=detname)
        if dis_to_sam==None:
            dis_to_sam=float(raw_input('please enter the detector distance'))
        if center==None or len(center)!=2:
            centerString=raw_input('please enter the coordinates of the beam center as c1,c2 or [c1,c2]:')
            center=[int(centerString.replace('[','').replace(']','').split(',')[0]),
                    int(centerString.replace('[','').replace(']','').split(',')[1])]
        self.__dict__[detname].addAzAv(phiBins=phiBins, qBin=qBin, center=center, dis_to_sam=dis_to_sam, eBeam=eBeam, azavName=name, Pplane=Pplane, userMask=userMask,tx=tx,ty=ty)

    def getAzAvs(self,detname=None):
      if detname is None:
        detname, img, avImage = self.getAvImage(detname=None)   
        if detname is None:
          return
      azintArray = [ self.__dict__[detname][key] for key in self.__dict__[detname].__dict__.keys() if isinstance(self.__dict__[detname][key], ab.azimuthalBinning) ]
      azintNames = [ key for key in self.__dict__[detname].__dict__.keys() if isinstance(self.__dict__[detname][key], ab.azimuthalBinning) ]
      return azintNames, azintArray

    def AzInt(self, detname=None, use_mask=False, use_mask_local=False, plotIt=False, azintName=None, data=None, imgName=None):
        avImage=None
        if data is not None:
            if detname is None:
                detname=raw_input('type the name detector alias')
            img=data
        else:
            detname, img, avImage = self.getAvImage(detname=detname, imgName=None)
        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
            img = (img*mask)
        if use_mask_local:
            if avImage is None:
                detname_dump, img_dump, avImage = self.getAvImage(detname=detname)
            if self.__dict__.has_key('_mask_'+avImage):
                mask = self.__dict__['_mask_'+avImage]
                img = (img*mask)
            else:
                print 'no local mask defined for ',avImage

        azIntNames,azIntegrations = self.getAzAvs(detname)
        if len(azIntegrations)>1:
            if azintName is None:
                print 'we have the following options: ',azIntNames
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
                plt.plot(self.__dict__[detname].__dict__[azintName+'_q'],azintValues,'o')
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
                print 'we have the following options: ',azIntNames
                azintName=raw_input('type the name of the Azimuthal integral to use:')
            azIntegs = [ iiaz for iName, iiaz in zip(azIntNames,azIntegrations) if iName==azintName]
            azInteg = azIntegs[0]
            azint = azintName
        else:
            azInteg = azIntegrations[0]
            azint = azIntNames[0]
        if azint=='':
            print 'did not find azimuthal integral asked for'
            return
        else:
            if ('_azint_'+azint) not in self.__dict__.keys():
                print 'did not find azint ',azint,', all keys are: ',self.__dict__.keys()
        try:
            azintValues = self.__dict__['_azint_'+azint]
            fig=plt.figure(figsize=(8,5))
            if len(azintValues.shape)==1:
                print azintValues.shape, self.__dict__[detname].__dict__[azint+'_q'].shape
                plt.plot(self.__dict__[detname].__dict__[azint+'_q'],azintValues,'o')
            elif len(azintValues.shape)==2:
                plt.imshow(azintValues,aspect='auto',interpolation='none')
                plt.colorbar()
        except:
            pass

    def AzInt_centerVar(self, detname=None, use_mask=False, center=None, data=None, varCenter=110., zoom=[-1,-1]):
        if data is not None:
            if detname is None:
                detname=raw_input('type the name detector alias')
            img=data
        else:
            detname, img, avImage = self.getAvImage(detname=detname)
        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
            img = (img*mask)

        if center==None or len(center)!=2:
            centerString=raw_input('please enter the coordinates of the beam center as c1,c2 or [c1,c2]:')
            center=[int(centerString.replace('[','').replace(']','').split(',')[0]),
                    int(centerString.replace('[','').replace(']','').split(',')[1])]
        self.addAzInt(detname=detname, phiBins=13, qBin=0.001, eBeam=9.5, center=[center[0],center[1]], dis_to_sam=1000., name='c00')
        self.addAzInt(detname=detname, phiBins=13, qBin=0.001, eBeam=9.5, center=[center[0]-varCenter,center[1]], dis_to_sam=1000., name='cm10')
        self.addAzInt(detname=detname, phiBins=13, qBin=0.001, eBeam=9.5, center=[center[0]+varCenter,center[1]], dis_to_sam=1000., name='cp10')
        self.addAzInt(detname=detname, phiBins=13, qBin=0.001, eBeam=9.5, center=[center[0],center[1]-varCenter], dis_to_sam=1000., name='c0m1')
        self.addAzInt(detname=detname, phiBins=13, qBin=0.001, eBeam=9.5, center=[center[0],center[1]+varCenter], dis_to_sam=1000., name='c0p1')
        self.AzInt(detname=detname, use_mask=use_mask, data=data, azintName='c00')
        self.AzInt(detname=detname, use_mask=use_mask, data=data, azintName='cm10')
        self.AzInt(detname=detname, use_mask=use_mask, data=data, azintName='cp10')
        self.AzInt(detname=detname, use_mask=use_mask, data=data, azintName='c0m1')
        self.AzInt(detname=detname, use_mask=use_mask, data=data, azintName='c0p1')

        try:
            azintValues_c00 = self.__dict__['_azint_c00']
            azintValues_cm10 = self.__dict__['_azint_cm10']
            azintValues_cp10 = self.__dict__['_azint_cp10']
            azintValues_c0m1 = self.__dict__['_azint_c0m1']
            azintValues_c0p1 = self.__dict__['_azint_c0p1']

            fig=plt.figure(figsize=(10,6))
            from matplotlib import gridspec
            ymin=0;ymax=azintValues_c00.shape[1]-1
            #plt.subplot2grid((3,3),(1,1)).set_title('center= %9.2f, %9.2f'%(center[0],center[1]))
            #plt.subplot2grid((3,3),(1,0)).set_title('center= %9.2f, %9.2f'%(center[0]-varCenter,center[1]))
            #plt.subplot2grid((3,3),(1,2)).set_title('center= %9.2f, %9.2f'%(center[0]+varCenter,center[1]))
            #plt.subplot2grid((3,3),(0,1)).set_title('center= %9.2f, %9.2f'%(center[0],center[1]-varCenter))
            #plt.subplot2grid((3,3),(2,1)).set_title('center= %9.2f, %9.2f'%(center[0],center[1]+varCenter))
            plt.subplot2grid((3,3),(1,1)).imshow(azintValues_c00[:,ymin:ymax],aspect='auto',interpolation='none')
            plt.subplot2grid((3,3),(1,0)).imshow(azintValues_cm10[:,ymin:ymax],aspect='auto',interpolation='none')
            plt.subplot2grid((3,3),(1,2)).imshow(azintValues_cp10[:,ymin:ymax],aspect='auto',interpolation='none')
            plt.subplot2grid((3,3),(0,1)).imshow(azintValues_c0m1[:,ymin:ymax],aspect='auto',interpolation='none')
            plt.subplot2grid((3,3),(2,1)).imshow(azintValues_c0p1[:,ymin:ymax],aspect='auto',interpolation='none')

            print zoom, (zoom[0]!=zoom[1])
            while zoom[0]!=zoom[1]:
                if zoom[1]<zoom[0] or zoom[0]<0 or zoom[1]<0:
                    yString=raw_input('please enter x-boundaries of the zoomed figure as c1,c2 or [c1,c2]:')
                    ymin=int(yString.replace('[','').replace(']','').split(',')[0])
                    ymax=int(yString.replace('[','').replace(']','').split(',')[1])
                else:
                    ymin=zoom[0]
                    ymax=zoom[1]
                #plt.subplot2grid((3,3),(1,1)).set_title('center= %9.2f, %9.2f'%(center[0],center[1]))
                #plt.subplot2grid((3,3),(1,0)).set_title('center= %9.2f, %9.2f'%(center[0]-varCenter,center[1]))
                #plt.subplot2grid((3,3),(1,2)).set_title('center= %9.2f, %9.2f'%(center[0]+varCenter,center[1]))
                #plt.subplot2grid((3,3),(0,1)).set_title('center= %9.2f, %9.2f'%(center[0],center[1]-varCenter))
                #plt.subplot2grid((3,3),(2,1)).set_title('center= %9.2f, %9.2f'%(center[0],center[1]+varCenter))
                plt.subplot2grid((3,3),(1,1)).imshow(azintValues_c00[:,ymin:ymax],aspect='auto',interpolation='none')
                plt.subplot2grid((3,3),(1,0)).imshow(azintValues_cm10[:,ymin:ymax],aspect='auto',interpolation='none')
                plt.subplot2grid((3,3),(1,2)).imshow(azintValues_cp10[:,ymin:ymax],aspect='auto',interpolation='none')
                plt.subplot2grid((3,3),(0,1)).imshow(azintValues_c0m1[:,ymin:ymax],aspect='auto',interpolation='none')
                plt.subplot2grid((3,3),(2,1)).imshow(azintValues_c0p1[:,ymin:ymax],aspect='auto',interpolation='none')
                if raw_input("done? (y/n):\n") in ["y","Y"]:
                    zoom=[-1,-1]
                else:
                    yString=raw_input('please enter x-boundaries of the zoomed figure as c1,c2 or [c1,c2] - was: ',zoom)
                    ymin=int(yString.replace('[','').replace(']','').split(',')[0])
                    ymax=int(yString.replace('[','').replace(']','').split(',')[1])
        except:
            pass

    def SelectRegionDroplet(self, detname=None, limits=[5,99.5]):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('AvImg')==0:
                if detname is not None and key.find(detname)>=0:
                    avImages.append(key)
                elif detname is None:
                    avImages.append(key)
        if len(avImages)==0:
            print 'please create the AvImage first!'
            return
        elif len(avImages)>1:
            print 'we have the following options: ',avImages
            avImage=raw_input('type the name of the AvImage to use:')
        else:
            avImage=avImages[0]
        detname = self._getDetName_from_AvImage(avImage)
        img = self.__dict__[avImage]
        print img.shape

        plotMax = np.percentile(img, limits[1])
        plotMin = np.percentile(img, limits[0])
        print 'plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax)

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
        iX = self.__dict__[detname+'_iX']
        iY = self.__dict__[detname+'_iY']
        extent=[x.min(), x.max(), y.min(), y.max()]

        fig=plt.figure(figsize=(10,6))
        from matplotlib import gridspec
        gs=gridspec.GridSpec(1,2,width_ratios=[1,1])
        
        mask=None
        mask_r_nda=None

        plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')

        print 'select two corners: '
        p =np.array(ginput(2))
        mask_roi=np.zeros_like(image)
        mask_roi[p[:,1].min():p[:,1].max(),p[:,0].min():p[:,0].max()]=1
        if needsGeo:
            mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(iX,iY)] )
        else:
            mask_r_nda = mask_roi

        if mask_r_nda is not None:
            #print 'created a mask....'
            if mask is None:
                mask = mask_r_nda.astype(bool).copy()
            else:
                mask = np.logical_or(mask,mask_r_nda)
        print 'masked: ',np.ones_like(x)[mask.astype(bool)].sum()

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
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            fname='%s-end.data'%self.run
            det.save_txtnda(dirname+fname,mask, fmt='%d',addmetad=True)
        elif raw_input("Save to local?\n") in ["y","Y"]:
            mask = (~(mask.astype(bool))).astype(int)
            np.savetxt('%s_mask_run%s.data'%(self.sda.expname,self.run),mask)
        return mask
         
    def plotCalib(self, detname='None',common_mode=0):
        if not detname in self.__dict__.keys() or self.__dict__[detname].common_mode!=common_mode:
            detname = self.addDetInfo(detname=detname, common_mode=common_mode)
            if detname == 'None':
                print 'need detector name as input! '
                return
        det=self.__dict__[detname].det
        rms = self.__dict__[detname+'_rms']
        pedestals = self.__dict__[detname+'_pedestals']
        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        np.seterr(divide='ignore', invalid='ignore')
        if det.dettype==26:
            hstPed = np.histogram(pedestals[0].flatten(), np.arange(0, pedestals.max()*1.05))
            hstRms = np.histogram(rms[0].flatten(), np.arange(0, rms.max()*1.05,0.05))
            hstPed1 = np.histogram(pedestals[1].flatten(), np.arange(0, pedestals.max()*1.05))
            hstRms1 = np.histogram(rms[1].flatten(), np.arange(0, rms.max()*1.05,0.05))
            hstPed2 = np.histogram(pedestals[2].flatten(), np.arange(0, pedestals.max()*1.05))
            hstRms2 = np.histogram(rms[2].flatten(), np.arange(0, rms.max()*1.05,0.05))
        else:
            hstPed = np.histogram(pedestals.flatten(), np.arange(0, pedestals.max()*1.05))
            hstRms = np.histogram(rms.flatten(), np.arange(0, rms.max()*1.05,0.05))

        if needsGeo:
            if det.dettype==26:
                pedImg = det.image(self.run,pedestals[0])
                rmsImg = det.image(self.run,rms[0])
                if pedImg is None and pedestals.shape[1]==1:
                    pedImg = pedestals[0].squeeze()
                    rmsImg = rms[0].squeeze()
            else:
                pedImg = det.image(self.run,pedestals)
                rmsImg = det.image(self.run,rms)
        else:
            pedImg = pedestals
            rmsImg = rms

        figDark=plt.figure(figsize=(11,6))
        gsPed=gridspec.GridSpec(2,2,width_ratios=[1,1])
        plt.subplot(gsPed[1]).plot(hstPed[1][:-1],np.log(hstPed[0]),'o')
        plt.subplot(gsPed[3]).plot(hstRms[1][:-1],np.log(hstRms[0]),'o')
        if det.dettype==26:
            plt.subplot(gsPed[1]).plot(hstPed1[1][:-1],np.log(hstPed1[0]),'o')
            plt.subplot(gsPed[1]).plot(hstPed2[1][:-1],np.log(hstPed2[0]),'o')
            plt.subplot(gsPed[3]).plot(hstRms1[1][:-1],np.log(hstRms1[0]),'o')
            plt.subplot(gsPed[3]).plot(hstRms2[1][:-1],np.log(hstRms2[0]),'o')

            im0 = plt.subplot(gsPed[0]).imshow(pedImg,clim=[np.percentile(pedestals[0],1),np.percentile(pedestals[0],99)])
            cbar0 = plt.colorbar(im0)
            im2 = plt.subplot(gsPed[2]).imshow(rmsImg,clim=[np.percentile(rms[0],1),np.percentile(rms[0],99)])
        else:
            im0 = plt.subplot(gsPed[0]).imshow(pedImg,clim=[np.percentile(pedestals,1),np.percentile(pedestals,99)])
            cbar0 = plt.colorbar(im0)
            im2 = plt.subplot(gsPed[2]).imshow(rmsImg,clim=[np.percentile(rms,1),np.percentile(rms,99)])
        cbar2 = plt.colorbar(im2)

    def calibHisto(self, detname='None', common_mode=0, printVal=[-1]):
        if not detname in self.__dict__.keys() or self.__dict__[detname].common_mode!=common_mode:
            detname = self.addDetInfo(detname=detname, common_mode=common_mode)
            if detname == 'None':
                print 'need detector name as input! '
                return
        det=self.__dict__[detname].det
        rms = self.__dict__[detname+'_rms']
        pedestals = self.__dict__[detname+'_pedestals']
        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True
        #now look at directory and get filenames of darks. Hmm. extra function?
        detNameStr = det.name.__str__()
        if detNameStr.find('Epix')>=0:
            detTypeStr='Epix100a::CalibV1'
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

        currRun=0
        icurrRun=0
        allPeds=[]
        allRms=[]
        allPedsImg=[]
        allRmsImg=[]
        pedRuns.sort()
        for ipedRun,pedRun in enumerate(pedRuns):
            if pedRun <= self.run and self.run-pedRun < self.run-currRun:
                currRun = pedRun
                icurrRun = ipedRun
            allPeds.append(det.pedestals(pedRun))
            allRms.append(det.rms(pedRun))
            if needsGeo:
                allPedsImg.append(allPeds[-1])
                allRmsImg.append(allRms[-1])
            else:
                allPedsImg.append(allPeds[-1])
                allRmsImg.append(allRms[-1])
            if len(printVal)<2:
                print 'getting pedestal from run ',pedRun
            elif len(printVal)>=4:
                print 'run %d, pixel cold/hot, low noise, high noise: %d / %d / %d / %d pixels'%(pedRun,(allPeds[-1]<printVal[0]).sum(), (allPeds[-1]>printVal[1]).sum(),(allRms[-1]<printVal[2]).sum(), (allRms[-1]>printVal[3]).sum())
            elif len(printVal)>=2:
                print 'run %d, pixel cold/hot: %d / %d pixels'%(pedRun,(allPeds[-1]<printVal[0]).sum(), (allPeds[-1]>printVal[1]).sum())

        ped2d=[] 
        rms2d=[] 
        print 'maxmin peds ',allPeds[icurrRun].shape, allPeds[icurrRun].max()  , allPeds[icurrRun].min() ,np.percentile(allPeds[icurrRun],0.5), np.percentile(allPeds[icurrRun],99.5)
        for thisPed,thisRms in zip(allPeds,allRms):
            #ped2d.append(np.histogram(thisPed.flatten(), np.arange(np.percentile(allPeds[icurrRun],0.5), np.percentile(allPeds[icurrRun],99.5)))[0],10)
            ped2d.append(np.histogram(thisPed.flatten(), np.arange(0, np.percentile(allPeds[icurrRun],99.5)*1.2,10))[0])
            rms2d.append(np.histogram(thisRms.flatten(), np.arange(0, allRms[icurrRun].max()*1.05,0.1))[0])
        ped2dNorm=[] 
        rms2dNorm=[] 
        for thisPed,thisRms in zip(ped2d,rms2d):
            ped2dNorm.append(np.array(thisPed).astype(float)/np.array(ped2d[icurrRun]).astype(float))
            rms2dNorm.append(np.array(thisRms).astype(float)/np.array(rms2d[icurrRun]).astype(float))
        
        ped2dNorm = np.array(ped2dNorm)
        rms2dNorm = np.array(rms2dNorm)
        print 'shapes:',rms2dNorm.shape, ped2dNorm.shape

        self.__dict__['pedHisto_'+detname] = np.array(ped2d)
        self.__dict__['rmsHisto_'+detname] = np.array(rms2d)

        figDark=plt.figure(figsize=(11,10))
        gsPed=gridspec.GridSpec(2,2)
        im0 = plt.subplot(gsPed[0]).imshow(ped2d,clim=[np.percentile(ped2d,1),np.percentile(ped2d,99)],interpolation='none',aspect='auto')
        cbar0 = plt.colorbar(im0)

        im01 = plt.subplot(gsPed[1]).imshow(ped2dNorm,clim=[np.percentile(ped2dNorm,1),np.percentile(ped2dNorm,99)],interpolation='none',aspect='auto')
        cbar01 = plt.colorbar(im01)

        im2 = plt.subplot(gsPed[2]).imshow(rms2d,clim=[np.percentile(rms2d,1),np.percentile(rms2d,99)],interpolation='none',aspect='auto')
        cbar2 = plt.colorbar(im2)

        im21 = plt.subplot(gsPed[3]).imshow(rms2dNorm,clim=[np.percentile(rms2dNorm,1),np.percentile(rms2dNorm,99)],interpolation='none',aspect='auto')
        cbar21 = plt.colorbar(im21)


    def compareCommonMode(self, detname='None',common_modes=[], numEvts=100, thresADU=0., thresRms=0.):
        if detname is 'None':
            detname = self.addDetInfo(detname=detname)
            if detname == 'None':
                print 'need detector name as input! '
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

        #print imgNames
        #for img in imgs:
        #    print img.shape

        fig=plt.figure(figsize=(15,10))
        gsCM=gridspec.GridSpec(len(common_modes),2)
        for icm,cm in enumerate(common_modes):
            if cm==0:
                lims=[np.percentile(imgs[icm*2],1),np.percentile(imgs[icm*2],99)]
            else:
                lims=[np.percentile(imgs[0],1),np.percentile(imgs[0],99)]
            limsStd=[np.percentile(imgs[1],1),np.percentile(imgs[1],99)]
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
        
        if (rank==0):
            print 'Variables to be read from xtc: ',myCube.targetVarsXtc

        detInData=[]
        for k in self.Keys():
            if k.alias()!='':
                detInData.append(k.alias())

        detNames=[]
        targetVarsXtc=[]
        for idet,det in enumerate(myCube.targetVarsXtc):
            if isinstance(det, dict):
                dName = det['source']
            else:
                dName = det
            #print 'DEBUG dName ',dName
            if dName in detInData:
                detNames.append(dName)
                targetVarsXtc.append(det)
            else:
                printR(rank, 'Detector with alias %s not in data '%det)
        myCube.targetVarsXtc = [ {'source':det, 'full':1} if isinstance(det, basestring) else det for det in targetVarsXtc]
        for det in myCube.targetVarsXtc:
            self.addDetInfo(det)

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
                fout.close()
            return

        printR(rank, 'now make big cube')
        #only run on rank=0 & broadcast.
        outFileName = dirname+'/Cube_'+self.sda.fname.split('/')[-1].replace('.h5','_%s.h5.inprogress'%cubeName)
        if (onoff==0):
            outFileName=outFileName.replace('.h5','_off.h5')
        elif (onoff==1):
            outFileName=outFileName.replace('.h5','_on.h5')
        printR(rank, 'now write outputfile to : %s'%outFileName)
        try:
            #print 'DEBUG: open file %s in rank %d '%(outFileName, rank)
            fout = h5py.File(outFileName, "w",driver='mpio',comm=MPI.COMM_WORLD)
            #print 'DEBUG: opened file %s in rank %d '%(outFileName, rank)
        except:
            try:
                fout2 = h5py.File('/tmp/%d.h5'%(int(np.random.rand()*1000)), "w",driver='mpio',comm=MPI.COMM_WORLD)
                printR(rank, 'could not open the desired file for MPI writing. Likely a permission issue:')
                fout2.close()
                try:
                    import pwd
                    print 'A: ',oct(os.stat(outFileName).st_mode)[-3:]
                    print 'B: ',pwd.getpwuid(os.stat(outFileName).st_uid).pw_name
                    print 'C: ',pwd.getpwuid(os.stat(outFileName).st_uid).pw_gecos
                    printR(rank, 'owner: %s (%s), permissions: %s '%(pwd.getpwuid(os.stat(outFileName).st_uid).pw_name,pwd.getpwuid(os.stat(outFileName).st_uid).pw_gecos,oct(os.stat(outFileName).st_mode)[-3:]))
                except:
                    printR(rank, 'failed at printing info about file')
            except:
                printR(rank, 'you will need an analysis release < ana.1.3.21 or >= ana-1.3.42 for the cube to work. Solution in progress...')
            printR(rank, 'we will save the small cubed data only and return')
            if rank==0:
                outFileName = dirname+'/Cube_'+self.sda.fname.split('/')[-1].replace('.h5','_%s.h5'%cubeName)
                if (onoff==0):
                    outFileName=outFileName.replace('.h5','_off.h5')
                elif (onoff==1):
                    outFileName=outFileName.replace('.h5','_on.h5')
                fout = h5py.File(outFileName, "w")

                print 'bin the data now....be patient'
                cubeData = self.sda.makeCubeData(cubeName,onoff=onoff)  
                print 'now write outputfile (only small data) to : ',outFileName
                for key in cubeData.variables:
                    addToHdf5(fout, key, cubeData[key].values)
                #ADD CONF STUFF HERE TOO
                fout.close()
            return

        #configuration for cube making
        cube, cubeName_onoff = self.sda.prepCubeData(cubeName)
        #compare the number of bins to mpi jobs.
        nBins=myCube.binBounds.shape[0]-1
        for key in cube.addBinVars:
            nBins*=cube.addBinVars[key].shape[0]-1
        if nBins<size:
            if 'random/random' not in self.sda.Keys('random'):
                myrandom=np.random.rand(self.sda.xrData.time.shape[0])
                self.sda.addVar('random/random',myrandom)
            if isinstance(size/nBins, int):
                cube.add_BinVar({'random/random':[0.,1.,int(size/nBins)]})
            else:
                cube.add_BinVar({'random/random':[0.,1.,int(size/nBins)+1]})

        sel = self.sda.Sels[cube.useFilter]
        selString=''
        for icut,cut in enumerate(sel.cuts):
            selString+=('Cut %i: %f < %s < %f\n'%(icut, cut[1], cut[0],cut[2]))
        dsetcnf = fout.create_dataset('cubeSelection', [1.], dtype='f')
        dsetcnf.attrs['cubeSelection'] = selString

        #back to original code.
        if offEventsCube>0:
            self.sda.getOffVar('fiducials','xon',nNbr=offEventsCube, mean=False)
            self.sda.getOffVar('event_time','xon',nNbr=offEventsCube, mean=False)
            addVars = [ 'offNbrs_event_time_xon_nNbr%02d'%offEventsCube,'offNbrs_fiducials_xon_nNbr%02d'%offEventsCube]
            self.sda.cubes[cubeName].addIdxVar(addVars)
        cubeData, eventIdxDict = self.sda.makeCubeData(cubeName, returnIdx=True,onoff=onoff)

        #add small data to hdf5
        for key in cubeData.variables:
            addToHdf5(fout, key, cubeData[key].values)

        runIdx = self.dsIdxRun
        numBin = np.array(cubeData.nEntries).flatten().shape[0]
        bins_per_job = numBin/size + int((numBin%size)/max(1,numBin%size))

        binID=np.ones(bins_per_job,dtype=int); binID*=-1
        detShapes=[]
        detArrays=[]
        detSArrays=[] #for error calculation
        detMArrays=[] #for error calculation
        detOffArrays=[]
        detIArrays=[]
        for thisdetName,thisdetDict in zip(detNames, (myCube.targetVarsXtc)):
            detShape = self.__dict__[thisdetName].ped.shape
            if self.__dict__[thisdetName].det.dettype==26:
                detShape = (self.__dict__[thisdetName].ped[0]).shape
            lS = list(detShape);lS.insert(0,bins_per_job);csShape=tuple(lS)
            detShapes.append(csShape)
            #print 'DEBUG detshape C: ',detShape, thisdetName, csShape, np.array(detArrays).shape
            det_arrayBin=np.zeros(csShape)
            detArrays.append(det_arrayBin)
            detSArrays.append(det_arrayBin)
            detMArrays.append(det_arrayBin)
            detOffArrays.append(det_arrayBin)
            detIArrays.append(det_arrayBin)
            if rank==0:
                print 'for detector %s assume shape: '%thisdetName, csShape, det_arrayBin.shape
            #print 'DEBUG -- for detector %s assume shape: '%thisdetName, csShape, det_arrayBin.shape

        ###
        #nominal cube
        ###
        for ib,fids,evttimes in itertools.izip(itertools.count(), eventIdxDict['fiducial'],eventIdxDict['evttime']):
            if not (ib>=(bins_per_job*rank) and ib < bins_per_job*(rank+1)):
                continue
            print 'bin: %d has %d events, will be treated in rank %d'%(ib, len(fids),rank)
            binID[ib%bins_per_job]=ib

            nEvts_bin=0
            for ievt,evtfid, evttime in itertools.izip(itertools.count(),fids,evttimes):
                if nEvtsPerBin>0 and nEvts_bin >= nEvtsPerBin-1:
                    break
                nEvts_bin=nEvts_bin+1

                evtt = psana.EventTime(int(evttime.values),int(evtfid.values))
                evt = runIdx.event(evtt)                #now loop over detectors in this event
                for thisdetName,thisdetDict,dArray,dMArray,dSArray,dIArray in zip(detNames, (myCube.targetVarsXtc), detArrays, detMArrays, detSArrays, detIArrays):
                    det = self.__dict__[thisdetName]
                    det.evt = dropObject()
                    det.getData(evt)
                    det.processDetector()
                    
                    thisDetDataDict=getUserData(det)
                    for key in thisDetDataDict.keys():
                        if not (key=='full' or key.find('ROI')>=0 or key.find('photon_img')>=0):
                            continue
                        if (key=='full' or key.find('ROI')>=0):
                            if thisdetDict.has_key('thresADU'):
                                thisDetDataDict[key][thisDetDataDict[key]<thisdetDict['thresADU']]=0
                            elif thisdetDict.has_key('thresRms'):
                                thisDetDataDict[key][thisDetDataDict[key]<thisdetDict['thresRms']*det.rms]=0
                            dArray[ib%bins_per_job]=dArray[ib%bins_per_job]+thisDetDataDict[key]
                        else:
                            dIArray[ib%bins_per_job]=dArray[ib%bins_per_job]+thisDetDataDict[key]

                        x = thisDetDataDict[key]
                        oldM = dMArray
                        dMArray = dMArray + (x-dMArray)/(ievt+1)
                        dSArray = dSArray + (x-dMArray)*(x-oldM)
                        
        ###
        #off events cube
        ###
        if offEventsCube>0:
            nEntries_off=[]
            for ib,fids,evttimes in itertools.izip(itertools.count(), eventIdxDict['offNbrs_fiducials_xon_nNbr%02d'%offEventsCube],eventIdxDict['offNbrs_event_time_xon_nNbr%02d'%offEventsCube]):
                if not (ib>=(bins_per_job*rank) and ib < bins_per_job*(rank+1)):
                    continue
                binID[ib%bins_per_job]=ib

                nEvts_bin=0
                try:
                    allfids = fids.values.flatten()
                    allevttimes = evttimes.values.flatten()
                except:
                    if fids==[] and evttimes==[]:
                        nEntries_off.append(0)
                        continue
                    else:
                        print 'something is wrong with this bin, got fiducials like: ',fids
                        break
                #get unique event time
                evttimes, evtIdx = np.unique(allevttimes, return_inverse=True)
                #count how often each of these appears to get right mean/sum later
                noff_evt = np.bincount(evtIdx)
                nEntries_off.append(np.array(noff_evt).sum())
                print 'bin: %d has %d off events with %d unique off events'%(ib, nEntries_off[-1], len(evttimes),)
                #get matching fiducials
                off_fids = np.bincount(evtIdx, allfids)
                off_fids = off_fids/noff_evt
                for ievt,evtfid, evttime, noff_evts in itertools.izip(itertools.count(),off_fids,evttimes,noff_evt):
                    if nEvtsPerBin>0 and nEvts_bin >= nEvtsPerBin:
                        break
                    nEvts_bin=nEvts_bin+noff_evts

                    #evtt = psana.EventTime(int(evttime.values),int(evtfid.values))
                    evtt = psana.EventTime(int(evttime),int(evtfid))
                    evt = runIdx.event(evtt)
                    #now loop over detectors in this event
                    for thisdetName,thisdetDict,dArray in zip(detNames, (myCube.targetVarsXtc), detOffArrays):
                        det = self.__dict__[thisdetName]
                        det.evt = dropObject()
                        det.getData(evt)
                        det.processDetector()
                        
                        thisDetDataDict=getUserData(det)
                        for key in thisDetDataDict.keys():
                            if not (key=='full' or key.find('ROI')>=0):
                                continue
                            if thisdetDict.has_key('thresADU'):
                                thisDetDataDict[key][thisDetDataDict[key]<thisdetDict['thresADU']]=0
                                dArray[ib%bins_per_job]=dArray[ib%bins_per_job]+thisDetDataDict[key]*noff_evts
                            elif thisdetDict.has_key('thresRms'):
                                thisDetDataDict[key][thisDetDataDict[key]<thisdetDict['thresRms']*det.rms]=0
                                dArray[ib%bins_per_job]=dArray[ib%bins_per_job]+thisDetDataDict[key]*noff_evts
                            else:
                                dArray[ib%bins_per_job]=dArray[ib%bins_per_job]+thisDetDataDict[key]*noff_evts        
            addToHdf5(fout, 'nEntries_off', np.array(nEntries_off))

        #loop over arrays & bins, make image & rebin if requested.
        for detName, dArray, dOArray, dSArray, dMArray, dIArray, detDict in zip(detNames, detArrays, detOffArrays, detSArrays,detMArrays, detIArrays, (myCube.targetVarsXtc)):
            det = self.__dict__[detName]
            imgArray=[]
            imgOArray=[]
            imgSArray=[]
            imgMArray=[]
            imgIArray=[]

            nEntries = cubeData['nEntries'].values
            if len(nEntries.shape)>1:
                print 'DEBUG reshape rank %d '%(rank),' nEntries shape ',nEntries.shape
                nEntries = np.array(nEntries).flatten()[rank*bins_per_job:(rank+1)*bins_per_job], ' -- ',nEntries
            #this goes wrong for bins w/o entries. Split this?
            for binData,binOData, binSData,binMData,binIData,nent in zip(dArray,dOArray, dSArray,dMArray,dIArray,nEntries):
                if nent!=0:
                    #now divide mean & std sums to be correct.
                    binSData = np.sqrt(binSData/(nent-1))
                    binMData = binMData/nent
                if  detDict.has_key('image'):
                    thisImg = det.det.image(self.run, binData)
                    thisImgO = det.det.image(self.run, binOData)
                    thisImgS = det.det.image(self.run, binSData)
                    thisImgM = det.det.image(self.run, binMData)
                    thisImgI = det.det.image(self.run, binIData)
                else:
                    thisImg = binData
                    thisImgO = binOData
                    thisImgS = binSData
                    thisImgM = binMData
                    thisImgI = binIData
                if  detDict.has_key('rebin'):
                    thisImg = rebin(thisImg, detDict['rebin'])
                    imgArray.append(thisImg)
                    thisImgO = rebin(thisImgO, detDict['rebin'])
                    imgOArray.append(thisImgO)
                    thisImgS = rebin(thisImgS, detDict['rebin'])
                    imgSArray.append(thisImgS)
                    thisImgM = rebin(thisImgS, detDict['rebin'])
                    imgMArray.append(thisImgM)
                else:
                    imgArray.append(thisImg)
                    imgOArray.append(thisImgO)
                    imgSArray.append(thisImgS)
                    imgMArray.append(thisImgM)
                    imgIArray.append(thisImgI)
                
            #write hdf5 file w/ mpi for big array
            arShape=(numBin,)
            for i in range(0,len(imgArray[0].shape)):
                arShape+=(imgArray[0].shape[i],)
#            if rank==0:
#                print 'print array shape before saving: ',arShape
            cubeBigData = fout.create_dataset('%s'%detName,arShape)
            if storeMeanStd:
                cubeBigSData = fout.create_dataset('%s_std'%detName,arShape)
                cubeBigMData = fout.create_dataset('%s_mean'%detName,arShape)
            if offEventsCube>0:
                cubeBigOData = fout.create_dataset('%s_off'%detName,arShape)
            if (np.array(imgIArray).sum()>0):
                cubeBigIData = fout.create_dataset('%s_photon'%detName,arShape)

            for iSlice,Slice,SliceO, SliceS,SliceM,SliceI in itertools.izip(itertools.count(),imgArray,imgOArray, imgSArray,imgMArray,imgIArray):
                if np.nansum(Slice)!=0:
                    cubeBigData[rank*bins_per_job+iSlice,:] = Slice
                    print 'bin %d (%d per job)  mean %g std %g'%(rank*bins_per_job+iSlice,iSlice,np.nanmean(cubeBigData[rank*bins_per_job+iSlice,:]), np.nanstd(cubeBigData[rank*bins_per_job+iSlice,:]))
                if offEventsCube>0 and np.nansum(SliceO)>0:
                    cubeBigOData[rank*bins_per_job+iSlice,:] = SliceO
                if storeMeanStd and np.nansum(SliceS)>0:
                    cubeBigSData[rank*bins_per_job+iSlice,:] = SliceS
                if storeMeanStd and np.nansum(SliceM)>0:
                    cubeBigMData[rank*bins_per_job+iSlice,:] = SliceM
                if np.nansum(SliceI)!=0 and (np.array(imgIArray).sum()>0):
                    cubeBigIData[rank*bins_per_job+iSlice,:] = SliceI

            if det.rms is not None:
                if not detDict.has_key('image'):
                    addToHdf5(fout, 'Cfg__'+detName+'__ped', det.ped)
                    addToHdf5(fout, 'Cfg__'+detName+'__rms', det.rms)
                    if det.gain is not None:
                        addToHdf5(fout, 'Cfg__'+detName+'__gain', det.gain)
                    addToHdf5(fout, 'Cfg__'+detName+'__mask', det.mask)
                    addToHdf5(fout, 'Cfg__'+detName+'__calib_mask', det.cmask)
                    if det.x is not None:
                        addToHdf5(fout, 'Cfg__'+detName+'__x', det.x)
                        addToHdf5(fout, 'Cfg__'+detName+'__y', det.y)
                else:
                    if det.det.dettype==26:
                        addToHdf5(fout, 'Cfg__'+detName+'__ped', det.det.image(self.run,det.ped[0]))
                        addToHdf5(fout, 'Cfg__'+detName+'__rms', det.det.image(self.run,det.rms[0]))
                        addToHdf5(fout, 'Cfg__'+detName+'__gain', det.det.image(self.run,det.gain[0]))
                    else:
                        addToHdf5(fout, 'Cfg__'+detName+'__ped', det.det.image(self.run,det.ped))
                        addToHdf5(fout, 'Cfg__'+detName+'__rms', det.det.image(self.run,det.rms))
                        addToHdf5(fout, 'Cfg__'+detName+'__gain', det.det.image(self.run,det.gain))
                    addToHdf5(fout, 'Cfg__'+detName+'__mask', det.det.image(self.run,det.mask))
                    addToHdf5(fout, 'Cfg__'+detName+'__calib_mask', det.det.image(self.run,det.cmask))
                    if det.x is not None:
                        addToHdf5(fout, 'Cfg__'+detName+'__x', det.x)
                        addToHdf5(fout, 'Cfg__'+detName+'__y', det.y)


        comm.Barrier()
        printR(rank, 'first,last img mean: %g %g '%(np.nanmean(fout['%s'%detName][0]),np.nanmean(fout['%s'%detName][-1])))
        fout.close()
        comm.Barrier()
        #print 'in rank now: ',rank
        if rank==0:
            print 'renaming file from %s to %s, remove random variable if applicable'%(outFileName,outFileName.replace('.inprogress',''))
            rename_reduceRandomVar(outFileName)
