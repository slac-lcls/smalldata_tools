"""
Created on Tue Dec  8 21:31:56 2015

@author: snelson
"""
import h5py
from os import path
from os import walk
import numpy as np
from scipy import interpolate
import time
import json
import subprocess
import socket
from scipy import sparse
import tables
from matplotlib import gridspec
from pylab import ginput
from matplotlib import pyplot as plt
from utilities import getBins as util_getBins
from utilities import dictToHdf5
from utilities import shapeFromKey_h5
from utilities import hist2d

#including xarray means that you have to unset DISPLAY when submitting stuff to batch
import xarray as xr

#works from ana-1.2.9 on
try:
    from pscache import client
except:
    pass

class photons(object):
    def __init__(self,h5file,detName='epix',photName='photon'):
        self._detName=detName
        self._photName=photName
        if self._photName[-1]!='_':
            self._photName+='_'
        self._h5 = h5file
        self.shape = h5file['UserDataCfg/'+detName+'_cmask'].value.shape

    def getKeys(self):
        keyList=[]
        for data in self.__dict__.keys():
            if not data[0]=='_':
                keyList.append(data)
        return keyList

    def fillPhotArrays(self):
        eh5_dir = self._h5[self._detName]
        for key in eh5_dir.keys():
            if not (key.find(self._photName)>=0):
                continue
            #likely a different photName set of variables. Need own fillPhotArray
            if key.replace(self._photName,'').find('_')>=0:
                continue
            h5Name = self._photName+key.replace(self._photName,'')
            #print 'h5name ',h5Name
            keyName=key.replace(self._photName,'')
            self.__dict__[keyName] = eh5_dir[h5Name].value
        
    def flattenPhotArray(self, filterArray=None):
        for key in self.getKeys():
            if filterArray is None:
                filterArray=np.ones(self.__dict__[key].shape[0])
            if filterArray.shape == self.__dict__[key].shape[0]:
                self[key] = self.__dict__[key][filterArray].flatten()

    def photImgEvt(self, iEvt):
        eh5_dir = self._h5[self._detName]
        data = eh5_dir[self._photName+'data'][iEvt]
        row = eh5_dir[self._photName+'row'][iEvt]
        col = eh5_dir[self._photName+'col'][iEvt]
        if max(row)>=self.shape[0] or max(col)>=self.shape[1]:
            print 'inconsistent shape ',self.shape, max(row), max(col)
        #print eh5_dir[self._photName+'data'][iEvt].shape
        #print eh5_dir[self._photName+'row'].shape
        #print eh5_dir[self._photName+'row'][iEvt].shape
        #print max(eh5_dir[self._photName+'row'][iEvt])
        #print 'inconsistent shape ',self.shape, max(row), max(col)
        return sparse.coo_matrix( (data, (row, col)),shape=self.shape ).todense()

    def photImg(self, filterArray=None):
        if 'pHist' not in self.__dict__.keys():
            self.fillPhotArrays()
        if filterArray is None:
            filterArray=np.ones(self.pHist.shape[0]).astype(bool)
        data = self.data[filterArray].flatten()
        data = data[data>0]
        row = self.row[filterArray].flatten()[data>0]
        col = self.col[filterArray].flatten()[data>0]
        img = sparse.coo_matrix( (data, (row, col)) ).todense()
        return img

    def photonHist(self, filterArray=None):
        if 'pHist' not in self.__dict__.keys():
            self.fillPhotArrays()
        if filterArray is None:
            filterArray=np.ones(self.pHist.shape[0]).astype(bool)
        if filterArray.shape[0] == self.pHist.shape[0]:
            pHist = self.pHist[filterArray].sum(axis=0)
        else:
            return
        return pHist

class droplets(object):
    def __init__(self,h5file,detName='epix',dropName='droplet'):
        self._detName=detName
        self._dropName=dropName
        if self._dropName[-1]!='_':
            self._dropName+='_'
        self._h5 = h5file
        self._h5dir = self._h5.get_node('/'+self._detName)

    def getKeys(self):
        keyList=[]
        for data in self.__dict__.keys():
            if not data[0]=='_':
                keyList.append(data)
        return keyList

    def fillDropArrays(self, only_XYADU=False):
        for node in self._h5dir._f_list_nodes():
            if not (node.name.find(self._dropName)>=0):
                continue
            #likely a different dropName set of variables. Need own fillDropArray
            if node.name.replace(self._dropName,'').find('_')>=0:
                continue
            h5Name = self._dropName+node.name.replace(self._dropName,'')
            keyName=node.name.replace(self._dropName,'')
            if not only_XYADU:
                print 'fill drop ',h5Name
                self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()
            elif (keyName=='X' or keyName=='Y' or keyName=='adu' or keyName=='npix'):
                print 'fill drop ',h5Name
                self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()

    def flattenDropArray(self, filterArray=None):
        if filterArray is None:
            if 'adau' in self.getKeys():
                filterArray = (self.__dict__['adu']>0)
            else:
                print 'did not find adu, will not flatten'
                return

        for key in self.getKeys():
            if filterArray.shape == self.__dict__[key].shape:
                self.__dict__[key] = self.__dict__[key][filterArray].flatten()

    def getDropPixels(self, ievt, debug=False):
        print 'will need to be reimplemented'
        return
        if not 'Pix' in self.__dict__.keys():
            nDroplets = self._h5[self._detName][self._dropName+'nDroplets'][ievt]
            Pix = self._h5[self._detName][self._dropName+'Pix'][ievt]
            sizeX = self._h5[self._detName][self._dropName+'bbox_x1'][ievt] - self._h5[self._detName][self._dropName+'bbox_x0'][ievt]
            sizeY = self._h5[self._detName][self._dropName+'bbox_y1'][ievt] - self._h5[self._detName][self._dropName+'bbox_y0'][ievt]
            adu = self._h5[self._detName][self._dropName+'adu'][ievt]
            dX =  self._h5[self._detName][self._dropName+'X'][ievt]
            dY =  self._h5[self._detName][self._dropName+'Y'][ievt]
        else:
            nDroplets = self.nDroplets[ievt]
            Pix = self.Pix[ievt]
            sizeX = self.bbox_x1[ievt] - self.bbox_x0[ievt]
            sizeY = self.bbox_y1[ievt] - self.bbox_y0[ievt]
            adu = self.adu[ievt]
            dX =  self.X[ievt]
            dY =  self.Y[ievt]
        sizes = (sizeX*sizeY)
        imgs=[]
        idxPix=0
        for iDrop in range(0,nDroplets):
            img= np.array(Pix[idxPix:(idxPix+sizes[iDrop])]).reshape(sizeX[iDrop],sizeY[iDrop])
            if debug:
                print 'adu ',adu[iDrop],img.sum()
            imgs.append(img)
            idxPix+=sizes[iDrop]
        ret_dict = {'images' : imgs}
        ret_dict['adu']=adu[:len(imgs)]
        ret_dict['X']=dX[:len(imgs)]
        ret_dict['Y']=dY[:len(imgs)]
        return ret_dict

    def getDropPixelsRoi(self, ievt, mask, debug=False):
        dropInfo = self.getDropPixelsRoi(ievt, debug=debug)
        imgsROI = []
        for img,x,y in zip(dropInfo['images'],dropInfo['X'],dropInfo['Y']):
            if mask(int(x), int(y))==1:
                imgsROI.append(img)
        return imgsROI

    def plotSpectrum(self, plotLog=True, aduRange=[]):
        if len(aduRange)==0:
            aduRange=[0,700]
        elif len(aduRange)==1:
            aduRange=[0,aduRange[0]]
        elif len(aduRange)==2:
            aduRange=[aduRange[0], aduRange[1]]
            
        hst = np.histogram(self.__dict__['adu'], np.arange(aduRange[0],aduRange[1]))
        if plotLog:
            plt.semilogy(hst[1][1:],hst[0],'o')
        else:
            plt.plot(hst[1][1:],hst[0],'o')
            
    def plotAduX(self, ADUrange=[120,180],maxLim=99.5):
        if len(self.__dict__['X'].shape)>1:
            self.flattenDropArray()
        plt.figure()
        hist2d(self.__dict__['X'],self.__dict__['adu'], numBins=[702,180], histLims=[0,702,ADUrange[0], ADUrange[1]],limits=[1,maxLim])

    def plotAduY(self, ADUrange=[120,180],maxLim=99.5):    
        if len(self.__dict__['Y'].shape)>1:
            self.flattenDropArray()
        plt.figure()
        plt.subplot(211)
        hist2d(self.__dict__['Y'][self.__dict__['X']<351],self.__dict__['adu'][self.__dict__['X']<351], numBins=[766,180], histLims=[0,766,ADUrange[0], ADUrange[1]],limits=[1,maxLim])
        plt.subplot(212)
        hist2d(self.__dict__['Y'][self.__dict__['X']>353],self.__dict__['adu'][self.__dict__['X']>353], numBins=[766,180], histLims=[0,766,ADUrange[0], ADUrange[1]],limits=[1,maxLim])

    def plotXY(self, ADUrange=[120,180], npix=0):            
        allX = self.__dict__['X'][self.__dict__['adu']>ADUrange[0]]
        allY = self.__dict__['Y'][self.__dict__['adu']>ADUrange[0]]
        alladu = self.__dict__['adu'][self.__dict__['adu']>ADUrange[0]]
        if npix!=0:
            allNpix=self.__dict__['npix'][self.__dict__['adu']>ADUrange[0]]
            
        if ADUrange[1]>ADUrange[0]:
            allX = allX[alladu<ADUrange[1]]
            allY = allY[alladu<ADUrange[1]]
            if npix!=0:
                allNpix = allNpix[alladu<ADUrange[1]]
            alladu = alladu[alladu<ADUrange[1]]
        if npix!=0:
            if npix>0:
                allX = allX[allNpix==npix] 
                allY = allY[allNpix==npix] 
                alladu = alladu[allNpix==npix] 
            else:
                allX = allX[allNpix>=abs(npix)] 
                allY = allY[allNpix>=abs(npix)] 
                alladu = alladu[allNpix>=abs(npix)] 

        plt.figure()
        ndrop_int = max(1,490000./allX.shape[0])
        hist2d(allX,allY, numBins=[int(702/ndropc_int),int(766/ndrop_int)])

    def aduSlices(self,axis='y', ADUrange=[0,-1], npix=0):
        allX = self.__dict__['X'][self.__dict__['adu']>ADUrange[0]]
        allY = self.__dict__['Y'][self.__dict__['adu']>ADUrange[0]]
        alladu = self.__dict__['adu'][self.__dict__['adu']>ADUrange[0]]
        if npix!=0:
            allNpix=self.__dict__['npix'][self.__dict__['adu']>ADUrange[0]]
        if ADUrange[1]>ADUrange[0]:
            allX = allX[alladu<ADUrange[1]]
            allY = allY[alladu<ADUrange[1]]
            if npix!=0:
                allNpix = allNpix[alladu<ADUrange[1]]
            alladu = alladu[alladu<ADUrange[1]]
        if npix!=0:
            if npix>0:
                allX = allX[allNpix==npix] 
                allY = allY[allNpix==npix] 
                alladu = alladu[allNpix==npix] 
            else:
                allX = allX[allNpix>=abs(npix)] 
                allY = allY[allNpix>=abs(npix)] 
                alladu = alladu[allNpix>=abs(npix)] 
                
        if axis=='y':
            nSlice=16
            allY+=(allX>351).astype(int)*768
            sliceSize=768*2./nSlice
            binVar=allY
        elif axis=='x':
            nSlice=14
            sliceSize=704./nSlice
            binVar=allX
            
        aduS=[]
        for i in range(0,nSlice):
            aduS.append([])
        for adu,bv in zip(alladu,binVar):
            aduS[int(bv/sliceSize)].append(adu)
        return aduS
 
 
class Cube(object):
    def __init__(self, binVar='', bins=[], cubeName=None, SelName=None):
        self.binVar = binVar
        self.bins = bins
        self.SelName = SelName

        nbins = len(bins)
        if nbins==3:
            if type(bins[2]) is int:
                nbins=len(np.linspace(min(bins[0],bins[1]),max(bins[0],bins[1]),bins[2]))
            else:
                nbins=len(np.arange(min(bins[0],bins[1]),max(bins[0],bins[1]),bins[2]))

        if cubeName is not None:
            self.cubeName = cubeName
        else:
            if binVar.find('/')>=0:
                self.cubeName = '%s_%i'%(binVar.replace('/','__'), nbins)
            else:
                self.cubeName = '%s_%i'%(binVar, nbins)
        self.targetVars=[]
        #convenience variables
        self.binBounds = None   #bin boundary array (could be same as bins, but does not need to)
        self.targetVarsXtc = [] #split out the detectors that are not in the smallData

    def addVar(self, tVar):
        if isinstance(tVar, basestring):
            self.targetVars.append(tVar)
        elif isinstance(tVar, list):
            for tv in tVar:
                self.targetVars.append(tv)

    def printCube(self, Sel):
        print 'cube: ',self.cubeName
        if len(self.bins)==3 and type(self.bins[2]) is int:
            print '%s binned from %g to %g in %i bins'%(self.binVar,self.bins[0],self.bins[1],self.bins[2])
        elif len(self.bins)==3:
            print '%s binned from %g to %g in %i bins'%(self.binVar,self.bins[0],self.bins[1],int((self.bins[1]-self.bins[0])/self.bins[2]))
        elif len(self.bins)==0:
            print 'will use scan steps for binning'
        elif len(self.bins)==2:
            print 'have a single bin from %g to %g'%(self.bins[0],self.bins[1])
        elif len(self.bins)==1:
            print 'use step size of %g, determine boundaries from data'%self.bins[0]
        else:
            print '%s binned from %g to %g in %i bins'%(self.binVar,min(self.bins),max(self.bins),len(self.bins))
        for icut,cut in enumerate(Sel.cuts):
            print 'Cut %i: %f < %s < %f'%(icut, cut[1], cut[0],cut[2])
        print 'we will bin: ',self.targetVars

    def readCube(self, cubeName):
        jsonCube = json.load(open('CubeSetup_'+cubeName+'.txt','r'))
        self.binVar = jsonCube[0]
        self.bins = jsonCube[1]
        self.targetVars = jsonCube[2]
        self.SelName = jsonCube[4]
        f.close()

    def writeCube(self, cubeName):
        f = open('CubeSetup_'+cubeName+'.txt','w')
        print 'write CubeSetup file for ',cubeName, ' to ','CubeSetup_'+cubeName+'.txt'
        if isinstance(self.bins,list):
            json.dump([self.binVar,self.bins,self.targetVars,self.Sels[self.SelName].cuts,self.SelName], f, indent=3)
        else:
            json.dump([self.binVar,self.bins.tolist(),self.targetVars,self.Sels[self.SelName].cuts,self.SelName], f, indent=3)
        f.close()

class Selection(object):
    def __init__(self):
        self.cuts=[]
        self._filter=None
    def _setFilter(self,newFilter):
        if isinstance(newFilter, np.ndarray):
            self._filter = newFilter
        elif isinstance(newFilter, list):
            self._filter = np.array(newFilter)
        else:
            print 'cannot set filter array with this ',newFilter
    def addCut(self, varName, varmin, varmax):
        self.cuts.append([varName, varmin, varmax])
        self._filter=None
    def removeCut(self, varName):
        for cut in self.cuts:
            if cut[0]==varName: self.cuts.remove(cut)
        self._filter=None
    def printCuts(self):
        for icut,cut in enumerate(self.cuts):
            print 'Cut %i: %f < %s < %f'%(icut, cut[1], cut[0],cut[2])
        if isinstance(self._filter, np.ndarray):
            print 'of %d events %d pass filter'%(self._filter.shape[0], self._filter.sum())
    def add(self, additionalSel):
        for cut in additionalSel.cuts:
            self.cuts.append(cut)
        self._filter=None

class SmallDataAna(object):
    def __init__(self, expname='', run=-1, dirname='', filename='',intable=None, liveList=None):
        self._fields={}
        self._live_fields=[]
        if isinstance(liveList, list) and intable=='redis':
            self._live_fields=liveList
        self.expname=expname
        self.run=run
        if len(expname)>3:
            self.hutch=self.expname[:3]
            if dirname=='':
                self.dirname='/reg/d/psdm/%s/%s/hdf5/smalldata'%(self.hutch,self.expname)
                #run 13 and past.
                if not path.isdir(self.dirname):
                    self.dirname='/reg/d/psdm/%s/%s/ftc'%(self.hutch,self.expname)
            else:
                self.dirname=dirname
        if filename == '':
            self.fname='%s/ldat_%s_Run%03d.h5'%(self.dirname,self.expname,self.run)
            if not path.isfile(self.fname):
                self.fname='%s/%s_Run%03d.h5'%(self.dirname,self.expname,self.run)
        else:
            self.fname='%s/%s'%(self.dirname,filename)
        print 'and now open in dir: ',self.dirname,' to open file ',self.fname

        self.Sels = {}
        self.cubes = {}
        self.jobIds=[]
        self._isRedis=False
        if run == -1 or (intable is not None and intable == 'redis'):
            self.fh5=client.ExptClient(expname, host='psdb3')
            self._isRedis=True
        elif intable is not None:
            if intable == 'redis':
                self.fh5=client.ExptClient(expname, host='psdb3')
                self._isRedis=True
            elif isinstance(intable, basestring) and path.isfile(intable):
                self.fh5=tables.open_file(self.fname,'r')
            else:
                print 'pass unknown input parameter or file cannot be found: ',intable
                return None
        elif path.isfile(self.fname):
            self.fh5=tables.open_file(self.fname,'r')
        else: #if path.isfile(self.fname):
            print 'could not find file: ',self.fname
            return None

        self.xrData = {}
        #keep keys in here. Start w/ Keys from original hdf5/table/REDIS
        self.Keys(printKeys=False, areaDet=False, cfgOnly=False, returnShape=True)
        #check that required live fields are actually present
        for key in self._live_fields:
            if not key in self._fields.keys():
                print 'cannot find required variable %s, will return None!'%key
                return None
        self.ttCorr, self.ttBaseStr = self._getTTstr()

        #keep an Xarray that will become bigger on request.
        #start w/ certain fields (all 1-d & ipm channels)?
        if 'event_time' in self._fields.keys():
            #cannot use getVar as getVar is using the presence of self._tStamp
            if self._isRedis:
                evttime = self.fh5.fetch_data(self.run,['event_time'])['event_time']
            else:
                evttime = self.fh5.get_node('/event_time').read()
            self._fields['event_time'][1]='inXr'
            evttime_sec = evttime >> 32
            evttime_nsec = evttime - (evttime_sec << 32)
            self._tStamp = np.array([np.datetime64(int(tsec*1e9+tnsec), 'ns') for tsec,tnsec in zip(evttime_sec,evttime_nsec)])
            #self._tStamp = np.datetime64(int(sec*1e9+nsec), 'ns')
            self.xrData = xr.DataArray(evttime, coords={'time': self._tStamp}, dims=('time'),name='event_time')
        elif not self._isRedis:
            timeData = self.fh5.get_node('/EvtID','time').read()
            if timeData is None:
                print 'could not find eventtime data ',self._fields.keys()
            #evt_id.time()[0] << 32 | evt_id.time()[1] 
            else:
                evttime = (timeData[:,0].astype(np.int64) << 32 | timeData[:,1])
                self._tStamp = np.array([np.datetime64(int(ttime[0]*1e9+ttime[1]), 'ns') for ttime in timeData])
                self.xrData = xr.DataArray(evttime, coords={'time': self._tStamp}, dims=('time'),name='EvtID__time') 
                self._fields['EvtID/time'][1]='inXr'
        else:
            print 'could not create xarray'
            return
        self._addXarray()
        #there won't be any xarray files of correct size when running "live"
        if not self._isRedis:
            self._readXarrayData()

    def __del__(self):
        try:
            from mpi4py import MPI
        except:
            self._writeNewData()
            return
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank==0:
            self._writeNewData()
        return

    def addToLive(self, liveKeys=[]):
        if isinstance(liveKeys, basestring):
            if liveKeys in self._fields.keys():
                self._live_fields.append(liveKeys)
        else:
            for key in liveKeys:
                if key in self._fields.keys():
                    self._live_fields.append(key)

    def _getXarrayDims(self,key,tleaf_name=None, setNevt=-1):
        coords=None
        dims=None
        try:
            if self._isRedis:
                dataShapeIn = self.fh5.keyinfo(run=self.run)[key][0]
                if setNevt<0:
                    setNevt = self._tStamp.shape[0]
                dataShape = [setNevt]
                for shp in dataShapeIn:
                    dataShape.append(shp)
                dataShape = tuple(dataShape)
                #dataShape = self.fh5.keyinfo(run=-1)[key]
                if tleaf_name is None:
                    tleaf_name = key
            else:                
                dataShape = self.fh5.get_node(key,tleaf_name).shape
                setNevt = dataShape[0]
        except:
            return np.array(0), coords, dims
        if len(dataShape)==1:
            coords={'time': self._tStamp[:setNevt]}
            dims=('time')
        elif len(dataShape)==2:
            if dataShape[1]==1:
                coords={'time': self._tStamp[:setNevt]}
                dims=('time')
            elif key=='/EvtID':
                return dataShape, coords, dims
            elif tleaf_name=='channels':
                dimArr = ['%d'%i for i in range(0,dataShape[1])]
                coords={'time': self._tStamp[:setNevt],'channels':dimArr}
                dims=('time','channels')
            elif tleaf_name.find('AngPos')>=0:
                dimArr = ['AngX','AngY','PosX','PosY']
                coords={'time': self._tStamp[:setNevt],'AngPos':dimArr}
                dims=('time','AngPos')
            elif tleaf_name.find('com')>=0:
                dimArr = ['axis0','axis1']
                coords={'time': self._tStamp[:setNevt],'axes':dimArr}
                dims=('time','axes')
            else: #currently save all 1-d data.
                dimArr = np.arange(0, dataShape[1])
                coords={'time': self._tStamp[:setNevt],'dim0':dimArr}
                dims=('time','dim0')
                #print '1-d data per event, not IPM-channels: ',key,tleaf_name, dataShape
        else:
            coords={'time': self._tStamp[:setNevt]}
            dims = ['time']
            for dim in range(len(dataShape)-1):
                thisDim = np.arange(0, dataShape[dim+1])
                dimStr = 'dim%d'%dim
                coords[dimStr] = thisDim
                dims.append(dimStr)
            #print 'more >-2-d data per event, fill on request ',key,tleaf_name, dataShape

        return dataShape, coords, dims

    def _updateXarray_fromREDIS(self):
        if not self._isRedis:
            return
        #redisInfo = self.fh5.keyinfo(run=-1)
        keys_for_xarray=[]
        if len(self._live_fields)>0:
            keys_for_xarray = self._live_fields
        else:
            redisInfo = self.fh5.keyinfo(run=self.run)
            for key in redisInfo.keys():
                #1-dim date
                if len(redisInfo[key][0])<=2:
                    keys_for_xarray.append(key)
        nEvts = self.fh5.fetch_data(run=self.run, keys=['event_time'])['event_time'].shape[0]
        nEvts_in_Xarray = -1
        try:
            nEvts_in_Xarray = self.xrData['event_time'].shape[0]
            #print 'DEBUG:---events--',nEvts_in_Xarray, nEvts
            if nEvts == nEvts_in_Xarray:
                #print 'DEBUG: no update needed'
                return
        except:
            pass

        #FIX ME: don't know how to deal with added datasets on update here. Bummer
        #KLUDGE ME: for now, just need to put all that code into the update function for the notebook
        keys_for_xarray.append('event_time')
        data_for_xarray = self.fh5.fetch_data(self.run,keys=keys_for_xarray)
        evttime = data_for_xarray['event_time']
        #print 'update from %d events to %d '%(nEvts_in_Xarray, evttime.shape[0])
        self._fields['event_time'][1]='inXr'
        evttime_sec = evttime >> 32
        evttime_nsec = evttime - (evttime_sec << 32)
        self._tStamp = np.array([np.datetime64(int(tsec*1e9+tnsec), 'ns') for tsec,tnsec in zip(evttime_sec,evttime_nsec)])
        #print 'DEBUG: ',self._tStamp.shape,data_for_xarray['event_time'].shape,' -- ',data_for_xarray['ipm2/sum'].shape
        minNevt=data_for_xarray[keys_for_xarray[0]].shape[0]
        for key in keys_for_xarray:
            if data_for_xarray[key].shape[0]<minNevt:
                #print 'DEBUG: mismatched data...',key
                minNevt=data_for_xarray[key].shape[0]
        self.xrData = xr.DataArray(evttime[:minNevt], coords={'time': self._tStamp[:minNevt]}, dims=('time'),name='event_time')
        for key in keys_for_xarray:
            if key == 'event_time':
                continue
            dataShape, coords, dims = self._getXarrayDims(key, setNevt=minNevt)
            tArrName = key.replace('/','__')
            self.xrData = xr.merge([self.xrData, xr.DataArray(data_for_xarray[key].squeeze()[:minNevt], coords=coords, dims=dims,name=tArrName) ])
        #now make xarray summary reflect new smaller xarray.
        for key in self._fields.keys():
            if key in keys_for_xarray:
                self._fields[key][1]='inXr'
            else:
                self._fields[key][1]='onDisk'

        return

    def _addXarray(self):
        #filling info from redis.
        if self._isRedis:
            self._updateXarray_fromREDIS()
            return
        #methods for h5-files.
        for node in self.fh5.root._f_list_nodes():
            key = node._v_pathname
            if not isinstance(node, tables.group.Group):
                fieldkey = key[1:]
                #if key != '/event_time':
                if self._fields[fieldkey][1]=='onDisk':
                    self.xrData = xr.merge([self.xrData, xr.DataArray(self.fh5.get_node(key).read(), coords={'time': self._tStamp}, dims=('time'),name=key[1:]) ])
                    self._fields[fieldkey][1]='inXr'
                continue
            for tleaf in node._f_list_nodes():
                if not isinstance(tleaf, tables.group.Group):         
                    fieldkey = key[1:]
                    tArrName = '%s__%s'%(fieldkey,tleaf.name)
                    fieldName = fieldkey+'/'+tleaf.name
                    if self.fh5.get_node(key,tleaf.name).shape[0]!=self._tStamp.shape[0]: 
                        continue
                    dataShape, coords,dims = self._getXarrayDims(key,tleaf.name)
                    #limit on dimensions here!
                    if len(dataShape)<=2 and coords is not None:
                        self.xrData = xr.merge([self.xrData, xr.DataArray(self.fh5.get_node(key,tleaf.name).read().squeeze(), coords=coords, dims=dims,name=tArrName) ])
                        self._fields[fieldName][1]='inXr'

    def addVar(self, name='newVar',data=[]):
        if not isinstance(data, np.ndarray):
            try:
                data = np.array(data)
            except:
                print 'data is not array and could not be cast to one either'
                return
        if data.shape[0] < self._tStamp.shape[0]:
            print 'only adding event based data is implemented so far'
            return
        if data.shape[0] > self._tStamp.shape[0]:
            print 'have more events, only attach ones matching time stamps'
            data=data[self._tStamp.shape[0]]

        if len(data.shape)==1:
            self.xrData = xr.merge([self.xrData, xr.DataArray(data, coords={'time': self._tStamp}, dims=('time'),name=name) ])
        else:
            coords={'time': self._tStamp}
            dims = ['time']
            dataShape = data.shape
            for dim in range(len(dataShape)-1):
                thisDim = np.arange(0, dataShape[dim+1])
                dimStr = 'dim%d'%dim
                coords[dimStr] = thisDim
                dims.append(dimStr)
            
            newArray = xr.DataArray(data, coords=coords, dims=dims,name=name)
            self.xrData = xr.merge([self.xrData, newArray])

        #if not self._isRedis and name[0]!='/': name='/'+name
        name = name.replace('__','/')
        if name not in self._fields.keys():
            #print 'add a new variable to Xarray: ',name
            self._fields[name]=[data.shape, 'inXr', 'mem']
        elif self._fields[name][2]=='main':
            #print 'add a variable from the main data to Xarray: ',name
            self._fields[name]=[data.shape, 'inXr', 'main']
        else:
            #print 'add a variable from an netcdf to Xarray: ',name
            self._fields[name]=[data.shape, 'inXr', 'xrfile']

    def _updateFromXarray(self):
        """
        this function looks for keys in the xarray that are NOT from the original files (e.g. created by xr.assign...)
        """
        for key in self.xrData.keys():
            fieldName = key.replace('__','/')
            if fieldName not in self._fields.keys():
                print 'added new data field %s to list',key
                self._fields[fieldName]={self.xrData[key].shape, 'inXr', 'mem'}

    def _writeNewData(self):
        """
        write newly created fields to netcdf fiels that can be loaded in future sessions
        """
        #do not write files for REDIS (as you are likely auto-updating...)
        if self._isRedis:
            return
        print 'save derived data to be loaded next time:'
        for key in self._fields.keys():
            if self._fields[key][2] == 'mem':
                print 'saving data for field: ',key, self._fields[key]
                data = self.getVar(key)
                #print 'DEBUG: shape ',data.shape
                if key[0]=='/': key = key[1:]
                if isinstance(data, xr.DataArray):
                    print 'get xarray for ',key
                    new_xrData = data
                elif isinstance(data, np.ndarray):
                    if len(data.shape)==1:
                        new_xrData = xr.DataArray(data, coords={'time': self._tStamp}, dims=('time'),name=key) 
                    else:
                        coords={'time': self._tStamp}
                        dims = ['time']
                        dataShape = data.shape
                        for dim in range(len(dataShape)-1):
                            thisDim = np.arange(0, dataShape[dim+1])
                            dimStr = 'dim%d'%dim
                            coords[dimStr] = thisDim
                            dims.append(dimStr)
                        new_xrData = xr.DataArray(data, coords=coords, dims=dims,name=name)
                else:
                    print 'was passed data which is neither xArray nor np. array. will not save ',key

                xrname='%s/xr_%s_%s_Run%03d.nc'%(self.dirname,key.replace('/','__'),self.expname,self.run)
                print 'data for %s is only in memory, write to file with name: %s '%(key,xrname)
                new_xrData.to_netcdf(xrname,engine='h5netcdf')

    def _readXarrayData(self):
        for (dirpath, dirnames, filenames) in walk(self.dirname):        
            for fname in filenames:
                if fname.find('xr')==0 and fname.find('Run%03d'%self.run)>=0:
                    add_xrDataSet = xr.open_dataset(self.dirname+'/'+fname,engine='h5netcdf')
                    key = fname.replace('xr_','').replace('%s_Run%03d.nc'%(self.expname,self.run),'')
                    if key[-1]=='_':key=key[:-1]
                    if len(add_xrDataSet[key].shape)>2:
                        continue
                    if (len(add_xrDataSet[key].shape)==2 and add_xrDataSet[key].shape[1]<10):
                        continue
                    key = key.replace('__','/')
                    self._fields[key]=[add_xrDataSet[key].shape, 'inXr', 'xrfile']
                    values = add_xrDataSet[key].values
                    self.addVar(key, values)
                    #need to print this dataset, otherwise this does not work. Why #DEBUG ME
                    print 'found filename %s, added data for key %s '%(fname, key), add_xrDataSet[key]

    #FIX ME: need to fix this! this will NOT work anymore....
    def setRun(self, run):
        self.run=run
        self.fname='%s/ldat_%s_Run%03d.h5'%(self.dirname,self.expname,self.run)
        if path.isfile(self.fname):
            self.fh5=tables.open_file(self.fname,'r')
            if self.hasKey('tt/ttCorr'):
                self.ttCorr = 'tt/ttCorr'
            elif self.hasKey('ttCorr/tt'):
                self.ttCorr = 'ttCorr/tt'
            else:
                self.ttCorr = None
            self.ttBaseStr = 'tt/'
            if not self.hasKey(self.ttBaseStr+'AMPL'):
                if self.hasKey('tt/XPP_TIMETOOL_AMPL'):
                    self.ttBaseStr = 'tt/XPP_TIMETOOL_'
                elif self.hasKey('tt/TIMETOOL_AMPL'):
                    self.ttBaseStr = 'tt/TIMETOOL_'
                elif self.hasKey('tt/TTSPEC_AMPL'):
                    self.ttBaseStr = 'tt/TTSPEC_'
            
    def setExperiment(self,expname):
        self.expname=expname
        self.fname='%s/ldat_%s_Run%03d.h5'%(self.dirname,self.expname,self.run)
        if path.isfile(self.fname):
            self.fh5=tables.open_file(self.fname,'r')
    def setSmallDataDir(self, dirname):
        self.dirname=dirname
        self.fname='%s/%s_Run%03d.h5'%(self.dirname,self.expname,self.run)
        if not path.isfile(self.fname):
            self.fname='%s/ldat_%s_Run%03d.h5'%(self.dirname,self.expname,self.run)
        if path.isfile(self.fname):
            self.fh5=tables.open_file(self.fname,'r')
    def Keys2d(self, inh5 = None, printKeys=False):
        return self.Keys(inh5 = inh5, printKeys=printKeys, areaDet=True)

    def Keys(self, name=None, inh5 = None, printKeys=False, areaDet=False, cfgOnly=False, returnShape=False):
        """
        return list of available keys, allowing for filter to only print subsets/
        returnShape=True will fill the initial _fields dict
        """
        keys = []
        keyShapes=[]
        keysFiltered = []
        keyShapesFiltered=[]

        if inh5 == None and self.fh5:
            fh5 = self.fh5
        else:
            fh5 = inh5        

        #DEBUG ME - do not go back to hdf5 or redis if fields are already loaded.
        if len(self._fields.keys())>0 and not returnShape:
            keys=self._fields.keys()
            for key in keys:
                keyShapes.append((0,))
        elif self._isRedis:
            #redisInfo = self.fh5.keyinfo(run=-1)
            redisInfo = self.fh5.keyinfo(run=self.run)
            for key in redisInfo.keys():
                keys.append(key)
                keyShapes.append(redisInfo[key][0])
        elif fh5:
            for node in fh5.root._f_list_nodes() :
                key = node._v_pathname
                if cfgOnly and not key.find('Cfg')>=0:
                    continue
                if areaDet and key.find('Cfg')>=0:
                    continue
                if isinstance(name, basestring) and key.find(name)<0:
                    continue
                thiskey=''
                if not isinstance(node, tables.group.Group):
                    thiskey='%s'%key
                    keys.append(thiskey)
                    keyShapes.append(shapeFromKey_h5(fh5, thiskey))
                else:
                    for tleaf in node._f_list_nodes():   
                        if isinstance(tleaf, tables.group.Group):
                            for leaf in tleaf._f_list_nodes():   
                                thiskey='%s/%s'%(tleaf._v_pathname,leaf.name)
                                keys.append(thiskey)
                                keyShapes.append(shapeFromKey_h5(fh5, thiskey))
                        else:
                            thiskey='%s/%s'%(key,tleaf.name)
                            keys.append(thiskey)
                            keyShapes.append(shapeFromKey_h5(fh5, thiskey))

        for thiskey,thiskeyshape in zip(keys,keyShapes):
            if isinstance(name, basestring) and thiskey.find(name)<0:
                continue
            if areaDet and len(thiskeyshape)<=2:
                continue
            if thiskey[0]=='/': thiskey=thiskey[1:]
            keysFiltered.append(thiskey)
            keyShapesFiltered.append(thiskeyshape)
            if printKeys:
                print thiskey
        if returnShape:
            for tkey, tshape in zip(keysFiltered, keyShapesFiltered):
                if tkey not in self._fields.keys():
                    self._fields[tkey] = [tshape, 'onDisk', 'main']

        return keysFiltered
        
    def printSelection(self, selName=None, brief=True):
        if selName is not None:
            if selName not in self.Sels:
                print 'could not find selection ',selName,', options are: ',
                for sel in self.Sels:
                    print sel
                return
            print self.Sels[selName].printCuts()
            return
        for sel in self.Sels:
            print sel
            if not brief:
                print self.Sels[sel].printCuts()
                print '--------------------------'

    def nEvts(self,printThis=False):
        if ('fiducials' in self._fields.keys()):
            nent = self._fields['fiducials'][0][0]
        elif ('EvtID/fid' in self._fields.keys()):
            nent = self._fields['EvtID/fid'][0][0]
        else:
            print 'could not find dataset with fiducials'
            nent=-1
        if printThis:
            print 'Number of events in %s is %i'%(self.fname, nent)
        return nent

    def printRunInfo(self):
        self.nEvts(printThis=True)
        isScan=False
        scanVar = self.getScanName()
        if scanVar is not None:
            isScan=True
        if isScan:
            nPoints=np.unique(self.getVar('scan/%s'%scanVar)).shape[0]
            print 'this run is a scan of %s with %d points'%(scanVar,nPoints)

    def hasKey(self, inkey):
        if inkey in self._fields.keys():
            return True

    def _getTTstr(self):
        """
        function to determine the string for the timetool variables in the desired run
        necessary as naming scheme evolved over time
        """
        ttCorr = None
        ttBaseStr = 'tt/'
        if 'tt/ttCorr' in self._fields.keys():
            ttCorr = 'tt/ttCorr'
        elif 'ttCorr/tt' in self._fields.keys():
            ttCorr = 'ttCorr/tt'
        if not ttBaseStr+'AMPL' in self._fields.keys():
            if 'tt/XPP_TIMETOOL_AMPL'  in self._fields.keys():
                ttBaseStr = 'tt/XPP_TIMETOOL_'
            elif 'tt/TIMETOOL_AMPL'  in self._fields.keys():
                ttBaseStr = 'tt/TIMETOOL_'
            elif 'tt/TTSPEC_AMPL'  in self._fields.keys():
                ttBaseStr = 'tt/TTSPEC_'
        return ttCorr, ttBaseStr

    def addCut(self, varName, cmin, cmax, SelName):
        if not self.Sels.has_key(SelName):
            self.Sels[SelName] = Selection()
        self.Sels[SelName].addCut(varName, cmin,cmax)
        Filter = self.getFilter(SelName=SelName)
        self.Sels[SelName]._setFilter(Filter)
    def removeCut(self, varName, SelName):
        if not self.Sels.has_key(SelName):
            print 'Selection with name %s does not exist, cannot remove cut'%SelName
            return
        self.Sels[SelName].removeCut(varName)
        Filter = self.getFilter(SelName=SelName)
        self.Sels[SelName]._setFilter(Filter)
    def printCuts(self, SelName):
        self.Sels[SelName].printCuts()
    def getFilterLaser(self, SelName, ignoreVar=[]):
        ignoreVar.append('lightStatus/laser')
        onFilter = self.getFilter(SelName=SelName, ignoreVar=ignoreVar)
        if self.ttCorr is not None:
            ignoreVar.append(self.ttCorr)
        if self.hasKey(self.ttBaseStr+'AMPL'):
            ignoreVar.append(self.ttBaseStr+'AMPL')
            ignoreVar.append(self.ttBaseStr+'FLTPOS')
            ignoreVar.append(self.ttBaseStr+'FLTPOS_PS')
            ignoreVar.append(self.ttBaseStr+'FLTPOSFWHM')
        offFilter = self.getFilter(SelName=SelName, ignoreVar=ignoreVar)
        #return [onFilter, offFilter]

        FilterOff = (offFilter&[self.getVar('lightStatus/laser')!=1])
        FilterOn = (onFilter&[self.getVar('lightStatus/laser')==1])
        return [FilterOn.squeeze() , FilterOff.squeeze()]
        
    def getFilter(self, SelName=None, ignoreVar=[]):
        try:
            total_filter=np.ones_like(self.getVar('EvtID/fid')).astype(bool)
        except:
            total_filter=np.ones_like(self.getVar('fiducials')).astype(bool)
        if SelName==None or SelName not in self.Sels.keys():
            return total_filter
        if self.Sels[SelName]._filter is not None and len(ignoreVar)==0:
            return self.Sels[SelName]._filter
        filters=[]
        for thiscut in self.Sels[SelName].cuts:
            if not thiscut[0] in ignoreVar:
                thisPlotvar=self.get1dVar(thiscut[0])
                filters.append(~np.isnan(thisPlotvar))
                if len(thiscut)==3:
                    filters.append((thisPlotvar > thiscut[1]) & (thisPlotvar < thiscut[2]))
                else:
                    filters.append(thisPlotvar != thiscut[1])
        for ft in filters:
            total_filter&=ft                
        return total_filter
        
    def saveFilter(self, baseName='boolArray',SelName=None, ignoreVar=[]):
        total_filter = self.getFilter(SelName=SelName, ignoreVar=ignoreVar)
        np.savetxt('%s_Run%03d.txt'%(baseName, self.run),total_filter.astype(bool),fmt='%5i')

    def getSelIdx(self, SelName):
        if self.hasKey('EvtID'):
            fids = self.getVar('/EvtID/fid')
            times = self.getVar('/EvtID/time')
        elif self.hasKey('/fiducials'):
            fids = self.getVar('/fiducials')
            times = self.getVar('/event_time')
        elif self.hasKey('fiducials'):
            fids = self.getVar('fiducials')
            times = self.getVar('event_time')
        Filter =  self.getFilter(SelName)
        selfids = [ (ifid,itime) for ifid,itime in zip(fids[Filter],times[Filter])]
        return selfids
        
    def getVar(self, plotvar, Filter=None, addToXarray=True):
        #get the signal variable
        if isinstance(plotvar,list):
            plotvar = plotvar[0]
            sigROI = plotvar[1]
        else:
            sigROI=[]

        if isinstance(Filter,basestring):
            Filter = self.getFilter(Filter)

        if self._fields[plotvar][1]=='inXr':
            fullData = self.xrData[plotvar.replace('/','__')]
            if Filter is None:
                return fullData.values
            else:
                return fullData[Filter].values

        #check if this variable has been added to xarray and needs to be added to fields
        if plotvar not in self._fields.keys():
            self._updateFromXarray()
        if plotvar not in self._fields.keys():
            print 'available keys are: ',self._fields.keys()
            print 'signal variable %s not in list'%(plotvar)
            return

        #FIX ME
        #if Filter picks > 50% of events, get all  data, add to xarray & return filtered data after
        #if only few events are picked, just get those events.
        try:
            if Filter is None:
                if not self._isRedis:
                    if len(plotvar.split('/'))>1:
                        vals = self.fh5.get_node('/'+'/'.join(plotvar.split('/')[:-1]),plotvar.split('/')[-1]).read()
                    else:
                        vals = self.fh5.get_node('/'+plotvar).read()
                else:
                    vals = self.fh5.fetch_data(run=self.run,keys=[plotvar])[plotvar]
                if vals.shape[0]==self._tStamp.shape[0]:
                    tArrName = plotvar.replace('/','__')
                    if addToXarray:
                        self.addVar(tArrName, vals)
            else:
                if not self._isRedis:
                    if len(plotvar.split('/'))>1:
                        #seems like this is not how it works.
                        #vals = self.fh5.get_node('/'+'/'.join(plotvar.split('/')[:-1]),plotvar.split('/')[-1]).__getitem__[Filter]
                        vals = self.fh5.get_node('/'+'/'.join(plotvar.split('/')[:-1]),plotvar.split('/')[-1]).read()[Filter]
                    else:
                        #seems like this is not how it works.
                        #vals = self.fh5.get_node('/'+plotvar).__getitem__(Filter)
                        vals = self.fh5.get_node('/'+plotvar).read()[Filter]
                else:
                    vals = self.fh5.fetch_data(run=self.run,keys=[plotvar])[plotvar][Filter]
            return vals.squeeze()
        except:
            print 'failed to get data for ',plotvar
            return

    def getRedVar(self, plotvar,threshold=-1e25):
        if isinstance(plotvar, list):
            sigROI=plotvar[1]
            plotvar=plotvar[0]
        vals = self.getVar(plotvar)
        if vals is None:
            return
        if threshold!=-1e25:
            vals = vals[vals<threshold]=0
        if len(vals.shape)>1:
            if sigROI!=[]:
                if len(vals.shape)==2:
                    if not isinstance(sigROI, list):
                        vals=vals[:,sigROI]
                    elif len(sigROI)>1:
                        vals=vals[:,sigROI[0]:sigROI[1]]
                    else:
                        vals=vals[:,sigROI[0]]
                elif len(vals.shape)==3:
                    if not isinstance(sigROI, list):
                        vals=vals[:,sigROI]
                    elif len(sigROI)==1:
                        vals=vals[:,sigROI[0]]
                    elif len(sigROI)==2:
                        vals=vals[:,sigROI[0]:sigROI[1],:]
                    else:
                        vals=vals[:,sigROI[0]:sigROI[1],sigROI[2]:sigROI[3]]
                else:
                    print 'this dimension is not yet implemented:',len(sig.shape)
        return vals

    def get1dVar(self, plotvar,threshold=-1e25):
        vals = self.getRedVar(plotvar,threshold)
        while len(vals.shape)>1:
            vals = vals.sum(axis=1)
        return vals

    #make delay another Xarray variable.
    def getDelay(self, use_ttCorr=True, addEnc=False):
        """
        function to get the xray-laser delay from the data
        usage:
        getDelay(): get the delay from lxt and/or encoder stage, add the timetool correction
        getDelay(use_ttCorr=False): get the delay from lxt and/or encoder stage, NO timetool correction
        getDelay(addEnc=True): get the delay from lxt, add encoder stage and timetool correction
        """
        ttCorrStr, ttBaseStr = self._getTTstr()
        if self.ttCorr is not None:
            ttCorr=self.getVar(self.ttCorr)
        if (np.nanstd(ttCorr)==0):
            ttCorr=self.getVar(self.ttBaseStr+'FLTPOS_PS')
        nomDelay=np.zeros_like(ttCorr)

        isDaqDelayScan=False
        scanVar = self.getScanName()
        if scanVar.find('lxt')>=0:
            isDaqDelayScan=True
            #print 'DEBUG: found that we have a delay scan'
            nomDelay=self.getVar('scan/%s'%scanVar)*1e12

        if not isDaqDelayScan:
            if self.hasKey('enc/lasDelay'):
                encVal = self.getVar('enc/lasDelay')
                #print 'DEBUG: encoder info',encVal.std()
                if encVal.std()>1e-9:
                    nomDelay=encVal
                    addEnc=False
                elif encVal.std()>1e-15:
                    nomDelay=encVal*1e12
                    addEnc=False
                elif self.hasKey('enc/ch0'):
                    encVal = self.getVar('enc/ch0')
                    if encVal.std()>1e-15 and encVal.std()<1e-9:
                        nomDelay=encVal*1e12
                        #now look at the EPICS PV if everything else has failed.
                    elif encVal.std()>1e-3:
                        nomDelay=encVal
                    else:
                        print 'strange encoder value for runs taken before encoder FEX....', encCal.std()
                else:
                    epics_delay = self.getVar('epics/lxt_ttc')
                    if epics_delay.std()!=0:
                        nomDelay = epics_delay

        if addEnc and self.hasKey('enc/lasDelay'):
            print 'required to add encoder value, did not find encoder!'
        if addEnc and self.hasKey('enc/lasDelay'):
            if self.getVar(fh5,'enc/lasDelay').std()>1e-6:
                nomDelay+=getVar(fh5,'enc/lasDelay').value

        if use_ttCorr:
            #print 'DEBUG adding ttcorr,nomdelay mean,std: ',ttCorr.mean(),nomDelay.mean(),ttCorr.std(),nomDelay.std()
            delay = (ttCorr+nomDelay)
        else:
            delay = nomDelay

        self.addVar('delay', delay)
        return delay

    def getPeak(self, plotvar, numBins=[100],setCuts=None, applyCuts=None, limits=[1,99,'p'],fig=None,asHist=False,sigROI=[]):
        hst = plotVar(plotvar, numBins=[100],setCuts=None, applyCuts=None, limits=[1,99,'p'],fig=None,asHist=False,sigROI=[])
        if len(hst)==2:
            return [hst[0].max(), hst[1][hst[0].argmax()]]
        else:
            print 'getPeak is not yet implemented for this type of data (need 1d histo)'

    def plotVar(self, plotvar, numBins=[100],setCuts=None, applyCuts=None, limits=[1,99,'p'],fig=None,asHist=False):
        if not isinstance(numBins, (list, tuple)):
            numBins = [numBins]
        if isinstance(plotvar, basestring) or (len(plotvar)==2 and (isinstance(plotvar[0], basestring) and not isinstance(plotvar[1], basestring))):
            if len(numBins)!=1:
                print 'bin# needs to be of same dimensions as plotvariables (1d)'
                return
            return self.plotVar1d(plotvar, numBins=numBins[0],setCuts=setCuts, applyCuts=applyCuts, limits=limits,fig=fig)
        elif len(plotvar)>2:
            print 'plotting of 3 variables is not defined yet'
            return
        if len(numBins)!=2:
            if len(numBins)==1:
                numBins=[numBins[0],numBins[0]]
            else:
                print 'bin# needs to be of same dimentions as plotvariables (2d)'
        return self.plotVar2d(plotvar, numBins=numBins,setCuts=setCuts, applyCuts=applyCuts, limits=limits,fig=fig,asHist=asHist)

    def plotVar1d(self, plotvar, numBins=100,setCuts=None, applyCuts=None, limits=[1,99,'p'],fig=None):
        if isinstance(plotvar,list):
            if not (self.hasKey(plotvar[0]) or plotvar[0]=='delay'): 
                print 'request variable %s not in littleData file'%plotvar
                return
        else:
            if not (self.hasKey(plotvar) or plotvar=='delay'): 
                print 'request variable %s not in littleData file'%plotvar
                return
        if plotvar=='delay':
            vals = self.getDelay(use_ttCorr=True)
        elif len(plotvar)==1 and plotvar.find('droplets')>=0:
            vals = self.getRedVar(plotvar)
        else:
            vals = self.get1dVar(plotvar)

        total_filter = np.ones(vals.shape[0]).astype(bool)
        if applyCuts is not None and self.Sels.has_key(applyCuts):
            total_filter =  self.getFilter(applyCuts, [plotvar])
        vals = vals[total_filter]

        if  len(plotvar)==1 and plotvar.find('droplets')>=0:
            vals = vals.flatten()[vals.flatten()>0]

        if limits[2]=='p':
            pmin = np.percentile(vals,limits[0])
            pmax = np.percentile(vals,limits[1])
            if np.isnan(pmin): pmin=np.nanmin(vals)
            if np.isnan(pmax): pmax=np.nanmax(vals)
        else:
            pmin=min(limits[0],limits[1])
            pmax=max(limits[0],limits[1])
        hst = np.histogram(vals[~np.isnan(vals)],np.linspace(pmin,pmax,numBins))
        print 'plot %s from %g to %g'%(plotvar,pmin,pmax)
        if fig is None:
            fig=plt.figure(figsize=(8,5))

        plt.plot(hst[1][:-1],hst[0],'o')
        plt.xlabel(plotvar)
        plt.ylabel('entries')
        if setCuts is not None and self.Sels.has_key(setCuts):
            p = np.array(ginput(2))
            p=[p[0][0],p[1][0]]
            self.Sels[setCuts].addCut(plotvar,min(p),max(p))
        return hst

    def plotVar2d(self, plotvars, setCuts=None, applyCuts=None, limits=[1,99,'p'], asHist=False,numBins=[100,100],fig=None):
        for plotvar in plotvars:
            if isinstance(plotvar,list):
                plotvar = plotvar[0]
            if not self.hasKey(plotvar) or plotvar == 'delay': 
                print 'request variable %s not in littleData file'%(plotvar)
                return
        vals=[]
        for plotvar in plotvars:
            if plotvar == 'delay':
                vals=self.getDelay(use_ttCorr=True)
            #elif len(plotvar)==1 and plotvar.find('droplets')>=0:
            #    vals = self.fh5[plotvar].value
            else:   
                vals.append(self.get1dVar(plotvar))
            if len(vals[-1].shape)!=1:
                print 'variable %s does not return a scalar, this is not yet implemented'%plotvar
        pmin0 = np.nanmin(vals[0]); pmin1 = np.nanmin(vals[1]);
        pmax0 = np.nanmax(vals[0]); pmax1 = np.nanmax(vals[1]);
        if limits[0]>0:
            if not np.isnan(np.percentile(vals[0],limits[0])):
                pmin0 = np.percentile(vals[0],limits[0])
            if not np.isnan(np.percentile(vals[1],limits[0])):
                pmin1 = np.percentile(vals[1],limits[0])
        if limits[1]<100:
            if not np.isnan(np.percentile(vals[0],limits[1])):
                pmax0 = np.percentile(vals[0],limits[1])
            if not np.isnan(np.percentile(vals[1],limits[1])):
                pmax1 = np.percentile(vals[1],limits[1])
        print 'plots %s from %g to %g and  %s from %g to %g '%(plotvars[0],pmin0,pmax0,plotvars[1],pmin1,pmax1)
        total_filter=np.ones(vals[0].shape[0]).astype(bool)
        filters=[]
        filters.append((vals[0] >= pmin0 ) & (vals[0] <= pmax0))
        filters.append((vals[1] >= pmin1 ) & (vals[1] <= pmax1))
        if applyCuts is not None and self.Sels.has_key(applyCuts):
            filters.append(self.getFilter(applyCuts,plotvars))
        for ft in filters:
            total_filter&=ft                

        print 'select ',total_filter.sum(),' of ',np.ones_like(total_filter).sum(),' events'
        
        if fig is None:
            fig=plt.figure(figsize=(8,5))

        if not asHist:
            plt.plot(vals[1][total_filter],vals[0][total_filter],'o',markersize=3)
            plt.xlabel(plotvars[1])
            plt.ylabel(plotvars[0])
        else:
            v0 = vals[0][total_filter]
            v1 = vals[1][total_filter]
            binEdges0 = np.linspace(np.nanmin(v0),np.nanmax(v0),numBins[0])
            binEdges1 = np.linspace(np.nanmin(v1),np.nanmax(v1),numBins[1])
            ind0 = np.digitize(v0, binEdges0)
            ind1 = np.digitize(v1, binEdges1)
            ind2d = np.ravel_multi_index((ind0, ind1),(binEdges0.shape[0]+1, binEdges1.shape[0]+1)) 
            iSig = np.bincount(ind2d, minlength=(binEdges0.shape[0]+1)*(binEdges1.shape[0]+1)).reshape(binEdges0.shape[0]+1, binEdges1.shape[0]+1) 
            extent=[binEdges1[1],binEdges1[-1],binEdges0[1],binEdges0[-1]]
            plt.imshow(iSig,aspect='auto', interpolation='none',origin='lower',extent=extent,clim=[np.percentile(iSig,limits[0]),np.percentile(iSig,limits[1])])
            plt.xlabel(plotvars[1])
            plt.ylabel(plotvars[0])
        if setCuts is not None and self.Sels.has_key(setCuts):
            p =np.array(ginput(2))
            p0=[p[0][1],p[1][1]]
            p1=[p[0][0],p[1][0]]
            self.Sels[setCuts].addCut(plotvars[0],min(p0),max(p0))
            self.Sels[setCuts].addCut(plotvars[1],min(p1),max(p1))
        if asHist:
            return iSig, extent
        else:
            return vals[0][total_filter], vals[1][total_filter]

    def getScanName(self):
        for key in self.Keys('scan'):
            if key.find('var')<0 and key.find('none')<0 and key.find('damage')<=0:
                return key.replace('/scan/','').replace('scan/','')

    def getScanValues(self,ttCorr=False,addEnc=False):
        #get the scan variable & time correct if desired
        scanOrg = self.getVar('scan/var0')
        scanVarName = self.getScanName()
        if scanVarName.find('lxt')>=0 or scanVarName=='':
            delays=self.getDelay(use_ttCorr=ttCorr,addEnc=addEnc)
            scan = delays
            if scanVarName == '':
                if ttCorr:
                    scanVarName='delay (tt corrected) [fs]'
                else:
                    scanVarName='delay [fs]'
        else:
            scan = scanOrg
        return scanVarName,scan

    def getScans(self, runs=[], ttCorr=False, sig='', i0='', numBins=100, applyCuts=None, offData=True):
        plotData=[]
        currRun=self.run
        for run in runs:
            self.setRun(run)
            plotData.append(self.getScan(ttCorr=ttCorr,sig=sig,i0=i0,numBins=numBins,applyCuts=applyCuts))
        self.setRun(currRun)
        if not plotData[0].has_key('scanOffPoints'):
            offData=False
        for ids,dataSet in enumerate(plotData):
            if ids==0:
                scanPoints = dataSet['scanPoints']
                scan = dataSet['scan']
            elif len(dataSet['scanPoints']) != len(scanPoints) or (dataSet['scanPoints']-scanPoints).sum()>0.:
                print 'scanPoints not the same for all runs'
                return
            else:
                scan += dataSet['scan']
        if offData:
            for ids,dataSet in enumerate(plotData):
                if ids==0:
                    scanOffPoints = dataSet['scanOffPoints']
                    scanOff = dataSet['scanOff']
                elif len(dataSet['scanOffPoints']) != len(scanPoints) or (dataSet['scanOffPoints']-scanPoints).sum()>0.:
                    print 'scanOffPoints not the same for all runs'
                    print scanOffPoints
                    print dataSet['scanOffPoints']
                    offData=False
                    break
                else:
                    scanOff += dataSet['scanOff']

            if offData:
                return {'scanVarName':dataSet['scanVarName'],'scanPoints':scanPoints,'scan':scan, 'scanOffPoints':scanOffPoints,'scanoff':scanoff}
        
        return {'scanVarName':dataSet['scanVarName'],'scanPoints':scanPoints,'scan':scan}

    def getScan(self, ttCorr=False, sig='', i0='', Bins=100, applyCuts=None):
        return self.plotScan(ttCorr=ttCorr, sig=sig, i0=i0, Bins=Bins, returnData=True, applyCuts=applyCuts, plotThis=False)

    def plotScan(self, ttCorr=False, sig='', i0='', Bins=100, plotDiff=True, plotOff=True, saveFig=False,saveData=False, returnData=False, applyCuts=None, fig=None, interpolation='', plotThis=True, addEnc=False, returnIdx=False, binVar=None):
        
        plotVar=''
        if sig!='':
            sigVal = self.get1dVar(sig)
            for sigp in sig:
                if isinstance(sigp,basestring):
                    plotVar+=sigp.replace('/','_')
                elif isinstance(sigp,list):
                    for bound in sigp:
                        plotVar+='-%g'%bound
        else:
            print 'could not get signal variable, please specify'
            return

        if i0!='':
            i0Val = self.get1dVar(i0)
            plotVar+='/'
            for i0p in i0:
                if isinstance(i0p,basestring):
                    plotVar+=i0p.replace('/','_')
                elif isinstance(i0p,list):
                    for bound in i0p:
                        plotVar+='-%g'%bound
        else:
            print 'please specify normalizing variable '
            return
        
        [FilterOn, FilterOff] = self.getFilterLaser(applyCuts)
        FilterOn = FilterOn & ~np.isnan(i0Val) & ~np.isnan(sigVal)
        FilterOff = FilterOff & ~np.isnan(i0Val) & ~np.isnan(sigVal)

        if binVar is not None:
            if binVar[0] != 'delay':
                binVal = self.get1dVar(binVar[0])
            else:
                binVal=self.getDelay(use_ttCorr=ttCorr,addEnc=addEnc)
                ttCorr = None; addEnc=None
            FilterOn = FilterOn & ~np.isnan(binVal)
            FilterOff = FilterOff & ~np.isnan(binVal)

        print 'from %i events keep %i (%i off events)'%(np.ones_like(i0Val).sum(),np.ones_like(i0Val[FilterOn]).sum(), np.ones_like(i0Val[FilterOff]).sum() )

        #get the scan variable & time correct if desired
        scanVarName,scan =  self.getScanValues(ttCorr, addEnc)
            
        # create energy bins for plot: here no need to bin!
        if (not ttCorr) and (not addEnc):
            scanPoints, scanOnIdx = np.unique(scan[FilterOn], return_inverse=True)
        else:
            if isinstance(Bins, int) or isinstance(Bins, float):
                scanUnique = np.unique(scan[FilterOn])                
                if isinstance(Bins,int):
                    scanPoints = np.linspace(scanUnique[0],scanUnique[-1],Bins)
                elif isinstance(Bins,float):
                    if (abs(scanUnique[0]-scanUnique[-1])/Bins) > 1e5:
                        print 'this are more than 100k bins! will not try....'
                        return
                    scanPoints = np.arange(scanUnique[0],scanUnique[-1],Bins)
            elif isinstance(Bins,list) or isinstance(Bins,np.ndarray):
                scanPoints = Bins
            else:
                print 'Bins: ',isinstance(Bins,list),' -- ',Bins
            scanOnIdx = np.digitize(scan[FilterOn], scanPoints)
            scanPoints = np.concatenate([scanPoints, [scanPoints[-1]+(scanPoints[1]-scanPoints[0])]],0)
        OffData=False
        if scan[FilterOff].sum()!=0:
            scanOffPoints, scanOffIdx = np.unique(scan[FilterOff], return_inverse=True)
            if len(scanOffPoints) > len(scanPoints):
                scanOffPoints = scanPoints.copy()
                
            scanOffIdx = np.digitize(scan[FilterOff], scanOffPoints)
            OffData = True

        if returnIdx:
            return scanOnIdx

        if binVar is not None:
            if len(binVar)==1:
                nbin=100
            else:
                nbin=binVar[1]
            if len(binVar)<3:
                min = np.percentile(binVal,1)
                max = np.percentile(binVal,99)
            else:
                min = binVar[2]
                max = binVar[3]
            if isinstance(nbin, int):
                binPoints = np.linspace(min,max,nbin)
            elif isinstance(nbin, float):
                binPoints = np.arange(min,max,nbin)
            binIdx = np.digitize(binVal[FilterOn], binPoints)

            indOn2d = np.ravel_multi_index((scanOnIdx, binIdx),(scanPoints.shape[0]+1, binPoints.shape[0]+1)) 

            # calculate the normalized intensity for each bin
            iNorm = np.bincount(indOn2d, weights=i0Val[FilterOn], minlength=(scanPoints.shape[0]+1)*(binPoints.shape[0]+1)).reshape(scanPoints.shape[0]+1, binPoints.shape[0]+1)    
            iSig = np.bincount(indOn2d, weights=sigVal[FilterOn], minlength=(scanPoints.shape[0]+1)*(binPoints.shape[0]+1)).reshape(scanPoints.shape[0]+1, binPoints.shape[0]+1)    
        else:
            iNorm = np.bincount(scanOnIdx, i0Val[FilterOn])
            iSig = np.bincount(scanOnIdx, sigVal[FilterOn])

        #print 'evts ',np.bincount(scanOnIdx)
        #print 'i0',iNorm
        #print 'sig',iSig
        scan = iSig/iNorm
        #scan = scan/np.mean(scan[1]) # normalize to 1 for first energy point?

        if OffData:
            #same for off shots
            iNormoff = np.bincount(scanOffIdx, i0Val[FilterOff])
            iSigoff = np.bincount(scanOffIdx, sigVal[FilterOff])
            scanoff = iSigoff/iNormoff
            if scanOffIdx.shape > scanOffPoints.shape:
                shapeDiff = abs(scanOnIdx.max()+1-scanPoints.shape[0])
                shapeDiff = abs(scanOffIdx.max()+1-scanOffPoints.shape[0])
                if scanOffIdx.min()>0:
                    scanoff=scanoff[shapeDiff:]
                else:
                    scanoff=scanoff[:-shapeDiff]
        if (not OffData):
            plotDiff = False
        #now save data if desired
        if OffData:
            if saveData:
                np.savetxt('Scan_Run%i.txt'%self.run, (scanPoints, scan, scanOffPoints,scanoff))
            returnDict= {'scanVarName':scanVarName,'scanPoints':scanPoints,'scan':scan, 'scanOffPoints':scanOffPoints,'scanOff':scanoff,'plotVarName':plotVar}
        else:
            if saveData:
                np.savetxt('Scan_Run%i.txt'%self.run, (scanPoints, scan))
            returnDict= {'scanVarName':scanVarName,'scanPoints':scanPoints,'scan':scan,'plotVarName':plotVar}
        if binVar is not None:
            returnDict['binPoints']=binPoints
        if plotThis:
            #print 'retData: ',returnDict
            self.plotScanDict(returnDict, plotDiff=plotDiff, interpolation=interpolation,fig=fig,plotOff=plotOff,saveFig=saveFig)
        return returnDict

    def plotScanDict(self, returnDict, plotDiff=True, fig=None, plotOff=True, interpolation='', saveFig=False):
        plotVarName = returnDict['plotVarName']
        scanVarName = returnDict['scanVarName']
        scanPoints = returnDict['scanPoints']
        scan = returnDict['scan']
        print 'plot ',plotVarName, scanVarName, ' shape ',scan.shape,' plot diff ',plotDiff

        if interpolation!='' and returnDict.has_key('scanOffPoints'):
            finter_off = interpolate.interp1d(returnDict['scanOffPoints'], returnDict['scanOff'],kind=interpolation)
            scanoff_interp = finter_off(scanPoints[:-1])
        if len(scan.shape)>1:
            if fig is None:
                fig=plt.figure(figsize=(10,5))
            extent = [scanPoints.min(), scanPoints.max(), returnDict['binPoints'].min(), returnDict['binPoints'].max()]
            plt.imshow(scan, interpolation='none', aspect='auto', clim=[np.nanpercentile(scan,1), np.nanpercentile(scan,98)],extent=extent, origin='lower')
            plt.xlabel(scanVarName)
        #elif plotDiff and interpolation!='' and returnDict.has_key('scanOffPoints'):
        elif plotDiff and returnDict.has_key('scanOffPoints'):
            if fig is None:
                fig=plt.figure(figsize=(10,10))
            gs=gridspec.GridSpec(2,1,width_ratios=[1])
            plt.subplot(gs[0]).set_xlabel(scanVarName)
            plt.subplot(gs[0]).set_ylabel(plotVarName)
            plt.subplot(gs[0]).plot(scanPoints, scan, 'ro--', markersize=5)
            if interpolation!='':
                plt.subplot(gs[0]).plot(scanPoints[:-1], scanoff_interp, 'o--', markersize=5,markerfacecolor='none',markeredgecolor='b')
            plt.subplot(gs[0]).plot(returnDict['scanOffPoints'], returnDict['scanOff'], 'ko--', markersize=5)
            plt.subplot(gs[0]).set_ylim(np.nanmin(scan)*0.95,np.nanmax(scan)*1.05)
            if interpolation!='':
                plt.subplot(gs[1]).plot(scanPoints[:-1], (scan[:-1]-scanoff_interp), 'bo--', markersize=5)
            else:
                plt.subplot(gs[1]).plot(scanPoints[:-1], (scan[:-1]-returnDict['scanOff'][:-1]), 'bo--', markersize=5)
            plt.subplot(gs[1]).set_xlabel(scanVarName)
            plt.subplot(gs[1]).set_ylabel(plotVarName)
        else:
            if fig is None:
                fig=plt.figure(figsize=(10,5))
            plt.plot(scanPoints, scan, 'ro--', markersize=5)
            plt.xlabel(scanVarName)
            plt.ylabel(plotVarName)
            if returnDict.has_key('scanOffPoints') and plotOff:
                if interpolation!='':
                    plt.plot(scanPoints[:-1], (scanoff_interp), 'o--', markersize=5,markerfacecolor='none',markeredgecolor='b')
                plt.plot(returnDict['scanOffPoints'], returnDict['scanOff'], 'ko--', markersize=5)
            plt.ylim(np.nanmin(scan)*0.95,np.nanmax(scan)*1.05)
            plt.xlabel(scanVarName)
            plt.ylabel(plotVarName)
        if saveFig:
            fig.savefig('Scan_Run%i.jpg'%self.run)
            
    def defPlots(self, applyCuts=None):
        scanVarName,scan =  self.getScanValues(True)
        total_filter = np.ones_like(scan).astype(bool)
        if applyCuts is not None and self.Sels.has_key(applyCuts):
            total_filter =  self.getFilter(applyCuts, [plotvar])

        fig=plt.figure(figsize=(10,6))
        plt.title('Standard Plots for Run %i'%self.run)
        
        gs=gridspec.GridSpec(2,2,width_ratios=[2,2])
        self.plotVar('ipm2/sum',fig=plt.subplot(gs[0]),applyCuts=applyCuts)
        self.plotVar(['ipm2/sum','ebeam/L3Energy'],fig=plt.subplot(gs[1]),asHist=True,applyCuts=applyCuts)
        if len(scan)<200:
            pmin=scan[0]
            pmax=scan[-1]
        else:
            pmin = np.percentile(scan[total_filter],0.1)
            pmax = np.percentile(scan[total_filter],99.9)
        values = scan[total_filter]
        values = values[~np.isnan(values)]
        hst = np.histogram(scan[total_filter],np.linspace(pmin,pmax,100))
        plt.subplot(gs[2]).plot(hst[1][:-1],hst[0],'o')
        #plt.subplot(gs[2]).xlabel(scanVarName)
        #plt.subplot(gs[2]).ylabel('entries')
        plt.xlabel(scanVarName)
        plt.ylabel('entries')
        
        if self.hasKey(self.ttBaseStr+'AMPL'):
            if self.ttCorr is not None and np.nanstd(self.getVar(self.ttCorr))>0:
                self.plotVar(self.ttCorr,fig=plt.subplot(gs[3]))
            else:
                self.plotVar(self.ttBaseStr+'FLTPOS_PS',fig=plt.subplot(gs[3]))

    #########################################################
    ###
    ### functions for easy cube creation
    ###
    #########################################################

    #cube might be better to be its own class
    def addCube(self, cubeName, binVar='', bins=[], SelName=''):    
        self.cubes[cubeName] = Cube(binVar, bins, cubeName=cubeName, SelName=SelName)
        
    def addToCube(self, cubeName, targetVariable):
        if cubeName in self.cubes.keys():
            self.cubes[cubeName].addVar(targetVariable)
 
    def getCube(self, cubeName):
        if cubeName in self.cubes.keys():
            return self.cubes[cubeName]
        
    def printCubes(self, printDetail=True):
        cubeNames=[]
        if len(self.cubes.keys())>0:
            print 'defined cubes:'
            for cubeName in self.cubes.keys():
                cube = self.cubes[cubeName]
                if printDetail:
                    cube.printCube(self.Sels[cube.SelName])
                else:
                    print cubeName
                cubeNames.append(cubeName)
        return cubeNames

    def prepCubeData(self, cubeName):
        cube = self.getCube(cubeName)
        if not (self.hasKey(cube.binVar) or cube.binVar == 'delay'):
            print 'selection variable not in littleData, cannot make the data for this cube'
            return None

        # create the bins for the cube
        if len(cube.bins)>3:
            Bins = cube.bins
        else:
            Bins = util_getBins(cube.bins, self.fh5)
        cube.binBounds = Bins

        #now look through targetVars & split out ones not in xarray/hdf5
        targetVarsLocal = []
        for tVar in cube.targetVars:
            if isinstance(tVar, basestring):
                if tVar not in self._fields.keys():
                    cube.targetVarsXtc.append(tVar)
                else:
                    targetVarsLocal.append(tVar)
            else:
                cube.targetVarsXtc.append(tVar)
        cube.targetVars = targetVarsLocal

        #now get the filter & create a new one taking bins & detector damage into account.
        orgSel = cube.SelName
        if orgSel.find(cube.cubeName)!=0:
            self.Sels['%s_%s'%(cube.cubeName,cube.SelName)] = Selection()
            self.Sels['%s_%s'%(cube.cubeName,cube.SelName)].add(self.Sels[orgSel])
            self.Sels['%s_%s'%(cube.cubeName,cube.SelName)].addCut(cube.binVar, min(Bins), max(Bins) )
            #add cuts with detector damage - if we have damage detector info.
            for txVar in targetVarsLocal:
                if txVar[0]=='/':txVar=txVar[1:] 
                if 'damage/%s'%txVar  in self._fields.keys(): 
                    newSel.addCut('damage/%s'%txVar.split('/')[0],0.5,1.5)
            for txVar in cube.targetVarsXtc:
                if isinstance(txVar, dict):
                    try:
                        txVar=txVar['source']
                    except:
                        continue
                if 'damage/%s'%txVar  in self._fields.keys(): 
                        newSel.addCut('damage/%s'%txVar,0.5,1.5)
            cube.SelName='%s_%s'%(cube.cubeName,cube.SelName)

        return cube

    def makeCubeData(self, cubeName, debug=False, toHdf5=False, replaceNan=False, onoff=2, returnIdx=False, addIdxVar=''):
        cube = self.prepCubeData(cubeName)
        if cube is None:
            return 
        Bins = cube.binBounds
        numBin=len(Bins)-1

        if cube.binVar == 'delay':
            binVar = self.getDelay()
        else:
            binVar = self.get1dVar(cube.binVar)

        if debug:
            cube.printCube(self.Sels[cube.SelName])
        cubeFilter = self.getFilter(cube.SelName)
        [cubeOn, cubeOff] = self.getFilterLaser(cube.SelName, ignoreVar=[])
        if onoff==1:
            cubeFilter = cubeOn
            cubeName = cubeName+'_laserOn'
        elif onoff==0:
            cubeFilter = cubeOff
            cubeName = cubeName+'_laserOff'

        binVar = binVar[cubeFilter]
        if binVar.shape[0] == 0:
            print 'did not select any event, quit now!'
            return

        if debug:
            print 'bin boundaries: ',Bins

        timeFiltered = self._tStamp[cubeFilter]
        newXr = xr.DataArray(np.ones(timeFiltered.shape[0]), coords={'time': timeFiltered}, dims=('time'),name='nEntries')
        newXr = xr.merge([newXr, xr.DataArray(binVar, coords={'time': timeFiltered}, dims=('time'),name='binVar') ])       
        for tVar in cube.targetVars:
            if not self.hasKey(tVar):                
                continue
            #print 'addvar: ',tVar,self.getVar(tVar,cubeFilter).shape
            filteredVar = self.getVar(tVar,cubeFilter).squeeze()
            tVar=tVar.replace('/','__')
            if len(filteredVar.shape)==1:
                newXr = xr.merge([newXr, xr.DataArray(filteredVar, coords={'time': timeFiltered}, dims=('time'),name=tVar) ])
            else:
                coords={'time': timeFiltered}
                dims = ['time']
                dataShape = filteredVar.shape
                for dim in range(len(dataShape)-1):
                    thisDim = np.arange(0, dataShape[dim+1])
                    dimStr = '%s_dim%d'%(tVar,dim)
                    coords[dimStr] = thisDim
                    dims.append(dimStr)
                newXr = xr.merge([newXr, xr.DataArray(filteredVar, coords=coords, dims=dims,name=tVar)])

        cubeData = newXr.groupby_bins('binVar',Bins,labels=(Bins[1:]+Bins[:-1])*0.5).sum(dim='time')                  
        #could add error using the std of the values.
        cubeDataErr = newXr.groupby_bins('binVar',Bins,labels=(Bins[1:]+Bins[:-1])*0.5).std(dim='time')
        for key in cubeDataErr:
            if key.replace('std_','').replace('__','/') in cube.targetVars:
                cubeDataErr.rename({key: 'std_%s'%key}, inplace=True)
        for key in cubeDataErr:
            if key not in cubeData.keys():
                cubeData = xr.merge([cubeData, cubeDataErr[key]])

        if toHdf5:
            fname = 'Cube_%s_Run%03d.nc'%(cubeName, self.run)
            cubeData.to_netcdf(fname,engine='h5netcdf')
        
        if not returnIdx:
            return cubeData

        fidVar='fiducials'
        evttVar='event_time'
        if '/fiducials' in self.Keys():
            fidVar='/fiducials'
            evttVar='/event_time'
        elif 'EvtID/fid' in self.Keys():
            fidVar='EvtID/fid'
            evttVar='EvtID/time'
        print 'we will use fiducials from here: ',fidVar

        evtIDXr = xr.DataArray(self.getVar(fidVar,cubeFilter), coords={'time': timeFiltered}, dims=('time'),name='fiducial')
        evtIDXr = xr.merge([evtIDXr,xr.DataArray(self.getVar(evttVar,cubeFilter), coords={'time': timeFiltered}, dims=('time'),name='evttime')])
        evtIDXr = xr.merge([evtIDXr, xr.DataArray(binVar, coords={'time': timeFiltered}, dims=('time'),name='binVar') ])       
        if addIdxVar!='':
            evtIDXr = xr.merge([evtIDXr, xr.DataArray(self.getVar(addIdxVar,cubeFilter), coords={'time': timeFiltered}, dims=('time'),name=addIdxVar) ])       
        #else:
        #    print 'could not find event idx in data'
        #    return cubeData,None
        cubeIdxData = evtIDXr.groupby_bins('binVar',Bins)   
        keys = cubeIdxData.groups.keys()
        keys.sort()

        fidArray=[]
        timeArray=[]
        addArray=[]
        for key in keys:
            fidArray.append(evtIDXr.fiducial[cubeIdxData.groups[key]])
            timeArray.append(evtIDXr.evttime[cubeIdxData.groups[key]])
            if addIdxVar!='':
                addArray.append(evtIDXr[addIdxVar][cubeIdxData.groups[key]])
        retDict={'keys': keys}
        retDict['fiducial']=fidArray
        retDict['evttime']=timeArray
        if addIdxVar!='':
            retDict[addIdxVar]=addArray
        return cubeData,retDict

    #CHECK ME: REWRITE ACCESS TO NON-SMALL data....
    def submitCube(self, cubeName, run=None, expname=None, image=False, rebin=-1):                                    
        if socket.gethostname().find('psana')<0:
            print 'we can only submit jobs from psana, you are on ',socket.gethostname()
            return

        if expname == None:
            expname = self.expname
        if run == None:
            run = self.run
            
        #check if cube file exists
        haveFile = True
        if not path.isfile('CubeSetup_'+cubeName+'.txt'):
            haveFile=False
            for cube in self.cubes:
                if cube.cubeName == cubeName:
                    if raw_input('this cube has not been written yet, do it?(y/n)') == 'y':
                        self.writeCubeSetup(cubeName)
                        haveFile=True
        if not haveFile:
            print 'cube %s has not been defined yet'%cubeName
            return

        if path.isfile('./cubeRun'):
            cmd = './cubeRun -r %i -e %s -c %s'%(run, expname, 'CubeSetup_'+cubeName+'.txt')
        else:
            cmd = '/reg/d/psdm/xpp/%s/res/littleData/xppmodules/scripts/cubeRun -r %i -e %s -c %s'%(expname, run, expname, 'CubeSetup_'+cubeName+'.txt')
        if image:
            cmd+=' -i'
        if rebin>0:
            cmd+=' -R %i'%rebin
        print 'command is: %s'%cmd
        sout = subprocess.check_output([cmd],shell=True)
        print 'and submission returned %s'%sout
        self.jobIds.append(sout.split('Job <')[1].split('> is')[0])

    def checkJobs(self):
        remainingIds=[]
        cmd = "bjobs | awk {'print $1'} | grep -v JOBID | grep -v psana"
        cout = subprocess.check_output([cmd],shell=True)
        print cout
        for jobid in self.jobIds:
            if cout.find(jobid)>=0:
                remainingIds.append(jobid)
            else:
                print 'Job %s is done'%jobid
        self.jobIds = remainingIds

    ##########################################################################
    ###
    ### functions for image treatment - starting w/ assembled 2-d image
    ###
    ##########################################################################

    def AvImage(self, detname='None', numEvts=100, nSkip=0, thresADU=0., thresRms=0.,applyCuts=None, mean=False, std=False):
        #look for detector
        if detname=='None':
            aliases=self.Keys2d()
            if len(aliases)<1:
                print 'no area detectors in littleData, all keys present are:'
                self.Keys(printKeys=True)
            if len(aliases)==1:
                detname = aliases[0]
            else:
                print 'detectors in event: \n',
                for alias in aliases:
                    print alias
                detname = raw_input("Select detector to select ROI of?:\n")
        print 'we are looking at ',detname

        #arrays useful for thresholding
        detsrc = detname.split('/')[0]
        roi = self.getVar('UserDataCfg/%s_bound'%(detname.replace('/','_')))
        try:
            rmsFull = self.getVar('UserDataCfg/%s_rms'%detsrc)
            maskFull = self.getVar('UserDataCfg/%s_mask'%detsrc)
            rms = rmsFull[roi[0,0]:roi[0,1], roi[1,0]:roi[1,1], roi[2,0]:roi[2,1]].squeeze()
            mask = maskFull[roi[0,0]:roi[0,1], roi[1,0]:roi[1,1], roi[2,0]:roi[2,1]].squeeze()
        except:
            rms=None
            mask = None

        #only events requested
        if applyCuts is not None:
            Filter = self.getFilter(SelName=applyCuts)
            dataAr = self.getVar(detname,Filter)
            dataAr = dataAr[nSkip:min(nSkip+numEvts, dataAr.shape[0])].squeeze()
        else:
            #now look at subarray
            dataAr = self.getVar(detname)[nSkip:min(nSkip+numEvts, self._tStamp.shape[0])].squeeze()

        #now apply threshold is requested:
        data='AvImg_'
        if std:
            thresDat = dataAr.mean(axis=0)
            data+='std_'
        elif mean:
            thresDat = dataAr.std(axis=0)
            data+='mean_'
        else:
            thresDat = np.zeros_like(dataAr[0])
            for shot in dataAr:
                if thresADU != 0:
                    shot[shot<abs(thresADU)]=0
                    #shot[shot>=abs(thresADU)]=1
                if thresRms > 0 and rms is not None:
                    shot[shot<thresRms*rms]=0
                    #shot[shot>=thresRms*rms]=1
                thresDat += shot

        if thresADU!=0:
            data+='thresADU%d_'%int(thresADU)
        if thresRms!=0:
            data+='thresRms%d_'%int(thresRms)
        data+=detname.replace('/','_')
        self.__dict__[data]=thresDat
        
    def getAvImage(self,detname=None, ROI=[]):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('AvImg')>=0:
                if detname is not None and key.find(detname)>=0:
                    avImages.append(key)
                elif detname is None:
                    avImages.append(key)
        if len(avImages)==0:
            print 'creating the AvImage first!'
            return
        elif len(avImages)>1:
            print 'we have the following options: ',avImages
            avImage=raw_input('type the name of the AvImage to use:')
        else:
            avImage=avImages[0]
        detname = avImage.replace('AvImg_','')
        img = self.__dict__[avImage]
        return img
        
    def plotAvImage(self,detname=None, ROI=[],limits=[5,99.5]):
        img = self.getAvImage(detname=detname, ROI=ROI)
        print img.shape

        plotMax = np.percentile(img, limits[1])
        plotMin = np.percentile(img, limits[0])
        print 'using the percentiles %g/%g as plot min/max: (%g, %g)'%(limits[0],limits[1],plotMin,plotMax)

        image = img
        
        fig=plt.figure(figsize=(10,6))
        if ROI!=[]:
            gs=gridspec.GridSpec(1,2,width_ratios=[2,1])        
            plt.subplot(gs[1]).imshow(img[ROI[0][0],ROI[1][0]:ROI[1][1],ROI[2][0]:ROI[2][1]],clim=[plotMin,plotMax],interpolation='None')
        else:
            gs=gridspec.GridSpec(1,2,width_ratios=[99,1])        
        plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')

        plt.show()
        
    def getPeakAvImage(self,detname=None, ROI=[]):
        img=self.getAvImage(detname=detname, ROI=ROI)
        return [img.max(), img.mean(axis=0).argmax(),img.mean(axis=1).argmax()]

    def FitCircle(self, detname=None, mask=None, method=None, thres=None):
        try:
            from utilities import fitCircle
        except:
            print 'could not import underlying code, this does not work yet'
            return
        print 'nearly done, but there is an issue in the image display and x/y coordinates that needs figuring out with faster x-respose.....'
        #return
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('AvImg')>=0:
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
        detname = avImage.split('_')[1]
        print 'detname: ',detname,avImage
        image = self.__dict__[avImage]
        if len(image.shape)!=2:
            print 'not a 2-d image! Will return. image %s has %d dims'%(avImage,len(image.shape))
            return

        plotMax = np.percentile(image, 99.5)
        plotMin = np.percentile(image, 5)
        print 'using the 5/99.5% as plot min/max: (',plotMin,',',plotMax,')'

        if mask:
            image = (image*mask)

        #get x & y array from data to get extent
        x = self.getVar('UserDataCfg/%s_x'%detname)
        y = self.getVar('UserDataCfg/%s_y'%detname)
        #get the ROI bounds
        if len(avImage.split('_'))>2:
            roiname = avImage.split('_')[2]
            ROI = self.getVar('UserDataCfg/%s_%s_bound'%(detname,roiname))
            if len(ROI)==2:
                x = x[ROI[0,0]:ROI[0,1], ROI[1,0]:ROI[1,1]].squeeze()
                y = y[ROI[0,0]:ROI[0,1], ROI[1,0]:ROI[1,1]].squeeze()
            elif len(ROI)==3:
                x = x[ROI[0,0]:ROI[0,1], ROI[1,0]:ROI[1,1], ROI[2,0]:ROI[2,1]].squeeze()
                y = y[ROI[0,0]:ROI[0,1], ROI[1,0]:ROI[1,1], ROI[2,0]:ROI[2,1]].squeeze()
        extent=[x.min(), x.max(), y.min(), y.max()]

        fig=plt.figure(figsize=(10,10))
        plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None',aspect='auto')

        if method == None:
            method=raw_input("Select circle points by mouse or threshold [m/t]?:\n")
        if method not in ["m","t"]:
            method=raw_input("Please select m or p (mouse/threshold) or we will return\n")
            if method not in ["m","t"]:
                return

        if method=="m":
            happy = False
            while not happy:
                points=ginput(n=0)
                parr=np.array(points)
                #res: xc, yc, R, residu
                res = fitCircle(parr[:,0],parr[:,1])
                #draw the circle.
                circle = plt.Circle((res[0],res[1]),res[2],color='b',fill=False)
                plt.gca().add_artist(circle)
                plt.plot([res[0],res[0]],[y.min(),y.max()],'r')
                plt.plot([x.min(),x.max()],[res[1],res[1]],'r')

                if raw_input("Happy with this fit:\n") in ["y","Y"]:
                    happy = True
                print 'x,y: ',res[0],res[1],' R ',res[2]
                print 'avs: ',parr[:,0].mean(),parr[:,1].mean()
        else:
            happy = False
            while not happy:
                if thres is None:
                    thres = raw_input("percentile in % of selected points min[,max]:\n")
                if thres.find(',')>=0:
                    thresMin=float(thres.split(',')[0])
                    thresMax=np.percentile(image, float(thres.split(',')[1]))
                else:
                    thresMin=float(thres.split(',')[0])
                    thresMax=np.nanmax(image)
                thresP = np.percentile(image, float(thresMin))
                print 'thresP',thresP
                imageThres=image.copy()
                imageThres[image>thresP]=1
                imageThres[image<thresP]=0
                imageThres[image>thresMax]=0
                fig=plt.figure(figsize=(5,5))
                plt.imshow(imageThres,clim=[-0.1,1.1],extent=extent,aspect='auto')
                if thres is None:
                    if raw_input("Happy with this threshold (y/n):\n") in ["y","Y"]:
                        happy=True
                else:
                    happy=True

            #res = fitCircle(x.flatten()[image.flatten()>thresP],y.flatten()[image.flatten()>thresP])
            res = fitCircle(x.flatten()[imageThres.flatten().astype(bool)],y.flatten()[imageThres.flatten().astype(bool)])
            print 'res',res
            print 'x,y av: ',(x.flatten()[imageThres.flatten().astype(bool)]).mean(),(y.flatten()[imageThres.flatten().astype(bool)].mean())
            circleM = plt.Circle((res[0],res[1]),res[2],color='b',fill=False)
            fig=plt.figure(figsize=(10,10))
            #will need to check of rotation necessary here???
            #plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
            plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None',aspect='auto')
            plt.gca().add_artist(circleM)
            plt.plot([res[0],res[0]],[y.min(),y.max()],'r')
            plt.plot([x.min(),x.max()],[res[1],res[1]],'r')
            print 'x,y: ',res[0],res[1],' R ',res[2]
    
            plt.show()

    def MakeMask(self, detname=None):
        print ' not yet implemented, exists in LittleDataAna_psana.py'

    def azimuthalBinning(self, detname=None):
        print ' not yet implemented, exists in LittleDataProduced.py, uses code in xppmodules/src. Not sure if good idea'

    ##########################################################################
    ###
    ### functions for droplet analysis
    ###
    ##########################################################################

    def DropletCube(self, applyCuts='', i0='ipm3/sum', rangeAdu=[], rangeX=[], rangeY=[], addName='', returnData=False, writeFile=False):
        data='DropletCube_Run%d_%s'%(self.run,addName)
        if applyCuts!='':
            data+='_'+applyCuts
        self.__dict__[data]=None

        #get basename of droplets
        dkey = [ key for key in self.Keys() if key.find('dropletsAdu')>=0]
        if len(dkey)==0:
            print 'did not find any droplets in this littleData file: ',self.fh5.filename
            return
        if len(dkey)>1:
            print 'we have the following options: ',dbkey
            basename=raw_input('type the name of the droplets to use:')
        else:
            basename = dkey[0].split('dropletsAdu')[0]

        #get filtered list of events
        i0_all = self.getVar(i0)
        if applyCuts is not None:
            Filter = self.getFilter(SelName=applyCuts)
        else:
            Filter = np.ones_like(i0_all)

            dataAr = self.getVar(detname,Filter)
            #dataAr = self.fh5[detname].value.squeeze()[Filter]
            dataAr = dataAr[nSkip:min(nSkip+numEvts, dataAr.shape[0])].squeeze()
        
        #get sum of i0
        i0_sum = i0_all[Filter].sum().astype(float)

        #get all droplets in selected events, ADU>0
        adu = self.getVar(basename+'dropletsAdu',Filter=Filter).flatten()
        x = self.getVar(basename+'dropletsX',Filter=Filter).flatten()[adu>0]
        y = self.getVar(basename+'dropletsY',Filter=Filter).flatten()[adu>0]
        adu=adu[adu>0]
        #adu = self.fh5[basename+'dropletsAdu'][Filter,:].flatten()[self.fh5[basename+'dropletsAdu'][Filter,:].flatten()>0]
        #x = self.fh5[basename+'dropletsX'][Filter,:].flatten()[self.fh5[basename+'dropletsAdu'][Filter,:].flatten()>0]
        #y = self.fh5[basename+'dropletsY'][Filter,:].flatten()[self.fh5[basename+'dropletsAdu'][Filter,:].flatten()>0]
        

        #make 3d histo of ADU, X, Y
        if rangeAdu==[]:
            rangeAdu=[np.percentile(adu,1), np.percentile(adu,99.9)]
        if rangeX==[]:
            rangeX=[np.percentile(x,1), np.percentile(x,99.9)]
        if rangeY==[]:
            rangeY=[np.percentile(y,1), np.percentile(y,99.9)]

        binAdu = np.arange(int(rangeAdu[0]), int(rangeAdu[1])).astype(int)
        binX = np.arange(int(rangeX[0]), int(rangeX[1])).astype(int)
        binY = np.arange(int(rangeY[0]), int(rangeY[1])).astype(int)
            
        indA = np.digitize(adu, binAdu)
        indX = np.digitize(x, binX)
        indY = np.digitize(y, binY)
        ind3d = np.ravel_multi_index((indA, indX, indY),(binAdu.shape[0]+1, binX.shape[0]+1, binY.shape[0]+1)) 
        cube = np.bincount(ind3d, minlength=(binAdu.shape[0]+1)*(binX.shape[0]+1)*(binY.shape[0]+1)).reshape(binAdu.shape[0]+1, binX.shape[0]+1, binY.shape[0]+1)

        returnDict= {'i0_sum':i0_sum, 'binAdu':binAdu.tolist(), 'binX':binX.tolist(), 'binY':binY.tolist(),'cube':cube.tolist()}
        self.__dict__[data]=returnDict        

        if writeFile:
            f = open(data+'.txt','w')
            print 'write DropletCube file for ',data, ' to ',data,'.txt'
            #indent does NOT work here...
            #json.dump(returnDict,f,indent=2)
            json.dump(returnDict,f)
            f.close()

        if returnData:
            return returnDict

