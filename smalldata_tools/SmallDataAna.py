"""
Created on Tue Dec  8 21:31:56 2015

@author: snelson
"""
from os import makedirs
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
from utilities import dictToHdf5, shapeFromKey_h5
from utilities import hist2d
from utilities import running_median_insort
from utilities import get_startOffIdx, get_offVar
from utilities import getBins as util_getBins
from utilities import printR
from utilities_plotting import plotImageBokeh, plotMarker, plotImage
import bokeh
import bokeh.plotting as bp
from bokeh.models import WheelZoomTool, BoxZoomTool, Div
from bokeh.models import PanTool, SaveTool, HoverTool, ResetTool
try:
    from bokeh.models import ResizeTool
except:
    pass
from bokeh.layouts import column
import sys
from epicsarchive import EpicsArchive

from pscache import client #works from ana-1.2.9 on

#including xarray means that you have to unset DISPLAY when submitting stuff to batch
import xarray as xr

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

 
class Cube(object):
    """ class to define binning operation for small/big data in a given run"""
    def __init__(self, binVar='', bins=[], cubeName=None, useFilter=None, addBinVars=None):
        self.binVar = binVar
        self.bins = bins
        self.useFilter = useFilter

        nbins = len(bins)
        if nbins==3:
            if type(bins[2]) is int:
                nbins=len(np.linspace(min(bins[0],bins[1]),max(bins[0],bins[1]),bins[2]))
            else:
                nbins=len(np.arange(min(bins[0],bins[1]),max(bins[0],bins[1]),bins[2]))

        self.add_BinVar(addBinVars)

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
        self.addIdxVars = []
        self.dropletProc = {} #for droplet treatment in cube

    def add_BinVar(self, addBinVars):
        """
        add extra dimensions to bin the data in
        parameters: addBinVars: dict or list
                    list: only 1 extra variable [varname, bin1, bin2, ....]
                    dict: {varname: bins}
        """
        if addBinVars is None:
            addBinVars={}
        if isinstance(addBinVars,list):
            if isinstance(addBinVars[1], list):
                inputBins = addBinVars[1]
            else:
                inputBins = addBinVars[1:]
            addbins = util_getBins(inputBins)
            if addbins is None:
                print 'Could not define bin boundaries for addBinVars:',inputBins
                print 'automatic detection only works for primary binning axis'
            addBinVars = { addBinVars[0]: addbins}
        elif isinstance(addBinVars,dict):
            for addVar in addBinVars.keys():
                addbins = util_getBins(addBinVars[addVar])
                if addbins is None:
                    print 'Could not define bin boundaries for addBinVars:',addBinVars[addVar]
                    print 'automatic detection only works for primary binning axis'
                addBinVars[addVar] = addbins
        else:
            print 'please pass addBinVar as list: [ varName, bins ]'

        self.addBinVars = addBinVars
        return

    def addVar(self, tVar):
        """
        add variables that should be binned
        parameters: tVar (either string of list of strings)

        variables are either variables present in the smallData hdf5 file
                  for droplet binning, specify if you want an image/bin or the X/Y/adu array
                  e.g. 'droplet:epix_2/droplet:image' 
                  or detector aliases for data present only in the xtc file
                  for the latter, the cube needs to be created using the 
                           SmallDataAna_psana class
                  example: epixDict = {'source':'epix_2','full':1,'image':1}
        """
        if isinstance(tVar, basestring):
            self.targetVars.append(tVar)
        elif isinstance(tVar, list):
            for tv in tVar:
                self.targetVars.append(tv)
        elif isinstance(tVar, dict):
            self.targetVars.append(tVar)

    def addIdxVar(self, tVar):
        """
        parameter: tVar as string or list of strings
        add variables that you would like to have returned as lists for each bin 
            rather than the sum. Will be returned as extra dictionary in addition 
            to the binned xArray
        """
        if isinstance(tVar, basestring):
            self.addIdxVars.append(tVar)
        elif isinstance(tVar, list):
            for tv in tVar:
                self.addIdxVars.append(tv)
        else:
            print 'addIdxVar need to be a string or list of strings, got: ',tVar

    def printCube(self, Sel=None):
        """
        parameter: name of filter/selection
        prints the setting of this cube.
        """
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
        if Sel is not None:
            for icut,cut in enumerate(Sel.cuts):
                print 'Cut %i: %f < %s < %f'%(icut, cut[1], cut[0],cut[2])
        print 'we will bin: ',self.targetVars


class Selection(object):
    """ class to define a subset of events using simple square cuts"""
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
        """
        add a variable to the selection
        parameters: varName, varmin, varmax
        varName: variable name as in smallData hdf5 file (string)
        varmin: minimum of variable varName that passes
        varmax: maximum of variable varName that passes
        """
        self.removeCut(varName)
        self.cuts.append([varName, varmin, varmax])
        self._filter=None
    def removeCut(self, varName):
        """
        remove a variable from the selection
        parameters: varName
        """
        for cut in self.cuts:
            if cut[0]==varName: self.cuts.remove(cut)
        self._filter=None
    def printCuts(self):
        """ print the currently defined list of square cuts for selection/filter"""
        for icut,cut in enumerate(self.cuts):
            print 'Cut %i: %f < %s < %f'%(icut, cut[1], cut[0],cut[2])
        if isinstance(self._filter, np.ndarray):
            print 'of %d events %d pass filter'%(self._filter.shape[0], self._filter.sum())
    def add(self, additionalSel):
        """ add all cuts defined in a different filter to this one"""
        for cut in additionalSel.cuts:
            self.cuts.append(cut)
        self._filter=None

class SmallDataAna(object):
    """ 
    class to deal with data in smallData hdf5 file. 
    Most functions assume standard variables to be present
    """
    def __init__(self, expname='', run=-1, dirname='', filename='',intable=None, liveList=None, plotWith='matplotlib'):
        self._fields={}
        self._live_fields=[]
        if isinstance(liveList, list) and intable=='redis':
            self._live_fields=liveList
        self.expname=expname
        self.run=run
        self.runLabel='Run%03d'%run
        self.plotWith=plotWith
        if len(expname)>3:
            self.hutch=self.expname[:3]
            self.plot_dirname='/reg/d/psdm/%s/%s/results/smalldata_plots/'%(self.hutch,self.expname)
            if dirname=='':
                self.dirname='/reg/d/psdm/%s/%s/hdf5/smalldata'%(self.hutch,self.expname)
                #run 13 and past.
                if not path.isdir(self.dirname):
                    self.dirname='/reg/d/psdm/%s/%s/ftc'%(self.hutch,self.expname)
                    self.plot_dirname='/reg/d/psdm/%s/%s/res/smalldata_plots/'%(self.hutch,self.expname)
            else:
                self.dirname=dirname
                self.plot_dirname = dirname+'/smalldata_plots'
            if not path.isdir(self.plot_dirname):
                makedirs(self.plot_dirname)

        if filename == '':
            self.fname='%s/%s_Run%03d.h5'%(self.dirname,self.expname,self.run)
            print 'setting up dirs:'
            if not path.isdir('/reg/d/psdm/%s/%s/results/'%(self.hutch,self.expname)):
                self.fname='%s/ldat_%s_Run%03d.h5'%(self.dirname,self.expname,self.run)
            if not path.isfile(self.fname):
                self.dirname='/reg/d/psdm/%s/%s/results/arphdf5'%(self.hutch,self.expname)
                self.fname='%s/%s_Run%03d.h5'%(self.dirname,self.expname,self.run)
        else:
            self.fname='%s/%s'%(self.dirname,filename)
        if intable != 'redis':
            print 'and now open in dir: ',self.dirname,' to open file ',self.fname

        self.Sels = {}
        self.cubes = {}
        self.jobIds=[]
        self._isRedis=False
        if run == -1 or (intable is not None and intable == 'redis'):
            self.fh5=client.ExptClient(expname, host='psdb3')
            self._isRedis=True
            plot_dirname='/reg/neh/operator/%sopr/experiments/%s/smalldata_plots/'%(self.hutch,self.expname)
            if not path.isdir(plot_dirname):
                try:
                    makedirs(plot_dirname)
                    self.plot_dirname = plot_dirname
                except:
                    pass
            try:
                import RegDB.experiment_info
                currRun=RegDB.experiment_info.experiment_runs(self.expname[:3].upper)[-1]['num']
                self.runLabel='Run%03d'%currRun
            except:
                try:
                    runs = self.fh5.runs()
                    lastRun=-1
                    while len(runs)>0:
                        thisRun = runs.pop()
                        lastRun = max(thisRun, lastRun)
                    self.runLabel='Run%03d'%int(lastRun)
                except:
                    self.runLabel='RunFromRedis'
        elif intable is not None:
            if intable == 'redis':
                self.fh5=client.ExptClient(expname, host='psdb3')
                self._isRedis=True
                plot_dirname='/reg/neh/operator/%sopr/experiments/%s/smalldata_plots/'%(self.hutch,self.expname)
                if not path.isdir(plot_dirname):
                    try:
                        makedirs(plot_dirname)
                        self.plot_dirname = plot_dirname
                    except:
                        pass
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
                if len(evttime)<1:
                    print 'read from REDIS, found no entries for run ',self.run,', cannot make SmallData object'
                    return None
                print 'read data for run %d from REDIS, have %d events'%(self.run,len(evttime))
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

        self._delay_ttCorr=None
        self._delay_addLxt=None
        self._delay_addEnc=None
        #define bokeh palette based on matplotlib job for consistency of matplotlib & bokeh image plots
        import matplotlib.cm as cm
        import matplotlib as pltm
        colormap =cm.get_cmap("jet")
        self.bokehpalette = [pltm.colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]
        self._epicsArchive=None

    def __del__(self):
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
        """ helper function to get the shape of data from hdf5 file"""
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
                dimArr = ['%02d'%i for i in range(0,dataShape[1])]
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
        #print 'DEBUG ',keys_for_xarray
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

###
# function to deal with extra variables added to smallData
###
    def _addXarray(self):
        """  add xarray object to main class instance as xrData """
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
                    try:
                        self.xrData = xr.merge([self.xrData, xr.DataArray(self.fh5.get_node(key).read(), coords={'time': self._tStamp}, dims=('time'),name=key[1:].replace('/','__')) ])
                        self._fields[fieldkey][1]='inXr'
                    except:
                        print 'failed to create dataset for: ',fieldkey, self._fields[fieldkey]
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
                        if len(dataShape)==2 and dataShape[1]>25 and (dataShape[0]*dataShape[1])>1e7:
                            print 'array for %s is too large (shape %d * %d), will load when asked'%(fieldName,dataShape[0], dataShape[1])
                            continue
                        self.xrData = xr.merge([self.xrData, xr.DataArray(self.fh5.get_node(key,tleaf.name).read().squeeze(), coords=coords, dims=dims,name=tArrName) ])
                        self._fields[fieldName][1]='inXr'

    def addVar(self, name='newVar',data=[]):
        """  add new variables to xrData. Keep track so that it will be saved on 
        destuction of main object"""
        if name.find('__')>=0 and name.replace('__','/') not in self._fields.keys():
            print 'Names of newly defined variables may not contain __, will replace with single _'
            name = name.replace('__','_')
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

        #if not self._isRedis and name[0]!='/': name='/'+name
        name = name.replace('__','/')
        if name not in self._fields.keys():
            #print 'DEBUG: add a new variable to Xarray: ',name
            self._fields[name]=[data.shape, 'inXr', 'mem']
        elif self._fields[name][2]=='main':
            #print 'add a variable from the main data to Xarray: ',name
            self._fields[name]=[data.shape, 'inXr', 'main']
        elif  self._fields[name][2]=='xrfile':
            #print 'add a variable from an netcdf to Xarray: ',name
            self._fields[name]=[data.shape, 'inXr', 'xrfile']
            try:
                testShape = self.xrData[name].data.shape
                if (self.xrData[name].data==data).sum()!=data.shape[0]:
                    self.xrData[name].data = data
                    #reset this to memory so that file will get overwritten
                    self._fields[name]=[data.shape, 'inXr', 'mem']
                return
            except:
                pass
        else:
            #print 'try to add a variable already present in xrData, only replace values!'
            self.xrData[name].data = data
            return

        #create new xrData to be merged
        name = name.replace('/','__')
        data = data.squeeze()
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
        write newly created fields to netcdf files that can be loaded in future sessions
        """
        #do not write files for REDIS (as you are likely auto-updating...)
        if self._isRedis:
            return
        print 'save derived data to be loaded next time:'
        for key in self._fields.keys():
            #delay is special: make sure it gets redefined on each new SmallDataAna creation 
            #to allow to use different definitions. getDelay will be called internally w/ default params 
            #so it needs to return the current defined variable first if applicable
            #saving does not help performace this is a cheap calculation.
            if key=='delay': continue 
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
        """ read data previously saved as netcdf files for this run"""
        for (dirpath, dirnames, filenames) in walk(self.dirname):        
            for fname in filenames:
                if fname.find('xr')==0 and fname.find('Run%03d'%self.run)>=0:
                    try:
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
                        add_xrDataSet.close()
                        print 'closed the xarray file '
                    except:
                        print 'failed at xr.open_dataset for: ',self.dirname+'/'+fname

###
# functions to add extra variables to smallData
###
    def addMedianVar(self, name='newVar', windowsize=31):
        """
        add a variable containing the median of the input variable to the data.
        parameters: name='newVar', windowsize=31
                    name: variable name as already present in data
                    windowsize: number of events prior to current used to calc median
        """
        dataOrg = self.getVar(name)
        dataNew=running_median_insort(dataOrg, windowsize=windowsize)
        medVarName = ('median_%s'%name).replace('/','_')
        print 'add variable named: ',medVarName
        self.addVar(medVarName, dataNew)

    def _getStartOffIdx(self, selName, nNbr=3, overWrite=False):
        """ function to get start of the indices for neighboring off events"""
        if '%s_offIdx_nNbr%02d'%(selName, nNbr) in self.Keys() and not overWrite:
            startOffIdx = self.getVar('%s_offIdx_nNbr%02d'%(selName,nNbr))
        else:
            tStamp = self.xrData.event_time.values
            filterOff = self.getFilter(selName.split('__')[0]+'__off')
            startOffIdx = get_startOffIdx(tStamp, filterOff, nNbr=nNbr)
            print 'add offIdx to data: ',('%s_offIdx_nNbr%02d'%(selName, nNbr))
            self.addVar('%s_offIdx_nNbr%02d'%(selName, nNbr), startOffIdx)
        return startOffIdx

    def getOffVar(self, varName, selName, nNbr=3, mean=True, returnMe=False):
        """
        function to add data from neighboring off events for each (on) event.
        parameters: varName, selName, nNbr=3, mean=True, returnMe=False)
           varName: variable that you want the off-event mean/sum for
           selName: selection used to define which events are good enough on/off events
           nNbrs: number of off-events neighboring the current on event
           mean: if True, return mean of variable for close off events
                 if False, ???
           returnMe: if False, add data to xrData; if True, also return
        """
        if '%s_offIdx_nNbr%02d'%(selName, nNbr) in self.Keys():
            startOffIdx = self.getVar('%s_offIdx_nNbr%02d'%(selName,nNbr))
        else:
            startOffIdx = self._getStartOffIdx(selName, nNbr=nNbr)
        filterOff = self.getFilter(selName.split('__')[0]+'__off')
        varArray = self.getVar(varName)
        varArrayOff =  get_offVar(varArray, filterOff, startOffIdx, nNbr=nNbr, mean=mean)
        if mean:
            self.addVar('offNbrsAv_%s_%s_nNbr%02d'%(varName.replace('/','_'),selName, nNbr), varArrayOff)
        else:
            self.addVar('offNbrs_%s_%s_nNbr%02d'%(varName.replace('/','_'),selName, nNbr), varArrayOff)
        if returnMe:
            return varArrayOff
            
    def Keys2d(self, inh5 = None, printKeys=False):
        """ return list of variables names in data that contain 2d data for each event """
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
                try:
                    keyShapes.append(self._fields[key][0])
                except:
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
            if areaDet and thiskey.find('Cfg')>=0:
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
        
    def nEvts(self,printThis=False):
        """ return number of events in the smallData hdf5 file for this event"""
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
        """ print information for this run: total number of events and if we had a scan"""
        self.nEvts(printThis=True)
        isScan=False
        scanVar = self.getScanName()
        if scanVar is not None and scanVar!='':
            isScan=True
        if isScan:
            nPoints=np.unique(self.getVar('scan/%s'%scanVar)).shape[0]
            print 'this run is a scan of %s with %d points'%(scanVar,nPoints)

    def hasKey(self, inkey):
        """ return boolean to reflect if a given variable is present in data
            parameter: inkey (variable name as string)"""
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

    def addCut(self, varName, varmin=0, varmax=0, useFilter=None):
        """
        add a variable to the selection
        parameters: varName, varmin, varmax
        varName: variable name as in smallData hdf5 file (string)
        varmin: minimum of variable varName that passes
        varmax: maximum of variable varName that passes
        useFilter: name of selection to add this cut to (required)
        """
        if useFilter is None:
            print 'you need to pass the name for the filter/selection you want to use'
            return
        if not self.Sels.has_key(useFilter):
            self.Sels[useFilter] = Selection()
        if varmin!=varmax:
            self.Sels[useFilter].addCut(varName, varmin, varmax)
        else:
            if isinstance(varName, basestring):
                if self.Sels.has_key(varName):
                    self.Sels[useFilter].add(self.Sels[varName])
                else:
                    print 'Cannot add cuts for filter %s as its not defined yet'%varName
            elif isinstance(varName, Selection):
                self.Sels[useFilter].add(varName)
            else:
                print 'need to pass selection to add as Selection or its name'
        Filter = self.getFilter(useFilter=useFilter)
        self.Sels[useFilter]._setFilter(Filter)
        #if self.Sels[useFilter]._filter is not None:
        #    print 'of %d  event, %d pass:'%(self.Sels[useFilter]._filter.shape[0],self.Sels[useFilter]._filter.sum())
    def removeCut(self, varName, useFilter):
        """
        remove a variable from the selection
        parameters: varName, useFilter
           varName: name of variable on which to not cut anymore
           useFilter: name of selection/filter
        """
        if not self.Sels.has_key(useFilter):
            print 'Selection with name %s does not exist, cannot remove cut'%useFilter
            return
        self.Sels[useFilter].removeCut(varName)
        Filter = self.getFilter(useFilter=useFilter)
        self.Sels[useFilter]._setFilter(Filter)

    def printSelections(self, selName=None, brief=True):
        """
        print the square cuts that currently define the selection/filter
        parameter: selName=None, brief=True
           if selName is None, all available selections are listed
           if brief is True, the cuts are not printed out
        """
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
    def printCuts(self, selName=None, brief=True):
        """
        alias to printSelections
        """
        self.printSelections(selName=selName, brief=brief)
    def getFilterLaser(self, useFilter, ignoreVar=[]):
        #moved functionality into getFilter.
        useFilterBase = useFilter.split('__')[0]
        return [self.getFilter(useFilter=useFilterBase+"__on", ignoreVar=ignoreVar).squeeze(),self.getFilter(useFilter=useFilterBase+"__off", ignoreVar=ignoreVar).squeeze()]
        
    def getFilter(self, useFilter=None, ignoreVar=[], debug=False):
        if 'fiducials' in self._fields.keys():
            total_filter=np.ones_like(self.getVar('fiducials')).astype(bool)
        elif 'EvtID/fid' in self._fields.keys():
            total_filter=np.ones_like(self.getVar('EvtID/fid')).astype(bool)
        else:
            print 'cannot find fiducials, cannot define a filter'
            return None

        useFilterBase = useFilter.split('__')[0]
        if useFilter==None or useFilterBase not in self.Sels.keys():
            if debug and useFilterBase not in self.Sels.keys():
                print 'Selection %s is not available, defined Selections are:'%useFilter
                self.printSelections()
            return total_filter.squeeze()

        #if useFilter ends in __off, drop on requirements
        LaserReq = -1
        if len(useFilter.split('__'))>1:
            if (useFilter.split('__')[1] == 'on' or useFilter.split('__')[1] == 'On'):
                LaserReq = 1
            if (useFilter.split('__')[1] == 'off' or useFilter.split('__')[1] == 'Off'):
                LaserReq = 0

        if self.Sels[useFilterBase]._filter is not None and len(ignoreVar)==0:
            if LaserReq == -1:
                return self.Sels[useFilterBase]._filter.squeeze()
            if LaserReq == 1:
                return (self.Sels[useFilterBase]._filter&[self.getVar('lightStatus/laser')==1]).squeeze()

        if LaserReq == 0:
            ignoreVar.append('lightStatus/laser')
            if self.ttCorr is not None:
                ignoreVar.append(self.ttCorr)
            if self.hasKey(self.ttBaseStr+'AMPL'):
                ignoreVar.append(self.ttBaseStr+'AMPL')
                ignoreVar.append(self.ttBaseStr+'FLTPOS')
                ignoreVar.append(self.ttBaseStr+'FLTPOS_PS')
                ignoreVar.append(self.ttBaseStr+'FLTPOSFWHM')

        filters=[]
        for thiscut in self.Sels[useFilterBase].cuts:
            if not thiscut[0] in ignoreVar:
                thisPlotvar=self.get1dVar(thiscut[0])
                filters.append(~np.isnan(thisPlotvar))
                if len(thiscut)==3:
                    filters.append((thisPlotvar > thiscut[1]) & (thisPlotvar < thiscut[2]))
                else:
                    filters.append(thisPlotvar != thiscut[1])
                if debug: 
                    total_filter_test = np.empty_like(total_filter)
                    total_filter_test[:] = total_filter
                    for ft in filters:
                        total_filter_test&=ft     
                        
                    print 'getFilter: Cut %f < %s < %f passes %d events of %d, total passes up to now: %d '%(thiscut[1], thiscut[0], thiscut[2], filters[-1].sum(), thisPlotvar.shape[0], total_filter_test.sum())

        for ft in filters:
            total_filter&=ft     
           
        if LaserReq == 2:
            return total_filter.squeeze()
        if LaserReq == 1:
            return (total_filter&[self.getVar('lightStatus/laser')==1]).squeeze()
        if LaserReq == 0:
            return (total_filter&[self.getVar('lightStatus/laser')==0]).squeeze()

        return total_filter
        
    def saveFilter(self, baseName='boolArray',useFilter=None, ignoreVar=[]):
        total_filter = self.getFilter(useFilter=useFilter, ignoreVar=ignoreVar)
        np.savetxt('%s_Run%03d.txt'%(baseName, self.run),total_filter.astype(bool),fmt='%5i')

    def getSelIdx(self, useFilter):
        if self.hasKey('EvtID'):
            fids = self.getVar('/EvtID/fid')
            times = self.getVar('/EvtID/time')
        elif self.hasKey('/fiducials'):
            fids = self.getVar('/fiducials')
            times = self.getVar('/event_time')
        elif self.hasKey('fiducials'):
            fids = self.getVar('fiducials')
            times = self.getVar('event_time')
        Filter =  self.getFilter(useFilter)
        selfids = [ (ifid,itime) for ifid,itime in zip(fids[Filter],times[Filter])]
        return selfids
        
    def getVar(self, plotvar, useFilter=None, addToXarray=True):
        #get the signal variable
        if isinstance(plotvar,list):
            plotvar = plotvar[0]
            sigROI = plotvar[1]
        else:
            sigROI=[]

        if isinstance(useFilter,basestring):
            Filter = self.getFilter(useFilter)
        else:
            Filter = useFilter

        if plotvar in self._fields.keys():
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
                        #this can fail when the data is too large. better debugging?
                        #print 'DEBUG B ',('/'+'/'.join(plotvar.split('/')[:-1]),plotvar.split('/')[-1])
                        #print 'node: ',self.fh5.get_node('/'+'/'.join(plotvar.split('/')[:-1]),plotvar.split('/')[-1])
                        vals = self.fh5.get_node('/'+'/'.join(plotvar.split('/')[:-1]),plotvar.split('/')[-1]).read()
                    else:
                        vals = self.fh5.get_node('/'+plotvar).read()
                else:
                    vals = self.fh5.fetch_data(run=self.run,keys=[plotvar])[plotvar]
                if vals.shape[0]==self._tStamp.shape[0]:
                    tArrName = plotvar.replace('/','__')
                    if addToXarray:
                        print 'add me to xarray...',plotvar
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
        sigROI=[]
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
            vals = np.nansum(vals,axis=1)
        return vals

    def setDelay(self, use_ttCorr=True, addEnc=False, addLxt=True, reset=False):
        delay=self.getDelay(use_ttCorr=use_ttCorr, addEnc=addEnc, addLxt=addLxt, reset=reset)

    #make delay another Xarray variable.
    def getDelay(self, use_ttCorr=True, addEnc=False, addLxt=True, reset=False):
        """
        function to get the xray-laser delay from the data
        usage:
        getDelay(): get the delay from lxt and/or encoder stage, add the timetool correction
        getDelay(use_ttCorr=False): get the delay from lxt and/or encoder stage, NO timetool correction
        getDelay(addEnc=True): get the delay from lxt, add encoder stage and timetool correction
        return previously define delay unless reset=True
        """
        if 'delay' in self._fields.keys() and not reset:
            return self.xrData['delay']
        
        self._delay_ttCorr=use_ttCorr
        self._delay_addLxt=addLxt
        self._delay_addEnc=addEnc

        ttCorrStr, ttBaseStr = self._getTTstr()
        ttCorr = np.zeros_like(self.xrData.fiducials)
        if self.ttCorr is not None:
            ttCorr=self.getVar(self.ttCorr)
        if (np.nanstd(ttCorr)==0):
            if (self.ttBaseStr+'FLTPOS_PS') in self._fields.keys():
                ttCorr=self.getVar(self.ttBaseStr+'FLTPOS_PS')
        nomDelay=np.zeros_like(ttCorr)
        if len(nomDelay.shape) == 0:
            return None

        isDaqDelayScan=False
        scanVar = self.getScanName()
        if scanVar.find('lxt')>=0:
            isDaqDelayScan=True
            #print 'DEBUG: found that we have a delay scan'
            nomDelay=self.getVar('scan/%s'%scanVar)*1e12

        if not isDaqDelayScan:
            if self.hasKey('enc/lasDelay'):
                encVal = self.getVar('enc/lasDelay')
                #print 'DEBUG: encoder info',np.nanstd(encVal)
                if np.nanstd(encVal)>1e-9:
                    nomDelay=encVal
                    addEnc=False
                elif np.nanstd(encVal)>1e-15:
                    nomDelay=encVal*1e12
                    addEnc=False
                elif self.hasKey('enc/ch0'):
                    encVal = self.getVar('enc/ch0')
                    if np.nanstd(encVal)>1e-15 and np.nanstd(encVal)<1e-9:
                        nomDelay=encVal*1e12
                        #now look at the EPICS PV if everything else has failed.
                    elif np.nanstd(encVal)>1e-3:
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
                nomDelay=nomDelay.copy()+self.getVar('enc/lasDelay')

        if addLxt and scanVar!='lxt':
            try:
                nomDelay=nomDelay.copy()+self.getVar('epics/lxt_ttc')*1e12
            except:
                pass

        if use_ttCorr:
            #print 'DEBUG adding ttcorr,nomdelay mean,std: ',ttCorr.mean(),nomDelay.mean(),ttCorr.std(),nomDelay.std()
            delay = (ttCorr+nomDelay)
        else:
            delay = nomDelay

        self.addVar('delay', delay)
        return delay

    def getPeak(self, plotvar, numBins=[100], useFilter=None, limits=[1,99,'p'],fig=None,asHist=False,sigROI=[]):
        hst = plotVar(plotvar, numBins=[100], useFilter=None, limits=[1,99,'p'],fig=None,asHist=False,sigROI=[])
        if len(hst)==2:
            return [hst[0].max(), hst[1][hst[0].argmax()]]
        else:
            print 'getPeak is not yet implemented for this type of data (need 1d histo)'

    def plotVar(self, plotvar, numBins=[100], useFilter=None, limits=[1,99,'p'],fig=None,asHist=False,plotWith=None, plotMultidimMean=False):
        if not isinstance(numBins, (list, tuple)):
            numBins = [numBins]
        if isinstance(plotvar, basestring) or (len(plotvar)==2 and (isinstance(plotvar[0], basestring) and not isinstance(plotvar[1], basestring))):
            if len(numBins)!=1:
                print 'bin# needs to be of same dimensions as plotvariables (1d)'
                return
            return self.plotVar1d(plotvar, numBins=numBins[0], useFilter=useFilter, limits=limits,fig=fig,plotWith=plotWith, plotMultidimMean=plotMultidimMean)
        elif len(plotvar)>2:
            print 'plotting of 3 variables is not defined yet'
            return
        if len(numBins)!=2:
            if len(numBins)==1:
                numBins=[numBins[0],numBins[0]]
            else:
                print 'bin# needs to be of same dimentions as plotvariables (2d)'
            if plotMultidimMean:
                print 'plotting multidim means of two variables is not implemented yet'
                return
        return self.plotVar2d(plotvar, numBins=numBins, useFilter=useFilter, limits=limits,fig=fig,asHist=asHist,plotWith=plotWith)

    def plotVar1d(self, plotvar, numBins=100, useFilter=None, limits=[1,99,'p'],fig=None, plotWith=None, plotMultidimMean=False):
        if plotWith is None:
            plotWith = self.plotWith

        if isinstance(plotvar,list):
            if not (self.hasKey(plotvar[0]) or plotvar[0]=='delay'): 
                print 'request variable %s not in smallData file'%plotvar
                return
        else:
            if not (self.hasKey(plotvar) or plotvar=='delay'): 
                print 'request variable %s not in smallData file'%plotvar
                return

        if plotvar=='delay':
            vals = self.getDelay()
        elif plotMultidimMean:
            vals = self.getRedVar(plotvar)
        elif len(plotvar)==1 and plotvar.find('droplets')>=0:
            vals = self.getRedVar(plotvar)
        else:
            vals = self.get1dVar(plotvar)

        total_filter = np.ones(vals.shape[0]).astype(bool)
        if useFilter is not None:
            total_filter =  self.getFilter(useFilter, [plotvar])
        vals = vals[total_filter]

        if  len(plotvar)==1 and plotvar.find('droplets')>=0:
            vals = vals.flatten()[vals.flatten()>0]

        if not plotMultidimMean:
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
        else:
            vals = np.nanmean(vals,axis=0)
            if len(vals.shape)==0:
                print 'mean of plotvar is ',vals
                return 

            if len(vals.shape)==1:
                hst = [vals, np.arange(0,vals.shape[0]+1)]
            elif len(vals.shape)==2:
                if plotWith.find('matplotlib')>=0:
                    plotImage(vals, ylim=[np.nanpercentile(vals,1),np.nanpercentile(vals,99)], xLabel='', yLabel='', plotWith=plotWith, plotTitle='%s for %s'%(plotvar, self.runLabel), fig=fig)
                    return
                elif plotWith.find('bokeh')>=0:
                    plotImage(vals, xLabel='', yLabel='', plotWith=plotWith, plotTitle='%s for %s'%(plotvar, self.runLabel))
                    return
                
            elif len(vals.shape)==3:
                print 'cannot plot 3-dim data'
                return

        if isinstance(plotvar, list):
            plotXlabel = plotvar[0]
        else:
            plotXlabel = plotvar
        plotMarker(hst[0], xData=hst[1][:-1], xLabel=plotXlabel, yLabel='entries', plotWith=plotWith, runLabel=self.runLabel, plotTitle="%s histogram for %s"%(plotvar, self.runLabel), plotDirname=self.plot_dirname)
        return hst

    def plotVar2d(self, plotvars, useFilter=None, limits=[1,99,'p'], asHist=False,numBins=[100,100],fig=None, plotWith=None):
        if plotWith is None:
            plotWith = self.plotWith

        for plotvar in plotvars:
            if isinstance(plotvar,list):
                plotvar = plotvar[0]
            if not self.hasKey(plotvar) or plotvar == 'delay': 
                print 'request variable %s not in smallData file'%(plotvar)
                return
        vals=[]
        for plotvar in plotvars:
            if plotvar == 'delay':
                vals=self.getDelay()
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
        if useFilter is not None and self.Sels.has_key(useFilter):
            filters.append(self.getFilter(useFilter,plotvars))
        for ft in filters:
            total_filter&=ft                

        print 'select ',total_filter.sum(),' of ',np.ones_like(total_filter).sum(),' events'
        if not asHist:
            msize=2
            if len(vals[1][total_filter])<100:
                msize=5
            elif len(vals[1][total_filter])<1000:
                msize=3
            plotMarker(vals[0][total_filter], xData=vals[1][total_filter], xLabel=plotvars[1], yLabel=plotvars[0], plotWith=plotWith, runLabel=self.runLabel, plotTitle="%s vs %s for %s"%(plotvars[0], plotvars[1], self.runLabel), plotDirname=self.plot_dirname, markersize=msize)
            return vals[0][total_filter], vals[1][total_filter]
        
        #asHist only
        v0 = vals[0][total_filter]
        v1 = vals[1][total_filter]
        binEdges0 = np.linspace(np.nanmin(v0),np.nanmax(v0),numBins[0])
        binEdges1 = np.linspace(np.nanmin(v1),np.nanmax(v1),numBins[1])
        ind0 = np.digitize(v0, binEdges0)
        ind1 = np.digitize(v1, binEdges1)
        ind2d = np.ravel_multi_index((ind0, ind1),(binEdges0.shape[0]+1, binEdges1.shape[0]+1)) 
        iSig = np.bincount(ind2d, minlength=(binEdges0.shape[0]+1)*(binEdges1.shape[0]+1)).reshape(binEdges0.shape[0]+1, binEdges1.shape[0]+1) 
        extent=[binEdges1[1],binEdges1[-1],binEdges0[1],binEdges0[-1]]

        plotImage(iSig, extent=extent, ylim=[np.nanpercentile(iSig,limits[0]),np.percentile(iSig,limits[1])], xLabel=plotvars[1], yLabel=plotvars[0], plotWith=plotWith, fig=fig)

        return iSig, extent

    def getScanName(self):
        for key in self.Keys('scan'):
            if key.find('var')<0 and key.find('none')<0 and key.find('damage')<0:
                return key.replace('/scan/','').replace('scan/','')
        return ''

    def getScanValues(self):
        #get the scan variable & time correct if desired
        scanVarName = self.getScanName()
        if scanVarName.find('lxt')>=0 or scanVarName=='':
            delays=self.getDelay()
            #CHECK ME: not sure why I required both mean&std to be==0 for not scan?
            if delays is None or delays.mean()==0 or delays.std()==0: 
                return '',[]
            scan = delays
            if scanVarName == '': 
                if self._delay_ttCorr:
                    scanVarName='delay (tt corrected) [fs]'
                else:
                    scanVarName='delay [fs]'
        else:
            try:
                scanOrg = self.getVar('scan/var0')
                scan = scanOrg
            except:
                scan=[]
        return scanVarName,scan


    def getBins(self,bindef=[], debug=False):
        Bins = util_getBins(bindef, debug)
        if Bins is not None:
            return Bins

        #have no input at all, assume we have unique values in scan. If not, return empty list 
        scanVarName, scan =  self.getScanValues()

        if len(bindef)==0:
            if scanVarName=='':
                print 'this run is no scan, will need bins as input, quit now'
                return []
            print 'no bins as input, we will use the scan variable %s '%scanVarName
            Bins = np.unique(scan)
            if scanVarName.find('lxt')>=0:
                Bins*=1e12
            if debug:
                print 'Bins: ',Bins
            return Bins

        #give a single number (binwidth or numBin)
        if len(bindef)==1:
            #this is a "normal" scan, use scanVar
            if scanVarName!='':
                Bins = np.unique(scan)
                if scanVarName.find('lxt')>=0:
                    Bins*=1e12
                valBound = [ min(Bins), max(Bins)]
                if type(bindef[0]) is int:
                    Bins=np.linspace(valBound[0],valBound[1],bindef[0],endpoint=True)
                else:
                    Bins=np.arange(valBound[0],valBound[1],bindef[0])
                    Bins=np.append(Bins,valBound[1])
                return Bins
        #free running scan...rely on getDelay
        minEnc = (int(min(scan)*10.))
        if minEnc<0:
            minEnc+=-1
        minEnc /= 10.
        maxEnc = (int(max(scan)*10.)+1)/10.
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

    def getScan(self, sig='', i0='', Bins=100, useFilter=None):
        return self.plotScan(sig=sig, i0=i0, Bins=Bins, returnData=True, useFilter=useFilter, plotThis=False)

    def plotScan(self, sig='', i0='', Bins=100, plotDiff=True, plotOff=True, saveFig=False,saveData=False, returnData=False, useFilter=None, fig=None, interpolation='', plotThis=True, returnIdx=False, binVar=None, plotWith=None):
        if plotWith is None:
            plotWith = self.plotWith

        plotVar=''
        if sig!='':
            sigVal = self.get1dVar(sig)
            for sigp in sig:
                if isinstance(sigp,basestring):
                    plotVar+=sigp.replace('/','__')
                elif isinstance(sigp,list):
                    for bound in sigp:
                        plotVar+='-%g'%bound
        else:
            print 'could not get signal variable %s, please specify'%plotVar.replace('__','/')
            return

        if i0!='':
            i0Val = self.get1dVar(i0)
            plotVar+='/'
            for i0p in i0:
                if isinstance(i0p,basestring):
                    plotVar+=i0p.replace('/','__')
                elif isinstance(i0p,list):
                    for bound in i0p:
                        plotVar+='-%g'%bound
        else:
            i0Val = np.ones_like(sigVal)
        
        [FilterOn, FilterOff] = self.getFilterLaser(useFilter)
        FilterOn = FilterOn & ~np.isnan(i0Val) & ~np.isnan(sigVal)
        FilterOff = FilterOff & ~np.isnan(i0Val) & ~np.isnan(sigVal)

        #get the binning variable here so that points where this is not good can be thrown out.
        if binVar is not None:
            if binVar[0] != 'delay':
                if isinstance(binVar, basestring): binVar=[binVar]
                binVal = self.get1dVar(binVar[0])
            else:
                binVal=self.getDelay()
            FilterOn = FilterOn & ~np.isnan(binVal)
            FilterOff = FilterOff & ~np.isnan(binVal)

        print 'from %i events keep %i (%i off events)'%(np.ones_like(i0Val).sum(),np.ones_like(i0Val[FilterOn]).sum(), np.ones_like(i0Val[FilterOff]).sum() )

        #get the scan variable & time correct if desired
        scanVarName,scan =  self.getScanValues()
            
        usedDigitize = 0
        # create energy bins for plot: here no need to bin!
        if (not self._delay_ttCorr) and (not self._delay_addEnc):
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
            usedDigitize = 1

        if returnIdx:
            return scanOnIdx

        #now do the same for laser off data
        OffData=False
        if scan[FilterOff].sum()!=0:
            scanPointsT, scanOnIdxT = np.unique(scan[FilterOn], return_inverse=True)
            scanOffPoints, scanOffIdx = np.unique(scan[FilterOff], return_inverse=True)

            #unique & digitize do not behave the same !!!!!
            if len(scanOffPoints) > len(scanPoints):
                scanOffPoints = scanPoints.copy()
                usedDigitize = 1
            if usedDigitize>0:
                scanOffIdx = np.digitize(scan[FilterOff], scanOffPoints)
                
            OffData = True

        #now get the binning information for second variable.
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
            iNorm = np.bincount(indOn2d, weights=i0Val[FilterOn], minlength=(scanPoints.shape[0]+1)*(binPoints.shape[0]+1)).reshape(scanPoints.shape[0]+1, binPoints.shape[0]+1).T    
            iSig = np.bincount(indOn2d, weights=sigVal[FilterOn], minlength=(scanPoints.shape[0]+1)*(binPoints.shape[0]+1)).reshape(scanPoints.shape[0]+1, binPoints.shape[0]+1).T    
        else:
            iNorm = np.bincount(scanOnIdx, i0Val[FilterOn], minlength=len(scanPoints)+1)
            iSig = np.bincount(scanOnIdx, sigVal[FilterOn], minlength=len(scanPoints)+1)

        if OffData:
            print '#scan points: -- 3 ',len(scanPointsT),len(scanOffPoints),len(scanPoints)
        else:
            print '#scan points: -- 3 ',len(scanPoints)
        scan = iSig/iNorm
        #scan = scan/np.mean(scan[1]) # normalize to 1 for first energy point?
        scan = scan[1:-1]
        scanPoints = (scanPoints[:-1]+scanPoints[1:])*0.5

        if OffData:
            #same for off shots
            iNormoff = np.bincount(scanOffIdx, i0Val[FilterOff], minlength=len(scanOffPoints)+1)
            iSigoff = np.bincount(scanOffIdx, sigVal[FilterOff], minlength=len(scanOffPoints)+1)
            scanoff = iSigoff/iNormoff
            scanoff = scanoff[1:-1]
            scanOffPoints = (scanOffPoints[:-1]+scanOffPoints[1:])*0.5
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
            returnDict['binVar']=binVar[0]
        if plotThis:
            #print 'retData: ',returnDict
            self.plotScanDict(returnDict, plotDiff=plotDiff, interpolation=interpolation,fig=fig,plotOff=plotOff,saveFig=saveFig, plotWith=plotWith)
        return returnDict

    def plotScanDict(self, returnDict, plotDiff=True, fig=None, plotOff=True, interpolation='', saveFig=False, plotWith=None):
        if plotWith is None:
            plotWith = self.plotWith

        plotVarName = returnDict['plotVarName']
        scanVarName = returnDict['scanVarName']
        scanPoints = returnDict['scanPoints']
        scan = returnDict['scan']
        print 'plot ',plotVarName, scanVarName, ' shape ',scan.shape,' plot diff ',plotDiff

        if interpolation!='' and returnDict.has_key('scanOffPoints'):
            finter_off = interpolate.interp1d(returnDict['scanOffPoints'], returnDict['scanOff'],kind=interpolation)
            scanoff_interp = finter_off(scanPoints[:-1])
        if plotWith=='matplotlib':
            if len(scan.shape)>1:
                extent = [scanPoints.min(), scanPoints.max(), returnDict['binPoints'].min(), returnDict['binPoints'].max()]
                plotImage(scan, ylim=[np.nanpercentile(scan,1), np.nanpercentile(scan,98)], xLabel=scanVarName, yLabel=returnDict['binVar'], plotWith=plotWith, plotTitle='%s  for %s'%(returnDict['plotVarName'], self.runLabel), extent=extent)
                return
            else:
                scanYvals=[scan]
                scanXvals=[scanPoints]
                markers = ['o']
                colors = ['red']
                if returnDict.has_key('scanOffPoints') and plotOff:
                    markers.append('o')
                    colors.append('black')
                    scanYvals.append(returnDict['scanOff'])
                    scanXvals.append(returnDict['scanOffPoints'])
                    if interpolation!='':
                        markers.append('D')
                        colors.append('black')
                        scanYvals.append(scanoff_interp)
                        scanXvals.append(scanPoints)

                if plotDiff and returnDict.has_key('scanOffPoints') and (interpolation!='' or len(scan)==len(returnDict['scanOff'])):
                    if fig is None:
                        fig=plt.figure(figsize=(10,10))
                    gs=gridspec.GridSpec(2,1,width_ratios=[1])
                    if interpolation!='':
                        ydata = (scan[:-1]-scanoff_interp)
                    else:
                        ydata = (scan[:-1]-returnDict['scanOff'][:-1])
                    plotMarker(scanYvals, xData=scanXvals, xLabel=scanVarName, yLabel=plotVarName, plotWith=plotWith, runLabel=self.runLabel, plotTitle="%s vs %s for %s"%(plotVarName, scanVarName, self.runLabel), plotDirname=self.plot_dirname, markersize=5, plotFilename=('Scan_%s_vs_%s'%(plotVarName, scanVarName)), line_dash='dashed', color=colors, marker=markers, ylim=[np.nanmin(scan)*0.95,np.nanmax(scan)*1.05], fig=gs[0])
                    plotMarker( ydata, xData=scanPoints[:-1], xLabel=scanVarName, yLabel='on-off '+plotVarName, plotWith=plotWith, runLabel=self.runLabel, plotTitle='', plotDirname=self.plot_dirname, markersize=5, plotFilename=('Scan_diff_%s_vs_%s'%(plotVarName, scanVarName)), line_dash='dashed', color=['blue'], fig=gs[1])
                else:
                    if saveFig:
                        plotWith='matplotlib_file'
                    plotMarker(scanYvals, xData=scanXvals, xLabel=scanVarName, yLabel=plotVarName, plotWith=plotWith, runLabel=self.runLabel, plotTitle="%s vs %s for %s"%(plotVarName, scanVarName, self.runLabel), plotDirname=self.plot_dirname, markersize=5, plotFilename=('Scan_%s_vs_%s'%(plotVarName, scanVarName)), line_dash='dashed', color=colors, marker=markers, ylim=[np.nanmin(scan)*0.95,np.nanmax(scan)*1.05])

        elif plotWith.find('bokeh')>=0:
            if len(scan.shape)>1:
                extent = [scanPoints.min(), scanPoints.max(), returnDict['binPoints'].min(), returnDict['binPoints'].max()]
                output_file='%s/%s_Scan_%s_%s.html'%(self.plot_dirname,self.runLabel, scanVarName.replace('/','_'), plotVarName.replace('/','_'))
                plotImage(scan, ylim=[np.nanpercentile(scan,1), np.nanpercentile(scan,98)], xLabel=scanVarName, yLabel=returnDict['binVar'], plotWith=plotWith, plotTitle='%s  for %s'%(returnDict['plotVarName'], self.runLabel), extent=extent,plotFilename=output_file)
                return
            else:
                scanYvals=[scan]
                scanXvals=[scanPoints]
                markers = ['o']
                colors = ['red']
                if returnDict.has_key('scanOffPoints') and plotOff:
                    markers.append('o')
                    colors.append('black')
                    scanYvals.append(returnDict['scanOff'])
                    scanXvals.append(returnDict['scanOffPoints'])
                    if interpolation!='':
                        markers.append('D')
                        colors.append('black')
                        scanYvals.append(scanoff_interp)
                        scanXvals.append(scanPoints)

                if plotDiff and returnDict.has_key('scanOffPoints') and (interpolation!='' or len(scan)==len(returnDict['scanOff'])):

                    p = plotMarker(scanYvals, xData=scanXvals, xLabel=scanVarName, yLabel=plotVarName, plotWith=plotWith, runLabel=self.runLabel, plotTitle="%s vs %s for %s"%(plotVarName, scanVarName, self.runLabel), plotDirname=self.plot_dirname, markersize=5, plotFilename=('Scan_%s_vs_%s'%(plotVarName, scanVarName)), line_dash='dashed', color=colors, marker=markers, ylim=[np.nanmin(scan)*0.95,np.nanmax(scan)*1.05], width_height=(750,350), fig='return')

                    if interpolation!='':
                        ydata = (scan[:-1]-scanoff_interp)
                    else:
                        ydata = (scan[:-1]-returnDict['scanOff'][:-1])
                    pdiff = plotMarker( ydata, xData=scanPoints[:-1], xLabel=scanVarName, yLabel='on-off '+plotVarName, plotWith=plotWith, runLabel=self.runLabel, plotTitle='', plotDirname=self.plot_dirname, markersize=5, plotFilename=('Scan_diff_%s_vs_%s'%(plotVarName, scanVarName)), line_dash='dashed', color=['blue'], fig='return', width_height=(750,400))

                    grid = bokeh.plotting.gridplot([[p], [pdiff]])
                    layout = column(Div(text='<h1>%s as a function of %s for %s</h1>'%(plotVarName, scanVarName, self.runLabel)),grid)

                    if plotWith=='bokeh_notebook':
                        bp.show(layout)
                    else:
                        bp.save(layout)

                else:
                    plotMarker(scanYvals, xData=scanXvals, xLabel=scanVarName, yLabel=plotVarName, plotWith=plotWith, runLabel=self.runLabel, plotTitle="%s vs %s for %s"%(plotVarName, scanVarName, self.runLabel), plotDirname=self.plot_dirname, markersize=5, plotFilename=('Scan_%s_vs_%s'%(plotVarName, scanVarName)), line_dash='dashed', color=colors, marker=markers, ylim=[np.nanmin(scan)*0.95,np.nanmax(scan)*1.05])

        elif plotWith != 'no_plot':
            print 'plotting using %s is not implemented yet, options are matplotlib, bokeh_notebook, bokeh_html or no_plot'

            
    def defPlots(self, useFilter=None):
        scanVarName,scan =  self.getScanValues(True)
        total_filter = np.ones_like(scan).astype(bool)
        if useFilter is not None and self.Sels.has_key(useFilter):
            total_filter =  self.getFilter(useFilter, [plotvar])

        fig=plt.figure(figsize=(10,6))
        plt.title('Standard Plots for Run %i'%self.run)
        
        gs=gridspec.GridSpec(2,2,width_ratios=[2,2])
        self.plotVar('ipm2/sum',fig=plt.subplot(gs[0]),useFilter=useFilter)
        self.plotVar(['ipm2/sum','ebeam/L3Energy'],fig=plt.subplot(gs[1]),asHist=True,useFilter=useFilter)
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
    def addCube(self, cubeName, binVar='', bins=[], useFilter='', addBinVars=None):    
        self.cubes[cubeName] = Cube(binVar, bins, cubeName=cubeName, useFilter=useFilter, addBinVars=addBinVars)
        
    def addToCube(self, cubeName, targetVariable, isIdxVar=False):
        if cubeName in self.cubes.keys():
            if not isIdxVar:
                self.cubes[cubeName].addVar(targetVariable)
            else:
                self.cubes[cubeName].addIdxVar(targetVariable)
        else:
            print 'could not add variable %s to cube %s as this cube was not found'%(targetVariable, cubeName)

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
                    cube.printCube(self.Sels[cube.useFilter])
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
        elif len(cube.bins)==0:
            scanValues = self.getScanValues()
            if scanValues[0]!='' and len(scanValues[1])>0:
                Bins = np.unique(scanValues[1])
                Bins = np.insert(Bins,0,Bins[0]-1e-5)
                Bins = np.append(Bins,Bins[-1]+1e-5)
        else:
            Bins = self.getBins(cube.bins)
        cube.binBounds = Bins

        #now look through targetVars & split out ones not in xarray/hdf5
        targetVarsLocal = []
        for tVar in cube.targetVars:
            if isinstance(tVar, basestring):
                if tVar.split(':')[0]=='droplet': #ex: droplet:epix_2/droplet:image
                    dropBaseName = tVar.split(':')[1]
                    dropDet = dropBaseName.split('/')[0]
                    dropFct = tVar.split(':')[2]
                    #add droplet source.
                    if dropBaseName not in cube.dropletProc.keys():
                        if (dropBaseName+'_X') in self._fields.keys():
                            cube.dropletProc[dropBaseName] = {}
                            cube.dropletProc[dropBaseName]['source'] = dropDet
                            cube.dropletProc[dropBaseName]['fct'] = [dropFct]
                            cube.dropletProc[dropBaseName]['vars'] = [dropBaseName+'_X', dropBaseName+'_Y', dropBaseName+'_adu', dropBaseName+'_npix']
                            #get shape for image.
                            cube.dropletProc[dropBaseName]['shape'] = ( np.nanmax(self.getVar('UserDataCfg/%s/iX'%dropDet)), 
                                                                        np.nanmax(self.getVar('UserDataCfg/%s/iY'%dropDet)))
                        else:
                            print 'could not find droplets with base name: ',dropBaseName
                    else:
                        cube.dropletProc[dropBaseName]['fct'].append(dropFct)                        

                elif tVar not in self._fields.keys():
                    cube.targetVarsXtc.append(tVar)
                else:
                    targetVarsLocal.append(tVar)
            else:
                cube.targetVarsXtc.append(tVar)
        cube.targetVars = targetVarsLocal

        #now get the filter & create a new one taking bins & detector damage into account.
        onoff=2
        useFilterBase = cube.useFilter.split('__')[0]
        if len(cube.useFilter.split('__'))>1:
            if cube.useFilter.split('__')[1]=='on': onoff=1
            elif cube.useFilter.split('__')[1]=='off': onoff=0
        if cube.useFilter.find(cube.cubeName)!=0 and '%s_%s'%(cube.cubeName,useFilterBase) not in self.Sels.keys():
            self.Sels['%s_%s'%(cube.cubeName,useFilterBase)] = Selection()
            self.Sels['%s_%s'%(cube.cubeName,useFilterBase)].add(self.Sels[useFilterBase])
            #if cube.binVar.find('lxt')>=0:
            #    Bins/=1e12
            self.Sels['%s_%s'%(cube.cubeName,useFilterBase)].addCut(cube.binVar, min(Bins), max(Bins) )
            for addVar in cube.addBinVars.keys():
                self.Sels['%s_%s'%(cube.cubeName,useFilterBase)].addCut(addVar, min(cube.addBinVars[addVar]), max(cube.addBinVars[addVar]) )
            #add cuts with detector damage - if we have damage detector info.
            for txVar in targetVarsLocal:
                if txVar[0]=='/':txVar=txVar[1:] 
                if 'damage/%s'%txVar  in self._fields.keys(): 
                    self.Sels['%s_%s'%(cube.cubeName,cube.useFilter)].addCut('damage/%s'%txVar.split('/')[0],0.5,1.5)
            for txVar in cube.targetVarsXtc:
                if isinstance(txVar, dict):
                    try:
                        txVar=txVar['source']
                    except:
                        continue
                if 'damage/%s'%txVar  in self._fields.keys(): 
                        self.Sels['%s_%s'%(cube.cubeName,cube.useFilter)].addCut('damage/%s'%txVar,0.5,1.5)
            for dropletBase in cube.dropletProc.keys():
                dropDet = cube.dropletProc[dropletBase]['source']
                if 'damage/%s'%dropDet in self._fields.keys(): 
                    self.Sels['%s_%s'%(cube.cubeName,cube.useFilter)].addCut('damage/%s'%dropDet,0.5,1.5)
            cube.useFilter='%s_%s'%(cube.cubeName,cube.useFilter)

        return cube, onoff

    def makeCubeData(self, cubeName, debug=False, toHdf5=None, replaceNan=False, onoff=2, returnIdx=False):
        cube, cubeName_onoff = self.prepCubeData(cubeName)
        if onoff == 2:
            onoff = cubeName_onoff

        if cube is None:
            return 
        Bins = cube.binBounds

        if cube.binVar == 'delay':
            binVar = self.getDelay()
        else:
            binVar = self.get1dVar(cube.binVar)
                
        #floating point error really mess this up otherwise....
        if abs(max(Bins))<1e-9:
            #XXX
            print 'Bins: ',Bins
            testIdx=np.digitize(binVar, Bins)
            print 'DEBUG: ',np.bincount(testIdx)
            testIdx=np.digitize(binVar*1e12, Bins*1e12)
            print 'DEBUG2: ',np.bincount(testIdx)
            
            binVar*=1e12
            Bins*=1e12

        if debug and rank==0:
            cube.printCube()
            #self.printSelections(self.Sels[cube.useFilter])
            self.printSelections(cube.useFilter)

        cubeFilter = self.getFilter(cube.useFilter)
        [cubeOn, cubeOff] = self.getFilterLaser(cube.useFilter, ignoreVar=[])

        if onoff==1:
            cubeFilter = cubeOn
            cubeName = cubeName+'_laserOn'
        elif onoff==0:
            cubeFilter = cubeOff
            cubeName = cubeName+'_laserOff'

        binVar = binVar[cubeFilter]
        if binVar.shape[0] == 0:
            printR(rank, 'did not select any event, quit now!')
            ft=self.getFilter(cube.useFilter, debug=True)
            return
        if debug and rank==0:
            print 'bin boundaries: ',Bins

        nTotBins=Bins.shape[0]-1
        binShp=[Bins.shape[0]-1]
        if len(cube.addBinVars.keys())>0:
            binIdx=[np.digitize(binVar, Bins)]
            if np.array(binIdx).min()<1:
                printR(rank, 'something went wrong in the setting of the cube selection, please fix me....')
                return
            binIdx = (np.array(binIdx)-1).tolist()
            for addVar in cube.addBinVars.keys():
                if addVar == 'delay':
                    addBinVar = self.getDelay()
                else:
                    addBinVar = self.get1dVar(addVar)
                addBinVar = addBinVar[cubeFilter]
                addBinIdx = np.digitize(addBinVar, cube.addBinVars[addVar])
                if np.array(addBinIdx).min()<1:
                    printR(rank, 'something went wrong in the setting of the cube selection for variable %s, please fix me....'%addVar)
                    return
                addBinIdx = (np.array(addBinIdx)-1).tolist()
                binShp.append(cube.addBinVars[addVar].shape[0]-1)
                nTotBins=nTotBins*(cube.addBinVars[addVar].shape[0]-1)
                binIdx.append(addBinIdx)

            #binShp & binIdx as tuple for use in ravel_multi_index. Append them to cube

            indMultiD = np.ravel_multi_index(tuple(binIdx),tuple(binShp))
            binVar = indMultiD
            orgBins = Bins
            Bins = np.arange(0,nTotBins+1)

        timeFiltered = self._tStamp[cubeFilter]
        newXr = xr.DataArray(np.ones(timeFiltered.shape[0]), coords={'time': timeFiltered}, dims=('time'),name='nEntries')
        newXr = xr.merge([newXr, xr.DataArray(binVar, coords={'time': timeFiltered}, dims=('time'),name='binVar') ])       
        
        for tVar in cube.targetVars:
            if not self.hasKey(tVar):                
                continue
            #printR(rank, 'addvar: ',tVar,self.getVar(tVar,cubeFilter).shape)
            filteredVar = self.getVar(tVar,cubeFilter).squeeze()
            tVar=tVar.replace('/','__')
            if len(filteredVar.shape)==1:
                #newXr = xr.merge([newXr, xr.DataArray(filteredVar, coords={'time': timeFiltered}, dims=('time'),name=tVar) ])
                newDar = xr.DataArray(filteredVar, coords={'time': timeFiltered}, dims=('time'),name=tVar)
            else:
                coords={'time': timeFiltered}
                dims = ['time']
                dataShape = filteredVar.shape
                for dim in range(len(dataShape)-1):
                    thisDim = np.arange(0, dataShape[dim+1])
                    dimStr = '%s_dim%d'%(tVar,dim)
                    coords[dimStr] = thisDim
                    dims.append(dimStr)
                newDar = xr.DataArray(filteredVar, coords=coords, dims=dims,name=tVar)
                #newXr = xr.merge([newXr, xr.DataArray(filteredVar, coords=coords, dims=dims,name=tVar)])
            newXr = xr.merge([newXr, newDar])

        #now we actually bin.
        cubeData = newXr.groupby_bins('binVar',Bins,labels=(Bins[1:]+Bins[:-1])*0.5,include_lowest=True, right=False).sum(dim='time')                  
        #could add error using the std of the values.
        cubeDataErr = newXr.groupby_bins('binVar',Bins,labels=(Bins[1:]+Bins[:-1])*0.5,include_lowest=True, right=False).std(dim='time')
            
        if len(cube.addBinVars.keys())>0:
            newXr = None
            for key in cubeData.variables:
                ##treat only actual data -- this is taken care of by using .variables instead.
                #if key in cubeData.dims:
                #    continue
                #get dimensions & coords for reshaped data
                dims = [cube.binVar]
                coords={cube.binVar: (orgBins[1:]+orgBins[:-1])*0.5}
                for addVar in cube.addBinVars.keys():
                    dims.append(addVar)
                    coords[addVar]=0.5*(cube.addBinVars[addVar][1:]+cube.addBinVars[addVar][:-1])
                for thisdim in cubeData[key].dims:
                    if thisdim=='binVar_bins':
                        continue
                    dims.append(thisdim)
                    coords[thisdim]=cubeData[key].coords[thisdim].data
                #do not reshape the coordinates!
                isCoord=True
                for thisdim in cubeData[key].dims:
                    if thisdim=='binVar_bins':
                        isCoord=False
                if isCoord: continue
                #now get data & create new data and new shape:
                dataShp = cubeData[key].shape
                newShp = tuple(np.append(np.array(binShp), np.array(dataShp)[1:]))
                data = cubeData[key].data.reshape(newShp)
                dataArray = xr.DataArray(data, coords=coords, dims=dims,name=key)
                if newXr is None:
                    newXr = dataArray
                else:
                    newXr = xr.merge([newXr, dataArray])

            for key in cubeDataErr.variables:
                ##treat only actual data
                #if key in cubeDataErr.dims:
                #    continue
                if key == 'nEntries' or key == 'binVar':
                    continue
                #get dimensions & coords for reshaped data
                dims = [cube.binVar]
                coords={cube.binVar: (orgBins[1:]+orgBins[:-1])*0.5}
                for addVar in cube.addBinVars.keys():
                    dims.append(addVar)
                    coords[addVar]=0.5*(cube.addBinVars[addVar][1:]+cube.addBinVars[addVar][:-1])
                for thisdim in cubeDataErr[key].dims:
                    if thisdim=='binVar_bins':
                        continue
                    dims.append(thisdim)
                    coords[thisdim]=cubeDataErr[key].coords[thisdim].data
                #do not reshape the coordinates!
                isCoord=True
                for thisdim in cubeDataErr[key].dims:
                    if thisdim=='binVar_bins':
                        isCoord=False
                if isCoord: continue
                #now get data & create new data and new shape:
                dataShp = cubeDataErr[key].shape
                newShp = tuple(np.append(np.array(binShp), np.array(dataShp)[1:]))
                data = cubeDataErr[key].data.reshape(newShp)
                newKey = ('std_%s'%key)
                dataArray = xr.DataArray(data, coords=coords, dims=dims,name=newKey)
                newXr = xr.merge([newXr, dataArray])
            cubeData = newXr
        else:
            for key in cubeDataErr.variables:
                if key.replace('std_','').replace('__','/') in cube.targetVars:
                    cubeDataErr.rename({key: 'std_%s'%key}, inplace=True)
            for key in cubeDataErr.variables:
                if key not in cubeData.variables:
                    cubeData = xr.merge([cubeData, cubeDataErr[key]])

        if not returnIdx and len(cube.addIdxVars)==0 and len(cube.dropletProc.keys())==0:
            if toHdf5 == 'h5netcdf':
                fname = '%s/Cube_%s_Run%03d_%s.nc'%(self.dirname,self.expname,self.run,cubeName)
                cubeData.to_netcdf(fname,engine='h5netcdf')
            elif toHdf5 == 'h5':
                fname = '%s/Cube_%s_Run%03d_%s.h5'%(self.dirname,self.expname,self.run,cubeName)
                h5Dict={}
                for key in cubeData.variables:
                    h5Dict[key] = cubeData[key].values
                    dictToHdf5(fname, h5Dict)
            return cubeData

        fidVar='fiducials'
        evttVar='event_time'
        if '/fiducials' in self.Keys():
            fidVar='/fiducials'
            evttVar='/event_time'
        elif 'EvtID/fid' in self.Keys():
            fidVar='EvtID/fid'
            evttVar='EvtID/time'
        printR(rank, 'we will use fiducials from here: %s'%fidVar)

        evtIDXr = xr.DataArray(self.getVar(fidVar,cubeFilter), coords={'time': timeFiltered}, dims=('time'),name='fiducial')
        evtIDXr = xr.merge([evtIDXr,xr.DataArray(self.getVar(evttVar,cubeFilter), coords={'time': timeFiltered}, dims=('time'),name='evttime')])
        evtIDXr = xr.merge([evtIDXr, xr.DataArray(binVar, coords={'time': timeFiltered}, dims=('time'),name='binVar') ])       

        for thisIdxVar in cube.addIdxVars:
            varData = self.getVar(thisIdxVar,cubeFilter)
            coords={'time': timeFiltered}
            dims = ['time']
            for dim in range(len(varData.shape)-1):
                thisDim = np.arange(0, varData.shape[dim+1])
                dimStr = 'dim%d'%dim
                coords[dimStr] = thisDim
                dims.append(dimStr)
            addArray = xr.DataArray(varData, coords=coords, dims=dims,name=thisIdxVar)
            evtIDXr = xr.merge([evtIDXr, addArray])
        for thisDropDictKey in cube.dropletProc.keys():
            for thisVar in cube.dropletProc[thisDropDictKey]['vars']:
                varData = self.getVar(thisVar,cubeFilter)
                coords={'time': timeFiltered}
                dims = ['time']
                for dim in range(len(varData.shape)-1):
                    thisDim = np.arange(0, varData.shape[dim+1])
                    dimStr = 'dim%d'%dim
                    coords[dimStr] = thisDim
                    dims.append(dimStr)
                addArray = xr.DataArray(varData, coords=coords, dims=dims,name=thisVar)
                evtIDXr = xr.merge([evtIDXr, addArray])
                
        cubeIdxData = evtIDXr.groupby_bins('binVar',Bins,labels=(Bins[1:]+Bins[:-1])*0.5,include_lowest=True, right=False)
        keys = cubeIdxData.groups.keys()
        keys.sort()

        fidArray=[]
        timeArray=[]
        addArray=[]
        for addIdxVar in cube.addIdxVars:
            addArray.append([])
        for key in (Bins[1:]+Bins[:-1])*0.5:
            if key in cubeIdxData.groups.keys():
                fidArray.append(evtIDXr.fiducial[cubeIdxData.groups[key]])
                timeArray.append(evtIDXr.evttime[cubeIdxData.groups[key]])
                for iv,thisIdxVar in enumerate(cube.addIdxVars):
                    addArray[iv].append(evtIDXr[thisIdxVar][cubeIdxData.groups[key]])
            else:
                fidArray.append([])
                timeArray.append([])
                for iv,thisIdxVar in enumerate(cube.addIdxVars):
                    addArray[iv].append([])
                
        for droplet in cube.dropletProc.keys():
            if 'image' in cube.dropletProc[droplet]['fct']:
                imageList=[]
                for key in (Bins[1:]+Bins[:-1])*0.5:
                    if key in cubeIdxData.groups.keys():
                        dropVar = cube.dropletProc[droplet]['vars']
                        npix = evtIDXr[dropVar[3]][cubeIdxData.groups[key]].data.flatten()
                        tx = evtIDXr[dropVar[0]][cubeIdxData.groups[key]].data.flatten()[npix>0]
                        ty = evtIDXr[dropVar[1]][cubeIdxData.groups[key]].data.flatten()[npix>0]
                        tadu = evtIDXr[dropVar[2]][cubeIdxData.groups[key]].data.flatten()[npix>0]
                        #make image
                        data = tadu
                        row = tx
                        col = ty
                        if max(row)>=cube.dropletProc[droplet]['shape'][0] or max(col)>=cube.dropletProc[droplet]['shape'][1]:
                            if rank==0:
                                print 'inconsistent shape ',self.shape, max(row), max(col)
                        imageAsMatrix = sparse.coo_matrix( (data, (row, col)),shape=cube.dropletProc[droplet]['shape']).todense()
                        imageList.append(np.asarray(imageAsMatrix))
                    else:
                        nanImage = np.zeros(cube.dropletProc[droplet]['shape'])
                        nanImage.fill(np.nan)
                        imageList.append(nanImage)
                        
                if len(cube.addBinVars.keys())>0:
                    dims = [cube.binVar]
                    coords={cube.binVar: (orgBins[1:]+orgBins[:-1])*0.5}
                    for addVar in cube.addBinVars.keys():
                        dims.append(addVar)
                        coords[addVar]=0.5*(cube.addBinVars[addVar][1:]+cube.addBinVars[addVar][:-1])
                    dataShape = np.array(imageList[0]).shape
                    for dim in range(len(dataShape)):
                        thisDim = np.arange(0, dataShape[dim])
                        dimStr = 'imgDim%d'%dim
                        coords[dimStr] = thisDim
                        dims.append(dimStr)

                    #now get data & create new data and new shape:
                    dataShp = np.array(imageList).shape
                    newShp = tuple(np.append(np.array(binShp), np.array(dataShp)[1:]))
                    data = np.array(imageList).reshape(newShp)
                    newKey = (droplet+'_image').replace('/','_')
                    dataArray = xr.DataArray(data, coords=coords, dims=dims,name=newKey)
                    cubeData = xr.merge([cubeData, dataArray])

            if 'array' in cube.dropletProc[droplet]['fct']:
                xArrayList=[]
                yArrayList=[]
                aduArrayList=[]
                nDrops=[]
                for key in (Bins[1:]+Bins[:-1])*0.5:
                    if key in cubeIdxData.groups.keys():
                        dropVar = cube.dropletProc[droplet]['vars']
                        npix = evtIDXr[dropVar[3]][cubeIdxData.groups[key]].data.flatten()
                        xArrayList.append(evtIDXr[dropVar[0]][cubeIdxData.groups[key]].data.flatten()[npix>0])
                        yArrayList.append(evtIDXr[dropVar[1]][cubeIdxData.groups[key]].data.flatten()[npix>0])
                        aduArrayList.append(evtIDXr[dropVar[2]][cubeIdxData.groups[key]].data.flatten()[npix>0])
                        nDrops.append((npix>0).sum())
                    else:
                        xArrayList.append([0])
                        yArrayList.append([0])
                        aduArrayList.append([0])
                #now filter zeros out and make or max size.
                xArray=np.zeros((len(xArrayList), np.nanmax(np.array(nDrops))))
                yArray=np.zeros((len(xArrayList), np.nanmax(np.array(nDrops))))
                aduArray=np.zeros((len(xArrayList), np.nanmax(np.array(nDrops))))
                for ia, (thisx, thisy, thisadu) in enumerate(zip(xArrayList, yArrayList, aduArrayList)):
                    xArray[ia,:len(thisx)] = thisx
                    yArray[ia,:len(thisx)] = thisy
                    aduArray[ia,:len(thisx)] = thisadu

                if len(cube.addBinVars.keys())>0:
                    dims = [cube.binVar]
                    coords={cube.binVar: (orgBins[1:]+orgBins[:-1])*0.5}
                    for addVar in cube.addBinVars.keys():
                        dims.append(addVar)
                        coords[addVar]=0.5*(cube.addBinVars[addVar][1:]+cube.addBinVars[addVar][:-1])                        
                    dimStr = 'dropDim%d'%dim
                    coords[dimStr] = np.arange(0, xArray[0].shape[0])
                    dims.append(dimStr)
                    #now get data & create new data and new shape:
                    newShp = tuple(np.append(np.array(binShp), xArray.shape[1:]))
                    xdata = np.array(xArray).reshape(newShp)
                    ydata = np.array(yArray).reshape(newShp)
                    adudata = np.array(aduArray).reshape(newShp)
                    dataArray = xr.DataArray(xdata, coords=coords, dims=dims,name=((droplet+'_array_X').replace('/','_')))
                    cubeData = xr.merge([cubeData, dataArray])
                    dataArray = xr.DataArray(ydata, coords=coords, dims=dims,name=((droplet+'_array_Y').replace('/','_')))
                    cubeData = xr.merge([cubeData, dataArray])
                    dataArray = xr.DataArray(adudata, coords=coords, dims=dims,name=((droplet+'_array_adu').replace('/','_')))
                    cubeData = xr.merge([cubeData, dataArray])

        h5Dict={}
        for key in cubeData.variables:
            h5Dict[key] = cubeData[key].values

        #write final file.
        if toHdf5 == 'h5netcdf':
            fname = '%s/Cube_%s_Run%03d_%s.nc'%(self.dirname,self.expname,self.run,cubeName)
            cubeData.to_netcdf(fname,engine='h5netcdf')
        elif toHdf5 == 'h5':
            fname = '%s/Cube_%s_Run%03d_%s.h5'%(self.dirname,self.expname,self.run,cubeName)
            for key in cubeData.variables:
                h5Dict[key] = cubeData[key].values
                dictToHdf5(fname, h5Dict)
            
        if not returnIdx and len(cube.addIdxVars)==0:
            return cubeData

        retDict={'keys': keys}
        retDict['fiducial']=fidArray
        retDict['evttime']=timeArray
        for iv,thisIdxVar in enumerate(cube.addIdxVars):
            retDict[thisIdxVar]=addArray[iv]

        return cubeData,retDict


    ##########################################################################
    ###
    ### functions for image treatment - starting w/ assembled 2-d image
    ###
    ##########################################################################

    def AvImage(self, detname='None', numEvts=100, nSkip=0, thresADU=0., thresRms=0.,useFilter=None, mean=False, std=False):
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
        if useFilter is not None:
            Filter = self.getFilter(useFilter=useFilter)
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
            from utilities_FitCenter import fitCircle
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
                circle = plt.Circle((res['xCen'],res['yCen']),res['R'],color='b',fill=False)
                plt.gca().add_artist(circle)
                plt.plot([res['xCen'],res['xCen']],[y.min(),y.max()],'r')
                plt.plot([x.min(),x.max()],[res['yCen'],res['yCen']],'r')

                if raw_input("Happy with this fit:\n") in ["y","Y"]:
                    happy = True
                print 'x,y: ',res['xCen'],res['yCen'],' R ',res['R']
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
            circleM = plt.Circle((res['xCen'],res['yCen']),res['R'],color='b',fill=False)
            fig=plt.figure(figsize=(10,10))
            #will need to check of rotation necessary here???
            #plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
            plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None',aspect='auto')
            plt.gca().add_artist(circleM)
            plt.plot([res['xCen'],res['xCen']],[y.min(),y.max()],'r')
            plt.plot([x.min(),x.max()],[res['yCen'],res['yCen']],'r')
            print 'x,y: ',res['xCen'],res['yCen'],' R ',res['R']
    
            plt.show()

    def MakeMask(self, detname=None):
        print ' not yet implemented, exists in SmallDataAna_psana.py'

    def azimuthalBinning(self, detname=None):
        print ' not yet implemented, exists in SmallDataProduced.py, uses code in xppmodules/src. Not sure if good idea'

    ##########################################################################
    ###
    ### functions for droplet analysis
    ###
    ##########################################################################

    def DropletCube(self, useFilter='', i0='ipm3/sum', rangeAdu=[], rangeX=[], rangeY=[], addName='', returnData=False, writeFile=False):
        data='DropletCube_Run%d_%s'%(self.run,addName)
        if useFilter!='':
            data+='_'+useFilter
        self.__dict__[data]=None

        #get basename of droplets
        dkey = [ key for key in self.Keys() if key.find('dropletsAdu')>=0]
        if len(dkey)==0:
            print 'did not find any droplets in this smallData file: ',self.fh5.filename
            return
        if len(dkey)>1:
            print 'we have the following options: ',dbkey
            basename=raw_input('type the name of the droplets to use:')
        else:
            basename = dkey[0].split('dropletsAdu')[0]

        #get filtered list of events
        i0_all = self.getVar(i0)
        if useFilter is not None:
            Filter = self.getFilter(useFilter=useFilter)
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

    def makeDataArray(self, varList=[], useFilter=None):
        if useFilter is None:
            arrayFilter = np.ones_like(self._tStamp).astype(bool)
        else:
            arrayFilter = self.getFilter(useFilter)
        timeFiltered = self._tStamp[arrayFilter]
        newXr = None
        #newXr = xr.DataArray(timeFiltered, coords={'time': timeFiltered}, dims=('time'),name='time')
        for tVar in varList:
            if not (self.hasKey(tVar)):
                continue
            #printR(rank, 'addvar: ',tVar,self.getVar(tVar,cubeFilter).shape)
            filteredVar = self.getVar(tVar,arrayFilter).squeeze()
            tVar=tVar.replace('/','__')
            if len(filteredVar.shape)==1:
                if newXr is None:
                    newXr = xr.DataArray(filteredVar, coords={'time': timeFiltered}, dims=('time'),name=tVar)
                else:
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
                if newXr is None:
                    newXr = xr.DataArray(filteredVar, coords=coords, dims=dims,name=tVar)
                else:
                    newXr = xr.merge([newXr, xr.DataArray(filteredVar, coords=coords, dims=dims,name=tVar)])

        return newXr


    def addEpicsArchiveVar(self, PVname, name=None, returnRaw=False):
        if self._epicsArchive is None:
            try:
                self._epicsArchive = EpicsArchive()
            except:
                print 'failed to create the EPICS archiver'
                return None

        #check if PV is present
        PVlist = self._epicsArchive.search_pvs(PVname, do_print=False)
        if len(PVlist)==0:
            print 'PV %s is not in archiver'%PVname
            return

        #get start&stop time to ask the archiver for.
        evtt = self.getVar('event_time')
        tStart = evtt.min()
        tStop = evtt.max()
        tStart_sec = tStart >> 32
        tStart_nsec = tStart - (tStart_sec << 32)
        tStop_sec = tStop >> 32
        tStop_nsec = tStop - (tStop_sec << 32)

        #now get the data.
        timePoints,valPoints = self._epicsArchive.get_points(PV=PVname, start=tStart_sec, end=tStop_sec, two_lists=True, raw=True)

        #add to data
        if name is None: name = PVname.replace(':','_')

        #make a dataset with values * time points
        tStamp_epics = np.array([np.datetime64(int(tsec), 's') for tsec in timePoints])
        da_epics = xr.DataArray(valPoints, coords={'time': tStamp_epics}, dims=('time'),name=name)

        #match with other dataset (event_time or fiducals: closest or interpolated values)
        #da.resample() #need to look more at docs to find out how that shoudl be used.
        binBoundaries = da_epics.time.copy()
        #bb = np.append(binBoundaries[0]-1, binBoundaries)
        #binBoundaries = np.append(bb, binBoundaries[-1]+1)
        binnedData = self.xrData.fiducials.time.groupby_bins('time',binBoundaries, labels=(binBoundaries[:-1].astype(int)),include_lowest=True, right=False)
        
        newArray=np.array([])
        for epicsTime, epicsValue in zip(np.array(da_epics.time), da_epics.data):
            try:
                nevents = len(binnedData.groups[epicsTime.astype(float)])
                newArray = np.append(newArray,np.ones(nevents)*epicsValue)
            except:
                if (epicsTime.astype(float)<binnedData.groups.keys()[0]):
                    nevents = len(self.xrData.fiducials.time.data[self.xrData.fiducials.time.data<epicsTime])
                    newArray = np.append(np.ones(nevents)*da_epics.data[0], newArray)
                elif (epicsTime.astype(float)>binnedData.groups.keys()[-1]):
                    nevents = len(self.xrData.fiducials.time.data[self.xrData.fiducials.time.data>epicsTime])
                    newArray = np.append(newArray, np.ones(nevents)*da_epics.data[-1])
                else:
                    print 'no group for time: ',epicsTime.astype(float)
                pass
                #print 'no data for key: ',epicsTime

        da_epics_new = xr.DataArray(newArray, coords={'time': self.xrData.fiducials.time}, dims=('time'),name=name)
        self.addVar(name, da_epics_new)
                
        if returnRaw:
            return da_epics
