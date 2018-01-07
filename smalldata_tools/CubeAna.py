import glob
from os import path
import tables
import numpy as np
from utilities_plotting import plotMarker, plotImage
import xarray as xr
#plot functions.
class CubeAna(object):
    def __init__(self, expname='', run=0, dirname='', cubeName='', cubeDict=None, plotWith='matplotlib', debug=False):
        print 'DEBUG INIT; ',expname, run
        self._fields={}
        self.expname=expname
        self._initrun=run
        self.runs=[run]
        self.runLabel='Run%03d'%run
        self.plotWith=plotWith
        self._inDict=cubeDict
        self.debug = debug
        if len(expname)>3:
            self.hutch=self.expname[:3]
            #do not look at ftc anymore. If that old data is used, 
            #pass the directory name in constructor
            if dirname=='':
                self.dirname='/reg/d/psdm/%s/%s/hdf5/smalldata'%(self.hutch,self.expname)
                self.plot_dirname='/reg/d/psdm/%s/%s/results/smalldata_plots/'%(self.hutch,self.expname)
            else:
                self.dirname=dirname
                self.plot_dirname = dirname+'/smalldata_plots'
            if not path.isdir(self.plot_dirname):
                makedirs(self.plot_dirname)

        self._fname=''
        allFiles = glob.glob('%s/Cube_%s_Run%03d_*.h5'%(self.dirname,self.expname,run))
        if cubeName != '':
            for thisFile in allFiles:
                if thisFile.find(cubeName)>=0:
                    self._fname = thisFile
                    break
            if self._fname == '':
                print 'Could not find a file with that name'
        if cubeName == '' or self._fname=='':
            cubeNames=[]
            for thisFile in allFiles:
                cubeNames.append(thisFile.replace('%s/Cube'%self.dirname,'').replace('_%s_Run%03d_'%(self.expname,run),'').replace('.h5',''))
            if len(cubeNames)==0 and self._inDict is None:
                print 'no cube files found and no cube dictionary passed in, return'
                return None
            if len(cubeNames)==1:
                self._fname = '%s/Cube_%s_Run%03d_%s.h5'%(self.dirname,self.expname,run, cubeNames[0])
                self._cubeName = cubeNames[0]
            else:
                print 'Options for Cube names are: ',cubeNames
                self._cubeName = raw_input('Select one of these options: \n ')
                self._fname = self._filenameFromRun(run)

        thisXrdata = None
        if self._inDict is not None:
            self._bins = thisDict['binVar_bins']
            thisXrdata = xr.DataArray(self._bins, coords={'bins': self._bins}, dms=('bins'), name='Bins')
            for dkey in thisDict.keys():
                if thisDict[dkey].shape[0] != bins.shape[0]:
                    print 'not same shape as bin variables'
                    continue
                dataShape, coords, dims = self._getXarrayDims(dkey, thisDict[dbkey].shape)
                thisXrdata = xr.merge([thisXrdata, xr.DataArray(thisDict[dkey], coords=coords, dims=dims,name=dkey) ])
        elif path.isfile(self._fname):
            thisXrdata = self._xrFromRun(run)
        else:
            print 'Could not find cubeFile %s nor was passed a dictionay'%self._fname
            thisXrdata = None
        self._cubeDict = { "Run%03d"%run: thisXrdata}
        
    def _getXarrayDims(self,key,dataShape):
        coords=None
        dims=None
        if len(dataShape)==1:
            coords={'bins': self._bins}
            dims=('bins')
        elif len(dataShape)==2:
            if dataShape[1]==1:
                coords={'bins': self._bins}
                dims=('bins')
            elif key.find('channels')>=0:
                dimArr = ['%d'%i for i in range(0,dataShape[1])]
                coords={'bins': self._bins,'channels':dimArr}
                dims=('bins','channels')
            elif key.find('AngPos')>=0:
                dimArr = ['AngX','AngY','PosX','PosY']
                coords={'bins': self._bins,'AngPos':dimArr}
                dims=('bins','AngPos')
            elif key.find('com')>=0:
                dimArr = ['axis0','axis1']
                coords={'bins': self._bins,'axes':dimArr}
                dims=('bins','axes')
            else: #currently save all 1-d data.
                dimArr = np.arange(0, dataShape[1])
                coords={'bins': self._bins,'dim0':dimArr}
                dims=('bins','dim0')
        else:
            coords={'bins': self._bins}
            dims = ['bins']
            for dim in range(len(dataShape)-1):
                thisDim = np.arange(0, dataShape[dim+1])
                dimStr = 'dim%d'%dim
                coords[dimStr] = thisDim
                dims.append(dimStr)

        return dataShape, coords, dims

    def _filenameFromRun(self, run):
        return '%s/Cube_%s_Run%03d_%s.h5'%(self.dirname,self.expname,run, self._cubeName)

    def _xrFromRun(self, run):
        fname = self._filenameFromRun(run)
        cubeTable=tables.open_file(fname,'r')
        self._bins = cubeTable.get_node('/binVar_bins').read()
        thisXrdata = xr.DataArray(self._bins, coords={'bins': self._bins}, dims=('bins'), name='Bins')
        self._detConfig=[]
        for node in cubeTable.root._f_list_nodes():
            key = node._v_pathname
            if key == '/binVar_bins':
                continue
            tArrName = key[1:]
            #make a dictionary with detector config data is present
            if key.find('Cfg')>=0:
                detname = key.split('_')[1]
                if not (detname in self._detConfig.keys()):
                    self._detConfig[detname]={}
                self._detConfig[detname][key.split('_')[2]] = cubeTable.get_node(key).read()
                continue
            dataShape = cubeTable.get_node(key).shape
            if dataShape[0]!=self._bins.shape[0]:
                if self.debug:
                    print 'data for %s is not of the right shape: '%(key), dataShape
                continue
            dataShape,coords, dims = self._getXarrayDims(key, dataShape)
            thisXrdata = xr.merge([thisXrdata, xr.DataArray(cubeTable.get_node(key).read().squeeze(), coords=coords, dims=dims,name=tArrName) ])
        return thisXrdata

    def Keys(self, name=None):
        keys=[]
        for key in self._cubeDict[0]:
            if name is not None and key.find(name)<0:
                continue
                keys.append(key)
        return keys

    def addDataFromRun(self, run):
        fname = self._filenameFromRun(run)
        if not path.isfile(fname):
            print 'Could not find file for run %d for cube %s, looked for: %s'%(run, self._cubeName, fname)
            return            
        thisXrdata = self._xrFromRun(run)
        self._cubeDict = { "Run%03d"%run: thisXrdata}

    # add support for defining ROIs and reducing as well as plotting as img
    def plotCube(self, run=None, sig=None, i0='nEntries', plotWith=None):
        if plotWith==None:
            plotWith=self.plotWith
        runKey = self._cubeDict.keys()[0]
        if run is not None and isinstance(run, int):
            runKey = 'Run%03d'%run
        elif  run is not None and isinstance(run, basestring):
            if run[:3]=='Run': 
                runKey = run
            else:
                runKey = 'Run%03d'%int(run)
                
        cubeDict = self._cubeDict[runKey]
        binVar = cubeDict['bins'].data
        if sig not in cubeDict.keys():
            print 'could not find sig, return'; return
            sigVar = cubeDict[sig].data
        if i0 is not None and i0 not in cubeDict.keys():
            print 'could not find i0, return'; return

        if i0 is None:
            i0Var = np.ones_like(sigVar)
            sigvar = '%s'%(sig)
        else:
            i0Var=cubeDict[i0].data
            sigvar = '%s/%s'%(sig,i0)
        sigVar =  np.divide(np.array(cubeDict[sig].data.T),np.array(i0Var))
        plotvar = 'binVar'

        if len(sigVar.shape)==1:
            plotMarker(sigVar, xData=[binVar], xLabel=plotvar, yLabel=sigvar, plotWith=plotWith, runLabel=runKey, plotTitle="%s vs %s for %s"%(sigvar, plotvar, runKey), plotDirname=self.plot_dirname)
        elif len(sigVar.shape)==2:
            extent=[binVar[0], binVar[-1],0,sigVar.shape[0]]
            plotImage(sigVar, xLabel=plotvar, plotWith=plotWith, runLabel=runKey, plotTitle="%s vs %s for %s"%(sigvar, plotvar, runKey), plotDirname=self.plot_dirname, extent=extent)
        else:
            print 'no code yet to plot 3-d data'
            return

        return None

