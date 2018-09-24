import glob
from os import path
import tables
import numpy as np
from utilities import rebinShape
from utilities_plotting import plotMarker, plotImage
from utilities_plotting import plot2d_from3d, plot3d_img_time
import xarray as xr
from bokeh.io import show
#plot functions: 3-d plotting, fix origin.
#example cube for xpptut - 2 runs or Tims data (azav,...)
#test function to add cubes & to sum runs & extend plot (azav as fct of run) 

import bokeh_utils
cmaps = bokeh_utils.get_all_mpl_palettes(allmaps=['jet','gray','cool','hot','seismic','CMRmap','nipy_spectral'])

class CubeAna(object):
    def __init__(self, expname='', run=0, dirname='', cubeName='', cubeDict=None, plotWith='matplotlib', debug=False):
        if debug: print 'DEBUG INIT; ',expname, run
        self._fields={}
        self.expname=expname
        self._initrun=run
        self.runs=[run]
        self.runLabel='Run%03d'%run
        self._plotWith=plotWith
        self._inDict=cubeDict
        self._debug = debug
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
                    self._cubeName = cubeName
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
            return None
        self._cubeDict = { "Run%03d"%run: thisXrdata}
        self._cubeSumDict = thisXrdata.copy(deep=True)
        self._cubeExtendDict = thisXrdata.copy(deep=True)
        self._variableDefs={} #dictionary of variable definitions: key is name, value is list of [varname, ROI]

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
            elif key.find('acq')>=0: 
                detname = key.split('__')[0]
                if (detname in self._detConfig.keys()) and ('wfx' in self._detConfig[detname].keys()) and (self._detConfig[detname]['wfx'].shape[0] == dataShape[1]):
                    dirArr =  self._detConfig[detname]['wfx']
                else:
                    dimArr = np.arange(0, dataShape[1])
                coords={'bins': self._bins,'time':dimArr}
                dims=('bins','time')                
            else: #currently save all 1-d data.
                dimArr = np.arange(0, dataShape[1])
                coords={'bins': self._bins,'dim0':dimArr}
                dims=('bins','dim0')
        else:
            coords={'bins': self._bins}
            dims = ['bins']
            detname = key.split('__')[0].replace('/','')
            if self._debug: print 'get dims&coords: ',key, key.split('__')[0], self._detConfig.keys()
            if (detname in self._detConfig.keys()) and ('x' in self._detConfig[detname].keys()) and ('y' in self._detConfig[detname].keys()):
                for idim,thisDimStr in enumerate( ['x','y']):
                    thisDim = np.linspace(self._detConfig[detname][thisDimStr].min(), self._detConfig[detname][thisDimStr].max(),dataShape[idim+1])
                    coords[thisDimStr] = thisDim
                    dims.append(thisDimStr)
            else:
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
        #first loop to fill _detConfig so that we can have meaningful dimensions where possible
        self._detConfig={}
        for node in cubeTable.root._f_list_nodes():
            key = node._v_pathname
            if self._debug: print 'DEBUG config: looking at key: ',key
            #make a dictionary with detector config data is present
            if key.find('Cfg')>=0:
                detname = key.split('_')[1]
                #if self._debug: print 'DEBUG config: ',self._detConfig
                if not (detname in self._detConfig.keys()):
                    self._detConfig[detname]={}
                if key.split('_')[2] not in self._detConfig[detname].keys():
                    self._detConfig[detname][key.split('_')[2]] = cubeTable.get_node(key).read()
        if self._debug: 
            for key in self._detConfig.keys():
                for kkey in self._detConfig[key].keys():
                    print 'DEBUG: config for %s is %s '%(key, kkey)
        #second loop for data.
        for node in cubeTable.root._f_list_nodes():
            key = node._v_pathname
            if self._debug: print 'DEBUG: looking at key for data: ',key
            if key == '/binVar_bins' or key.find('Cfg')>=0:
                continue
            if isinstance(node, tables.group.Group): #skip for now, should become loop over keys of group
                continue
            tArrName = key[1:]
            dataShape = cubeTable.get_node(key).shape
            if dataShape[0]!=self._bins.shape[0]:
                if self._debug:
                    print 'data for %s is not of the right shape: '%(key), dataShape
                continue
            dataShape,coords, dims = self._getXarrayDims(key, dataShape)
            thisXrdata = xr.merge([thisXrdata, xr.DataArray(cubeTable.get_node(key).read().squeeze(), coords=coords, dims=dims,name=tArrName) ])
            #print 'data has been added for ',key, dataShape
        return thisXrdata

    def _reduceData(self, inArray, sigROI=None):
        if sigROI is None or inArray is None:
            return inArray
        if isinstance(sigROI, np.ndarray):
            sigROI = sigROI.flatten().tolist()
        if len(inArray.shape)==1:
            return inArray
        if len(inArray.shape)==2:
            if not isinstance(sigROI, list):
                inArray=inArray[:,sigROI]
            elif len(sigROI)>1:
                inArray=inArray[:,sigROI[0]:sigROI[1]]
            else:
                inArray=inArray[:,sigROI[0]]
        elif len(inArray.shape)==3:
            if not isinstance(sigROI, list):
                inArray=inArray[:,sigROI]
            elif len(sigROI)==1:
                inArray=inArray[:,sigROI[0]]
            elif len(sigROI)==2:
                inArray=inArray[:,sigROI[0]:sigROI[1],:]
            else:
                inArray=inArray[:,sigROI[0]:sigROI[1],sigROI[2]:sigROI[3]]
        else:
            print 'this dimension is not yet implemented:',len(sig.shape)
        return inArray
        
    def setPlottingTool(self, plotWith):
        self._plotWith=plotWith

    def setDebug(self, debug):
        self._debug=bool(debug)

    def Keys(self, name=None, img=False):
        keys=[]
        for key in self._cubeDict[self._cubeDict.keys()[0]].variables:
            if name is not None and key.find(name)<0:
                continue
            if img and len(self._cubeDict[self._cubeDict.keys()[0]][key].shape)<3:
                continue
            keys.append(key)
        return keys

    def addDataFromRun(self, run):
        fname = self._filenameFromRun(run)
        if not path.isfile(fname):
            print 'Could not find file for run %d for cube %s, looked for: %s'%(run, self._cubeName, fname)
            return            
        thisXrdata = self._xrFromRun(run)
        self._cubeDict["Run%03d"%run] = thisXrdata
        for k in self._cubeSumDict.variables:
            if k == 'bins' or k == 'Bins': continue
            if self._cubeSumDict[k].shape[0]==self._cubeSumDict['Bins'].shape[0]:
                self._cubeSumDict[k].data = self._cubeSumDict[k].data + thisXrdata[k].data

    def extendData(self):
        if self._cubeExtendDict is not None and 'runs' in self._cubeExtendDict.dims and self._cubeExtendDict.runs.shape[0] == len(self._cubeDict.keys()):
            return
        
        extendRun=[]
        runKeys = self._cubeDict.keys()
        runKeys.sort()
        for runKey in runKeys:
            run = int(runKey.replace('Run',''))
            extendRun.append(run)

        extendKeys=[]
        extendCoords=[]
        extendDims=[]
        extendShapes=[]
        extendData={}
        allDims=[]
        for runKey in runKeys:
            thisCube = self._cubeDict[runKey]
            for ik,key in enumerate(thisCube.keys()):
                if key not in extendKeys:
                    extendKeys.append(key)
                    extendDims.append(thisCube[key].dims)
                    for d in thisCube[key].dims:
                        if d not in allDims:
                            allDims.append(d)
                    extendShapes.append(thisCube[key].shape)
                    extendData[ik] = thisCube[key].data
                else:
                    extendData[ik] = np.append(extendData[ik], thisCube[key].data)

        for ik,exShp in enumerate(extendShapes):
            extendShapes[ik] =  (len(runKeys),)+exShp
            extendDims[ik] = ('runs',)+extendDims[ik]
            #print 'ik coords ',ik,extendCoords[ik]
            #DEBUG/FIX ME: how do I add additional coordinates
            #make coords from dims &add run using dict.
            extendCoords[ik]={'runs':extendRun}
            for d in extendDims[ik]:
                extendCoords[ik][d] = thisCube[ik][d].data
        #now make new dictionary
        thisXrdata = xr.DataArray(self._bins, coords={'bins': self._bins}, dims=('bins'), name='Bins')
        for ik,k in enumerate(extendKeys):
            print 'reshaping ',extendData[ik].shape, extendShapes[ik]
        for ik,k in enumerate(extendKeys):
            #reject if k isdim
            if k in allDims:
                continue;
            print 'reshaping ',k,extendData[ik].shape, extendShapes[ik]
            thisXrdata = xr.merge([thisXrdata, xr.DataArray(extendData[ik].reshape(extendShapes[ik]), coords=extendCoords[ik], dims=extendDims[ik],name=k) ])
        self._cubeExtendDict = thisXrdata


    # add support for defining ROIs and reducing as well as plotting as img
    def plotCube(self, run=None, sig=None, i0='nEntries', plotWith=None, addData=None, nPoints=None, plotLog=False):
        if plotWith==None:
            plotWith=self._plotWith
        runKey = None
        if run is not None and isinstance(run, int):
            runKey = 'Run%03d'%run
        elif  run is not None and isinstance(run, basestring):
            if run[:3]=='Run': 
                runKey = run
            else:
                runKey = 'Run%03d'%int(run)
        if runKey is None:
            cubeDict = self._cubeSumDict
        else:
            cubeDict = self._cubeDict[runKey]
        binVar = cubeDict['bins'].data
        plotvar = 'binVar'

        if i0 is not None and i0 not in cubeDict.variables:
            print 'could not find i0, return'; return

        sigROI = None
        aimDim=0
        if isinstance(sig, list):
            if len(sig)>2:
                aimDim=sig[2]
            sigROI = sig[1]
            sig = sig[0]
        if sig not in cubeDict.variables:
            print 'could not find sig, return'; return

        if i0 is None:
            i0Var = np.ones_like(sigVar)
            sigvar = '%s'%(sig)
        else:
            i0Var=cubeDict[i0].data
            sigvar = '%s/%s'%(sig,i0)

        sigVar = cubeDict[sig].data
        sigVar = self._reduceData(sigVar, sigROI)   
        if aimDim>0:
            while len(sigVar.shape)>aimDim:
                sigVar = np.nanmean(sigVar,axis=1)
        if len(sigVar.shape)==1 and nPoints is not None and nPoints < sigVar.shape[0]:
            i0Var = rebinShape(i0Var,(nPoints,))
            sigVar = rebinShape(sigVar,(nPoints,))
            binVar = rebinShape(binVar,(nPoints,))
        sigVar = np.divide(np.array(sigVar.T),np.array(i0Var))

        #if addVar were dict of form {varname:ROI} woth ROI=None->no ROI
        #loop & get list of points to be passed to plotMarker.

        if plotWith is 'returnData':
            return sigVar
        if len(sigVar.shape)==1:
            if addData is None:
                if nPoints is not None and nPoints < sigVar.shape[0]:
                    sigVar = rebinShape(sigVar,(nPoints,))
                    binVar = rebinShape(binVar,(nPoints,))
                plotMarker(sigVar, xData=[binVar], xLabel=plotvar, yLabel=sigvar, plotWith=plotWith, runLabel=runKey, plotTitle="%s vs %s for %s"%(sigvar, plotvar, runKey), plotDirname=self.plot_dirname, line_dash='dashed')
            else:
                print 'plot2Data NOT SURE WHAT THIS WAS SUPPOSED TO ACTUALLY DO'
                plotMarker([sigVar, addData], xData=[binVar, binVar], xLabel=plotvar, yLabel=sigvar, plotWith=plotWith, runLabel=runKey, plotTitle="%s vs %s for %s"%(sigvar, plotvar, runKey), plotDirname=self.plot_dirname, line_dash='dashed')
        elif len(sigVar.shape)==2:
            extent=[binVar[0], binVar[-1],0,sigVar.shape[0]]
            plotImage(sigVar, xLabel=plotvar, plotWith=plotWith, runLabel=runKey, plotTitle="%s vs %s for %s"%(sigvar, plotvar, runKey), plotDirname=self.plot_dirname, extent=extent, plotLog=plotLog)
        else:
            print 'no code yet to plot 3-d data'
            return

        return None

    #change this: do not use what typically contains ROI to pick slice, add explicit value.
    def plotCubeImage(self, run=None, sig=None, plotWith=None, plotLog=False, plot3d=False, sigIdx=None):
        if plotWith==None:
            plotWith=self._plotWith
        runKey = None
        if run is not None and isinstance(run, int):
            runKey = 'Run%03d'%run
        elif  run is not None and isinstance(run, basestring):
            if run[:3]=='Run': 
                runKey = run
            else:
                runKey = 'Run%03d'%int(run)
        if runKey is None:
            cubeDict = self._cubeSumDict
        else:
            cubeDict = self._cubeDict[runKey]

        if sig is None:
            if len(self.Keys(img=True))>1:
                print 'we have the following cubed images: ',
                for k in self.Keys(img=True):
                    print k
                sig = raw_input('Please select one:')
            elif len(self.Keys(img=True))==1:
                sig = self.Keys(img=True)[0]
            else:
                print 'have no 2d image, return'
                return
        
        if isinstance(sig, list):
            if len(sig)>1:
                sigROI = sig[1] 
            else:
                sigROI=None
            sig = sig[0]
        else:
            sigROI=None
        if sig not in cubeDict.variables:
            print 'could not find sig, return'; return

        if plot3d:
            if plotWith=='matplotlib':
                print 'only supported with bokeh'
                return
            bins = self._cubeSumDict['bins'].data
            data2plot = cubeDict[sig].data
            data2plot = self._reduceData(data2plot, sigROI).copy()
            if sigIdx is not None:
                data2plot=data2plot[bins.shape[0]/2-sigIdx/2:bins.shape[0]/2+sigIdx/2]
                bins=bins[bins.shape[0]/2-sigIdx/2:bins.shape[0]/2+sigIdx/2]
            layout,im,p=plot2d_from3d(data2plot=data2plot,cmaps=cmaps,coord=bins.tolist(),init_plot=bins[bins.shape[0]/2])
            show(layout)
        else:
            if sigIdx is not None:
                data2plot=cubeDict[sig].data[sigIdx]
            else:
                data2plot=cubeDict[sig].data.sum(axis=0)

            data2plot = (self._reduceData(np.array([data2plot]), sigROI).copy()).squeeze()

            extent=[0,data2plot.shape[0],0,data2plot.shape[1]]
            plotImage(data2plot, xLabel='', plotWith=plotWith, runLabel=runKey, plotTitle="image for %s for run %s"%(sig, runKey), plotDirname=self.plot_dirname, extent=extent, plotLog=plotLog)

    #change this: do not use what typically contains ROI to pick slice, add explicit value.
    def inspectCube(self, run=None, sig=None, sigIdx=None, plotWith=None):
        if plotWith==None:
            plotWith=self._plotWith
            if plotWith=='matplotlib':
                print 'only supported with bokeh'
                return
            
        runKey = None
        if run is not None and isinstance(run, int):
            runKey = 'Run%03d'%run
        elif  run is not None and isinstance(run, basestring):
            if run[:3]=='Run': 
                runKey = run
            else:
                runKey = 'Run%03d'%int(run)
        if runKey is None:
            cubeDict = self._cubeSumDict
        else:
            cubeDict = self._cubeDict[runKey]

        if sig is None:
            if len(self.Keys(img=True))>1:
                print 'we have the following cubed images: ',
                for k in self.Keys(img=True):
                    print k
                sig = raw_input('Please select one:')
            elif len(self.Keys(img=True))==1:
                sig = self.Keys(img=True)[0]
            else:
                print 'have no 2d image, return'
                return
        
        if isinstance(sig, list):
            sigROI = sig[1] #this is not used yet...
            sig = sig[0]
        else:
            sigROI=None
        if sig not in cubeDict.variables:
            print 'could not find sig, return'; return

        bins = self._cubeSumDict['bins'].data
        data2plot = cubeDict[sig].data
        data2plot = self._reduceData(data2plot, sigROI).copy() #apply ROI
        #replace nan & inf for JavaScript
        data2plot[np.isinf(data2plot)]=np.nan 
        data2plot[np.isneginf(data2plot)]=np.nan 
        data2plot[np.isnan(data2plot)]=0

        if sigIdx is not None:
            if isinstance(sigIdx, int):
                data2plot=data2plot[bins.shape[0]/2-sigIdx/2:bins.shape[0]/2+sigIdx/2]
                bins=bins[bins.shape[0]/2-sigIdx/2:bins.shape[0]/2+sigIdx/2]
            elif len(np.array(sigIdx)):
                data2plot=data2plot[sigIdx[0]:sigIdx[1]]
                bins=bins[sigIdx[0]:sigIdx[1]]
                
        layout,im,p,p1d=plot3d_img_time(data2plot=data2plot,cmaps=cmaps,coord=bins.tolist(),init_plot=bins[bins.shape[0]/2])
        print 'dataIdx: ',data2plot.shape,' bins: ',bins
        show(layout)
        print 'done'
        
    def getReducedCube(self, aimShape=None, aimSize=4000000, runKey = None):
        print 'make reduced data'
        if runKey is None:
            cubeDict = self._cubeSumDict
        else:
            cubeDict = self._cubeDict[runKey]

        maxSize=0
        maxVar=''
        for thisVar in cubeDict.variables:
            if cubeDict[thisVar].size> maxSize:
                maxSize = cubeDict[thisVar].size
            maxVar = thisVar

        if maxSize <= aimSize:
            return cubeDict
            
        #too big, figure out how/if to rebin the bin array rather than image    
        if aimShape is None or aimShape[0]>cubeDict[thisVar].shape[0]:
            aimShape = (cubeDict[maxVar].shape[0],)

        binVar = cubeDict['nEntries'].bins
        nentriesVar = cubeDict['nEntries'].data
        if aimShape[0]<cubeDict[thisVar].shape[0]:
            binVar = rebinShape(cubeDict['bins'].data, (aimShape[0],))
            nentriesVar = rebinShape(cubeDict['nEntries'].data, (aimShape[0],))
        newXr = xr.DataArray(nentriesVar, coords={'bins': binVar}, dims=('bins'),name='nEntries')

        for thisVar in cubeDict.variables:
            locShape = aimShape
            coords={'bins': binVar}
            dims = ['bins']
            if thisVar=='nEntries' or thisVar=='bins':
                continue
            if cubeDict[thisVar].shape[0]!=cubeDict['nEntries'].shape[0]:
                if self._debug: 'skip here as likely coordinate: ',thisVar
                continue
            thisVarData = cubeDict[thisVar].data
            dataShape = thisVarData.shape
            if cubeDict[thisVar].size<=aimSize and locShape[0]>=cubeDict[thisVar].shape[0]:
                newXr = xr.merge([newXr, cubeDict[thisVar]]) #nothing changes.
            else: #need reshaping.
                if self._debug: print 'we need to reshape ',thisVar
                if len(cubeDict[thisVar].shape)==2:
                    locShape=(locShape[0], cubeDict[thisVar].shape[1])
                elif len(cubeDict[thisVar].shape)==3:
                    locShape=(locShape[0], cubeDict[thisVar].shape[1], cubeDict[thisVar].shape[2])
                elif len(cubeDict[thisVar].shape)==4:
                    locShape=(locShape[0], cubeDict[thisVar].shape[1], cubeDict[thisVar].shape[2], cubeDict[thisVar].shape[3])
                #check if I need to make image smaller.
                rebinFac=1.
                if cubeDict[thisVar].size*locShape[0]/cubeDict[thisVar].shape[0]>aimSize:
                    if self._debug:
                        print 'make image smaller. deal with coordinates!'
                    imgSize = maxSize/cubeDict[thisVar].shape[0]
                    rebinFac = np.sqrt(1.*aimSize/locShape[0]/imgSize)
                    locShape = (locShape[0], int(cubeDict[thisVar].shape[1]*rebinFac), int(cubeDict[thisVar].shape[2]*rebinFac))                
                if self._debug: print 'rebin: ',cubeDict[thisVar].data.shape, '***',locShape
                thisVarData = rebinShape(cubeDict[thisVar].data, locShape)
                #
                if len(cubeDict[thisVar].shape)>1:
                    dataShape = thisVarData.shape

                    for ic,imgCoord in enumerate(cubeDict[thisVar].coords.keys()):
                        if imgCoord=='bins': continue
                        thisDim = rebinShape(cubeDict[imgCoord].data, (locShape[ic],))
                        dimStr = imgCoord

                    #for dim in range(len(dataShape)-1):
                    #    thisDim = np.arange(0, dataShape[dim+1])
                    #    dimStr = '%s_dim%d'%(thisVar,dim)
                        coords[dimStr] = thisDim
                        dims.append(dimStr)
                    newXr = xr.merge([newXr, xr.DataArray(thisVarData, coords=coords, dims=dims,name=thisVar)])
                else:
                    newXr = xr.merge([newXr, xr.DataArray(thisVarData, coords={'bins': binVar}, dims=('bins'),name=thisVar) ])       

        return newXr


    def plotReducedCube(self, sig=None, i0='nEntries', inCube=None, aimShape=None, aimSize=4000000, runKey = None):
        if inCube is None:
            inCube = self.getReducedCube(aimShape=aimShape, aimSize=aimSize, runKey=runKey)
            
        #now get possible images.
        if sig is not None:
            if not sig in inCube.keys():
                print 'Variable %s is not in cube, check for all options.'%sig
                sig = None
        if sig is None:
            imgList=[]
            for key in inCube.variables:
                if len(inCube[key].shape)>2:
                    for k in self.Keys(img=True):
                        imgList.append(k)
            if len(imgList)==1:
                sig = imgList[0]
            else:
                print 'we have the following cubed images: ',
                for k in imgList: print k
                sig = raw_input('Please select one:')
            if not sig in inCube.variables:
                print 'Well, this did not work. Likely we had no images or you made a typo. We give up for now'
                return
        
        #get i0
        sigData = inCube[sig].data
        i0Data = inCube[i0].data
        data2plot = np.divide(sigData.T,i0Data).T
        #replace nan & inf for JavaScript
        data2plot[np.isinf(data2plot)]=np.nan 
        data2plot[np.isneginf(data2plot)]=np.nan 
        data2plot[np.isnan(data2plot)]=0

        #get bins
        bins = inCube[sig].bins.data
            
        #plot stuff.
        layout,im,p,p1d=plot3d_img_time(data2plot=data2plot,cmaps=cmaps,coord=bins.tolist(),init_plot=bins[bins.shape[0]/2])
        show(layout)
