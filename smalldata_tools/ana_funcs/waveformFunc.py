import tables
import numpy as np
import itertools
import scipy
from scipy import signal as scipy_signal
from smalldata_tools.DetObjectFunc import DetObjectFunc
from smalldata_tools.ana_funcs.roi_rebin import spectrumFunc
from smalldata_tools.utilities import templateArray as utility_templateArray
from smalldata_tools.utilities_waveforms import hsdBaselineFourierEliminate
from smalldata_tools.utilities_waveforms import hitFinder_CFD

#find the left-most peak (or so.....)
class getCMPeakFunc(DetObjectFunc):
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','cmPeak')
        super(getCMPeakFunc, self).__init__(**kwargs)
        #need to get bin range
        #need to get ROI ? ROI w/ ROI func before!
        self.nPeak =  kwargs.get('nPeak',None)
        self.minPeakNum =  kwargs.get('minPeakNum',None)

    def setFromFunc(self, parentFunc=None):
        super(getCMPeakFunc, self).setFromFunc()
        if parentFunc is not None and isinstance(parentFunc, spectrumFunc):
            print('DEBUG: set bins for finding the peak....')
            self.bins = parentFunc.bins

    def process(self, data):
        retDict={}
        
        peakIdx = scipy_signal.find_peaks_cwt(data, widths=[3])
        maxIdx = [ x for _,x in sorted(zip(data[peakIdx],peakIdx), key=lambda pair: pair[0], reverse=True)]
        maxIdxFilter = [ x for x in maxIdx if data[max(0, x-2): min(x+3, data.shape[0])].sum()>self.minPeakNum ]

        if len(maxIdxFilter)==0:
            retDict['peak']=0
            retDict['peakCom']=0
            retDict['zeroPeak']=0
            return retDict

        leftPeak = min(maxIdxFilter[:self.nPeak])
        com=0.
        vals=0.

        for ibin in np.arange(leftPeak-2, leftPeak+3,1):
            com += data[ibin] * ibin #all in binnum coordinates.
            vals += data[ibin]

        retDict['peak']=leftPeak
        try:
            retDict['peakCom']=com/vals
        except:
            pass
        retDict['zeroPeak']=data[leftPeak]

        return retDict

# ## Peak Shape and position from sum of two templates.

# Each peak is the sum of two templates, shifted by a bin, the
# floating point part of the main peak position corresponds to the
# ratio of the two peaks. The figure below allows you to see the
# effect of the floating point position (which is also needed for the
# optimization to work as written).


class templateFitFunc(DetObjectFunc):
    def __init__(self, **kwargs):
        self.nPeaks = kwargs.get('nPeaks',1)
        self._name = kwargs.get('name','peak_%d'%self.nPeaks)
        super(templateFitFunc, self).__init__(**kwargs)
        self.template = kwargs.get('template',None)
        if isinstance(self.template,list):
            self.template = np.array(template)
        elif isinstance(self.template, basestring):    #why is this hardcoded?        
            templateFile = '/reg/d/psdm/xpp/xpptut15/results/smallDataAnalysis/SingleAcqirisPeak_notebook.h5'
            try:
                peakTemplate = tables.open_file(templateFile).root.singlePicked_alt
                self.template = peakTemplate[60:250].copy()
            except:
                print('could not find requested waveform in file: ',self.template)
        if self.template is None:
            print('We do not have a template waveform, will return None')
            return None

        #normalize waveform - set maximum to 1
        self.template = self.template/self.template.max()
        
        self.saturationFraction = kwargs.get('saturationFraction',0.98)
        self.nMax = kwargs.get('nMax',2) #allowed pixels above saturation fraction w/o clipping applied
        self._fitMethod = kwargs.get('fitMethod','pah_trf') #pah/sn _ trf/dogbox/lm
        self.saveParams = ['success','x', 'cost', 'fun']
        self._debug = kwargs.get('debug',False)
        self.fitShape = kwargs.get('fitShape',None)
        self.invert = kwargs.get('invert',False)
        self.baseline = kwargs.get('baseline',None)

    def templateArray(self, args, templateShape):
        return utility_templateArray(args, self.template, self.nPeaks, templateShape)

    def findPars(self, trace):
        """
        The initial estimate step finds the largest derivitive as peak position and uses the maximum value after the peak as guess for the amplitude.
        """
        difftrace = np.diff(trace)
        peakIdx = scipy_signal.find_peaks_cwt(difftrace, widths=[3])
        maxIdx = [ x for _,x in sorted(zip(difftrace[peakIdx],peakIdx), key=lambda pair: pair[0], reverse=True)]
        peakPos=[];peakVal=[]
        for iPk in range(self.nPeaks):
            idx=maxIdx[iPk]
            while difftrace[idx]>0:
                idx+=1
            peakPos.append(idx)
            peakVal.append(trace[idx])
        #sort in peak position.
        pPos = [ pos for pos,_ in sorted(zip(peakPos,peakVal))]
        pPk = [ pk for _,pk in sorted(zip(peakPos,peakVal))]
        for pVal in pPk:
            pPos.append(pVal)
        return pPos

    ##points > max contribute 0 to optimization var delta.
    def clippedDelta(self, par, trace, maxTrace):  ##PAH
        delta = self.templateArray(par, trace.shape[0])-trace
        clip = [trace>maxTrace]
        delta[clip] = np.clip(delta[clip], -666666666, 0)
        return delta

    #def fitTemplateLeastsq(self, trace, debug=False):
    def process(self, trace):
        trace=trace.squeeze()
        if self.fitShape is None:
            self.fitShape=trace.shape[0]
        elif self.fitShape>trace.shape[0]:
            trace=np.append(np.zeros(int(self.fitShape>trace.shape[0])), trace)
        elif self.fitShape<trace.shape[0]:
            print('templateFitFunc: truncate the input trace!', trace.shape, self.fitShape)
            trace=trace[:self.fitShape]
        if self.invert:
            trace *= -1.
        ret_dict={}
        if self.baseline is not None:
            try:
                traceBase = np.nanmedian(trace[self.baseline[0]:self.baseline[1]])
                trace = trace-traceBase
                ret_dict['base']=traceBase
            except:
                pass

        if trace.ndim>1:
            print('input data is not a waveform: ', trace.shape)
            return ret_dict
        args0 = self.findPars(trace)
        ret_dict['initialGuess']=np.array(args0)

        if self._debug: print('DEBUG: initial parameters for leastsq fit:', args0)
        maxTrace = np.nanmax(trace)*self.saturationFraction
        ret_dict['maxTrace']=maxTrace
        ret_dict['nmaxTrace']=(trace>maxTrace).sum()
        if self._debug: print('maxtrace: ',maxTrace,' n high pix ',(trace>maxTrace).sum() )
        resObj=None
        if self._fitMethod.split('_')[0] == 'sn':
            #my way to fit this.
            if (trace>maxTrace).sum() > self.nMax:
                maskSaturation=[trace<maxTrace]
            else:
                maskSaturation=np.ones_like(trace).astype(bool)
            errorfunction = lambda p: self.templateArray(p, trace.shape[0])[maskSaturation]-trace[maskSaturation]
            if self._debug: print('func is defined - silke ',self._fitMethod.split('_')[1])
            if self._fitMethod.split('_')[1]=='old':
                resObj, success = scipy.optimize.leastsq(errorfunction, args0)
            else:
                resObj = scipy.optimize.least_squares(errorfunction, args0, method=self._fitMethod.split('_')[1]) 
        elif self._fitMethod.split('_')[0] == 'pah':
            #philips way to fit this.
            minFunc = lambda p: self.clippedDelta(p, trace, maxTrace)
            resObj = scipy.optimize.least_squares(minFunc, args0)#, method=self._fitMethod('_')[1]) #chokes if I pass a method.
        else:
            print('this fit method is not defined', self._fitMethod)
        if resObj is not None:
            if isinstance(resObj, np.ndarray):
                ret_dict['fit_params']=resObj
            else:
                for param in self.saveParams:
                    try:
                        if param=='success':
                            ret_dict[param]=[int(getattr(resObj, param))]
                        elif param=='fun':
                            if isinstance(getattr(resObj, param), list):
                                #print 'getattr. fun ',getattr(resObj, param)
                                fun=np.array(getattr(resObj, param))
                            else:
                                fun=getattr(resObj, param)
                            if fun.shape[0]<self.fitShape:
                                fun=np.append(fun, np.array([0]*(self.fitShape-fun.shape[0])))
                            elif len(fun)>self.fitShape:
                                fun=fun[:self.fitShape]
                            ret_dict[param]=fun
                        else:
                            ret_dict[param]=getattr(resObj, param)
                    except:
                        pass

        return ret_dict

class hsdsplitFunc(DetObjectFunc):
    """
    function to rebin input data to new shape
    shape: desired shape of array
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','full')
        super(hsdsplitFunc, self).__init__(**kwargs)
        self.writeHsd = kwargs.get('writeHsd', True)
        self.hsdName = kwargs.get('hsdName', {})
        #print('init hsdname: ', hsdName)
        #self.hsdName = {}

    def setFromDet(self, det):
        super(hsdsplitFunc, self).setFromDet(det)
        self._cidx = det.cidx
        self._wfxlen = None
        self._det=det

    def process(self, data):
        pass_dict={}
        if self._wfxlen is None:
            self._wfxlen = self._det.wfxlen
            self._wfx = self._det.wfx
            #split the wfx array into separate pieces
            startidx=0
            for cidx,wfxlen in zip(self._cidx, self._wfxlen):
                timesName = 'times_%d'%cidx
                if ('hsd_%d'%cidx) in self.hsdName:
                    timesName = 'times_%s'%(self.hsdName['hsd_%d'%cidx])
                setattr(self, timesName, self._wfx[startidx:(startidx+wfxlen)])
                startidx+=wfxlen
                
        startidx=0
        for cidx,wfxlen in zip(self._cidx, self._wfxlen):
            cName = 'hsd_%d'%cidx
            if ('hsd_%d'%cidx) in self.hsdName:
                cName = self.hsdName['hsd_%d'%cidx]
            pass_dict[cName] = data[startidx:(startidx+wfxlen)]
            startidx+=wfxlen

        #save full waveforms if desired
        ret_dict = {}
        if self.writeHsd:
            ret_dict = { k:pass_dict[k] for k in pass_dict}

        for cidx,wfxlen in zip(self._cidx, self._wfxlen):
            tName = 'times_hsd_%d'%cidx
            tOrgName = 'times_%d'%cidx
            if ('hsd_%d'%cidx) in self.hsdName:
                tName = 'times_%s'%self.hsdName['hsd_%d'%cidx]
                tOrgName = 'times_%s'%self.hsdName['hsd_%d'%cidx]

            pass_dict[tName] = getattr(self, tOrgName)

        self.dat = pass_dict
        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]
        return ret_dict

class hsdROIFunc(DetObjectFunc):
    """
    function to rebin input data to new shape
    shape: desired shape of array
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','hsd_0__ROI')
        super(hsdROIFunc, self).__init__(**kwargs)
        ## later add different methods.
        self._hsdname = kwargs.get('name','hsd_0__ROI').split('__')[0]
        self.ROI = kwargs.get('ROI', [0,-1])
        self.writeArea = kwargs.get('writeArea', False)
        self.calcPars = kwargs.get('calcPars', True)            

    def setFromDet(self, det):
        super(hsdROIFunc, self).setFromDet(det)
        self._det=det

    def process(self, data):
        ret_dict={}
        pass_dict={}

        kdata = None
        if self._hsdname in data:
            kdata = data[self._hsdname]
        if kdata is not None:
            pass_dict['data'] = kdata[self.ROI[0]:self.ROI[1]]
        if self.writeArea:
            ret_dict = {'wf':pass_dict[k] for k in pass_dict}
        if self.calcPars:
            for k in pass_dict:
                ret_dict['sum'] = np.nansum(pass_dict[k])
                #ret_dict['%s_com'%(k)] = np.nansum(pass_dict[k]*np.arange(pass_dict[k].shape[0]))/np.nansum(pass_dict[k])

        self.dat = pass_dict
        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]
        return ret_dict


class hsdBaselineCorrectFunc(DetObjectFunc):
    """
    function to rebin input data to new shape
    shape: desired shape of array
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','baseline')
        #print('init baselinefunc')
        super(hsdBaselineCorrectFunc, self).__init__(**kwargs)
        ## later add different methods.
        self._method = kwargs.get('method','Fourier')

    def setFromDet(self, det):
        super(hsdBaselineCorrectFunc, self).setFromDet(det)
        self._det=det

    def process(self, data):
        ret_dict={}
        pass_dict={}
        if self._method == 'Fourier':
            if isinstance(data, np.ndarray):
                ret_dict['wf'] = hsdBaselineFourierEliminate(wf, self._det.wfx)
            elif isinstance(data, dict):
                for k in data:
                    if k.find('times')>=0: continue
                    if 'times_%s'%k in data:
                        ret_dict[k] = hsdBaselineFourierEliminate(data[k], data['times_%s'%k])
                        pass_dict['times_%s'%k] = data['times_%s'%k]
                    else:
                        ret_dict[k] = hsdBaselineFourierEliminate(data[k], self._det.wfx)
                        pass_dict['times_%s'%k] = dataself._det.wfx
                    pass_dict[k]=ret_dict[k]

        self.dat = pass_dict
        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]
        return ret_dict

class hitFinderCFDFunc(DetObjectFunc):
    """
    function to rebin input data to new shape
    shape: desired shape of array
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','hitFinderCFD')
        super(hitFinderCFDFunc, self).__init__(**kwargs)
        self.convFilterLength = kwargs.get('convFilterLength', 35)
        self.CFDOffset = kwargs.get('CFDOffset', 25)
        self.inverseMultiplier = kwargs.get('inverseMultiplier',-0.75)
        self.threshold = kwargs.get('threshold', 4)
        self.nmax_hits = kwargs.get('nmax_hits', 100)

    def setFromDet(self, det):
        super(hitFinderCFDFunc, self).setFromDet(det)
        self._det=det

    def hitFinder(self, trace):
        hitIndices = hitFinder_CFD(trace, self.convFilterLength, self.CFDOffset, 
                                   self.inverseMultiplier, self.threshold)
        nHits = len(hitIndices)
        if nHits>self.nmax_hits: 
            hitIndices = hitIndices[:self.nmax_hits]
        elif len(hitIndices)<self.nmax_hits: 
            hitIndices.extend([0]*(self.nmax_hits-nHits))
        return nHits, np.array(hitIndices)

    def process(self, data):
        ret_dict={}
        pass_dict={}
        #
        if isinstance(data, np.ndarray):
            nHits, hitIndices = self.hitFinder(data)
            ret_dict['nHits'] = nHits
            ret_dict['hitIndices'] = hitIndices

        elif isinstance(data, dict):
            for k in data:
                if k.find('times')>=0: continue
                nHits, hitIndices = self.hitFinder(data[k])
                ret_dict['%s_nHits'%k] = nHits
                ret_dict['%s_hitIndices'%k] = hitIndices

        #self.dat = pass_dict
        #subfuncResults = self.processFuncs()
        #for k in subfuncResults:
        #    for kk in subfuncResults[k]:
        #        ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]
        return ret_dict

class fimSumFunc(DetObjectFunc):
    """
    function to rebin input data to new shape
    shape: desired shape of array
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','fimSum')
        super(fimSumFunc, self).__init__(**kwargs)
        ## later add different methods.
        self._bkgROI = kwargs.get('bkgROI',None)
        if self._bkgROI is not None:
            if isinstance(self._bkgROI, slice):
                self.bkgROI=[self._bkgROI.start, self._bkgROI.stop]
            else:
                self.bkgROI = self._bkgROI.copy()
                self._bkgROI=slice(self.bkgROI[0], self.bkgROI[1])
        self._sigROI = kwargs.get('sigROI',[0,-1])
        if isinstance(self._sigROI, slice):
            self.sigROI=[self._sigROI.start, self._sigROI.stop]
        else:
            self.sigROI = self._sigROI.copy()
            self._sigROI=slice(self.sigROI[0], self.sigROI[1])

    def process(self, data):
        ret_dict={}
        bkgData=np.zeros(data.shape[0])
        if self._bkgROI is not None:
            bkgData = data[:,self._bkgROI].mean(axis=1)
        data = data-(bkgData.reshape((bkgData.shape[0],1)))
        dataS = data[:,self._sigROI]
        ret_dict['sum'] = data[:,self._sigROI].sum(axis=1)

        #self.dat = ret_dict
        #subfuncResults = self.processFuncs()
        #for k in subfuncResults:
        #    for kk in subfuncResults[k]:
        #        ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]
        return ret_dict

