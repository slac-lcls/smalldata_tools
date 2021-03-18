import tables
import numpy as np
import scipy
from scipy import signal as scipy_signal
from smalldata_tools.DetObject import DetObjectFunc
from smalldata_tools.roi_rebin import spectrumFunc
from smalldata_tools.utilities import templateArray as utility_templateArray

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

