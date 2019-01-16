import tables
import numpy as np
import scipy
from scipy import signal as scipy_signal
import holoviews as hv
hv.extension('bokeh')

# ## Peak Shape and position from sum of two templates.

# Each peak is the sum of two templates, shifted by a bin, the
# floating point part of the main peak position corresponds to the
# ratio of the two peaks. The figure below allows you to see the
# effect of the floating point position (which is also needed for the
# optimization to work as written).


class templatePeakFit:
    def __init__(self, templateWaveform, nPeaks=1, saturationFraction=0.98, nMax=2, fitMethod='pah_trf', name=''):
        self.peakTemplateWf = None
        if isinstance(templateWaveform,np.ndarray):
            self.peakTemplateWf = templateWaveform
        elif isinstance(templateWaveform,list):
            self.peakTemplateWf = np.array(templateWaveform)
        elif isinstance(templateWaveform, basestring):            
            templateWaveform = '/reg/d/psdm/xpp/xpptut15/results/smallDataAnalysis/SingleAcqirisPeak_notebook.h5'
            try:
                peakTemplate = tables.open_file(templateFile).root.singlePicked_alt
                self.peakTemplateWf = peakTemplate[60:250].copy()
            except:
                print('could not find requested waveform in file: ',templateWaveform)
        if self.peakTemplateWf is None:
            print('We do not have a template waveform, will return None')
            return None

        #normalize waveform - set maximum to 1
        self.peakTemplateWf = self.peakTemplateWf/self.peakTemplateWf.max()
        
        self.nPeaks = nPeaks
        if name!='':
            self.name=name
        else:
            self.name='peak_%d'%self.nPeaks
        self.saturationFraction=saturationFraction
        self.nMax = nMax #allowed pixels above saturation fraction w/o clipping applied
        self.fitMethod = fitMethod #pah/sn _ trf/dogbox/lm
        self.saveParams = ['success','x', 'cost', 'fun']
        self.templateShapeMax=8000#Huh? figure out what 8000 was again, think shape of acqiris trace.

    def templateArray(self, args, templateShape):
        template =self.peakTemplateWf#[10:110]
        templateMaxPos = np.argmax(template)
        templateSum = np.zeros(self.templateShapeMax)
        for i in range(self.nPeaks):
            if args[i] < 0:
                print("nPeaks %d, nonsensical args[%d], bailing" %(self.nPeaks, i), args)
                return np.zeros(self.templateShapeMax)
            if (args[i]>templateMaxPos):
                templatePk = np.append(np.zeros(args[i]-templateMaxPos), template)
            else:
                templatePk = template[templateMaxPos-args[i]:]
            if (templateShape-templatePk.shape[0])>0:
                templatePk = np.append(templatePk, np.zeros(templateShape-templatePk.shape[0]))
            elif (templateShape-templatePk.shape[0])<0:
                templatePk = templatePk[:templateShape]
            templatePkp = np.append(np.array([0]), templatePk[:-1])
            frac1 = args[i+self.nPeaks]-int(args[i+self.nPeaks])
            templatep = templatePk*(1.-frac1)+templatePkp*frac1
            ##        if args[3]==0:
            ##            return template1*args[i+self.nPeaks]
            ##    print(args[1], )
            try:
                templateSum += templatep*args[i+self.nPeaks]
            except:
                "something unknown went wrong, peak %d, bailing" %(i)
                return np.zeros(self.templateShapeMax)
        return templateSum


    ##what to do with hardcoded numbers?
    ##make them members of the code.
    #def templatePlot(peakpos, templateShape=150):
    #    args = [peakpos, min(peakpos-250,100), min(peakpos-250,150), min(peakpos-250,200), 1., 1., 1., 1.]
    #    template = self.templateArray(args, templateShape, self.nPeaks)
    #    c1 = hv.Curve(template,('time','time'),('V','V')).options(width=300)
    #    c2 = hv.Curve(template,('Time','Time'),('Volt','Volt')).redim.range(Time=(17,32),Volt=(0.3,1.)).options(width=250)
    #    return c1+c2

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

    def fitTemplateLeastsq(self, trace, debug=False):
        ret_dict={}
        args0 = self.findPars(trace)
        ret_dict['initialGuess']=args0
        if debug: print('DEBUG: initial parameters for leastsq fit:', args0)
        maxTrace = np.nanmax(trace)*self.saturationFraction
        ret_dict['maxTrace']=maxTrace
        ret_dict['nmaxTrace']=(trace>maxTrace).sum()
        if debug: print('maxtrace: ',maxTrace,' n high pix ',(trace>maxTrace).sum() )
        resObj=None
        if self.fitMethod.split('_')[0] == 'sn':
            #my way to fit this.
            if (trace>maxTrace).sum() > self.nMax:
                maskSaturation=[trace<maxTrace]
            else:
                maskSaturation=np.ones_like(trace).astype(bool)
            errorfunction = lambda p: self.templateArray(p, trace.shape[0])[maskSaturation]-trace[maskSaturation]
            print('func is defined - silke ',self.fitMethod.split('_')[1])
            if self.fitMethod.split('_')[1]=='old':
                resObj, success = scipy.optimize.leastsq(errorfunction, args0)
            else:
                resObj = scipy.optimize.least_squares(errorfunction, args0, method=self.fitMethod.split('_')[1]) 
        elif self.fitMethod.split('_')[0] == 'pah':
            #philips way to fit this.
            minFunc = lambda p: self.clippedDelta(p, trace, maxTrace)
            resObj = scipy.optimize.least_squares(minFunc, args0)#, method=self.fitMethod('_')[1]) #chokes if I pass a method.
        else:
            print('this fit method is not defined', self.fitMethod)
        if resObj is not None:
            if isinstance(resObj, np.ndarray):
                ret_dict['fit_params']=resObj
            else:
                for param in self.saveParams:
                    try:
                        if param=='success':
                            ret_dict[param]=[int(getattr(resObj, param))]
                        elif param=='fun':
                            fun=list(getattr(resObj, param))
                            if len(fun)<self.templateShapeMax:
                                fun.append(np.zeros(self.templateShapeMax-len(fun)))
                            elif len(fun)>self.templateShapeMax:
                                fun=fun[:self.templateShapeMax]
                            ret_dict[param]=fun
                        else:
                            ret_dict[param]=getattr(resObj, param)
                    except:
                        pass
        return ret_dict

