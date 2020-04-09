import tables
import numpy as np
import scipy
from scipy import signal as scipy_signal
import holoviews as hv
hv.extension('bokeh')

nPeaksToFit = 4

templateFile = '/reg/d/psdm/xpp/xpptut15/results/smallDataAnalysis/SingleAcqirisPeak_notebook.h5'
peakTemplate = tables.open_file(templateFile).root.singlePicked_alt

peakMain = peakTemplate[60:250].copy()
peakMain = peakMain/peakMain.max()

## Produce the template for an n-peak waveform

def templateArray(args, templateShape, nPeaks=4):
    template = peakMain#[10:110]
    templateMaxPos = np.argmax(template)
    templateSum = np.zeros(8000)
    for i in range(nPeaks):
        if args[i] < 0:
            print "nPeaks %d, nonsensical args[%d], bailing" %(nPeaks, i), args
            return np.zeros(8000)
        if (args[i]>templateMaxPos):
            templatePk = np.append(np.zeros(args[i]-templateMaxPos), template)
        else:
            templatePk = template[templateMaxPos-args[i]:]
        if (templateShape-templatePk.shape[0])>0:
            templatePk = np.append(templatePk, np.zeros(templateShape-templatePk.shape[0]))
        elif (templateShape-templatePk.shape[0])<0:
            templatePk = templatePk[:templateShape]
        templatePkp = np.append(np.array([0]), templatePk[:-1])
        frac1 = args[i+nPeaks]-int(args[i+nPeaks])
        templatep = templatePk*(1.-frac1)+templatePkp*frac1
##        if args[3]==0:
##            return template1*args[i+nPeaks]
##    print args[1], 
        try:
            templateSum += templatep*args[i+nPeaks]
        except:
            "something unknown went wrong, peak %d, bailing" %(i)
            return np.zeros(8000)
    return templateSum


# ## Peak Shape and position from sum of two templates.

# Each peak is the sum of two templates, shifted by a bin, the
# floating point part of the main peak position corresponds to the
# ratio of the two peaks. The figure below allows you to see the
# effect of the floating point position (which is also needed for the
# optimization to work as written).


def templatePlot(peakpos, templateShape=150):
    args = [peakpos, min(peakpos-250,100), min(peakpos-250,150), min(peakpos-250,200), 1., 1., 1., 1.]
    template = templateArray(args, templateShape, nPeaksToFit)
    c1 = hv.Curve(template,('time','time'),('V','V')).options(width=300)
    c2 = hv.Curve(template,('Time','Time'),('Volt','Volt')).redim.range(Time=(17,32),Volt=(0.3,1.)).options(width=250)
    return c1+c2

##print hv.Dimension('peakpos',range=(19.75,21.25))

##dmap = hv.DynamicMap(templatePlot, kdims=[hv.Dimension('peakpos',range=(19.75,21.25))])
##dmap


# ## Functions for acquiris trace template fitting
# The initial estimate step finds the largest derivitive as peak position and uses the maximum value after the peak as guess for the amplitude.

def findPars(trace, nPeak=4):
    difftrace = np.diff(trace)
    peakIdx = scipy_signal.find_peaks_cwt(difftrace, widths=[3])
    maxIdx = [ x for _,x in sorted(zip(difftrace[peakIdx],peakIdx), key=lambda pair: pair[0], reverse=True)]
    peakPos=[];peakVal=[]
    for iPk in range(nPeak):
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

##no longer quite exclude points > max from calc
def clippedDelta(par, trace, maxTrace):  ##PAH
    delta = templateArray(par, trace.shape[0], nPeaksToFit)-trace
    clip = [trace>maxTrace]
    delta[clip] = np.clip(delta[clip], -666666666, 0)
    return delta

def fitTemplateLeastsq(trace, saturationFrac=0.98, debug=False):
    args0 = findPars(trace, nPeaksToFit)
    print "args0:", args0
    maxTrace = np.nanmax(trace)*saturationFrac
    if debug: print 'maxtrace: ',maxTrace,' n high pix ',(trace>maxTrace).sum() 
    if (trace>maxTrace).sum() > 2:
        maskSaturation=[trace<maxTrace]
    else:
        maskSaturation=np.ones_like(trace).astype(bool)
    errorfunction = lambda p: templateArray(p, trace.shape[0], nPeaksToFit)[maskSaturation]-trace[maskSaturation]
    minFunc = lambda p: clippedDelta(p, trace, maxTrace)
    resObj = scipy.optimize.least_squares(minFunc, args0)
##    boundLow = [args0[0]-1, args0[1]-1, 0, 0]
##    boundHigh = [args0[0]+1, args0[1]+1, np.inf, np.inf]
##    resObj = scipy.optimize.least_squares(errorfunction, args0)  ## calls Silke's cost function
##    resObj = scipy.optimize.least_squares(minFunc, args0, bounds=(boundLow, boundHigh))  ## good way to crash
    return resObj

expname='xcsx31116'
dirname='/reg/d/psdm/%s/%s/scratch/prod_July18/' %(expname[0:3],expname)

## PAH stuff mostly developed using run 106
#run=101
#run=102
run = 103
run = 102
run = 106
##run = 108
fname='%s_Run%d.h5'%(expname, run)
Dat = tables.open_file('%s/%s'%(dirname, fname)).root
acqTraces = Dat.acq01.ROI.read()
framesA = Dat.uxi.frameA.read()
framesB = Dat.uxi.frameB.read()
l3Energy = Dat.ebeam.L3_energy.read()

# # Fit some example curves.

# In[7]:


xRange=(3100,3350)
#print findPars(trace1)
dimx=('xt','Time (trace1)')
dimy=('yt','Volt (trace1)')
dimx2=('xt2','Time (trace2)')
dimy2=('yt2','Volt (trace2)')
blob = []
residMax = 0
costMax = 0
indR = 666
indC = 666
for aT in range(0, 6, 2):
##for aT in range(0, len(acqTraces), 2):
    print "trace", aT
    trace = acqTraces[aT].copy()
    trace[50:] += acqTraces[aT+1][:-50]
    fitObj = fitTemplateLeastsq(trace)
    fitPars = fitObj.x
    cost = fitObj.cost
    estCurve = templateArray(findPars(trace, nPeaksToFit), trace.shape[0], nPeaksToFit)
    fitCurve = templateArray(fitPars, trace.shape[0], nPeaksToFit)
##    resid = fitObj.fun[3125:3350]  ## this was really chi**2/dof
    if cost>costMax:
            costMax = cost
            indC = aT
    blob.append((trace, fitPars, estCurve, fitCurve, cost))
    hv.Curve(trace,dimx,dimy,label='data').redim.range(xt=xRange).options(width=400) *    hv.Curve(estCurve,dimx,dimy,label='estimate').redim.range(xt=xRange) *    hv.Curve(fitCurve,dimx,dimy,label='fit').redim.range(xt=xRange) +    hv.Image(framesA[indC]) + hv.Image(framesB[indC])

    print aT,
print

plot(range(3100,3350), trace[3100:3350])
plot(range(3100,3350), fitCurve[3100:3350])
