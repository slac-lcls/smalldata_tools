#!/usr/bin/env python
import psana
import numpy as np
import holoviews as hv
hv.extension('bokeh')
import panel as pn
pn.extension()
import os
import argparse
import sys
import logging
try:
    basestring
except NameError:
    basestring = str
fpath=os.path.dirname(os.path.abspath(__file__))
fpathup = '/'.join(fpath.split('/')[:-1])
try:
    fpath = os.environ.get('MYDIR', fpathup).replace('/arp_scripts','')
except:
    fpath = fpathup
sys.path.append(fpath)
from smalldata_tools.utilities import rebin

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run', type=str, default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--postElog', help='post plot to elog', action='store_true', default=True)
parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
args = parser.parse_args()
logger.debug('Args to be used for pedestal plots: {0}'.format(args))

save_elog = args.postElog
expname = args.experiment
run = int(args.run)

detImgMaxSize = 500
det_name = 'Epix10kaQuad2'
dsname = 'exp={}:run={}:smd'.format(expname,run)
ds = psana.DataSource(dsname)
det = psana.Detector(det_name)

peds = det.pedestals(run)
try:
    peds_pre = det.pedestals(run-1)
except:
    peds_pre = None
noise = det.rms(run)
xcoords, ycoords = det.coords_xy(run)

xDim = hv.Dimension(('x','x in micron'))
xrange=(np.nanmin(xcoords), np.nanmax(xcoords))
yDim = hv.Dimension(('y','y in micron'))
yrange=(np.nanmin(ycoords), np.nanmax(ycoords))

def ped_rms_histograms(nCycles, peds, noise, diff):
    min5Ped=1e6
    max95Ped=0
    max95Noise=0
    min5Diff=1e6
    max95Diff=0
    for i in range(nCycles):
        if nCycles>1:
            thisPed = peds[i]
            thisNoise = noise[i]
            thisDiff = diff[i]
        else:
            thisPed = peds
            thisNoise = noise
            thisDiff = diff

        if np.nanpercentile(thisNoise,95) > max95Noise:
            max95Noise = np.nanpercentile(thisNoise,95)
        if np.nanpercentile(thisPed,95) > max95Ped:
            max95Ped = np.nanpercentile(thisPed,95)
        if np.nanpercentile(thisPed,5) < min5Ped:
            min5Ped = np.nanpercentile(thisPed,5)
        if np.nanpercentile(thisDiff,98) > max95Diff:
            max95Diff = np.nanpercentile(thisDiff,98)
        if np.nanpercentile(thisDiff,2) < min5Diff:
            min5Diff = np.nanpercentile(thisDiff,2)

    max95Noise*=1.25
    max95Ped*=1.25
    min5Ped*=0.9
    max95Diff*=1.25
    min5Diff*=0.8
    
    pedBins = np.linspace(min5Ped, max95Ped, 200)
    pedHistograms = []
    pedBinDim = hv.Dimension(('ped_bin','pedestal in ADU'), range=(min5Ped, max95Ped))
    noiseBins = np.linspace(0, max95Noise, 200)
    noiseHistograms = []
    noiseBinDim = hv.Dimension(('rms_bin','noise in ADU'), range=(0, max95Noise))
    diffBins = np.linspace(min5Diff, max95Diff, 200)
    diffHistograms = []
    diffBinDim = hv.Dimension(('diff_bin','diff in ADU'), range=(min5Diff, max95Diff))

    noiseMax=0
    pedMax=0
    diffMax=0
    for i in range(nCycles-2):
        if nCycles>1:
            thisNoise = noise[i]
            thisPed = peds[i]
            thisDiff = diff[i]
        else:
            thisPed = peds
            thisNoise = noise
            thisDiff = diff

        pedHistograms.append(np.histogram(thisPed.flatten(), pedBins))
        noiseHistograms.append(np.histogram(thisNoise.flatten(), noiseBins))
        diffHistograms.append(np.histogram(thisDiff.flatten(), diffBins))
        if pedMax < np.nanmax(pedHistograms[-1][0]):
            pedMax = np.nanmax(pedHistograms[-1][0])
        if noiseMax < np.nanmax(noiseHistograms[-1][0]):
            noiseMax = np.nanmax(noiseHistograms[-1][0])
        if diffMax < np.nanmax(diffHistograms[-1][0]):
            diffMax = np.nanmax(diffHistograms[-1][0])
    pedMax*=1.1 
    noiseMax*=1.1 
    diffMax*=1.1
    evtsPDim = hv.Dimension(('evtsP','N pixels / pedestal'), range=(0,pedMax))
    evtsDim = hv.Dimension(('evtsN','N pixels / noise'), range=(0,noiseMax))
    evtsDDim = hv.Dimension(('evtsD','N pixels / diff'), range=(0,diffMax))

    pedHists = []
    noiseHists = []
    diffHists = []
    i=0
    for pedH,noiseH,diffH in zip(pedHistograms, noiseHistograms, diffHistograms):
        pedHists.append(hv.Points((0.5*(pedBins[1:]+pedBins[:-1]), pedH[0]), label='Cycle %d'%i, 
                                    kdims=[pedBinDim, evtsPDim]) *
                         hv.Curve((0.5*(pedBins[1:]+pedBins[:-1]), pedH[0])))
        noiseHists.append(hv.Points((0.5*(noiseBins[1:]+noiseBins[:-1]), noiseH[0]), label='Cycle %d'%i, 
                                    kdims=[noiseBinDim, evtsDim]) *
                         hv.Curve((0.5*(noiseBins[1:]+noiseBins[:-1]), noiseH[0])))
        diffHists.append(hv.Points((0.5*(diffBins[1:]+diffBins[:-1]), diffH[0]), label='Cycle %d'%i, 
                                    kdims=[diffBinDim, evtsDDim]) *
                         hv.Curve((0.5*(diffBins[1:]+diffBins[:-1]), diffH[0])))
        i+=1
        
    return pedHists, noiseHists, diffHists

pedImgs=[]
rmsImgs=[]
diffImgs=[]
nCycles=1
if len(peds.shape)>3:
    nCycles=peds.shape[0]
for i in range(nCycles):
    if nCycles>1:
        thisPed = peds[i]
        thisNoise = noise[i]
        tpedDim = hv.Dimension(('ped_%d'%i,'pedestal in ADU, cycle %d'%i), range=(np.nanpercentile(thisPed,0.1), 
                                                            np.nanpercentile(thisPed,99.9)))
        trmsDim = hv.Dimension(('rms_%d'%i,'noise in ADU, cycle %d'%i), range=(np.nanpercentile(thisNoise,0.1), 
                                                         np.nanpercentile(thisNoise,99.9)))
        if peds_pre is not None:
            thisDiff = peds[i]-peds_pre[i]
            tdiffDim = hv.Dimension(('ped_%d'%i,'delta pedestal in ADU, cycle %d'%i), range=(np.nanpercentile(thisDiff,0.1), 
                                                            np.nanpercentile(thisDiff,99.9)))
    else:
        thisPed = peds
        thisNoise = noise
        tpedDim = hv.Dimension(('ped','pedestal in ADU'), range=(np.nanpercentile(thisPed,0.1), 
                                                            np.nanpercentile(thisPed,99.9)))
        trmsDim = hv.Dimension(('rms','noise in ADU'), range=(np.nanpercentile(thisNoise,0.1), 
                                                         np.nanpercentile(thisNoise,99.9)))
        if peds_pre is not None:
            thisDiff = peds-peds_pre
            tdiffDim = hv.Dimension(('ped_%d'%i,'delta pedestal in ADU'), range=(np.nanpercentile(thisDiff,0.1), 
                                                            np.nanpercentile(thisDiff,99.9)))

    pedImg = det.image(run,thisPed)    
    rmsImg = det.image(run,thisNoise)       
    diffImg = det.image(run,thisDiff)       
    if max(pedImg.shape[0], pedImg.shape[1])>detImgMaxSize:
        rebinFactor = float(detImgMaxSize)/max(pedImg.shape[0],pedImg.shape[1])
        pedImg = rebin(pedImg, [int(pedImg.shape[0]*rebinFactor), int(pedImg.shape[1]*rebinFactor)])
        rmsImg = rebin(rmsImg, [int(rmsImg.shape[0]*rebinFactor), int(rmsImg.shape[1]*rebinFactor)])
        diffImg = rebin(diffImg, [int(diffImg.shape[0]*rebinFactor), int(diffImg.shape[1]*rebinFactor)])
                        
    #print(pedImg.shape)
    pedImgs.append(hv.Image(pedImg, bounds=(xrange[0],yrange[0], xrange[1], yrange[1]),
         kdims=[xDim, yDim], vdims=[tpedDim], label='Pedestal, Cycle %d'%i).
         options(colorbar=True, aspect='equal',cmap='rainbow'))
    rmsImgs.append(hv.Image(rmsImg, bounds=(xrange[0],yrange[0], xrange[1], yrange[1]),
         kdims=[xDim, yDim], vdims=[trmsDim], label='Rms, Cycle %d'%i).
         options(colorbar=True, aspect='equal',cmap='rainbow'))
    diffImgs.append(hv.Image(diffImg, bounds=(xrange[0],yrange[0], xrange[1], yrange[1]),
         kdims=[xDim, yDim], vdims=[tdiffDim], label='Diff current-prevous, Cycle %d'%i).
         options(colorbar=True, aspect='equal',cmap='rainbow'))


gspec = pn.GridSpec(sizing_mode='stretch_both', max_width=500, name='Detector Images')
iwidth=3
iheight=3
gspec[0,0:(iwidth*3)] = pn.Row('# Pedestals&RMS - Run %04d'%run)
for i in range(nCycles):
    gspec[(i*iwidth+1):((i+1)*iwidth+1),(iheight*0):(iheight*1)] = pn.Column(pedImgs[i])
    gspec[(i*iwidth+1):((i+1)*iwidth+1),(iheight*1):(iheight*2)] = pn.Column(rmsImgs[i])
    gspec[(i*iwidth+1):((i+1)*iwidth+1),(iheight*2):(iheight*3)] = pn.Column(diffImgs[i])
#gspec

pedHists, noiseHists, diffHists = ped_rms_histograms(nCycles, peds, noise, (peds-peds_pre))

gspecH = pn.GridSpec(sizing_mode='stretch_width', max_width=500, name='Histogram')
gspecH[0,0:8] = pn.Row('# Pedestals&RMS Histograms - Run %04d'%run)
gspecH[1:4,0:8] = pn.Column(hv.Overlay(pedHists))
gspecH[4:7,0:8] = pn.Column(hv.Overlay(noiseHists))
gspecH[7:10,0:8] = pn.Column(hv.Overlay(diffHists))

tabs = pn.Tabs(gspecH)
tabs.append(gspec)

elogDir = '/reg/d/psdm/%s/%s/stats/summary/Pedestals_Run%03d'% (expname[0:3],expname,run)
if save_elog:
    import os
    if not os.path.isdir(elogDir):
        os.makedirs(elogDir)
    print('Made Directory to save data:', elogDir) 
    tabs.save(('%s/report.html'%elogDir))

