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
import requests
import socket
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
parser.add_argument('--pedImgs', help='make images of pedestals', action='store_true', default=True)
parser.add_argument('--pedDiffImgs', help='make images of first 10 subtracted images ', action='store_true', default=False)
parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
args = parser.parse_args()
logger.debug('Args to be used for pedestal plots: {0}'.format(args))

save_elog = args.postElog
make_ped_imgs = args.pedImgs
make_ped_data_imgs = args.pedImgs
expname = args.experiment
run = int(args.run)


def postRunTable(runtable_data):
    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(args.experiment)
    print('URL:',ws_url)
    user=args.experiment[:3]+'opr'
    with open('/cds/home/opr/%s/forElogPost.txt'%user,'r') as reader:
        answer = reader.readline()
    r = requests.post(ws_url, params={"run_num": args.run}, json=runtable_data, \
                      auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
    #we might need to use this for non=current expetiments. Currently does not work in ARP
    #krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    #r = requests.post(ws_url, headers=krbheaders, params={"run_num": args.run}, json=runtable_data)
    print(r)

def ped_rms_histograms(nCycles, peds, noise, diff, alias=''):
    min5Ped=1e6
    max95Ped=-1e6
    max95Noise=-1e6

    for i in range(nCycles):
        if nCycles>1:
            thisPed = peds[i]
            thisNoise = noise[i]
        else:
            thisPed = peds
            thisNoise = noise

        if np.nanpercentile(thisNoise,95) > max95Noise:
            max95Noise = np.nanpercentile(thisNoise,95)
        if np.nanpercentile(thisPed,95) > max95Ped:
            max95Ped = np.nanpercentile(thisPed,95)
        if np.nanpercentile(thisPed,5) < min5Ped:
            min5Ped = np.nanpercentile(thisPed,5)

    #for the UXI, most pixels read back as 0. Widen the range to look for decent min/max values.
    if max95Ped == min5Ped:
        for i in range(nCycles):
            if nCycles>1:
                thisPed = peds[i]
                thisNoise = noise[i]
            else:
                thisPed = peds
                thisNoise = noise
            if np.nanpercentile(thisNoise,99.9) > max95Noise:
                max95Noise = np.nanpercentile(thisNoise,99.9)
            if np.nanpercentile(thisPed,99.9) > max95Ped:
                max95Ped = np.nanpercentile(thisPed,99.9)
            if np.nanpercentile(thisPed,0.1) < min5Ped:
                min5Ped = np.nanpercentile(thisPed,0.1)

    max95Noise*=1.25
    max95Ped*=1.25
    min5Ped*=0.9
    
    pedBins = np.linspace(min5Ped, max95Ped, 200)
    pedHistograms = []
    pedBinDim = hv.Dimension(('ped_bin%s'%alias,'pedestal in ADU'), range=(min5Ped, max95Ped))
    noiseBins = np.linspace(0, max95Noise, 200)
    noiseHistograms = []
    noiseBinDim = hv.Dimension(('rms_bin%s'%alias,'noise in ADU'), range=(0, max95Noise))

    noiseMax=0
    pedMax=0
    for i in range(max(1, nCycles)):
        if nCycles>1:
            thisNoise = noise[i]
            thisPed = peds[i]
        else:
            thisPed = peds
            thisNoise = noise

        pedHistograms.append(np.histogram(thisPed.flatten(), pedBins))
        noiseHistograms.append(np.histogram(thisNoise.flatten(), noiseBins))
        if pedMax < np.nanmax(pedHistograms[-1][0]):
            pedMax = np.nanmax(pedHistograms[-1][0])
        if noiseMax < np.nanmax(noiseHistograms[-1][0]):
            noiseMax = np.nanmax(noiseHistograms[-1][0])
    pedMax*=1.1 
    noiseMax*=1.1 
    evtsPDim = hv.Dimension(('evtsP%s'%alias,'N pixels / pedestal'), range=(0,pedMax))
    evtsDim = hv.Dimension(('evtsN%s'%alias,'N pixels / noise'), range=(0,noiseMax))

    pedHists = []
    noiseHists = []
    i=0
    for pedH,noiseH in zip(pedHistograms, noiseHistograms):
        pedHists.append(hv.Points((0.5*(pedBins[1:]+pedBins[:-1]), pedH[0]), label='Cycle %d'%i, 
                                    kdims=[pedBinDim, evtsPDim]) *
                         hv.Curve((0.5*(pedBins[1:]+pedBins[:-1]), pedH[0])))
        noiseHists.append(hv.Points((0.5*(noiseBins[1:]+noiseBins[:-1]), noiseH[0]), label='Cycle %d'%i, 
                                    kdims=[noiseBinDim, evtsDim]) *
                         hv.Curve((0.5*(noiseBins[1:]+noiseBins[:-1]), noiseH[0])))
        i+=1

    if diff is None:
        return pedHists, noiseHists, None

    min5Diff=1e6
    max95Diff=-1e6

    for i in range(nCycles):
        if nCycles>1:
            thisDiff = diff[i]
        else:
            thisDiff = diff

        if np.nanpercentile(thisDiff,98) > max95Diff:
            max95Diff = np.nanpercentile(thisDiff,98)
        if np.nanpercentile(thisDiff,2) < min5Diff:
            min5Diff = np.nanpercentile(thisDiff,2)

    if max95Diff >= 0:
        max95Diff*=1.25
    else:
        max95Diff*=0.8
    if min5Diff >=0:
        min5Diff*=0.8
    else:
        min5Diff*=1.25
    
    diffBins = np.linspace(min5Diff, max95Diff, 200)
    diffHistograms = []
    diffBinDim = hv.Dimension(('diff_bin%s'%alias,'diff in ADU'), range=(min5Diff, max95Diff))

    diffMax=0
    for i in range(max(1, nCycles)):
        if nCycles>1:
            thisDiff = diff[i]
        else:
            thisDiff = diff

        diffHistograms.append(np.histogram(thisDiff.flatten(), diffBins))
        if diffMax < np.nanmax(diffHistograms[-1][0]):
            diffMax = np.nanmax(diffHistograms[-1][0])
    diffMax*=1.1
    evtsDDim = hv.Dimension(('evtsD%s'%alias,'N pixels / diff'), range=(0,diffMax))

    diffHists = []
    i=0
    for diffH in diffHistograms:
        diffHists.append(hv.Points((0.5*(diffBins[1:]+diffBins[:-1]), diffH[0]), label='Cycle %d'%i, 
                                    kdims=[diffBinDim, evtsDDim]) *
                         hv.Curve((0.5*(diffBins[1:]+diffBins[:-1]), diffH[0])))
        i+=1
        

    return pedHists, noiseHists, diffHists

def plotPedImgs(nCycles, det, run, peds, noise, peds_pre = None, detImgMaxSize=500, plotInfo=None, isLCLS2=False):
    pedImgs=[]
    rmsImgs=[]
    diffImgs=[]

    for i in range(nCycles):
        if nCycles>1:
            thisPed = peds[i]
            thisNoise = noise[i]
            tpedDim = hv.Dimension(('ped_%d'%i,'pedestal in ADU, cycle %d'%i), range=(np.nanpercentile(thisPed,0.1), 
                                                            np.nanpercentile(thisPed,99.9)))
            trmsDim = hv.Dimension(('rms_%d'%i,'noise in ADU, cycle %d'%i), range=(np.nanpercentile(thisNoise,0.1), 
                                                         np.nanpercentile(thisNoise,99.9)))
            if peds_pre is not None:
                try:
                    thisDiff = peds[i]-peds_pre[i]
                    tdiffDim = hv.Dimension(('ped_%d'%i,'delta pedestal in ADU, cycle %d'%i), range=(np.nanpercentile(thisDiff,0.1), 
                                                            np.nanpercentile(thisDiff,99.9)))
                except:
                    peds_pre = None
        else:
            thisPed = peds
            thisNoise = noise
            tpedDim = hv.Dimension(('ped','pedestal in ADU'), range=(np.nanpercentile(thisPed,0.1), 
                                                            np.nanpercentile(thisPed,99.9)))
            trmsDim = hv.Dimension(('rms','noise in ADU'), range=(np.nanpercentile(thisNoise,0.1), 
                                                         np.nanpercentile(thisNoise,99.9)))
            if peds_pre is not None:
                try:
                    thisDiff = peds-peds_pre
                    tdiffDim = hv.Dimension(('ped_%d'%i,'delta pedestal in ADU'), range=(np.nanpercentile(thisDiff,0.1), 
                                                                                         np.nanpercentile(thisDiff,99.9)))
                except:
                    peds_pre = None

        if not isLCLS2:
            pedImg = det.image(run,thisPed)    
            rmsImg = det.image(run,thisNoise)       
            if peds_pre is not None:
                diffImg = det.image(run,thisDiff)       
        else:
            pedImg = det.raw.image(run,thisPed)    
            rmsImg = det.raw.image(run,thisNoise)       
            if peds_pre is not None:
                diffImg = det.raw.image(run,thisDiff)       

        if pedImg is None:
            pedImg = thisPed
            rmsImg = thisNoise
            if peds_pre is not None:
                diffImg = thisDiff       
        if max(pedImg.shape[0], pedImg.shape[1])>detImgMaxSize:
            rebinFactor = float(detImgMaxSize)/max(pedImg.shape[0],pedImg.shape[1])
            pedImg = rebin(pedImg, [int(pedImg.shape[0]*rebinFactor), int(pedImg.shape[1]*rebinFactor)])
            rmsImg = rebin(rmsImg, [int(rmsImg.shape[0]*rebinFactor), int(rmsImg.shape[1]*rebinFactor)])
            if peds_pre is not None:
                diffImg = rebin(diffImg, [int(diffImg.shape[0]*rebinFactor), int(diffImg.shape[1]*rebinFactor)])
                        
        if plotInfo is not None:
            xrange, yrange, xDim, yDim = plotInfo
        if xrange is not None: xrange=[0,pedImg.shape[0]]
        if yrange is not None: yrange=[0,pedImg.shape[1]]
        pedImgs.append(hv.Image(pedImg, bounds=(xrange[0],yrange[0], xrange[1], yrange[1]),
            kdims=[xDim, yDim], vdims=[tpedDim], label='Pedestal, Cycle %d'%i).
                options(colorbar=True, aspect='equal',cmap='rainbow'))
        rmsImgs.append(hv.Image(rmsImg, bounds=(xrange[0],yrange[0], xrange[1], yrange[1]),
            kdims=[xDim, yDim], vdims=[trmsDim], label='Rms, Cycle %d'%i).
            options(colorbar=True, aspect='equal',cmap='rainbow'))
        if peds_pre is not None:
            diffImgs.append(hv.Image(diffImg, bounds=(xrange[0],yrange[0], xrange[1], yrange[1]),
                kdims=[xDim, yDim], vdims=[tdiffDim], label='Diff current-prevous, Cycle %d'%i).
                options(colorbar=True, aspect='equal',cmap='rainbow'))

    return pedImgs, rmsImgs, diffImgs

def plotDataImgs(expname, run, det_name, nCycles, plotInfo=None):
    import smalldata_tools.SmallDataAna_psana as sda
    anaps = sda.SmallDataAna_psana(expname,run, plotWith=None)

    gspecI = pn.GridSpec(sizing_mode='stretch_both', max_width=900, name='Pedestal data - subtracted')
    iwidth=3
    iheight=3
    gspecI[0,0:(iwidth*3)] = pn.Row('# Pedestal data - subtracted - Run %04d'%run)

    print(expname, run, det_name)
    for i in range(min(5,nCycles)):
        common_mode=None
        if det_name.find('Epix')>=0:
            common_mode=80
        anaps.AvImage(det_name, common_mode=common_mode, nSkip=1200*i, numEvts=10)

        retp = anaps.plotAvImage(returnIt=True, plotWith=None)
        imageDim = hv.Dimension(('ped_subtr','in keV'), range=(np.nanpercentile(retp,0.1), 
                                                                    np.nanpercentile(retp,99.9)))
        if plotInfo is not None:
            xrange, yrange, xDim, yDim = plotInfo
        if xrange is not None: xrange=[0,retp.shape[0]]
        if yrange is not None: yrange=[0,retp.shape[1]]
        timg = hv.Image(retp, bounds=(xrange[0],yrange[0], xrange[1], yrange[1]),
                        kdims=[xDim, yDim], vdims=[imageDim], label='Image, Cycle %d'%i).options(colorbar=True, aspect='equal',cmap='rainbow')
        row = 1+iheight*int(i*0.5)
        col1 = iwidth*(i%2)
        col2 = iwidth*(i%2+1)
        gspecI[row, col1:col2] = pn.Column(timg)

    return gspecI

def allPlots(det_name, run, make_ped_imgs=False, make_ped_data_imgs=False, tabs=None, 
             detImgMaxSize=400, isLCLS2=False):
    print('Working on plots for ',det_name)
    if not isLCLS2:
        det = psana.Detector(det_name)
        peds = det.pedestals(run)
        try:
            peds_pre = det.pedestals(run-1)
        except:
            peds_pre = None
        noise = det.rms(run)
        runnum=run
    else:    
        det = run.Detector(det_name)
        runnum = run.runnum
        detrawid = det.raw._uniqueid
        from psana.pscalib.calib.MDBWebUtils import calib_constants
        peds = calib_constants(detrawid, exp=expname, ctype='pedestals', run=runnum)[0]
        try:
            peds_pre = calib_constants(detrawid, exp=expname, ctype='pedestals', run=runnum-1)[0]
        except:
            peds_pre = None
        noise = calib_constants(detrawid, exp=expname, ctype='pixel_rms', run=runnum)[0]

        #snelson - debug....this call here is necessary. Not 100% sure why...
        evt = next(run.events())
        print('call raw image...',det.raw.image(evt, peds[0]).shape)
        #snelson end debug

    xDim = hv.Dimension(('x','x in micron'))
    yDim = hv.Dimension(('y','y in micron'))
    try:
        xcoords, ycoords = det.coords_xy(run)
        xrange=(np.nanmin(xcoords), np.nanmax(xcoords))
        yrange=(np.nanmin(ycoords), np.nanmax(ycoords))
    except:
        if len(noise.shape)==2:
            xmax=noise.shape[0]
            ymax=noise.shape[1]
        else:
            xmax=noise[0].shape[0]
            ymax=noise[0].shape[1]
        xrange=(0, xmax)
        yrange=(0, ymax)
        xDim = hv.Dimension(('x','x in pixel'))
        yDim = hv.Dimension(('y','y in pixel'))
        
    plotInfo = [xrange, yrange, xDim, yDim]
    
    nCycles=1
    if len(peds.shape)>=3 and det_name.find('CsPad')<0:
        nCycles=peds.shape[0]
    if nCycles > 5: nCycles = 5

    if peds_pre is not None:
        diffPeds = (peds-peds_pre)
    else:
        diffPeds = None
    pedHists, noiseHists, diffHists = ped_rms_histograms(nCycles, peds, noise, diffPeds, det_name)
    gspecH = pn.GridSpec(sizing_mode='stretch_width', max_width=500, name='Histogram - %s'%det_name)
    gspecH[0,0:8] = pn.Row('# Pedestals&RMS Histograms - Run %04d'%(runnum))
    gspecH[1:4,0:8] = pn.Column(hv.Overlay(pedHists))
    gspecH[4:7,0:8] = pn.Column(hv.Overlay(noiseHists))
    if diffHists is not None:
        gspecH[7:10,0:8] = pn.Column(hv.Overlay(diffHists))
    if tabs is None:
        tabs = pn.Tabs(gspecH)
        #this is for debugging.....
        #return tabs
    else:
        tabs.append(gspecH)
        
    if make_ped_imgs:
        if nCycles== 1:
            pedImgs, rmsImgs, diffImgs = plotPedImgs(nCycles, det, runnum, peds, noise, peds_pre, detImgMaxSize=detImgMaxSize,
                                                     plotInfo=plotInfo, isLCLS2=isLCLS2)
        else:
            pedImgs, rmsImgs, diffImgs = plotPedImgs(nCycles, det, runnum, peds, noise, detImgMaxSize=detImgMaxSize,
                                                plotInfo=plotInfo, isLCLS2=isLCLS2)
   
        gspec = pn.GridSpec(sizing_mode='stretch_both', max_width=1000, name='Det Imgs - %s'%det_name)
        iwidth=3
        iheight=3
        gspec[0,0:(iwidth*3)] = pn.Row('# Pedestals&RMS - Run %04d'%(runnum))
        if nCycles == 1:
            gspec[1:(1*iheight)+1,0:iwidth] = pn.Column(pedImgs[0])
            gspec[(1*iheight)+1:(2*iheight)+1,0:iwidth] = pn.Column(rmsImgs[0])
            if len(diffImgs)==len(pedImgs):
                gspec[(2*iheight)+1:(3*iheight)+1,0:iwidth] = pn.Column(diffImgs[0])
        else:
            for i in range(nCycles):
                gspec[(i*iheight+1):((i+1)*iheight+1),(iwidth*0):(iwidth*1)] = pn.Column(pedImgs[i])
                gspec[(i*iheight+1):((i+1)*iheight+1),(iwidth*1):(iwidth*2)] = pn.Column(rmsImgs[i])
                if len(diffImgs)==len(pedImgs):
                    gspec[(i*iheight+1):((i+1)*iheight+1),(iwidth*2):(iwidth*3)] = pn.Column(diffImgs[i])
        tabs.append(gspec)

    if make_ped_data_imgs:
        gspecI = plotDataImgs(det.env.experiment(), runnum, det.name.__str__(), nCycles, plotInfo=plotInfo)
        tabs.append(gspecI)
    
    return tabs

def getPixelStatus(det_name, run):
    det = psana.Detector(det_name)
    return det.status(run)

def countBadPixels(det_name, run):
    status = getPixelStatus(det_name, run)
    status_sum = status.sum(axis=0)
    return (status_sum>0).sum()

def plotPedestals(expname='mfxc00118', run=364, save_elog=False, make_ped_imgs=False, make_ped_data_imgs=False,
                 detImgMaxSize=400):
    isLCLS2=False
    if expname[:3] in ['tmo','rix','ued']: isLCLS2=True
    if not isLCLS2:
        ds_name = 'exp={}:run={}:smd'.format(expname,run)
        hostname = socket.gethostname()
        if hostname.find('drp')>=0:
            ds_name += ':dir=/cds/data/drpsrcf/%s/%s/xtc'%(expname[0:3],expname)
        ds = psana.DataSource(ds_name)

        det_names = [dn[0] for dn in psana.DetNames() if dn[0].find('Jungfrau')>=0 or dn[0].find('Epix')>=0 or dn[0].find('Cspad')>=0 or dn[0].find('Uxi')>=0]
        aliases = [dn[1] for dn in psana.DetNames() if dn[0].find('Jungfrau')>=0 or dn[0].find('Epix')>=0 or dn[0].find('Cspad')>=0 or dn[0].find('Uxi')>=0]
        runnum=run
    else:
        ds = psana.DataSource(exp=expname, run=run)
        thisrun = next(ds.runs())
        det_names = [dn for dn in thisrun.detnames if dn.find('epix')>=0]
        aliases = [dn for dn in thisrun.detnames if dn.find('epix')>=0]
        runnum=run
        run=thisrun

    tabs = None
    runTableData = {}
    for det_name, alias in zip(det_names, aliases):
        #print(det_name, alias)
        if alias == '':
            tabs = allPlots(det_name, run, make_ped_imgs=make_ped_imgs, 
                            make_ped_data_imgs=make_ped_data_imgs, tabs=tabs, detImgMaxSize=detImgMaxSize, isLCLS2=isLCLS2)
        else:
            tabs = allPlots(alias, run, make_ped_imgs=make_ped_imgs, 
                            make_ped_data_imgs=make_ped_data_imgs, tabs=tabs, detImgMaxSize=detImgMaxSize, isLCLS2=isLCLS2)
        runTableData[f'{det_name}_n_bad_pixels'] = countBadPixels(det_name, run)
      
    if save_elog:
        elogDir = '/reg/d/psdm/%s/%s/stats/summary/Pedestals/Pedestals_Run%03d'%(expname[0:3],expname,runnum)

        import os
        if not os.path.isdir(elogDir):
            os.makedirs(elogDir)
        print('Made Directory to save data:', elogDir)
        tabs.save(('%s/report.html'%elogDir))

        postRunTable(runTableData)
        
    return tabs


tabs = plotPedestals(expname,run, make_ped_imgs=True, make_ped_data_imgs=False, save_elog=save_elog)
