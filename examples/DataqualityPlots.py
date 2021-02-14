#!/usr/bin/env python

import panel as pn
import h5py
import os
import argparse
import logging
import requests
import numpy as np
from requests.auth import HTTPBasicAuth
import holoviews as hv
from holoviews import dim
hv.extension('bokeh')
pn.extension()
import sys
fpath=os.path.dirname(os.path.abspath(__file__))
fpathup = '/'.join(fpath.split('/')[:-1])
try:
    fpath = os.environ.get('MYDIR', fpathup).replace('/arp_scripts','')
except:
    fpath = fpathup
sys.path.append(fpath)
from smalldata_tools.SmallDataAna_psana import SmallDataAna_psana as sdaps
from smalldata_tools.utilities import rebin

def postRunTable(runtable_data):
    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(args.experiment)
    print('URL:',ws_url)
    user=args.experiment[:3]+'opr'
    with open('/cds/home/opr/%s/forElogPost.txt'%user,'r') as reader:
        answer = reader.readline()
    r = requests.post(ws_url, params={"run_num": args.run}, json=runtable_data, auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
    #we might need to use this for non=current expetiments. Currently does not work in ARP
    #krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    #r = requests.post(ws_url, headers=krbheaders, params={"run_num": args.run}, json=runtable_data)
    print(r)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run', type=str, default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--stn', help='hutch station', type=int, default=0)
parser.add_argument('--nevents', help='number of events', type=int, default=1e9)
parser.add_argument('--directory', help='directory to read files from (def <exp>/hdf5/smalldata)', default=None)
parser.add_argument('--postElog', help='post plot to elog', action='store_true', default=True)
parser.add_argument('--scatterVar', help='variable for scattering signal if present', default=None)
parser.add_argument('--inclDropped', help='variable for scattering signal if present', action='store_true', default=False)
parser.add_argument('--postStats', help='post summary numbers to run tables', action='store_true', default=False)
#parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-kerb/lgbk/")
parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
args = parser.parse_args()
logger.debug('Args to be used for data quality plots: {0}'.format(args))


## Loading the data (single run)
save_elog = args.postElog
detImgMaxSize = 500 #max dimension of image.
expname = args.experiment
run = int(args.run)

if (int(os.environ.get('RUN_NUM', '-1')) > 0):
    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Dataquality Plots: </b>", "value": "Started"}])

if args.directory is not None:
    anaps = sdaps(expname,run, dirname=args.directory)
else:
    anaps = sdaps(expname,run)

ana = anaps.sda
#ana = sda(expname,run)

## Defining initial selection (laser-on events)
iniFilter='initial'
ana.addCut('lightStatus/xray',0.5,1.5,iniFilter)
ana.addCut('lightStatus/laser',0.5,1.5,iniFilter)

#effectively, keep all events
if args.inclDropped:
    ana.addCut('lightStatus/xray',-0.5,1.5,iniFilter)
    ana.addCut('lightStatus/laser',-0.5,1.5,iniFilter)
    
#droppled sthots.
ana.addCut('lightStatus/xray',-0.5,0.5,'xoff')
n162 = ana.getVar('evr/code_162').sum()
nOff = ana.getFilter('xoff').sum()
nEvents = ana.getFilter('xoff').shape[0]

#data to be posted to the run table if so requested.
runtable_data = {"N dropped Shots":nOff,
                 "N BYKIK 162":n162}

### Get data & define axis title&ranges.

if expname.find('xcs')>=0:
    ipmUpDim = hv.Dimension(('ipm4/sum','ipm4 Sum'))
    ipmDownDim = hv.Dimension(('ipm5/sum','ipm5 Sum'))
elif expname.find('xpp')>=0:
    ipmUpDim = hv.Dimension(('ipm2/sum','ipm2 Sum'))
    ipmDownDim = hv.Dimension(('ipm3/sum','ipm3 Sum'))
else:
    ipmUpDim = None
    ipmDownDim = None
eventTimeDim = hv.Dimension(('eventTimeR','relative event time'))
l3eDim = hv.Dimension(('l3e','L3 Energy'))

if args.scatterVar is not None:
    scatterDim = hv.Dimension((args.scatterVar,'scatter intensity'))
else:
    scatterDim = None

#timing vars.
nevtsLxtDim = hv.Dimension(('neventslxt','N events / lxt'))
lxtDim = hv.Dimension(('epics/lxt','lxt'))

ipmUpVar = ana.getVar(ipmUpDim.name,useFilter=iniFilter)
ipmDownVar = ana.getVar(ipmDownDim.name,useFilter=iniFilter)
ipmUpVar_percentiles = np.nanpercentile(ipmUpVar,[25,50,75])
runtable_data["%s 25 percentile"%(ipmUpDim.name)]=ipmUpVar[0]
runtable_data["%s 50 percentile"%(ipmUpDim.name)]=ipmUpVar[1]
runtable_data["%s 75 percentile"%(ipmUpDim.name)]=ipmUpVar[2]
runtable_data["%s 25 percentile"%(ipmDownDim.name)]=ipmDownVar[0]
runtable_data["%s 50 percentile"%(ipmDownDim.name)]=ipmDownVar[1]
runtable_data["%s 75 percentile"%(ipmDownDim.name)]=ipmDownVar[2]

eventTimeRaw = ana.getVar('event_time',useFilter=iniFilter)
if ana.getFilter(iniFilter).sum()==0:
    print('No events pass Filter, quit')
    if args.postStats:
        print('posting to the run tables - ipm values.')
        postRunTable(runtable_data)
    sys.exit()

eventTime = (eventTimeRaw>>32).astype(float)+((eventTimeRaw<<32)>>32).astype(float)*1e-9
eventTimeR = eventTime-eventTime[0]
if 'ebeam/L3_energy' in ana.Keys():
    l3eVar = ana.getVar('ebeam/L3_energy',useFilter=iniFilter)
else:
    l3eVar = None

eventTimeRMed = [np.nanmedian(eventTimeR[i*120:i*120+120]) for i in range(int(eventTimeR.shape[0]/120))]
ipmUpMed =  [np.nanmedian(ipmUpVar[i*120:i*120+120]) for i in range(int(eventTimeR.shape[0]/120))]
ipmDownMed =  [np.nanmedian(ipmDownVar[i*120:i*120+120]) for i in range(int(eventTimeR.shape[0]/120))]

if args.scatterVar is not None:
    scatterData = None
    try:
        scatterData = ana.getVar(args.scatterVar,useFilter=iniFilter)
        if len(scatterVar.shape)>1:
            scatterData = np.nanmean(scatterVar,axis=1)
    except:
        scatterData = None

### Scan Variable
nevtsDim = hv.Dimension(('nevents','N events / scan point'))
scanVar = ana.getScanName()
scanDim = None
isStepScan=False
if scanVar is not None and scanVar!='':
    isStepScan=True
if isStepScan:
    if isinstance(scanVar, basestring):
        scanVar=[scanVar]
    nPoints=[]
    for thisScanVar in scanVar:

        if scanDim == None:
            scanDim = hv.Dimension(('scan/%s'%scanVar[0],'%s'%scanVar[0]))
        nPoints.append(np.unique(self.getVar('scan/%s'%thisScanVar)).shape[0])
    if len(scanVar)==1:
        print('this run is a scan of %s with %d points'%(scanVar[0],nPoints[0]))

    try:
        stepVar = ana.getVar('scan/stepVar',useFilter=iniFilter)
        scanNsteps = np.bincount(stepVar)
        if args.scatterVar is not None and scatterData is not None:
            scanVarBins = np.bincount(stepVar,weights=scatterData)
        else:
            scanVarBins = None
    except:
        scanVarBins = None

### Fast delay stage 

lxt_fast_his = None
try:
    if 'enc/lasDelay' in ana.Keys():
        lxt_fast = ana.getVar('enc/lasDelay',useFilter=iniFilter)
        print(np.nanstd(lxt_fast))
        if lxt_fast is not None and np.nanstd(lxt_fast)<1e-4:
            lxt_fast_his = np.histogram(lxt_fast, np.linspace(np.nanpercentile(lxt_fast,1), np.nanpercentile(lxt_fast,99),100))
except:
    pass


#plots.
ipmUpTime = hv.HexTiles((eventTimeR[ipmUpVar<np.nanpercentile(ipmUpVar,99)],
                         ipmUpVar[ipmUpVar<np.nanpercentile(ipmUpVar,99)]),
                        kdims=[eventTimeDim, ipmUpDim]).\
                        opts(cmap='Blues')
ipmUpTimeMed = hv.Points((eventTimeRMed, ipmUpMed), kdims=[eventTimeDim,ipmUpDim],label=ipmUpDim.label).\
    options(color='r')
ipmDownTimeMed = hv.Points((eventTimeRMed, ipmDownMed), kdims=[eventTimeDim,ipmUpDim],label=ipmDownDim.label).\
    options(color='m')
    
ipmTimeLayout = ipmUpTime*ipmUpTimeMed*ipmDownTimeMed

if l3eVar is not None:
    treeSel = (l3eVar>np.nanpercentile(l3eVar,1))
    treePlot = hv.HexTiles((l3eVar[treeSel],ipmUpVar[treeSel]),kdims=[l3eDim,ipmUpDim])
else:
    treePlot = None

ipmPlot = hv.HexTiles((ipmUpVar, ipmDownVar), kdims=[ipmUpDim, ipmDownDim])
ipmLayout = ipmPlot.hist(dimension=[ipmUpDim.name,ipmDownDim.name])
if args.scatterVar is not None and scatterData is not None:
    ipmscatterPlot = hv.HexTiles((scatterData, ipmDownVar), kdims=[scatterDim, ipmDownDim])

stepPlot = None
if isStepScan:
    try:
        stepPlot = hv.Points((scanVarBins/scanNsteps, scanNsteps), kdims=[scanDim,nevtsDim])
    except:
        print('Failed to make stepPlot', np.nanmax(stepVar), scanNsteps.shape, scanVarBins.shape)
        print('DataPoints: ',scanVarBins/scanNsteps)
        print('Dimensions:', scanDim, nevtsDim)
if lxt_fast_his is not None:
    lxtPlot = hv.Points( (0.5*(lxt_fast_his[1][:-1]+lxt_fast_his[1][1:]), lxt_fast_his[0]), \
                             kdims=[lxtDim,nevtsLxtDim])
else:
    lxtPlot = None

gspec = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name='Data Quality - Run %d'%run)
gspec[0:2,0:8] = pn.Column(ipmTimeLayout)
gspec[2:5,0:4] = pn.Column(ipmLayout)
if treePlot is not None:
    gspec[2:5,4:8] = pn.Column(treePlot)

gspecS = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name='Scan&Scatter')
maxRow=0
if args.scatterVar is not None and scatterData is not None:
    gspecS[maxRow:4,maxRow:8] = pn.Column(ipmscatterPlot)
    maxRow=4
if stepPlot is not None:
    gspecS[maxRow,0:8] = pn.Row('## Scan Variable')
    gspecS[maxRow+1:maxRow+3,0:8] = pn.Column(stepPlot)
    maxRow+=3
if lxtPlot is not None:
    gspecS[maxRow,0:8] = pn.Row('## Laser - xray Timing')
    gspecS[maxRow+1:maxRow+3,0:8] = pn.Column(lxtPlot)
    maxRow+=3
if maxRow == 0: gspecS = None



if (int(os.environ.get('RUN_NUM', '-1')) > 0):
    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Dataquality Plots: </b>", "value": "Detector Plots"}])
# Detector stuff. 
detImgs=[]
detGrids=[]
for detImgName in ana.Keys('Sums'):
    image = ana.fh5.get_node('/%s'%detImgName)
    if max(image.shape[0], image.shape[1])>detImgMaxSize:
        rebinFactor = float(detImgMaxSize)/max(image.shape[0],image.shape[1])
        imageR = rebin(image, [int(image.shape[0]*rebinFactor), int(image.shape[1]*rebinFactor)])/ana.getVar('fiducials').shape[0]
    else:
        imageR = image/ana.getVar('fiducials').shape[0]
    #imgArrays.append(imageR/ana.getVar('fiducials').shape[0])
    imgDim = hv.Dimension(('image',detImgName.replace('Sums/','').replace('_calib_img',' Mean Image')),
                                    range=(np.nanpercentile(imageR,1), np.nanpercentile(imageR,99.)))
    detImgs.append(hv.Image(imageR, vdims=[imgDim], name=imgDim.label).options(colorbar=True, cmap='rainbow'))
        
    detGrid = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name=detImgName.replace('Sums/',''))
    detGrid[0,0] = pn.Row(detImgs[-1])
    detGrids.append(detGrid)

if nOff>100:
    for detImgName in ana.Keys('Sums'):
        detName = detImgName.replace('_calib_img','').replace('Sums/','')
        try:
            anaps.AvImage(detName, useFilter=offFilter, numEvts=min(1000, nOff))
        except:
            print('failed to get off shot data for detector %s'%detName)
            continue
        avData = anaps.getAvImage(detName)[1]
        try:
            image = anaps.__dict__[detName].det.image(run, avData)
        except:
            print('failed to make image for detector %s'%detName)
            continue
        if max(image.shape[0], image.shape[1])>detImgMaxSize:
            rebinFactor = float(detImgMaxSize)/max(image.shape[0],image.shape[1])
            imageR = rebin(image, [int(image.shape[0]*rebinFactor), int(image.shape[1]*rebinFactor)])
        else:
            imageR = image
        imgOffDim = hv.Dimension(('image_off',detImgName.replace('Sums/','').replace('_calib_img',' Mean Image Off')),
                                    range=(np.nanpercentile(imageR,1), np.nanpercentile(imageR,99.)))
        detImgs.append(hv.Image(imageR, vdims=[imgOffDim], name=imgOffDim.label).options(colorbar=True, cmap='rainbow'))
        
        detGrid = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name='%s, dropped shots'%detName)
        detGrid[0,0] = pn.Row(detImgs[-1])
        detGrids.append(detGrid)       

tabs = pn.Tabs(gspec)
if gspecS is not None:
    tabs.append(gspecS)
for detGrid in detGrids:
    tabs.append(detGrid)

elogDir = '/reg/d/psdm/%s/%s/stats/summary/DataQuality_Run%03d'%\
                (expname[0:3],expname,run)

if (int(os.environ.get('RUN_NUM', '-1')) > 0):
    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Dataquality Plots: </b>", "value": "Processed"}])
if save_elog:
    import os
    if not os.path.isdir(elogDir):
        os.makedirs(elogDir)
    print('Made Directory to save data:', elogDir)
    #gspec.save(('%s/report.html'%elogDir)) 
    tabs.save(('%s/report.html'%elogDir))

    if (int(os.environ.get('RUN_NUM', '-1')) > 0):
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Dataquality Plots: </b>", "value": "Posted"}])

if args.postStats:
    print('posting to the run tables - ipm values.')
    postRunTable(runtable_data)
