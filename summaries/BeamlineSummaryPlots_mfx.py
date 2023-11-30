#!/usr/bin/env python
########################################################################
## load tools and packages required for data handling and plotting
########################################################################
import sys
import os
import argparse
import logging
import socket
import requests
from requests.auth import HTTPBasicAuth
from typing import Optional, Tuple

import h5py
import numpy as np
import panel as pn
import holoviews as hv
from holoviews import dim
hv.extension('bokeh')
pn.extension()

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
from smalldata_tools.SmallDataAna_psana import SmallDataAna_psana as sdaps
from smalldata_tools.utilities import image_from_dxy
from smalldata_tools.utilities import rebin
from summaries.elog_utilities import *

def postDetectorDamageMsg(
        detectors: list,
        exp: str,
        run: int,
        *,
        tag: str = "SUMMARY_EVENT_DAMAGE",
        title: str = "EVENT DAMAGE INFO",
        smd_dir: Optional[str] = None,
        post_thresh: float = 0.1
) -> None:
    """Post detector event damage info for a specified run to the eLog.

    Parameters
    ----------
    detectors (list[str]) Names of detectors to report on.
    exp (str) Experiment name.
    run (int) Run number. Usually the current run.
    tag (str) Optional. Tag for the event damage summary posts.
    title (str) Optional. Title for event damage summary posts.
    smd_dir (str) Optional. Alternative directory for smalldata HDF5 files.
    post_thresh (float) Optional. Damage threshold (as a percentage)
        required to post a message to eLog. At least 1 detector must pass
        the threshold to post to eLog. Only detectors passing the threshold
        will be included.
    """
    if smd_dir:
        anaps = sdaps(exp, run, dirname=smd_dir)
    else:
        anaps = sdaps(exp, run)

    ana = anaps.sda

    table_header: str = (
        "<thead><tr><th colspan=\"3\">"
        f"<center>{title}</center>"
        "</th></tr></thead>"
    )
    table_body: str = (
        "<tbody><tr>"
        "<td><b><center>Detector</center></b></td>"
        "<td><b><center>Missing/Damaged Events</center></b></td>"
        "<td><b><center>Percentage Missing/Damaged</center></b></td></tr>"
    )

    post_msg: bool = False

    for det_name in detectors:
        damage_var: np.ndarray = ana.getVar(f"damage/{det_name}")
        damage: int = len(damage_var[damage_var == 0])
        dmg_percent: float = damage/len(damage_var)
        if dmg_percent > post_thresh:
            post_msg = True
            det_entry: str = (
                f"<tr><td><center>{det_name}</center></td>"
                f"<td><center>{damage}</center></td>"
                f"<td><center>{dmg_percent:.2%}</center></td></tr>"
            )
            table_body += det_entry
    table_body += "</tbody>"
    msg: str = f"<table border=\"1\">{table_header}{table_body}</table>"
    if post_msg:
        postElogMsg(exp=exp, msg=msg, tag=tag, title=title)

## function that chops the 64 bit time integer into something a bit more useful
def evtt2Rt(event_time):
    evtt0 = event_time>>32
    evtt1 = (event_time<<32)>>32
    evtt_sec = evtt0.astype(float)
    evtt_ns = evtt1.astype(float)*1e-9
    Rt = evtt_sec + evtt_ns
    Rt = Rt-Rt[0]
    return Rt

def postRunTable(runtable_data):
    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(args.experiment)
    print('URL:',ws_url)
    user=args.experiment[:3]+'opr'
    elogPostFile='/cds/home/opr/%s/forElogPost.txt'%user
    hostname=socket.gethostname()
    if hostname.find('sdf')>=0:
        elogPostFile='/sdf/group/lcls/ds/tools/forElogPost.txt'
    with open(elogPostFile,'r') as reader:
        answer = reader.readline()
    r = requests.post(ws_url, params={"run_num": args.run}, json=runtable_data, \
                      auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
    #we might need to use this for non=current expetiments. Currently does not work in ARP
    #krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    #r = requests.post(ws_url, headers=krbheaders, params={"run_num": args.run}, json=runtable_data)
    print(r)

def makeRunTableData(ana, ipmUpDim, ipmDownDim, Filter, scanName):
    n162 = ana.getVar('evr/code_162').sum()
    ana.addCut('evr/code_162',-0.5,0.5,'xon')
    ana.addCut('evr/code_137',0.5,1.5,'xon')
    nOff = ana.getFilter('xon').shape[0]-ana.getFilter('xon').sum()
    #data to be posted to the run table if so requested.
    runtable_data = {"N dropped Shots":int(nOff),
                     "N BYKIK 162":int(n162)}
    if scanName != '':
        runtable_data['scanName'] = scanName

    ipmUpVar = ana.getVar(ipmUpDim.name,useFilter=Filter)
    ipmDownVar = ana.getVar(ipmDownDim.name,useFilter=Filter)
    ipmUpP = np.nanpercentile(ipmUpVar,[25,50,75])
    ipmDownP = np.nanpercentile(ipmDownVar,[25,50,75])
    runtable_data["%s_1qt"%(ipmUpDim.name.replace('/','__'))]=ipmUpP[0]
    runtable_data["%s_med"%(ipmUpDim.name.replace('/','__'))]=ipmUpP[1]
    runtable_data["%s_3qt"%(ipmUpDim.name.replace('/','__'))]=ipmUpP[2]
    runtable_data["%s_1qt"%(ipmDownDim.name.replace('/','__'))]=ipmDownP[0]
    runtable_data["%s_med"%(ipmDownDim.name.replace('/','__'))]=ipmDownP[1]
    runtable_data["%s_3qt"%(ipmDownDim.name.replace('/','__'))]=ipmDownP[2]
    
    return runtable_data

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run', type=str, default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--stn', help='hutch station', type=int, default=0)
parser.add_argument('--nevents', help='number of events', type=int, default=1e9)
parser.add_argument('--directory', help='directory to read files from (def <exp>/hdf5/smalldata)', default=None)
parser.add_argument('--postElog', help='post plot to elog', action='store_true')
parser.add_argument('--postStats', help='post summary numbers to run tables', action='store_true')
#parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-kerb/lgbk/")
parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
args = parser.parse_args()
logger.debug('Args to be used for data quality plots: {0}'.format(args))


##############################################
## Setup Global parameters and run numbers ###
##############################################
save_elog = args.postElog
detImgMaxSize = 500 #max dimension of image.
expname = args.experiment
run = int(args.run)

if (int(os.environ.get('RUN_NUM', '-1')) > 0):
    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>BeamlineSummary Plots </b>", "value": "Started"}])

######################################
### load data for the chosen run  ####
######################################
#get the ana & anaps objects (from smalldata_tools
if args.directory is not None:
    anaps = sdaps(expname,run, dirname=args.directory)
else:
    anaps = sdaps(expname,run)

ana = anaps.sda #

## Defining initial selection (laser-on events)
iniFilter='initial'
ana.addCut('lightStatus/xray',0.5,1.5,iniFilter)
ana.addCut('lightStatus/laser',0.5,1.5,iniFilter)

### Get data & define axis title&ranges.

ipmUpDim = hv.Dimension(('ipm_dg1/sum','ipm_dg1 Sum'))
ipmDownDim = hv.Dimension(('ipm_dg2/sum','ipm_dg2 Sum'))

#xes1Dim = hv.Dimension(("ePix100_1/var_full_sparse"))

rayonixDim = hv.Dimension(('Rayonix/ROI_0_sum','Rayonix intensity'))
eventTimeDim = hv.Dimension(('eventTimeR','relative event time'))
#l3eDim = hv.Dimension(('l3e','L3 Energy'))

scanVar = ana.getScanName()
try:
    scanDim = hv.Dimension(('scan/%s'%scanVar,'%s'%scanVar))
except:
    scanDim = None
nevtsDim = hv.Dimension(('nevents','N events / scan point'))
nevtsLxtDim = hv.Dimension(('neventslxt','N events / lxt'))

#timing vars.
lxtDim = hv.Dimension(('epics/lxt','lxt'))

ipmUpVar = ana.getVar(ipmUpDim.name,useFilter=iniFilter)
ipmDownVar = ana.getVar(ipmDownDim.name,useFilter=iniFilter)
stepVar = ana.getVar('scan/varStep',useFilter=iniFilter)
#l3eVar = ana.getVar('ebeam/L3_energy',useFilter=iniFilter)
eventTimeRaw = ana.getVar('event_time',useFilter=iniFilter)
eventTime = (eventTimeRaw>>32).astype(float)+((eventTimeRaw<<32)>>32).astype(float)*1e-9
eventTimeR = eventTime-eventTime[0]

eventTimeRMed = [np.nanmedian(eventTimeR[i*120:i*120+120]) for i in range(int(eventTimeR.shape[0]/120))]
ipmUpMed =  [np.nanmedian(ipmUpVar[i*120:i*120+120]) for i in range(int(eventTimeR.shape[0]/120))]
ipmDownMed =  [np.nanmedian(ipmDownVar[i*120:i*120+120]) for i in range(int(eventTimeR.shape[0]/120))]

try:
    azav = ana.getVar('epix10k2M/azav_azav',useFilter=iniFilter)
    azav_sum = np.nanmean(azav, axis=0)
    azav_peak = np.argmax(azav_sum)
    if len(azav.shape)>2:
        azav = np.nanmean(azav, axis=1)
    scatterVar = np.nanmean(azav[:,max(0,azav_peak-50):min(azav.shape[1],azav_peak+50)], axis=1)
    if len(scatterVar.shape)>1:
        scatterVar = np.nanmean(scatterVar,axis=1)
except:
    scatterVar = None

### Scan Variable

try:
    isStepScan = np.nanmax(stepVar)>0
    scanVarBins = np.bincount(stepVar,weights=scatterVar)
    scanNsteps = np.bincount(stepVar)
except:
    isStepScan = False

### Fast delay stage

lxt_fast_his = None
try:
    lxt_fast = ana.getVar('enc/lasDelay',useFilter=iniFilter)
    print(np.nanstd(lxt_fast))
    if lxt_fast is not None and np.nanstd(lxt_fast)<1e-4:
        lxt_fast_his = np.histogram(lxt_fast, np.linspace(np.nanpercentile(lxt_fast,1), np.nanpercentile(lxt_fast,99),100))
except:
    pass

#droppled sthots.
ana.addCut('lightStatus/xray',-0.5,0.5,'off')
ana.addCut('evr/code_137',-0.5,0.5,'hxroff')
if ana.getFilter('hxroff').sum() >  ana.getFilter('off').sum():
    offFilter = 'hxroff'
else:
    offFilter = 'off'
nOff = ana.getFilter(offFilter).sum()

#########################################
# INSERT DATA QUALITY PLOTS HERE
#########################################
ipmUpTime = hv.HexTiles((eventTimeR[ipmUpVar<np.nanpercentile(ipmUpVar,99)],
                         ipmUpVar[ipmUpVar<np.nanpercentile(ipmUpVar,99)]),
                        kdims=[eventTimeDim, ipmUpDim]).\
                        opts(cmap='Blues')
ipmUpTimeMed = hv.Points((eventTimeRMed, ipmUpMed), kdims=[eventTimeDim,ipmUpDim],label=ipmUpDim.label).\
    options(color='r')
ipmDownTimeMed = hv.Points((eventTimeRMed, ipmDownMed), kdims=[eventTimeDim,ipmUpDim],label=ipmDownDim.label).\
    options(color='m')

ipmTimeLayout = ipmUpTime*ipmUpTimeMed*ipmDownTimeMed


ipmPlot = hv.HexTiles((ipmUpVar, ipmDownVar), kdims=[ipmUpDim, ipmDownDim])
ipmLayout = ipmPlot.hist(dimension=[ipmUpDim.name,ipmDownDim.name])

stepPlot = None

if lxt_fast_his is not None:
    lxtPlot = hv.Points( (0.5*(lxt_fast_his[1][:-1]+lxt_fast_his[1][1:]), lxt_fast_his[0]), \
                             kdims=[lxtDim,nevtsLxtDim])
else:
    lxtPlot = None

gspec = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name='Data Quality - Run %d'%run)
gspec[0:2,0:8] = pn.Column(ipmTimeLayout)
gspec[2:5,0:4] = pn.Column(ipmLayout)

detNames: list = ["epix_1", "Rayonix", "FEE_SPEC0", "EBeam"] # for tracking events missing data in one of these categories
try:
    plots = []
    from holoviews.operation.timeseries import rolling
    for detname in detNames:
        damageDim = hv.Dimension((f'damage/{detname}', f'{detname} Present'))
        damageVar = ana.getVar(damageDim.name, useFilter=iniFilter)
        damagePlot = rolling(
            hv.Curve(
                damageVar,
                vdims=[damageDim],
                label=f'{detname}'
            ).opts(axiswise=True, color=hv.Palette('Spectral')),
            rolling_window=10
        )
        plots.append(damagePlot)
    multiDamagePlot = hv.Overlay(plots).opts(
        xlabel='Event (rolling average of 10)',
        ylabel='Present (Yes/No)',
        title='Missing/Damaged Data'
    )
    gspec[2:5, 4:8] = multiDamagePlot

except Exception as e:
    pass

#########################################
# Detector Images
#########################################
detNames = ["epix_1", "Rayonix"]

detImgs=[]
detGrids=[]
for detImgName in ana.Keys('Sums'):
    image = ana.fh5.get_node('/%s'%detImgName).read()
    if len(image.shape)>2:
        if detImgName.find('135')<0:
            detName = detImgName.replace('Sums/','').replace('_calib','')
            ix = ana.fh5.get_node('/UserDataCfg/%s/ix'%detName).read()
            iy = ana.fh5.get_node('/UserDataCfg/%s/iy'%detName).read()
            image = image_from_dxy(image, ix, iy)
        else:
            #somehow the epix10k135 has the wrong shape....
            image = image[0]
            #image = image.squeeze()
    if max(image.shape[0], image.shape[1])>detImgMaxSize:
        rebinFactor = float(detImgMaxSize)/max(image.shape[0],image.shape[1])
        imageR = rebin(image, [int(image.shape[0]*rebinFactor), int(image.shape[1]*rebinFactor)])/(ana.getVar('fiducials').shape[0])
    else:
        imageR = image/(ana.getVar('fiducials').shape[0])
    #imgArrays.append(imageR/ana.getVar('fiducials').shape[0])
    imgDim = hv.Dimension(('image',detImgName.replace('Sums/','').replace('_calib_img',' Mean Image')),
                                    range=(np.nanpercentile(imageR,1), np.nanpercentile(imageR,99.)))
    detImgs.append(hv.Image(imageR, vdims=[imgDim], name=imgDim.label).options(colorbar=True, cmap='rainbow'))
        
    detGrid = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name=detImgName.replace('Sums/',''))
    detGrid[0,0] = pn.Row(detImgs[-1])
    detGrids.append(detGrid)

if nOff>100:
    for detImgName in ana.Keys('Sums'):
        detName = detImgName.replace('_calib','').replace('_img','').replace('Sums/','')
        try:
            common_mode=0
            if detName.find('epix10k'): common_mode=80
            anaps.AvImage(detName, useFilter=offFilter, numEvts=min(1000, nOff), common_mode=common_mode)
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

########################
# XES Droplet Plots
#######################

def processDropletData(
        ana, #: smalldata_tools.SmallDataAna.SmallDataAna,
        det_name: str,
        adr_thr: Tuple[int,int] = (2, 15),
        rot_angle: Optional[float] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract image and spectra from droplet data for a detector.
    
    Constructs an image from the 2D histogram of droplet center-of-masses.
    The projections onto the spatial and spectral axes are then computed.
    A number of transformations can be applied to the image prior to computing
    the projections.
    
    Parameters
    ----------
    ana (smalldata_tools.SmallDataAna.SmallDataAna) ana object for pulling data
        from smalldata file.
    det_name (str) The detector to get droplet data for.
    adu_thr (tuple) 2-tuple of lower and upper threshold to apply to the data.
    rot_angle (float) Rotation to apply to the image prior to projections. In
        degrees. By default, unused.
    """
    col = ana.getVar(f"{det}/var_full_sparse/col")
    row = ana.getVar(f"{det}/var_full_sparse/row")
    adu = ana.getVar(f"{det}/var_full_sparse/data")
    n_droplets = len(adu)

    # Additional thresholding.
    adu_thr = [2, 15]
    indices = np.where((adu >= adu_thr[0]) & (adu <= adu_thr[1]))

    dimg, xedges, yedges = np.histogram2d(
        col[indices],
        row[indices],
        bins=[np.arange(800),np.arange(769)],
        weights=adu[indices]
    )

    if rot_angle:
        from scipy.ndimage import rotate
        dimg = rotate(dimg, angle=rot_angle)

    # Calculate projections
    spatial_proj = np.sum(dimg, axis=0)
    xes_proj = np.sum(dimg, axis=1)

    return dimg, spatial_proj, xes_proj


xesPlots = []
for det in ["epix_1", "epix_2"]:# detNames:
    try:
        dimg, spatial_proj, xes_proj = processDropletData(ana=ana, det_name=det)
        #dimg_flip = np.fliplr(dimg)

        xes_grid = pn.GridSpec(max_width=700, name=f"XES - {det}")
        img_dim = hv.Dimension(("Image", "Image"))
        spatial_dim = hv.Dimension(("Spatial Projection", "Spatial Projection"))
        xes_dim = hv.Dimension(("XES Spectrum", "XES Spectrum"))
        xes_grid[:2, :] = hv.Image(
            dimg,
            bounds=(0, 0, dimg.shape[0], dimg.shape[1]),
            vdims=[img_dim]).options(colorbar=True,
            clim=(np.nanpercentile(dimg,1),np.nanpercentile(dimg,99.))
        ).opts(xlabel="Spatial Axis", ylabel="Spectral Axis")
        spatial_plot = hv.Curve(spatial_proj, vdims=[spatial_dim])
        xes_plot = hv.Curve(xes_proj, vdims=[xes_dim])
        xes_grid[2, :] = spatial_plot.opts(
            axiswise=True,
            xlabel="Pixel",
            ylabel="Sum",
            title="Spatial Projection"
        )
        xes_grid[3, :] = xes_plot.opts(
            axiswise=True,
            xlabel="Energy (Pixel)",
            ylabel="I",
            title="XES Spectrum"
        )

        xesPlots.append(xes_grid)
    except Exception as e:
        print(e)

###########
# FEE Plots
###########
feeGrid = pn.GridSpec(max_width=700, name="FEE Stats")
try:
    feeDamageDim = hv.Dimension(('damage/FEE_SPEC0','FEE Present'))
    feeDamageVar = ana.getVar(feeDamageDim.name,useFilter=iniFilter)
    feeDamagePlot = hv.Curve(feeDamageVar, vdims=[feeDamageDim]).opts(
            axiswise=True,
            xlabel="Event",
            ylabel="FEE Present",
            title="Damage/FEE"
        )
    feeGrid[:2, :2] = feeDamagePlot
except Exception as e:
    pass

try:
    feeSpecDim = hv.Dimension(('feeBld/hproj','FEE Spec (mean)'))
    feeSpecVar = ana.getVar(feeSpecDim.name,useFilter=iniFilter)
    feeSpecPlot = hv.Curve(
        np.mean(feeSpecVar, axis=0),
        vdims=[feeSpecDim]
    ).opts(
            axiswise=True,
            xlabel="Energy (Pixel)",
            ylabel="I",
            title="FEE Spectrum"
        )
    feeGrid[:2, 2:4] = feeSpecPlot

    x = np.arange(feeSpecVar[0].shape[0])
    feeCOMDim = hv.Dimension(('ProcessedFee', 'FEE COM'))
    feeCOMVar = np.sum(feeSpecVar*x, axis=1)/np.sum(feeSpecVar, axis=1)
    eBeamDim = hv.Dimension(('ebeam/photon_energy', 'eBeam Energy'))
    eBeamVar = ana.getVar(eBeamDim.name,useFilter=iniFilter)
    energyCorrPlot = hv.HexTiles(
        (eBeamVar, feeCOMVar),
        kdims=[eBeamDim,feeCOMDim]
    ).opts(
            axiswise=True,
            xlabel="eBeam Energy",
            ylabel="FEE COM",
            title="FEE vs eBeam"
        )
    feeGrid[2:4, :2] = energyCorrPlot
except Exception as e:
    pass

########################
# Tabs construction and finish
########################

tabs = pn.Tabs(gspec)

for xes_grid in xesPlots:
    tabs.append(xes_grid)

for detGrid in detGrids:
    tabs.append(detGrid)

if (int(os.environ.get('RUN_NUM', '-1')) > 0):
    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>BeamlineSummary Plots </b>", "value": "Done"}])

elogDir = f"/sdf/data/lcls/ds/{expname[:3]}/{expname}/stats/summary/BeamlineSummary/BeamlineSummary_Run{run:04d}"

if save_elog:
    import os
    if not os.path.isdir(elogDir):
        os.makedirs(elogDir)
    print('Made Directory to save data:', elogDir)
    tabs.save(('%s/report.html'%elogDir))
    postDetectorDamageMsg(
        detectors=detNames,
        exp=expname,
        run=run,
        title=f"EVENT DAMAGE INFO- r{run:04d}",
        smd_dir=args.directory,
        post_thresh=0.1 # Percentage threshold to post to eLog
    )

    if (int(os.environ.get('RUN_NUM', '-1')) > 0):
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>BeamlineSummary Plots </b>", "value": "Posted"}])

if args.postStats:
    if scanVar == '':
        encDelay = ana.getVar('enc/lasDelay')
        delta_encDelay = np.nanmax(encDelay)-np.nanmin(encDelay)
        if delta_encDelay > 0.5:
            scanVar='delay'
    elif scanVar.find('lxt'):
        scanVar='delay'
    runtable_data = makeRunTableData(ana, ipmUpDim, ipmDownDim, iniFilter, scanVar)
    postRunTable(runtable_data)
