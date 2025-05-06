#!/usr/bin/env python
########################################################################
## load tools and packages required for data handling and plotting
########################################################################

import tables
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
from holoviews import opts

hv.extension("bokeh")
pn.extension()
import sys
#sys.path.append('/sdf/home/d/dgarratt/chemRIXS/Pyrazine_Dev/Modules')
#from chemRIXSAnalysis_V3 import *
fpath = os.path.dirname(os.path.abspath(__file__))
fpathup = "/".join(fpath.split("/")[:-1])
try:
    fpath = os.environ.get("MYDIR", fpathup).replace("/arp_scripts", "")
except:
    fpath = fpathup
sys.path.append(fpath)
from smalldata_tools.utilities import rebin
from smalldata_tools.utilities import evtt2Rt
from smalldata_tools.utilities import postRunTable


def click_policy(plot, element):
    plot.state.legend.click_policy = "hide"

def makeRunTableData(dat):
    xray = dat.lightStatus.xray.read()
    laser = dat.lightStatus.laser.read()
    on = xray.astype(bool) & laser.astype(bool)
    off = xray.astype(bool) & ~laser.astype(bool)
    dropped = ~xray.astype(bool)
    total = xray.shape[0]

    runtable_data = {
        "N dropped Shots": int(dropped.sum()),
        "N total": int(total),
        "N laser on": int(on.sum()),
        "N laser off": int(off.sum()),
    }

    return runtable_data

def process_ANDOR_Basic(andor, bg_roi_x = [0,100],bg_roi_y=[0,100],\
                        absolute_theshold = False,threshold = False,photon_count = False,\
                        threshold_min=4, threshold_max=4000):
    # 1D processing
    if len(andor.shape) ==2:
        andor_nbg = subtract_bg(andor,bg_roi_x)
        andor_proc = andor_nbg.copy()
        if threshold:
            astd,amean = andor_nbg[:,bg_roi_x[0]:bg_roi_x[1]].std(),andor_nbg[:,bg_roi_x[0]:bg_roi_x[1]].mean()
            andor_proc[andor_proc>threshold_max] = 0
            if absolute_threshold:
                andor_proc[andor_proc<(threshold_min)] = 0
            else:
                andor_proc[andor_proc<(threshold_min*astd)] = 0
            if photon_count:
                andor_proc[andor_proc>threshold_min] = 1
    elif len(andor.shape)==3:
        # Shape of full image is (2048 x 512)
        bgs = andor[:,bg_roi_y[0]:bg_roi_y[1],bg_roi_x[0]:bg_roi_x[1]].mean(axis = (1,2))
        andor_nbg = andor - bgs[:,np.newaxis,np.newaxis]
        andor_proc = andor_nbg.copy()
        if threshold:
            astd = andor_nbg[:,bg_roi_y[0]:bg_roi_y[1],bg_roi_x[0]:bg_roi_x[1]].std()
            amean = andor_nbg[:,bg_roi_y[0]:bg_roi_y[1],bg_roi_x[0]:bg_roi_x[1]].mean()
            andor_proc[andor_proc>threshold_max] = 0
            if absolute_threshold:
                andor_proc[andor_proc<(threshold_min)] = 0
            else:
                andor_proc[andor_proc<(threshold_min*astd)] = 0
            if photon_count:
                andor_proc[andor_proc>threshold_min] = 1
    return andor_proc
def sum_ANDOR(andor_ims):
    if len(andor_ims.shape) ==2:
        return np.sum(andor_ims,axis = 1)
    elif len(andor_ims.shape) ==3:
        return np.sum(andor_ims,axis = (1,2))


# Data loading
def load_int_data(run,fh,prefix = '', load_ttfex = True, evc_laser = 272,int_detector_key = 'intg/andor_vls',roll = True):
    keys_andor = [('intg/andor_dir',f'{prefix}full_area','dir'),
            ('intg/andor_vls',f'{prefix}full_area','vls'),
            ('intg/andor_vls',f'{prefix}hitfinder','vls_hitfound'),
            ('intg/andor_norm',f'{prefix}full_area','norm'),
           ]
    keys_int = [(int_detector_key ,f'{prefix}det_rix_fim1_sum_full_area','fim1_wf'),
            (int_detector_key ,f'{prefix}det_rix_fim0_sum_full_area','fim0_wf'),
            (int_detector_key ,f'{prefix}det_crix_w8_sum_full_area','apd_wf'),
           # (int_detector_key ,f'{prefix}hsd_sum_full_hsd_1__ROI_wf','hsd1_wf'),
            (int_detector_key ,f'{prefix}det_rix_fim1_sum_wfintegrate','fim1_sum'),
            (int_detector_key ,f'{prefix}det_rix_fim0_sum_wfintegrate','fim0_sum'),
            (int_detector_key ,f'{prefix}det_crix_w8_sum_wfintegrate','apd_sum'),
            (int_detector_key ,f'{prefix}mono_encoder_sum_interpolated_value','mono_encoder'),
            (int_detector_key ,f'{prefix}mono_encoder_sum_raw_value','mono_encoder_lowrate'),
            (int_detector_key ,f'{prefix}mono_hrencoder_sum_value','mono_hrencoder'),
            (int_detector_key ,f'{prefix}c_piranha_sum_full_area','pir'),
            (int_detector_key ,f'{prefix}timing_sum_eventcodes','evc'),
            (int_detector_key ,f'{prefix}timing_sum_destination','dest')
               ]
    # These should not be normalized to count
    keys_ts = [(int_detector_key ,f'{prefix}timestamp_min','timestamp_min'),
            (int_detector_key ,f'{prefix}timestamp_max','timestamp_max'),
            (int_detector_key ,f'{prefix}timing_sum_timestamp','timestamp_norm')
               ]
    if load_ttfex:
        keys_int = keys_int + [(int_detector_key ,f'{prefix}tt_sum_fltpos','pos'),
            (int_detector_key ,f'{prefix}tt_sum_fltposfwhm','fwhm'),
            (int_detector_key ,f'{prefix}tt_sum_ampl','amp')
           ]
    
    data = {}
    # Load timestamps 
    data['timestamp'] = np.array(fh['intg/timestamp'])
    # summed timestamps
    for key in keys_ts:
        #data[key[2]] = np.array(fh[key[0]][key[1]])
        try:
            data[key[2]] = np.array(fh[key[0]][key[1]]) 
        except:
            print('missing key {0}'.format(key))
    
    # Load andor images
    for key in keys_andor:
        try:
            data[key[2]] = np.array(fh[key[0]][key[1]]) 
        except:
                print('missing key {0}'.format(key))

    # Normalize these parameters to the number of single shots per image
    data['count'] = np.array(fh[int_detector_key][f'{prefix}count'])
    for key in keys_int:
        #data[key[2]] = np.array(fh[key[0]][key[1]], chunks=chunks) 
        try:
            data[key[2]] = np.array(fh[key[0]][key[1]])  
            # print('Key: {0} shape: %s'.format(key) %(data[key[2]].shape))
            if len(data[key[2]].shape)==1:
                data[key[2]] = data[key[2]]/data['count']
    
            elif len(data[key[2]].shape)==2:
                data[key[2]] = data[key[2]]/data['count'][:,np.newaxis]
            elif len(data[key[2]].shape)==3:
                data[key[2]] = data[key[2]]/data['count'][:,np.newaxis,np.newaxis]
        except:
            print('missing key {0}'.format(key))

    try:
        scan_key = list(fh['scan'].keys())[0]
        data['scan_key'] = scan_key
        data['step'] = np.squeeze(np.array(fh[int_detector_key][f'{prefix}scan_sum_step_value']))/data['count']
        data['x'] = np.squeeze(np.array(fh[int_detector_key][f'{prefix}scan_sum_{scan_key}']))/data['count']
    except:
        print('Error loading integrated scan variables')

    # Generate a xray on and laser on 
    data['laser_on'] = data['evc'][:,evc_laser]>0.5
    data['xray_on'] = data['dest']==4
    data['dest'] = data['dest']

    idx = np.argsort(data['timestamp'])
    # Sort data by timestamp
    for k in data.keys():
        if k != 'scan_key':
            data[k] = data[k][idx]
    # Roll the andor images to account for andor trigger delay
    if roll:
        for key in keys_andor:
            try:
                data[key[2]] = np.roll(data[key[2]],1,axis = 0)
            except:
                print('missing key {0}'.format(key))
        # cropp all data 
        for key, value in data.items():
            #print(value.shape)
            data[key] = value[1:]
            #print(data[key].shape)
            

    # Load epics_archiver variables
    epics_data = load_epics(fh)
    for k in epics_data.keys():
        data[f'epics_{k}']=epics_data[k]
    
    return data

def load_epics(fh):
    '''Loads the epics archiver variables from h5 file.
    Note the shape varies
    '''
    epics_keys = [('epics_archiver', 'CRIX_VLS_CAM_MMS_PITCH.RBV', 'VLS_camera_pitch'),
                  ('epics_archiver', 'CRIX_VLS_MMS_GP.RBV', 'VLS_grating_pitch_mm'),
                  ('epics_archiver', 'CRIX_VLS_MMS_MP.RBV', 'VLS_mirror_pitch_mm'),
                  ('epics_archiver', 'LM2K2_COM_MP2_DLY1.RBV', 'ATM_stage_mm'),
                  ('epics_archiver', 'SP1K1_MONO_MMS_G_PI.RBV', 'MONO_grating_pitch'),
                  ('epics_archiver', 'SP1K1_MONO_MMS_M_PI.RBV', 'MONO_premirror_pitch')
           ]
    dat = {}
    for i, k in enumerate(epics_keys):
        try:
            dat[k[2]] = fh[k[0]][k[1]][...][:,1]
        except: print(f'Missing epics key {k[1]}')
    return dat

def load_smd_ROIS(fh):
    roi_dict = {}
    keys = [
        ('det_crix_w8','wfintegrate__wfintegrate_bkg_roi','apd_bg'),
        ('det_crix_w8','wfintegrate__wfintegrate_sig_roi','apd_sig'),
        ('det_rix_fim0','wfintegrate__wfintegrate_bkg_roi','fim0_bg'),
        ('det_rix_fim0','wfintegrate__wfintegrate_sig_roi','fim0_sig'),
        ('det_rix_fim1','wfintegrate__wfintegrate_bkg_roi','fim1_bg'),
        ('det_rix_fim1','wfintegrate__wfintegrate_sig_roi','fim1_sig')
    ]
    for k in keys:
        try:
            roi_dict[k[2]] = fh['UserDataCfg'][k[0]][k[1]][...]
        except:
            print(f'Missing key {k}')
    return roi_dict

def subtract_bg(wf,bg_roi=[0,50],invert = False):
    bg = np.mean(wf[...,bg_roi[0]:bg_roi[1]],axis= -1)
    shape = wf.shape
    if invert:
        wf_nbg = -(wf - bg[...,np.newaxis])
    else:
        wf_nbg = wf - bg[...,np.newaxis]
    return wf_nbg


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument(
    "--run", help="run", type=str, default=os.environ.get("RUN_NUM", "")
)
parser.add_argument(
    "--experiment",
    help="experiment name",
    type=str,
    default=os.environ.get("EXPERIMENT", ""),
)
parser.add_argument("--stn", help="hutch station", type=int, default=0)
parser.add_argument("--nevents", help="number of events", type=int, default=1e9)
parser.add_argument(
    "--directory",
    help="directory to read files from (def <exp>/hdf5/smalldata)",
    default=None,
)
parser.add_argument(
    "--postElog", help="post plot to elog", action="store_true", default=True
)
parser.add_argument(
    "--postStats",
    help="post summary numbers to run tables",
    action="store_true",
    default=False,
)
# parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-kerb/lgbk/")
parser.add_argument("--url", default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
args = parser.parse_args()
logger.debug("Args to be used for data quality plots: {0}".format(args))


##############################################
## Setup Global parameters and run numbers ###
##############################################
save_elog = args.postElog
expname = args.experiment
run = int(args.run)

if int(os.environ.get("RUN_NUM", "-1")) > 0:
    requests.post(
        os.environ["JID_UPDATE_COUNTERS"],
        json=[{"key": "<b>BeamlineSummary Plots </b>", "value": "Started"}],
    )

######################################
### load data for the chosen run  ####
######################################
dir_roi_x = [0,-1]
dir_roi_y = [0,-1]
dir_bg_roi_x = [50,100]
dir_bg_roi_y = [50,100]

# x is 2048 pixel axis 
vls_roi_x = [0,-1]
vls_roi_y = [0,-1]
vls_bg_roi_x = [50,100]
vls_bg_roi_y = [50,100]
vls_threshold = 10
absolute_threshold = True
vls_threshold_max = 500

fim_bg_roi = np.array([0,80])
fim0_roi = [100,113]
fim1_roi = [122,128]
neig_fim = 2
idx_fim1 = np.array([4,5,6,7])
idx_fim0 = np.array([4,5,6,7])

apd_bg_roi = np.array([0,30])
apd_roi = [72,95]
neig_apd = 2
idx_apd = np.array([5,6,7])

hsd_roi = [1700,2100]
hsd_bg_roi = [0,1000]

key_norm = 'i0_fim'
# Shot Filtering
i0_thresh = 500
# Time tool parameters
px0 = 150

# Shot Filtering
i0_thresh = 500

# Time tool parameters
px0 = 150
b = 50

fname = f'/sdf/data/lcls/ds/{expname[:3]}/{expname}/hdf5/smalldata/{expname}_Run{run:04d}.h5'
fh = h5py.File(fname, 'r')
dat = load_int_data(run, fh,int_detector_key='intg/andor_vls',roll=True)
rois = load_smd_ROIS(fh)
###################################################################################################################################################
###################################################################################################################################################
print('Processing detectors')
# Process VLS 
if 'vls' in dat.keys():
    dat['vls_proc'] = process_ANDOR_Basic(dat['vls'],[vls_bg_roi_x[0],vls_bg_roi_x[1]],[vls_bg_roi_y[0],vls_bg_roi_y[1]],\
                                          threshold_min = vls_threshold,threshold_max=vls_threshold_max)
    dat['vls_nbg'] = process_ANDOR_Basic(dat['vls'],[vls_bg_roi_x[0],vls_bg_roi_x[1]],[vls_bg_roi_y[0],vls_bg_roi_y[1]],\
                                         threshold = False)
    if len(dat['vls'].shape)==3:
        dat['vls_lineout'] = np.sum(dat['vls_proc'],axis = -1)
    else:
        dat['vls_lineout'] = dat['vls_proc']

    dat['pfy'] = sum_ANDOR(dat['vls_lineout'][:,vls_roi_x[0]:vls_roi_x[1]])
    vls_av = np.nanmean(dat['vls_proc'],axis = 0)
    
if 'dir' in dat.keys():
    # Process direct andor
    dat['dir_proc'] = process_ANDOR_Basic(dat['dir'],[dir_bg_roi_x[0],dir_bg_roi_x[1]],[dir_bg_roi_y[0],dir_bg_roi_y[1]])
    dat['iT'] = sum_ANDOR(dat['dir_proc'])
    dir_av = np.nanmean(dat['dir_proc'],axis = 0)
    
# Subtract background and invert FIM traces
dat['fim0_wf_proc'] = subtract_bg(dat['fim0_wf'],fim_bg_roi,invert = True)
dat['fim1_wf_proc'] = subtract_bg(dat['fim1_wf'],fim_bg_roi,invert = True)
fims_concat = np.concatenate((dat['fim0_wf_proc'][:,4:,:],dat['fim1_wf_proc'][:,4:,:]),axis = 1)
# Combine fim detectors
# Subtract background and invert FIM traces
dat['i0_fim0'] = np.sum(dat['fim0_wf_proc'][:,idx_fim0,fim0_roi[0]:fim0_roi[1]],axis = (1,2))
dat['i0_fim1'] = np.sum(dat['fim1_wf_proc'][:,idx_fim1,fim1_roi[0]:fim1_roi[1]],axis = (1,2))
dat['i0_fim'] = dat['i0_fim0'] + dat['i0_fim1']
# APD processing
dat['apd_wf_proc'] = subtract_bg(dat['apd_wf'],apd_bg_roi,invert = True)
dat['iF_sum'] = np.sum(dat['apd_wf_proc'][:,idx_apd,apd_roi[0]:apd_roi[1]],axis = (1,2))
#dat['iF_SVD_all'] = process_FIM_SVD(dat['apd_wf_proc'][:,idx_apd,:],neig_apd)
#dat['iF_SVD'] = np.nanmean(dat['iF_SVD_all'],axis = 1)
#
if 'hsd1_wf' in dat.keys():
    dat['hsd_wf_proc'] = subtract_bg(dat['hsd1_wf'],hsd_bg_roi,invert = True)
    dat['iF_sum_hsd'] =   np.sum(dat['hsd_wf_proc'][:,hsd_roi[0]:hsd_roi[1]],axis = (1))
#    dat['iF_SVD_hsd'] =  np.nanmean(process_FIM_SVD(dat['hsd_wf_proc'][:,np.newaxis,hsd_roi[0]:hsd_roi[1]],neig_apd),axis = 1)
    hsd_av = np.nanmean(dat['hsd_wf_proc'],axis = 0)

fim0_av = np.nanmean(dat['fim0_wf_proc'],axis = 0) 
fim1_av = np.nanmean(dat['fim1_wf_proc'],axis = 0) 
apd_av = np.nanmean(dat['apd_wf_proc'],axis = 0)

count, ncount = np.unique(dat['count'],return_counts = True)
expected_count = count[np.argmax(ncount)]
print(count[ncount>2],ncount[ncount>2],f'Expected count {expected_count}')
count_msk = dat['count']==expected_count
# Single shot normalization 
if 'vls' in dat.keys():
    dat['vls_lineout_norm'] = dat['vls_lineout']/dat[key_norm][:,np.newaxis]
    dat['pfy_norm'] = dat['pfy']/dat[key_norm][:]
if 'dir' in dat.keys():
    dat['iT_norm'] = dat['iT']/dat[key_norm][:]
dat['iF_sum_norm'] = dat['iF_sum']/dat[key_norm][:]
#dat['iF_SVD_norm'] = dat['iF_SVD']/dat[key_norm][:]
if 'hsd1_wf' in dat.keys():
    dat['iF_sum_hsd_norm'] = dat['iF_sum_hsd']/dat[key_norm][:]
#    dat['iF_SVD_hsd_norm'] = dat['iF_SVD_hsd']/dat[key_norm][:]
###################################################################################################################################################
def click_policy(plot, element):
    plot.state.legend.click_policy = "hide"
######################################
### load data for the chosen run  ####
######################################
fname = f'/sdf/data/lcls/ds/{expname[:3]}/{expname}/hdf5/smalldata/{expname}_Run{run:04d}.h5'

laser = dat['laser_on']
on = laser.astype(bool)
off = ~laser.astype(bool)


fim0Dim = hv.Dimension(("fim0", "fim0"))
fim1Dim = hv.Dimension(("fim1", "fim1"))
apdDim = hv.Dimension(("apd", "apd"))
fimDim = hv.Dimension(("i0", "i0"))
apdPx =  hv.Dimension(("w8", "px"))
FIMPx =  hv.Dimension(("w8", "px"))
laserDim = hv.Dimension(("laser_on", "laser_on"))
xrayDim = hv.Dimension(("Destination", "Average beam dest"))
pixelDim = hv.Dimension(("pixel", "pixel"))

dirDim = hv.Dimension(("DIR", "andor dir"))
vlsDim = hv.Dimension(("VLS", "andor VLS"))
pfyDim = hv.Dimension(("DIR sum", "pfy"))
iTDim = hv.Dimension(("VLS sum", "iT"))
iFDim = hv.Dimension(("APD sum", "iF"))
eventTimeDim = hv.Dimension(("eventTimeR", "relative event time"))
ttposDim = hv.Dimension(("ttpos", "tt fitted position"), range=(0, 1000))
ttamplDim = hv.Dimension(("ttampl", "tt amplitude"))
# fim0All = -1.*dat.det_rix_fim0.full_fimSum_fimSum.read().sum(axis=1)

i0Var = dat['i0_fim']
i0Dim = fimDim

eventTimeRaw = dat['timestamp']
eventTime = (eventTimeRaw >> 32).astype(float) + ((eventTimeRaw << 32) >> 32).astype(
    float
) * 1e-9
eventTimeR = eventTime - eventTime[0]

# X-ray parameters
i0Time = hv.Points((eventTimeR, i0Var), kdims=[eventTimeDim, fimDim], label=fimDim.label).options(color="r")
laserTime = hv.Curve((eventTimeR, laser), kdims=[eventTimeDim, laserDim], label=laserDim.label).options(color="r")
xrayTime = hv.Curve((eventTimeR, dat['dest']), kdims=[eventTimeDim, xrayDim], label=xrayDim.label).options(color="r")

gspecS = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="X-ray Summary "
)
gspecS[0:2, 0:8] = pn.Column(i0Time)
gspecS[2:4, 0:8] = pn.Column(laserTime)
gspecS[4:6, 0:8] = pn.Column(xrayTime)
gspecS
###################################################################################################################################################
pixelDim = hv.Dimension(("pixel", "pixel"))

# DIR Processing
if 'dir' in dat.keys():
    dir_proc_mean = np.mean(dat['dir_proc'], axis=0)
    vl_dir = hv.VLines(dir_roi_x)
    if dir_proc_mean.ndim == 2:
        dir_heatmap = hv.Image(dir_proc_mean).opts(
            cmap='viridis', colorbar=True, width=400, height=400, title='DIR Full Mean'
        )
        dir_lineout = (hv.Curve(np.nanmean(dat['dir_proc'], axis=(0, 2)),pixelDim,dirDim)*vl_dir).opts(
            opts.Curve(xlabel='Pixel', ylabel='Value', title='VLS Lineout', width=400, height=400),
            opts.VLines(color = 'red'))
    else:
        dir_lineout = (hv.Curve(dir_proc_mean,pixelDim,dirDim)*vl_dir).opts(
            opts.Curve(xlabel='Pixel', ylabel='Value', title='DIR Lineout', width=400, height=400),
            opts.VLines(color = 'red'))
        dir_heatmap = None
        
    dirCorr = hv.HexTiles((i0Var,dat['iT']), kdims=[i0Dim,iTDim]).opts(title = 'I0 DIR correlation')
    iTTime = hv.Points((eventTimeR, dat['iT']), kdims=[eventTimeDim, iTDim], label='DIR Sum vs time').options(color="r")

    gspecDIR = pn.GridSpec(sizing_mode="stretch_both", max_width=700, name="Transmission Summary ")
    gspecDIR[0:2, 0:2] = pn.Column(dir_lineout)
    gspecDIR[0:2, 2:4] = pn.Column(dirCorr)
    gspecDIR[2:4, 0:2] = pn.Column(iTTime)
    if dir_heatmap is not None:
        gspecDIR[2:4, 2:4] = pn.Column(dir_heatmap)
else:
    gspecDIR = None
    
###################################################################################################################################################
# VLS Processing
if 'vls' in dat.keys():
    vls_proc_mean = np.mean(dat['vls_proc'], axis=0)
    vl_vls = hv.VLines(vls_roi_x)
    if vls_proc_mean.ndim == 2:
        vls_heatmap = hv.Image(vls_proc_mean).opts(
            cmap='viridis', colorbar=True, width=400, height=400, title='VLS Full Mean'
        )
        vls_lineout = (hv.Curve((np.nanmean(dat['vls_proc'], axis=(0, 2))),pixelDim,vlsDim)*vl_vls).opts(
            opts.Curve(xlabel='Pixel', ylabel='Value', title='VLS Lineout', width=400, height=400),
            opts.VLines(color = 'red'))
        
    else:
        vls_lineout = (hv.Curve(vls_proc_mean,pixelDim,vlsDim)*vl_vls).opts(
            opts.Curve(xlabel='Pixel', ylabel='Value', title='VLS Lineout', width=400, height=400),
            opts.VLines(color = 'red'))
        
        vls_heatmap = None

    vlsCorr = hv.HexTiles((i0Var,dat['pfy']), kdims=[i0Dim,vlsDim]).opts(title = 'I0 VLS correlation')
    
    vls_count_bins = np.arange(0,100)
    # Histograms
    vls_hist_data, vls_hist_edges = np.histogram(dat['vls_proc'].flatten(), bins=vls_count_bins)
    vls_nbg_hist_data, vls_nbg_hist_edges = np.histogram(dat['vls_nbg'].flatten(), bins=vls_count_bins)
    vls_hist_centers = (vls_hist_edges[:-1] + vls_hist_edges[1:]) / 2
    vls_nbg_hist_centers = (vls_nbg_hist_edges[:-1] + vls_nbg_hist_edges[1:]) / 2
    vlsCountDim = hv.Dimension(("VLS Counts", "VLS Counts"))
    countDim = hv.Dimension(("Counts", "Counts"))
    vls_hist1 = hv.Curve((vls_hist_centers, vls_hist_data), kdims=[vlsCountDim,countDim]).opts(color='blue', title="Processed Histogram",alpha = 0.5,logy=True)
    vls_hist2 = hv.Curve((vls_nbg_hist_centers, vls_nbg_hist_data), kdims=[vlsCountDim,countDim]).opts(color='blue', title="Raw Histogram",alpha = 0.5,logy=True)
    vls_hist = (vls_hist1 * vls_hist2).opts(title="VLS Histograms", legend_position='left')

    pfyTime = hv.Points((eventTimeR, dat['pfy']), kdims=[eventTimeDim, pfyDim], label='VLS sum Vs Time').options(color="r")

    gspecVLS = pn.GridSpec(sizing_mode="stretch_both", max_width=700, name="VLS Summary ")
    
    maxRow = 0
    gspecVLS[0:2, 0:2] = pn.Column(vls_lineout)
    gspecVLS[0:2, 2:4] = pn.Column(vlsCorr)
    gspecVLS[2:4, 0:2] = pn.Column(vls_hist)
    gspecVLS[2:4, 2:4] = pn.Column(pfyTime)
    if vls_heatmap is not None:
        gspecVLS[4:6, 0:4] = pn.Column(vls_heatmap)
else:
    gspecVLS = None


###################################################################################################################################################

apd_curves = [hv.Curve(wf, label=f'Channel {idx_apd[i]}') for i, wf in enumerate(apd_av[idx_apd,:])]
apd_plot = hv.Overlay(apd_curves,kdims=[apdPx, apdDim]).opts(title="APD", width=400, height=400)


apdCorr = hv.HexTiles((i0Var,dat['iF_sum']), kdims=[i0Dim,iFDim]).opts(title = 'I0 APD correlation')
iFTime = hv.Points((eventTimeR, dat['iF_sum']), kdims=[eventTimeDim, iFDim], label=apdDim.label).options(color="r")

gspecAPD = pn.GridSpec(sizing_mode="stretch_both", max_width=700, name="TFY Summary ")
gspecAPD[0:2, 0:4] = pn.Column(apd_plot)
gspecAPD[2:4, 0:2] = pn.Column(apdCorr)
gspecAPD[2:4, 2:4] = pn.Column(iFTime)

###################################################################################################################################################
# FIM0 Curves
fim0_curves = [hv.Curve(wf, FIMPx,fim0Dim, label=f'Channel {i}') for i, wf in enumerate(fim0_av[:])]

gspecFim0 = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="Summed FIM traces - Run %d" % run
)

for k, plot in enumerate(fim0_curves):
    irow = int(int(k) / 3)
    icol = int(k) % 3
    print(irow,icol)
    gspecFim0[irow * 3 : (irow + 1) * 3, icol * 3 : (icol + 1) * 3] = pn.Column(
        plot
    )
gspecFim0

# FIM1 Curves
fim1_curves = [hv.Curve(wf, label=f'Channel {i}') for i, wf in enumerate(fim1_av[:])]

gspecFim1 = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="Summed FIM traces - Run %d" % run
)

for k, plot in enumerate(fim1_curves):
    irow = int(int(k) / 3)
    icol = int(k) % 3
    print(irow,icol)
    gspecFim1[irow * 3 : (irow + 1) * 3, icol * 3 : (icol + 1) * 3] = pn.Column(
        plot
    )

# APD Curves
apd_curves = [hv.Curve(wf,apdPx,apdDim, label=f'Channel {i}') for i, wf in enumerate(apd_av[:])]

gspecW8 = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="Summed APD traces - Run %d" % run
)

for k, plot in enumerate(apd_curves):
    irow = int(int(k) / 3)
    icol = int(k) % 3
    print(irow,icol)
    gspecW8[irow * 3 : (irow + 1) * 3, icol * 3 : (icol + 1) * 3] = pn.Column(
        plot
    )
tabs = pn.Tabs(gspecS)
if gspecDIR is not None:
    tabs.append(gspecDIR)
if gspecVLS is not None:
    tabs.append(gspecVLS)# if maxRow>0:
tabs.append(gspecAPD)
tabs.append(gspecW8)
tabs.append(gspecFim0)
tabs.append(gspecFim1)

##################################
## save the html file
##################################

if save_elog:
    from summaries.summary_utils import prepareHtmlReport

    pageTitleFormat = "BeamlineSummary/BeamlineSummary_Run{run:04d}"
    prepareHtmlReport(tabs, expname, run, pageTitleFormat)

    if int(os.environ.get("RUN_NUM", "-1")) > 0:
        requests.post(
            os.environ["JID_UPDATE_COUNTERS"],
            json=[{"key": "<b>BeamlineSummary Plots </b>", "value": "Posted"}],
        )

if args.postStats:
    print("posting to the run tables - ipm values.")
    runtable_data = makeRunTableData(dat)
    postRunTable(runtable_data)
