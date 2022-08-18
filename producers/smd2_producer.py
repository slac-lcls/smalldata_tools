#!/usr/bin/env python

import numpy as np
import psana
import time
from datetime import datetime
begin_job_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
import time
start_job = time.time()
import argparse
import socket
import os
import logging 
import requests
import sys
from glob import glob
from requests.auth import HTTPBasicAuth
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

########################################################## 
##
## User Input start --> 
##
########################################################## 
##########################################################
# functions for run dependant parameters
# 
# See smalldata_producer_template.py for more analysis functions
##########################################################

# 1) REGIONS OF INTEREST
def getROIs(run):
    """ Set parameter for ROI analysis. Set writeArea to True to write the full ROI in the h5 file.
    See roi_rebin.py for more info
    """
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}

    #dict for HSD is simple: values is list of ROIs, no options.
    ret_dict['andor_dir'] = []
    ret_dict['andor_vls'] = []
    ret_dict['atmopal'] = []
    ret_dict['rix_fim0'] = []
    ret_dict['rix_fim1'] = []
    ret_dict['rix_fim2'] = []

    ROI_area_dict_FullWrite={'thresADU': None}
    ROI_area_dict_FullWrite['writeArea']=True
    ROI_area_dict_FullWrite['calcPars']=False
    ROI_area_dict_FullWrite['ROI']=None

    #save the full ANDOR.
    ret_dict['andor_dir'].append(ROI_area_dict_FullWrite)
    ret_dict['andor_vls'].append(ROI_area_dict_FullWrite)
    #and the FIM waveforms.
    ret_dict['rix_fim0'].append(ROI_area_dict_FullWrite)
    ret_dict['rix_fim1'].append(ROI_area_dict_FullWrite)
    ret_dict['rix_fim2'].append(ROI_area_dict_FullWrite)

    #ideally, the ROIs for all detectors (if used) are set here in a run-dependent way
    # the ROIs are saved in the hdf5 file so that we could figure out what was used later
    # if we have a run-dependent function here, we can simply reprocess all runs w/o any further editing

    #if the settings for any ROIs will change during the experiment, please create a setup like: 
    # if run <= 5:
    #   <settings1>
    # elif run <= 22:
    #   <settings2>

    ##adding ROIs to the ANDOR: this is basically two sums, meant to be signal&background.
    #ret_dict['andor_dir'].append({'name': 'ROI_sig', 'ROI': [ [0,1],[900, 1300] ]})
    #ret_dict['andor_dir'].append({'name': 'ROI_bkg', 'ROI': [ [0,1],[50, 600] ]})

    ROI =  [ [[840,900], [0,1023]] ] 

    ROI_dict={'thresADU': None}
    ROI_dict['calcPars']=False
    ROI_dict['writeArea']=False #save only projection (defined later)
    ROI_dict['ROI']=ROI
    ROI_dict['name']='ROI0'
    ret_dict['atmopal'].append(ROI_dict)
       
    #this is currently saving the full traces for input detectors.

    #ROIs of the HSDs are set here: we are setting the ranges here
    #[0,-1] is the full waveform (115k pixels). 
    #Names will be as: full_ROI_hsd_X_Y
    #   X is hsd# (0-4)
    #   Y is the ROI # - usually 0 or higher, depending on the list length
    hsd_dict = {}
    #hsd_dict['hsd_0']=[[25000,31500]]
    if not args.nohsd:
        #hsd_dict['hsd_0']=[[0,-1],[25000,31500]]
        hsd_dict['hsd_0']=[0,-1]
        hsd_dict['hsd_1']=[0,-1]
        hsd_dict['hsd_2']=[0,-1]
        ret_dict['hsd'] = hsd_dict

    return ret_dict


# DEFINE DETECTOR AND ADD ANALYSIS FUNCTIONS
def define_dets(run):
    detnames = ['andor_dir', 'andor_vls', 'hsd','rix_fim0' ,'rix_fim1' ,'rix_fim2', 'atmopal']
    dets = []
    # Load DetObjectFunc parameters (if defined)
    try:
        ROIs = getROIs(run)
    except:
        ROIs = []

    for detname in detnames:
        havedet = (detname in thisrun.detnames)
        # Common mode
        common_mode=0
        if havedet:
            if detname.find('fim')>=0:
                det = DetObject_lcls2(detname , thisrun, common_mode=common_mode, name='det_%s'%detname)
            else:
                det = DetObject_lcls2(detname , thisrun, common_mode=common_mode)
            # Analysis functions
            # ROIs:
            
            #HSD need special treatment due to their data structure.
            if detname.find('hsd')>=0:# and not args.nohsd:
                hsdsplit = hsdsplitFunc(writeHsd=False)
                if detname in ROIs:
                    for sdetname in ROIs[detname]:
                        funcname='%s__%s'%(sdetname, 'ROI')
                        RF = hsdROIFunc(name='%s__%s'%(sdetname, 'ROI'),writeArea=True, ROI = ROIs[detname][sdetname])
                        hsdsplit.addFunc(RF)
                det.addFunc(hsdsplit)

            elif detname in ROIs:
                for iROI,ROI in enumerate(ROIs[detname]):
                    thisROIFunc = ROIFunc(**ROI)
                    # this is to treat the FIM data. Should move settings up.
                    # adding projections down here. Should be rewritten.
                    if detname.find('fim')>=0:
                        fimFunc = fimSumFunc(sigROI=slice(105,135),bkgROI=slice(0,50))
                        thisROIFunc.addFunc(fimFunc)
                    elif detname.find('atmopal')>=0:
                        projFunc = projectionFunc(axis=0, thresADU=None)
                        thisROIFunc.addFunc(projFunc)
                    det.addFunc(thisROIFunc)

            det.storeSum(sumAlgo='calib')
            dets.append(det)
    return dets
    


##########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw'] 
epicsPV = []
#fix timetool calibration if necessary
#ttCalib=[0.,2.,0.]
ttCalib=[]
#ttCalib=[1.860828, -0.002950]
#decide which analog input to save & give them nice names
#aioParams=[[1],['laser']]
aioParams=[]
########################################################## 
##
## <-- User Input end
##
##########################################################

##########################################################
# Custom exception handler to make job abort if a single rank fails.
# Avoid jobs hanging forever and report actual error message to the log file.
import traceback as tb

def global_except_hook(exctype, value, exc_traceback):
    tb.print_exception(exctype, value, exc_traceback)
    sys.stderr.write("except_hook. Calling MPI_Abort().\n")
    sys.stdout.flush() # Command to flush the output - stdout
    sys.stderr.flush() # Command to flush the output - stderr
    # Note: mpi4py must be imported inside exception handler, not globally.     
    import mpi4py.MPI
    mpi4py.MPI.COMM_WORLD.Abort(1)
    sys.__excepthook__(exctype, value, exc_traceback)
    return

sys.excepthook = global_except_hook
##########################################################


# General Workflow
# This is meant for arp which means we will always have an exp and run
# Check if this is a current experiment
# If it is current, check in ffb for xtc data, if not there, default to psdm

fpath=os.path.dirname(os.path.abspath(__file__))
fpathup = '/'.join(fpath.split('/')[:-1])
sys.path.append(fpathup)
if rank==0: print(fpathup)

from smalldata_tools.utilities import printMsg, checkDet
from smalldata_tools.SmallDataUtils import setParameter, getUserData, getUserEnvData
from smalldata_tools.SmallDataUtils import defaultDetectors, detData
from smalldata_tools.SmallDataDefaultDetector import lcls2_epicsDetector, genlcls2Detector
from smalldata_tools.DetObject_lcls2 import DetObject_lcls2
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, imageFunc
from smalldata_tools.ana_funcs.sparsifyFunc import sparsifyFunc
from smalldata_tools.ana_funcs.waveformFunc import getCMPeakFunc, templateFitFunc
from smalldata_tools.ana_funcs.waveformFunc import hsdsplitFunc, hsdBaselineCorrectFunc
from smalldata_tools.ana_funcs.waveformFunc import hitFinderCFDFunc, hsdROIFunc
from smalldata_tools.ana_funcs.waveformFunc import fimSumFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning
from smalldata_tools.ana_funcs.smd_svd import svdFit

# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
HUTCHES = [
	'TMO',
	'RIX',
	'AMO',
	'SXR',
	'XPP',
	'XCS',
	'MFX',
	'CXI',
	'MEC',
	'DIA'
]

FFB_BASE = '/cds/data/drpsrcf'
PSDM_BASE = '/reg/d/psdm'
SD_EXT = '/hdf5/smalldata'

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run', type=str, default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--stn', help='hutch station', type=int, default=0)
parser.add_argument('--nevents', help='number of events', type=int, default=1e9)
parser.add_argument('--directory', help='directory for output files (def <exp>/hdf5/smalldata)')
parser.add_argument('--gather_interval', help='gather interval', type=int, default=100)
parser.add_argument('--norecorder', help='ignore recorder streams', action='store_true', default=False)
parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
parser.add_argument('--epicsAll', help='store all epics PVs', action='store_true', default=False)
parser.add_argument('--full', help='store all data (please think before usig this)', action='store_true', default=False)
parser.add_argument('--default', help='store only minimal data', action='store_true', default=False)
parser.add_argument('--image', help='save everything as image (use with care)', action='store_true', default=False)
parser.add_argument('--tiff', help='save all images also as single tiff (use with even more care)', action='store_true', default=False)
parser.add_argument("--postRuntable", help="postTrigger for seconday jobs", action='store_true', default=False)
parser.add_argument("--wait", help="wait for a file to appear", action='store_true', default=False)
parser.add_argument("--rawFim", help="save raw Fim data", action='store_true', default=False)
parser.add_argument("--nohsd", help="dont save HSD data", action='store_true', default=False)
parser.add_argument("--nosum", help="dont save sums", action='store_true', default=False)
args = parser.parse_args()

logger.debug('Args to be used for small data run: {0}'.format(args))

###### Helper Functions ##########

def get_xtc_files(base, hutch, run):
	"""File all xtc files for given experiment and run"""
	run_format = ''.join(['r', run.zfill(4)])
	data_dir = ''.join([base, '/', hutch.lower(), '/', exp, '/xtc'])
	xtc_files = glob(''.join([data_dir, '/', '*', '-', run_format, '*']))

	return xtc_files

def get_sd_file(write_dir, exp, hutch):
    """Generate directory to write to, create file name"""
    if write_dir is None:
        if useFFB:
            write_dir = ''.join([FFB_BASE, '/', hutch.lower(), '/', exp, '/scratch', SD_EXT])
        else:
            write_dir = ''.join([PSDM_BASE, '/', hutch.lower(), '/', exp, SD_EXT])
    if args.default:
        if useFFB:
            write_dir = write_dir.replace('hdf5','hdf5_def')
        else:
            write_dir = write_dir.replace('hdf5','scratch')
    h5_f_name = ''.join([write_dir, '/', exp, '_Run', run.zfill(4), '.h5'])
    if not os.path.isdir(write_dir):
        logger.info('{0} does not exist, creating directory'.format(write_dir))
        try:
            os.mkdir(write_dir)
        except OSError as e:
            logger.info('Unable to make directory {0} for output, exiting: {1}'.format(write_dir, e))
            sys.exit()
    logger.debug('Will write small data file to {0}'.format(h5_f_name))
    return h5_f_name

##### START SCRIPT ########

# Define hostname
hostname = socket.gethostname()

# Parse hutch name from experiment and check it's a valid hutch
exp = args.experiment
run = args.run
station = args.stn
logger.debug('Analyzing data for EXP:{0} - RUN:{1}'.format(args.experiment, args.run))

begin_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')

hutch = exp[:3].upper()
if hutch not in HUTCHES:
	logger.debug('Could not find {0} in list of available hutches'.format(hutch))
	sys.exit()	

xtc_files = []
# If experiment matches, check for files in ffb
useFFB=False
#with the new FFB, no need to check both on & offline as systems are independant.
if hostname.find('drp')>=0:
    nFiles=0
    logger.debug('On FFB')
    waitFilesStart=datetime.now()
    while nFiles==0:
        xtc_files = get_xtc_files(FFB_BASE, hutch, run)
        print (xtc_files)
        nFiles = len(xtc_files)
        if nFiles == 0:
            if not args.wait:
                print('We have no xtc files for run %s in %s in the FFB system, we will quit')
                sys.exit()
            else:
                print('We have no xtc files for run %s in %s in the FFB system, we will wait for 10 second and check again.'%(run,exp))
                time.sleep(10)
    waitFilesEnd=datetime.now()
    print('Files appeared after %s seconds'%(str(waitFilesEnd-waitFilesStart)))
    useFFB = True

# If not a current experiment or files in ffb, look in psdm
else:
    logger.debug('Not on FFB, use offline system')
    xtc_files = get_xtc_files(PSDM_BASE, hutch, run)
    if len(xtc_files)==0:
        print('We have no xtc files for run %s in %s in the offline system'%(run,exp))
        sys.exit()

print('deifne output')
# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.directory, exp, hutch)
print(h5_f_name)
#if args.default:
#    if useFFB:
#        h5_f_name = h5_f_name.replace('hdf5','hdf5_def')
#    else:
#        h5_f_name = h5_f_name.replace('hdf5','scratch')

# Define data source name and generate data source object, don't understand all conditions yet
os.environ['PS_SRV_NODES']='1'

if rank==0: print('Opening the data source:')
try:
    if useFFB:
        xtcdir = '/cds/data/drpsrcf/%s/%s/xtc'%(exp[0:3],exp)
        if args.nevents<1e9:
            ds = psana.DataSource(exp=exp, run=int(run), dir=xtcdir, max_events=args.nevents, live=True)
        else:
            ds = psana.DataSource(exp=exp, run=int(run), dir=xtcdir, live=True)
    else:
        if args.nevents<1e9:
            ds = psana.DataSource(exp=exp, run=int(run), max_events=args.nevents)
        else:
            ds = psana.DataSource(exp=exp, run=int(run))
except Exception as e:
    logger.info('Could not instantiate DataSource with {0},{1}: {2}'.format(exp, run, e))
    sys.exit()
#LCLS-2: need to get run to get detectors....
if rank==0: print('Opened the data source, now get run')
thisrun = next(ds.runs())

# Generate smalldata object
print('Opening the h5file %s, gathering at %d'%(h5_f_name,args.gather_interval))
small_data = ds.smalldata(filename=h5_f_name, batch_size=args.gather_interval)
print('smalldata file has been created on rank %d'%rank)

# Not sure why, but here
if rank is 0:
    logger.info('psana conda environment is {0}'.format(os.environ['CONDA_DEFAULT_ENV']))
########################################################## 
##
## Setting up the default detectors
##
########################################################## 
default_dets = defaultDetectors(hutch.lower(), thisrun)
#
# add stuff here to save all EPICS PVs.
#
if args.full or args.epicsAll:
    epicsPV = [k[0] for k in thisrun.epicsinfo]
    if len(epicsPV)>0:
        try:
            logger.info('epicsStore names for epicsAll', epicsPV)
        except:
            pass
        logger.info('adding all epicsPVs....')
        default_dets.append(lcls2_epicsDetector(PVlist=epicsPV, name='epicsAll', run=thisrun))
elif len(epicsPV)>0:
    default_dets.append(lcls2_epicsDetector(PVlist=epicsPV, name='epicsUser', run=thisrun))

default_det_aliases = [det.name for det in default_dets]
#default_det_aliases = [det.name for det in default_dets if det.name.find('fim')<0]
#print(default_det_aliases)
if args.rawFim:
    for fim in default_det_aliases:
        if fim.find('fim')>=0: #are you really a FIM?
            default_dets.append(genlcls2Detector(fim, run=thisrun, h5name='%s_raw'%fim))

dets = []
if not args.default:
    #try-except as not every rank seems to know thisdet....
    try:
        dets = define_dets(args.run)
        print([d for d in dets])
    except:
        pass

det_presence={}
if args.full:
    try:
        aliases = [ dn for dn in thisrun.detnames ]
        vetoDets = [ 'epicsinfo'] #at least for run 339 of rixx43518

        for alias in aliases:
            det_presence[alias]=1
            if alias in default_det_aliases: continue
            if alias in vetoDets: continue
            try:
                thisDet = DetObject_lcls2(alias , thisrun)
                if alias.find('hsd')>=0 and not args.nohsd:
                    hsdsplit = hsdsplitFunc()
                    thisDet.addFunc(hsdsplit)
                else:
                    fullROI = ROIFunc(writeArea=True)
                    thisDet.addFunc(fullROI)
                print('adding detector for %s'%alias)
                dets.append(thisDet)
            except:
                pass

    except:
        pass

event_iter = thisrun.events()

evt_num=-1 #set this to default until I have a useable rank for printing updates...
if rank==0: print('And now the event loop....')
for evt_num, evt in enumerate(event_iter):

    det_data = detData(default_dets, evt)
    if det_data is not None:
        small_data.event(evt, det_data)

    #detector data using DetObject 
    userDict = {}
    for det in dets:
        try:
            #this should be a plain dict. Really.
            det.getData(evt)
            det.processFuncs()
            userDict[det._name]=getUserData(det)
            #print('userdata ',det)
            try:
                envData=getUserEnvData(det)
                if len(envData.keys())>0:
                    userDict[det._name+'_env']=envData
            except:
                pass
            det.processSums()
            #print(userDict[det._name])
        except:
            pass
            
    #hits = findHits(hsd.evt.dat)
    small_data.event(evt,userDict)


    #the ARP will pass run & exp via the enviroment, if I see that info, the post updates
    if ((evt_num<10) or(evt_num<100 and (evt_num%10)==0) or (evt_num<1000 and evt_num%100==0) or (evt_num%1000==0)):
        if (int(os.environ.get('RUN_NUM', '-1')) > 0):
            if not ((evt_num<1000 and evt_num%100==0) or (evt_num%1000==0)): continue
            try:
                if size == 1:
                    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Current Event</b>", "value": evt_num+1}])
                elif size > 2 and rank == 2:
                    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Current Event / rank </b>", "value": evt_num+1}])
            except:
                if rank==0: print('Processed evt %d'%evt_num)
        else:
            if rank==0: print('Processed evt %d'%evt_num)

print('Sums:')
sumDict={'Sums': {}}
for det in dets:
    for key in det.storeSum().keys():
        try:
            sumData=small_data.sum(det.storeSum()[key])
            sumDict['Sums']['%s_%s'%(det._name, key)]=sumData
        except:
            print('Problem with data sum for %s and key %s'%(det._name,key))
if len(sumDict['Sums'].keys())>0 and small_data.summary:
    small_data.save_summary(sumDict)

userDataCfg={}
for det in default_dets:
    #make a list of configs not to be saved as lists of strings don't work in ps-4.2.5
    noConfigSave = ['scan']
    if det.name not in noConfigSave:
        userDataCfg[det.name] = det.params_as_dict()
for det in dets:
    try:
        userDataCfg[det._name] = det.params_as_dict()
    except:
        userDataCfg[det.name] = det.params_as_dict()
Config={'UserDataCfg':userDataCfg}
#if rank==0: print(Config)
if small_data.summary:
    small_data.save_summary(Config) # this only works w/ 1 rank!

end_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
end_job = time.time()
prod_time = (end_job-start_job)/60
if rank==0:
    print('########## JOB TIME: {:03f} minutes ###########'.format(prod_time))
logger.info('rank {0} on {1} is finished'.format(rank, hostname))

#finishing up here....
try:
    small_data.save()
except:
    small_data.done()
    pass

if (int(os.environ.get('RUN_NUM', '-1')) > 0):
    if size > 2 and rank == 2:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": "~ %d cores * %d evts"%(size,evt_num)}])
    else:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": evt_num}])
logger.info('Saved all small data')

if args.postRuntable and rank>0:
            
    print('posting to the run tables.')
    locStr=''
    if useFFB:
        locStr='_ffb'
    try: 
        runtable_data = {"Prod%s_end"%locStr:end_prod_time,
                         "Prod%s_start"%locStr:begin_prod_time,
                         "Prod%s_jobstart"%locStr:begin_job_time,
                         "Prod%s_ncores"%locStr:size}
    except:
        runtable_data = {"Prod%s_end"%locStr:end_prod_time,
                         "Prod%s_start"%locStr:begin_prod_time,
                         "Prod%s_jobstart"%locStr:begin_job_time}
    if args.default:
        runtable_data["SmallData%s"%locStr]="default"
    else:
        runtable_data["SmallData%s"%locStr]="done"
    time.sleep(5)
    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(args.experiment)
    print('URL:',ws_url)
    #krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    #r = requests.post(ws_url, headers=krbheaders, params={"run_num": args.run}, json=runtable_data)
    user=(args.experiment[:3]+'opr').replace('dia','mcc')
    with open('/cds/home/opr/%s/forElogPost.txt'%user,'r') as reader:
        answer = reader.readline()
    r = requests.post(ws_url, params={"run_num": args.run}, json=runtable_data, auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
    print(r)
    if det_presence!={}:
        rp = requests.post(ws_url, params={"run_num": args.run}, json=det_presence, auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
        print(rp)

# Debug stuff
# How to implement barrier from ds?
# if ds.rank == 0:
#    time.sleep(60)#ideally should use barrier or so to make sure all cores have fnished.
#     if 'temp' in h5_f_name: # only for hanging job investigation (to be deleted later)
#         h5_f_name_2 = get_sd_file(None, exp, hutch)
#         os.rename(h5_f_name, h5_f_name_2)
#         logger.debug('Move file from temp directory')
