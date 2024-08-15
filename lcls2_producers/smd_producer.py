#!/usr/bin/env python

import numpy as np
import psana
import time
from datetime import datetime
begin_job_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
start_job = time.time()
import argparse
import socket
import os
import logging 
import requests
import sys
from glob import glob
from requests.auth import HTTPBasicAuth
from pathlib import Path
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

if rank==0:
    print(f"MPI size: {size}")

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

## list detector names you'd like to save event by event
#save_def_dets={'save':['rix_fim2']}
# list detector names you'd like to NOT save event by event
save_def_dets={}
#save_def_dets['veto']=['lightStatus','rix_fim2','hsd','rix_fim2_raw']
#save_def_dets['save':['rix_fim0','rix_fim1' ,'crix_w8']}

#detnames = ['c_atmopal', 'hsd', 'rix_fim0', 'rix_fim1', 'crix_w8', 'c_piranha']
detnames = ['hsd', 'rix_fim0', 'rix_fim1', 'crix_w8', 'c_piranha']
intdetnames = ['andor_dir', 'andor_vls', 'andor_norm']


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
    ret_dict['andor_norm'] = []
    ret_dict['c_atmopal'] = []
    ret_dict['manta'] = []
    ret_dict['c_piranha'] = []
    ret_dict['rix_fim0'] = []
    ret_dict['rix_fim1'] = []
    ret_dict['crix_w8'] = []

    ROI_area_dict_FullWrite = {'thresADU': None}
    ROI_area_dict_FullWrite['writeArea'] = True
    ROI_area_dict_FullWrite['calcPars'] = False
    ROI_area_dict_FullWrite['ROI'] = None

    #save the full ANDOR.
    ret_dict['andor_dir'].append(ROI_area_dict_FullWrite)
    ret_dict['andor_vls'].append(ROI_area_dict_FullWrite)
    ret_dict['andor_norm'].append(ROI_area_dict_FullWrite)
    ret_dict['c_piranha'].append(ROI_area_dict_FullWrite)
    #and the FIM waveforms.
    ret_dict['rix_fim0'].append(ROI_area_dict_FullWrite)
    ret_dict['rix_fim1'].append(ROI_area_dict_FullWrite)
    ret_dict['crix_w8'].append(ROI_area_dict_FullWrite)

    # Ideally, the ROIs for all detectors (if used) are set here in a run-dependent way
    # the ROIs are saved in the hdf5 file so that we could figure out what was used later
    # if we have a run-dependent function here, we can simply reprocess all runs w/o any further editing

    #if the settings for any ROIs will change during the experiment, please create a setup like:
    # if run <= 5:
    #   <settings1>
    # elif run <= 22:
    #   <settings2>
    ROI =  [ [[0,1023], [0,1023]] ]
    ROI_dict={'thresADU': None}
    ROI_dict['calcPars'] = False
    ROI_dict['writeArea'] = True #save only projection (defined later)
    ROI_dict['ROI'] = ROI
    ROI_dict['name'] = 'ROI0'
    ret_dict['c_atmopal'].append(ROI_dict)


   # ROI =  [ [[0,500], [0,2048]] ]
   # ROI_dict={'thresADU': None}
   # ROI_dict['calcPars']=False
   # ROI_dict['writeArea']=False #save only projection (defined later)
   # ROI_dict['ROI']=ROI
   # ROI_dict['name']='ROI0'
   # ret_dict['manta'].append(ROI_dict)
    ret_dict['manta'].append(ROI_area_dict_FullWrite)

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


# DEFINE INTEGRATING DETECTOR AND ADD ANALYSIS FUNCTIONS
def define_intdets(run):
    #intdetnames = []
    #I have to figure this one out...
    #list of keys from the event based loop! Think what to do if we want raw data here, but not later.
    #normdetnames = {'rix_fim0':'xxx'}
#    intdets = ['andor_dir','andor_vls','andor_norm']
    intdets = []

    try:
        ROIs = getROIs(run)
    except:
        ROIs = []

    common_mode=0
    for detname in intdetnames:
        det = DetObject(detname , thisrun, common_mode=common_mode)
        if det is None: continue
        if detname in ROIs:
            for iROI,ROI in enumerate(ROIs[detname]):
                thisROIFunc = ROIFunc(**ROI)
                # this is to treat the FIM data. Should move settings up.
                # adding projections down here. Should be rewritten.
                det.addFunc(thisROIFunc)
        intdets.append(det)
    return intdets

# DEFINE DETECTOR AND ADD ANALYSIS FUNCTIONS
def define_dets(run):
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
            if detname.find('fim')>=0 or detname.find('w8')>=0:
                det = DetObject(detname , thisrun, common_mode=common_mode, name='det_%s'%detname)
            else:
                det = DetObject(detname , thisrun, common_mode=common_mode)
            logger.debug(f'Instantiated det {detname}: {det}')
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
                    if detname.find('fim0')>=0:
                        fimFunc = fimSumFunc(sigROI=slice(102,108), bkgROI=slice(0,80))
                        thisROIFunc.addFunc(fimFunc)
                    if detname.find('fim1')>=0:
                        fimFunc = fimSumFunc(sigROI=slice(116,122), bkgROI=slice(0,80))
                        thisROIFunc.addFunc(fimFunc)
                    elif detname.find('crix_w8')>=0:
                        fimFunc = fimSumFunc(sigROI=slice(69,76), bkgROI=slice(0,50))
                        thisROIFunc.addFunc(fimFunc)
                  #  elif detname.find('c_atmopal')>=0:
                  #      projFunc = projectionFunc(axis=0, thresADU=None)
                  #      thisROIFunc.addFunc(projFunc)
                    elif detname.find('manta')>=0:
                        projFunc = projectionFunc(axis=1, thresADU=None)
                        thisROIFunc.addFunc(projFunc)
                    det.addFunc(thisROIFunc)

            det.storeSum(sumAlgo='calib')
            logger.debug(f'Rank {rank} Add det {detname}: {det}')
            dets.append(det)

    return dets




##########################################################
# run independent parameters
##########################################################
# These lists are either PV names, aliases, or tuples with both.
#epicsPV = ['las_fs14_controller_time']
#epicsOncePV = ['m0c0_vset', ('TMO:PRO2:MPOD:01:M2:C3:VoltageMeasure', 'MyAlias'),
#               'IM4K4:PPM:SPM:VOLT_RBV', "FOO:BAR:BAZ", ("X:Y:Z", "MCBTest"), "A:B:C"]
#epicsOncePV = [('GDET:FEE1:241:ENRC', "MyTest"), 'GDET:FEE1:242:ENRC', "FOO:BAR:BAZ"]
epicsPV = []
epicsOncePV = []
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
from smalldata_tools.common.detector_base import detData, getUserData, getUserEnvData
from smalldata_tools.lcls2.default_detectors import detOnceData
from smalldata_tools.lcls2.hutch_default import defaultDetectors
from smalldata_tools.lcls2.default_detectors import epicsDetector, genericDetector
from smalldata_tools.lcls2.DetObject import DetObject

from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, imageFunc
from smalldata_tools.ana_funcs.sparsifyFunc import sparsifyFunc
from smalldata_tools.ana_funcs.waveformFunc import getCMPeakFunc, templateFitFunc
from smalldata_tools.ana_funcs.waveformFunc import hsdsplitFunc, hsdBaselineCorrectFunc
from smalldata_tools.ana_funcs.waveformFunc import hitFinderCFDFunc, hsdROIFunc
from smalldata_tools.ana_funcs.waveformFunc import fimSumFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning


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

S3DF_BASE = Path('/sdf/data/lcls/ds/')
FFB_BASE = Path('/cds/data/drpsrcf/')
PSANA_BASE = Path('/cds/data/psdm/')
PSDM_BASE = Path(os.environ.get('SIT_PSDM_DATA', S3DF_BASE))
SD_EXT = Path('./hdf5/smalldata/')
logger.debug(f"PSDM_BASE={PSDM_BASE}")

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
parser.add_argument('--postRuntable', help="postTrigger for seconday jobs", action='store_true', default=False)
parser.add_argument('--wait', help="wait for a file to appear", action='store_true', default=False)
parser.add_argument('--rawFim', help="save raw Fim data", action='store_true', default=False)
parser.add_argument('--nohsd', help="dont save HSD data", action='store_true', default=False)
parser.add_argument('--nosum', help="dont save sums", action='store_true', default=False)
parser.add_argument('--noarch', help="dont use archiver data", action='store_true', default=False)
parser.add_argument('--intg', help="use integrating detector psana mode", action='store_true', default=False)
args = parser.parse_args()

logger.debug('Args to be used for small data run: {0}'.format(args))

###### Helper Functions ##########
def get_xtc_files(base, exp, run):
    """File all xtc files for given experiment and run"""
    run_format = ''.join(['r', run.zfill(4)])
    data_dir = Path(base) / exp[:3] / exp / 'xtc'
    xtc_files = list(data_dir.glob(f'*{run_format}*'))
    if rank==0: logger.info(f'xtc file list: {xtc_files}')
    return xtc_files

def get_sd_file(write_dir, exp, hutch):
    """Generate directory to write to, create file name"""
    if write_dir is None:
        if useFFB and not onS3DF: # when on a drp node
            write_dir = FFB_BASE / hutch.lower() / exp / '/scratch' / SD_EXT
        elif onPSANA: # when on old psana system
            write_dir = PSANA_BASE / hutch.lower() / exp, SD_EXT
        elif onS3DF: # S3DF should now be the default
            write_dir = S3DF_BASE / hutch.lower() / exp / SD_EXT
        else:
            print('get_sd_file problem. Please fix.')
    logger.debug(f'hdf5 directory: {write_dir}')

    write_dir = Path(write_dir)
    h5_f_name = write_dir / f'{exp}_Run{run.zfill(4)}.h5'
    if not write_dir.exists():
        logger.info(f'{write_dir} does not exist, creating directory now.')
        try:
            write_dir.mkdir(parents=True)
        except (PermissionError, FileNotFoundError) as e:
            logger.info(f'Unable to make directory {write_dir} for output' \
                        f'exiting on error: {e}')
            sys.exit()
    if rank==0:
        logger.info('Will write small data file to {0}'.format(h5_f_name))
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
useFFB = False
onS3DF = False
onPSANA = False

if hostname.find('sdf')>=0:
    logger.debug('On S3DF')
    onS3DF = True
    if 'ffb' in PSDM_BASE.as_posix():
        useFFB = True
        # wait for files to appear
        nFiles = 0
        n_wait = 0
        max_wait = 20 # 10s wait per cycle.
        waitFilesStart=datetime.now()
        while nFiles == 0:
            if n_wait > max_wait:
                print(f"Waited {str(n_wait*10)}s, still no files available. " \
                       "Giving up, please check dss nodes and data movers. " \
                       "Exiting now.")
                sys.exit()
            xtc_files = get_xtc_files(PSDM_BASE, exp, run)
            nFiles = len(xtc_files)
            if nFiles == 0:
                print(f"We have no xtc files for run {run} in {exp} in the FFB system, " \
                      "we will wait for 10 second and check again.")
                n_wait+=1
                time.sleep(10)
        waitFilesEnd = datetime.now()
        print(f"Files appeared after {str(waitFilesEnd-waitFilesStart)} seconds")

    xtc_files = get_xtc_files(PSDM_BASE, exp, run)
    if len(xtc_files)==0:
        print(f'We have no xtc files for run {run} in {exp} in the offline system. Exit now.')
        sys.exit()

elif hostname.find('drp')>=0:
    nFiles=0
    logger.debug('On FFB')
    waitFilesStart=datetime.now()
    while nFiles==0:
        xtc_files = get_xtc_files(FFB_BASE, hutch, run)
        nFiles = len(xtc_files)
        if nFiles == 0:
            if not args.wait:
                print("We have no xtc files for run %s in %s in the FFB system,"\
                      "Quitting now.")
                sys.exit()
            else:
                print("We have no xtc files for run %s in %s in the FFB system," \
                      "we will wait for 10 second and check again."%(run,exp))
                time.sleep(10)
    waitFilesEnd = datetime.now()
    print('Files appeared after %s seconds'%(str(waitFilesEnd-waitFilesStart)))
    useFFB = True


# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.directory, exp, hutch)

# Define data source name and generate data source object, don't understand all conditions yet
#os.environ['PS_SRV_NODES']='1'
#os.environ['PS_SMD_N_EVENTS']='1'

datasourcedict={'exp':exp, 'run':int(run)}

if rank==0: print('Opening the data source:')
if args.nevents<1e9:
    datasourcedict['max_events']=args.nevents
    if useFFB:
        datasourcedict['live']=True
        if not onS3DF:
            datasourcedict['dir']=f':dir=/cds/data/drpsrcf/{exp[0:3]}/{exp}/xtc'
else:
    if useFFB:
        datasourcedict['live']=True
        if not onS3DF:
            datasourcedict['dir']=f':dir=/cds/data/drpsrcf/{exp[0:3]}/{exp}/xtc'

if len(intdetnames)>0 and args.intg:
    datasourcedict['intg_det'] = intdetnames[0]  # problem is we have more than 1 int det here?
    datasourcedict['batch_size']=1
    os.environ['PS_SMD_N_EVENTS']='1'
ds = psana.DataSource(**datasourcedict)

if ds.unique_user_rank():
    print("#### DATASOURCE AND PSANA VAR INFO ####")
    print(f"Instantiated data source with arguments: {datasourcedict}")
    print(f"MPI size: {size}")
    print(f"PS_EB_NODES={os.environ.get('PS_EB_NODES')}")
    print(f"PS_SRV_NODES={os.environ.get('PS_SRV_NODES')}")
    print(f"PS_SMD_N_EVENTS={os.environ.get('PS_SMD_N_EVENTS')}") # defaults to 1000
    print(f"DS batchsize: {ds.batch_size}")
    print("#### END DATASOURCE AND PSANA VAR INFO ####")

    print('Opened the data source, now get run')

#LCLS-2: need to get run to get detectors....
thisrun = next(ds.runs())

# Generate smalldata object
if ds.unique_user_rank():
    print('Opening the h5file %s, gathering at %d'%(h5_f_name,args.gather_interval))
small_data = ds.smalldata(filename=h5_f_name, batch_size=args.gather_interval)
if ds.unique_user_rank():
    print('smalldata file has been created on rank %d'%rank)

# Not sure why, but here
if ds.unique_user_rank():
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
        default_dets.append(epicsDetector(PVlist=epicsPV, name='epicsAll', run=thisrun))
elif len(epicsPV)>0:
    default_dets.append(epicsDetector(PVlist=epicsPV, name='epicsUser', run=thisrun))

if len(epicsOncePV)>0:
    EODet = epicsDetector(PVlist=epicsOncePV, name='epicsOnce', run=thisrun)
else:
    EODet = None
EODetData = {'epicsOnce': {}}
EODetTS   = None

default_det_aliases = [det.name for det in default_dets]

if args.rawFim:
    for fim in default_det_aliases:
        if fim.find('fim')>=0: # are you really a FIM?
            default_dets.append(genericDetector(fim, run=thisrun, h5name='%s_raw'%fim))

dets = []
intdets = []
if not args.default:
    if not ds.is_srv():
        print(f"This run: {thisrun}")
        dets = define_dets(args.run)

        if args.intg:
            intdets = define_intdets(args.run)
logger.debug(f'Rank {rank} dets: {dets}')

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
                thisDet = DetObject(alias , thisrun)
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
if ds.unique_user_rank() == 0: print('And now the event loop....')

normdict={}
for det in intdets:
    normdict[det._name]={'count':0}
    normdict[det._name]['timestamp_min']=0
    normdict[det._name]['timestamp_max']=0

for evt_num, evt in enumerate(event_iter):
    det_data = detData(default_dets, evt)

    # If we don't have the epics once data, try to get it!
    if EODet is not None and EODetData['epicsOnce'] == {}:
        EODetData = detData([EODet], evt)
        EODetTS   = evt._seconds + 631152000 # Convert to linux time.

    #detector data using DetObject
    userDict = {}
    for det in dets:
        # try:
        if True:
            det.getData(evt)
            det.processFuncs()
            userDict[det._name] = getUserData(det)
            #print('userdata ',det)
            try:
                envData=getUserEnvData(det)
                if len(envData.keys()) > 0:
                    userDict[det._name+'_env'] = envData
            except:
                pass
            det.processSums()
            #print(userDict[det._name])
        # except Exception as e:
        #     print(f'Failed analyzing det {det}')
        #     print(e)
        #     pass

    #sum default data & user Data into single dict.
    det_data.update(userDict)
    #print('det_data w/ user',det_data)

    #timing - inhibit counts - collect in default data?
    #save
    if len(intdets)>0:
        userDictInt = {}
        #userDict has keys we could sum & remove!
        #normdict['inhibit']+=det_data['timing']['inhibit'][2]
        for det in intdets:
            normdict[det._name]['count']+=1
            #for now, sum up all default data....
            for k,v in det_data.items():
                if isinstance(v, dict):
                    for kk,vv in v.items():
                        sumkey=k+'_sum_'+kk
                        if k not in normdict[det._name]:
                            normdict[det._name][k]={}
                        if sumkey in normdict[det._name][k].keys():
                            normdict[det._name][k][sumkey]+=np.array(vv)
                        else:
                            normdict[det._name][k][sumkey]=np.array(vv)
                else:
                    sumkey=k+'_sum'
                    if sumkey in normdict[det._name]:
                        normdict[det._name][sumkey]+=v
                    else:
                        normdict[det._name][sumkey]=v
            normdict[det._name]['timestamp_max']=max(normdict[det._name]['timestamp_max'], evt.timestamp)
            if normdict[det._name]['timestamp_min'] == 0: normdict[det._name]['timestamp_min']=evt.timestamp
            else: normdict[det._name]['timestamp_min']=min(normdict[det._name]['timestamp_min'], evt.timestamp)
            try:
            #if True:
                det.getData(evt)
                if det.evt.dat is not None: # do I need that or would the try-except work?
                    det.processFuncs()
                    userDictInt[det._name] = {}
                    tmpdict=getUserData(det)
                    for k,v in tmpdict.items():
                        userDictInt[det._name]['unaligned_'+k] = v

                    try:
                        envData=getUserEnvData(det)
                        if len(envData.keys())>0:
                            userDictInt[det._name+'_env'] = envData
                    except:
                        pass

                    #save data in integrating det dictionary & reset norm dictionary
                    for k,v in normdict[det._name].items():
                        if isinstance(v, dict):
                            for kk,vv in v.items():
                                userDictInt[det._name]['unaligned_norm_'+kk] = vv
                                normdict[det._name][k][kk] = vv*0 #may not work for arrays....
                        else:
                            userDictInt[det._name]['unaligned_norm_'+k] = v
                            normdict[det._name][k] = v*0 #may not work for arrays....
                    #print(userDictInt)
                    small_data.event(evt, userDictInt)
            except:
                print(f"Bad int_det processing on evt {evt_num}")
                pass

    #hits = findHits(hsd.evt.dat)
    #store event-based data
    if det_data is not None:
        #remove data fields from the save_def_dets list
        if 'veto' in save_def_dets:
            for k in save_def_dets['veto']:
                v = det_data.pop(k, None)
            #for k,v in det_data.items():
            #    if k not in save_def_dets['veto']:
            #        save_det_data[k]=v
        if 'save' in save_def_dets:
            save_det_data={}
            for k,v in det_data.items():
                if k in save_def_dets['save']:
                    save_det_data[k]=v
            det_data = save_det_data
        #save what was selected to be saved.
        #print('SAVE ',det_data)
        small_data.event(evt, det_data)

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
            if size > 2 and rank == 2: print('Processed evt %d'%evt_num)

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
if EODet is not None:
    EODetData = detOnceData(EODet, EODetData, EODetTS, args.noarch)
    if EODetData['epicsOnce'] != {}:
        userDataCfg[EODet.name] = EODet.params_as_dict()
        Config={'UserDataCfg':userDataCfg}
        Config.update(EODetData)
    else:
        Config={'UserDataCfg':userDataCfg}
else:
    Config={'UserDataCfg':userDataCfg}
#if rank==0: print(Config)
if small_data.summary:
    small_data.save_summary(Config) # this only works w/ 1 rank!

end_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
end_job = time.time()
prod_time = (end_job-start_job)/60
if ds.unique_user_rank():
    print('########## JOB TIME: {:03f} minutes ###########'.format(prod_time))
logger.debug('rank {0} on {1} is finished'.format(rank, hostname))

#finishing up here....
small_data.done()

#if args.sort:
#    if ds.unique_user_rank():
#        import subprocess
#        cmd = ['timestamp_sort_h5', h5_f_name, h5_f_name]


#if (int(os.environ.get('RUN_NUM', '-1')) > 0):
if os.environ.get('ARP_JOB_ID', None) is not None:
    if size > 2 and rank == 2:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": "~ %d cores * %d evts"%(size,evt_num)}])
    else:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": evt_num}])
logger.debug('Saved all small data')

if args.postRuntable and ds.unique_user_rank():
    print('Posting to the run tables.')
    locStr=''
    if useFFB and not onS3DF:
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
    print('URL: ', ws_url)
    #krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    #r = requests.post(ws_url, headers=krbheaders, params={"run_num": args.run}, json=runtable_data)
    user=(args.experiment[:3]+'opr').replace('dia','mcc')
    if os.environ.get("ARP_LOCATION", None) == "S3DF":
        with open('/sdf/group/lcls/ds/tools/forElogPost.txt') as reader: 
            answer = reader.readline()
    else:
        with open('/cds/home/opr/%s/forElogPost.txt'%user,'r') as reader:
            answer = reader.readline()

    r = requests.post(ws_url, params={"run_num": args.run}, json=runtable_data, auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
    print(r)
    if det_presence!={}:
        rp = requests.post(ws_url, params={"run_num": args.run}, json=det_presence, auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
        print(rp)

