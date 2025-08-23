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
from PIL import Image
from requests.auth import HTTPBasicAuth
from pathlib import Path
from importlib import import_module

import warnings
import tables
warnings.filterwarnings("ignore", category=tables.NaturalNameWarning)

from mpi4py import MPI

COMM = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

logger = logging.getLogger(__name__)
log_format = "[ %(asctime)s | %(levelname)-3s | %(filename)s] %(message)s"
logging.basicConfig(format=log_format)
logger.setLevel(logging.INFO)  # Set level here instead of in the basic config so other loggers are  not affected

def isDropped(def_data):
    if def_data['lightStatus']['xray'] == 0: 
        return True
    return False


# DEFINE DETECTOR AND ADD ANALYSIS FUNCTIONS
def define_dets(run):

    # Load DetObjectFunc parameters (if defined)
    # Assumes that the config file with function parameters definition
    # has been imported under "config"

    rois_args = {}  
    mectt_args = {}
    azint_args = {}
    azintpyfai_args = {}
    drop_args = {}
    d2p_args = {}
    detsum_args = {}
    
    # Get the functions arguments from the production config
    if "getROIs" in dir(config):
        roi_args = config.getROIs(run)
    if "getmectt" in dir(config):
        mectt_args = config.get_mectt(run)
    if "getAzIntParams" in dir(config):
        azint_args = config.getAzIntParams(run)
    if "getAzIntPyFAIParams" in dir(config):
        azintpyfai_args = config.getAzIntPyFAIParams(run)
    if "get_droplet" in dir(config):
        drop_args = config.get_droplet(run)
    if "get_droplet2photon" in dir(config):
        d2p_args = config.get_droplet2photon(run)
    if "getDetSums" in dir(config):
        detsum_args = config.getDetSums(run)
    if "getCompressionParameters" in dir(config):
        compress_args = config.getCompressionParameters(run)

    detnames = config.detnames
                    
    dets = []
    
    # Define detectors and their associated DetObjectFuncs
    for detname in detnames:
        havedet = checkDet(ds.env(), detname)
        compression = False
        # Common mode
        if havedet:
            if detname=='': 
                # change here to specify common mode for detname if desired. Else default is used
                common_mode=0
            else:
                common_mode=None
            det = DetObject(detname ,ds.env(), int(run), common_mode=common_mode)


            # Analysis functions
            # Get compression func
            if detname in compress_args:
                compress_decompress = CompressDecompress(compress_args[detname])
                compression = True

            # ROIs:
            if detname in roi_args:
                for ROI in roi_args[detname]:
                    pjs = ROI.pop('pj',[])
                    rfunc = ROIFunc(**ROI)
                    for pj in pjs:
                        rfunc.addFunc(projectionFunc(**pj))

                    # If we want to test compression, let's add it at the beginning of the chain
                    if compression:
                        compress_decompress.addFunc(rfunc)
                        det.addFunc(compress_decompress)
                    else:
                        det.addFunc(rfunc)
            # Azimuthal binning
            if detname in azint_args:
                for azint in azint_args[detname]:
                    det.addFunc(azimuthalBinning(**azint))
            # Azimuthal binning - PyFAI
            if detname in azintpyfai_args:
                for az_pyfai in azintpyfai_args[detname]:
                    det.addFunc(azav_pyfai(**az_pyfai))
            # Droplet algo
            if detname in drop_args:
                for tdrop in drop_args[detname]:
                    if 'nData' in tdrop:
                        nData = tdrop.pop('nData')
                    else:
                        nData = None
                    func = dropletFunc(**tdrop)
                    func.addFunc(sparsifyFunc(nData=nData))
                    det.addFunc(func)
            # Droplet to photons
            if detname in d2p_args:
                for td2p in d2p_args[detname]:
                    if 'nData' in td2p:
                        nData = td2p.pop('nData')
                    else:
                        nData = None
                    # getp droplet dict
                    droplet_dict = td2p['droplet']
                    #get droplet2Photon dict
                    d2p_dict = td2p['d2p']
                    dropfunc = dropletFunc(**droplet_dict)
                    drop2phot = droplet2Photons(**d2p_dict)
                    sparsify = sparsifyFunc(nData=nData)
                    drop2phot.addFunc(sparsify)
                    dropfunc.addFunc(drop2phot)
                    det.addFunc(dropfunc)
            # MEC timetoolAzimuthal binning
            if detname in mectt_args:
                for ttt in mectt_args[detname]:
                    logger.info('make mectt func ')
                    det.addFunc(mecttFunc(**ttt))
            # summed images
            if detname in detsum_args:
                for thissum in detsum_args[detname]:
                    det.storeSum(sumAlgo=thissum)

            dets.append(det)
    return dets

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
logger.info(f'\n{fpathup}')

from smalldata_tools.utilities import printMsg, checkDet
from smalldata_tools.common.detector_base import getUserData, getUserEnvData
from smalldata_tools.lcls1.default_detectors import detData, detOnceData
from smalldata_tools.lcls1.default_detectors import ttRawDetector
from smalldata_tools.lcls1.default_detectors import epicsDetector, eorbitsDetector
from smalldata_tools.lcls1.default_detectors import bmmonDetector, ipmDetector
from smalldata_tools.lcls1.default_detectors import encoderDetector, adcDetector
from smalldata_tools.lcls1.default_detectors import xtcavDetector
from smalldata_tools.lcls1.hutch_default import defaultDetectors
from smalldata_tools.lcls1.DetObject import DetObject
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, imageFunc
from smalldata_tools.ana_funcs.mecttFunc import mecttFunc
from smalldata_tools.ana_funcs.sparsifyFunc import sparsifyFunc
from smalldata_tools.ana_funcs.waveformFunc import getCMPeakFunc, templateFitFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.droplet2Photons import droplet2Photons
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning
from smalldata_tools.ana_funcs.azav_pyfai import azav_pyfai
# from smalldata_tools.ana_funcs.smd_svd import svdFit
from smalldata_tools.ana_funcs.correlations.smd_autocorr import Autocorrelation
from smalldata_tools.ana_funcs.compression import CompressDecompress

# Constants
HUTCHES = [
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
if rank == 0:
    logger.info(f"PSDM_BASE={PSDM_BASE}")

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument('--run', 
                    help='run', 
                    type=str, 
                    default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--experiment', 
                    help='experiment name', 
                    type=str, 
                    default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--stn', 
                    help='hutch station', 
                    type=int, 
                    default=0)
parser.add_argument('--nevents', 
                    help='number of events', 
                    type=int, 
                    default=1e9)
parser.add_argument('--directory',
                    help='directory for output files (def <exp>/hdf5/smalldata)')
parser.add_argument('--gather_interval', 
                    help='gather interval', 
                    type=int, 
                    default=25)
parser.add_argument('--norecorder', 
                    help='ignore recorder streams', 
                    action='store_true', 
                    default=False)
parser.add_argument('--url',
                    default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
parser.add_argument('--epicsAll', 
                    help='store all epics PVs', 
                    action='store_true', 
                    default=False)
parser.add_argument('--full', 
                    help='store all data (please think before using this)',
                    action='store_true', 
                    default=False)
parser.add_argument('--fullSum', help='store sums for all area detectors', 
                    action='store_true', 
                    default=False)
parser.add_argument('--default', 
                    help='store only minimal data', 
                    action='store_true',
                    default=False)
parser.add_argument('--image', 
                    help='save everything as image (use with care)',
                    action='store_true', 
                    default=False)
parser.add_argument('--tiff', 
                    help='save all images as single tiff (use with extreme care)',
                    action='store_true', 
                    default=False)
parser.add_argument('--centerpix', 
                    help='do not mask center pixels for epix10k detectors.', 
                    action='store_true',
                    default=False)
parser.add_argument("--postRuntable",
                    help="postTrigger for seconday jobs", 
                    action='store_true',
                    default=False)
parser.add_argument("--wait", 
                    help="wait for a file to appear",
                    action='store_true', 
                    default=False)
parser.add_argument("--xtcav",
                    help="add xtcav processing",
                    action='store_true', 
                    default=False)
parser.add_argument("--noarch", 
                    help="dont use archiver data", 
                    action='store_true',
                    default=False)
parser.add_argument("--ttRaw",
                    help="add timetool projections",
                    action='store_true')
parser.add_argument("--config",
                    help="special configuration")
parser.add_argument("--outfilename",
                    help="special output filename")
args = parser.parse_args()
logger.debug('Args to be used for small data run: {0}'.format(args))

###### Helper Functions ##########

def get_xtc_files(base, exp, run):
    """File all xtc files for given experiment and run"""
    run_format = ''.join(['r', run.zfill(4)])
    data_dir = Path(base) / exp[:3] / exp / 'xtc'
    xtc_files = list(data_dir.glob(f'*{run_format}*'))
    if rank == 0:
        logger.info(f'xtc file list: {xtc_files}')
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
            logger.error('get_sd_file problem. Please fix.')
    logger.debug(f'hdf5 directory: {write_dir}')

    write_dir = Path(write_dir)
    if args.outfilename is not None:
        h5_f_name = write_dir / f'{args.outfilename}.h5'
    else:
        h5_f_name = write_dir / f'{exp}_Run{run.zfill(4)}.h5'
    
    if not write_dir.exists():
        if rank == 0:
            logger.info(f'{write_dir} does not exist, creating directory now.')
        try:
            write_dir.mkdir(parents=True)
        except (PermissionError, FileNotFoundError) as e:
            logger.info(f'Unable to make directory {write_dir} for output' \
                        f'exiting on error: {e}')
            sys.exit()
    if rank == 0:
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

if args.config is None:
    prod_cfg = f"prod_config_{hutch.lower()}"
else:
    prod_cfg = args.config
if prod_cfg.find('/')>=0:
    cfg = prod_cfg.split('/')[-1]
    cfgdir = prod_cfg.replace(cfg,'')
    sys.path.append(cfgdir)
    prod_cfg = cfg
if rank == 0:
    logger.info(f"Producer cfg file: <{prod_cfg}>.")
config = import_module(prod_cfg)

# Figure out where we are and where to look for data
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
                raise RuntimeError("Waited {str(n_wait*10)}s, still no files available. Giving up.")
            xtc_files = get_xtc_files(PSDM_BASE, exp, run)
            nFiles = len(xtc_files)
            if nFiles == 0:
                if rank == 0 :
                    logger.info(f"We have no xtc files for run {run} in {exp} in the FFB system, " \
                          "we will wait for 10 second and check again.")
                n_wait+=1
                time.sleep(10)
        waitFilesEnd = datetime.now()
        if rank == 0:
            logger.info(f"Files appeared after {str(waitFilesEnd-waitFilesStart)} seconds")

    xtc_files = get_xtc_files(PSDM_BASE, exp, run)
    if len(xtc_files)==0:
        raise RuntimeError(f'We have no xtc files for run {run} in {exp} in the offline system.')

elif hostname.find('drp')>=0:
    nFiles=0
    logger.debug('On FFB')
    waitFilesStart=datetime.now()
    while nFiles==0:
        xtc_files = get_xtc_files(FFB_BASE, hutch, run)
        nFiles = len(xtc_files)
        if nFiles == 0:
            if not args.wait:
                if rank == 0:
                    logger.error("We have no xtc files for run %s in %s in the FFB system,"\
                          "Quitting now.")
                sys.exit()
            else:
                if rank == 0:
                    logger.info("We have no xtc files for run %s in %s in the FFB system," \
                          "we will wait for 10 second and check again."%(run,exp))
                time.sleep(10)
    waitFilesEnd = datetime.now()
    if rank == 0:
        logger.info('Files appeared after %s seconds'%(str(waitFilesEnd-waitFilesStart)))
    useFFB = True

# If not a current experiment or files in ffb, look in psdm
else:
    logger.debug('Not on FFB or S3DF, use old offline system')
    xtc_files = get_xtc_files(PSDM_BASE, hutch, run)
    if len(xtc_files)==0:
        logger.error('We have no xtc files for run %s in %s in the offline system'%(run,exp))
        sys.exit()

# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.directory, exp, hutch)
#if args.default:
#    if useFFB:
#        h5_f_name = h5_f_name.replace('hdf5','hdf5_def')
#    else:
#        h5_f_name = h5_f_name.replace('hdf5','scratch')

# Define data source name and generate data source object
ds_name = f'exp={exp}:run={run}:smd'
if args.norecorder:
        ds_name += ':stream=0-79'
if useFFB:
        ds_name += ':live'
        if not onS3DF:
            ds_name += f':dir=/cds/data/drpsrcf/{exp[0:3]}/{exp}/xtc'
        psana.setOption('PSXtcInput.XtcInputModule.liveTimeout', 300)
        # bigger timeout so that the live mode does not fail

logger.debug(f'DataSource name: {ds_name}')
try:
    ds = psana.MPIDataSource(ds_name)
except Exception as e:
    logger.debug('Could not instantiate MPIDataSource with {0}: {1}'.format(ds_name, e))
    sys.exit()

# Generate smalldata object
small_data = ds.small_data(h5_f_name, gather_interval=args.gather_interval)


########################################################## 
##
## Setting up the default detectors
##
########################################################## 

start_setup_dets = time.time()

default_dets = defaultDetectors(hutch.lower(), env=ds.env())
if args.xtcav and not args.norecorder:
    #default_dets.append(xtcavDetector('xtcav','xtcav',method='COM'))
    default_dets.append(xtcavDetector('xtcav','xtcav'))
#adding raw timetool traces:
if args.ttRaw:
    try:
        ttRawDet = ttRawDetector(env=ds.env())
        ttRawDet.setPars({'beamOff':[-137], 'refitData': False})
        default_dets.append(ttRawDet)
    except:
        pass

#
# add stuff here to save all EPICS PVs.
#
# is someone has provided a list, save in epicsUser
if len(config.epicsPV) > 0:
    default_dets.append(epicsDetector(PVlist=config.epicsPV, name='epicsUser'))
# make a list of all PVs in data
logger.debug('epicsStore names', ds.env().epicsStore().pvNames())
if args.experiment.find('dia') >= 0:
    epicsPVlist = ds.env().epicsStore().pvNames()
else:
    epicsPVlist = ds.env().epicsStore().aliases()
if (args.full or args.epicsAll) and len(epicsPVlist)>0:
    logger.debug('adding all epicsPVs....')
    default_dets.append(epicsDetector(PVlist=epicsPV, name='epicsAll'))
#save specified list of PVs once/run, not nothing has been passed, save all.
if len(config.epicsOncePV) > 0:
    EODet = epicsDetector(PVlist=epicsOncePV, name='epicsOnce')
elif len(epicsPVlist) > 0:
    EODet = epicsDetector(PVlist=epicsPVlist, name='epicsOnce')
else:
    EODet = None

default_dets.append(eorbitsDetector())
default_det_aliases = [det.name for det in default_dets]

if not args.default:
    dets = define_dets(args.run)
else:
    dets = []

det_presence={}
if args.full or args.fullSum:
    aliases = []
    for dn in psana.DetNames():
        if dn[1]!='':
            aliases.append(dn[1])
        else:
            aliases.append(dn[0])

    for alias in aliases:
        det_presence[alias]=1
        if alias in default_det_aliases: continue
        if alias=='FEEGasDetEnergy': continue #done by mpidatasource
        if alias=='PhaseCavity':     continue #done by mpidatasource
        if alias=='EBeam':           continue #done by mpidatasource
        if alias.find('evr')>=0:     continue #done by mpidatasource
        if alias=='ControlData':     continue #done by my code
        #done in standard default detectors.
        if alias=='HX2-SB1-BMMON' and args.experiment.find('xpp')>=0:  continue
        if alias=='XPP-SB2-BMMON' and args.experiment.find('xpp')>=0:  continue
        if alias=='XPP-SB3-BMMON' and args.experiment.find('xpp')>=0:  continue
        if alias=='XPP-USB-ENCODER-02' and args.experiment.find('xpp')>=0:  continue
        if alias=='XCS-SB1-BMMON' and args.experiment.find('xcs')>=0:  continue
        if alias=='XCS-SB2-BMMON' and args.experiment.find('xcs')>=0:  continue
        if alias=='CXI-DG2-BMMON' and args.experiment.find('cxi')>=0:  continue
        if alias=='CXI-DG3-BMMON' and args.experiment.find('cxi')>=0:  continue
        if alias=='MEC-XT2-BMMON-02' and args.experiment.find('mec')>=0:  continue
        if alias=='MEC-XT2-BMMON-03' and args.experiment.find('mec')>=0:  continue

        if args.full:
            if alias.find('BMMON')>=0:
                print('append a bmmon',alias)
                default_dets.append(bmmonDetector(alias))
                continue
            elif alias.find('IPM')>=0 or alias.find('Ipm')>0:
                default_dets.append(ipmDetector(alias, savePos=True))
                continue
            elif alias.find('DIO')>=0 or alias.find('Imb')>0:
                default_dets.append(ipmDetector(alias, savePos=False))
                continue
            elif alias.find('USB')>=0:
                default_dets.append(encoderDetector(alias))
                continue
            elif alias.find('adc')>=0:
                default_dets.append(adcDetector(alias))
                continue

        try:
            thisDet = DetObject(alias, ds.env(), int(run), name=alias, maskCentral=(not args.centerpix))
            hasGeom=False
            for keyword in ['cs','Cs','epix','Epix','jungfrau','Jungfrau']:
                if alias.find(keyword)>=0 and args.image: hasGeom=True
            if hasGeom:
                fullROI=ROIFunc()
                fullROI.addFunc(imageFunc(coords=['x','y']))
            else:    
                fullROI = ROIFunc(writeArea=True)
            if args.full:
                thisDet.addFunc(fullROI)
            if args.fullSum:
                thisDet.storeSum(sumAlgo='calib_dropped')
                #thisDet.storeSum(sumAlgo='calib_dropped_square')
                thisDet.storeSum(sumAlgo='calib_thresADU1')
                thisDet.storeSum(sumAlgo='calib_max')

            dets.append(thisDet)
        except:
           pass


# save detector config data
userDataCfg={}
for det in default_dets:
    if det.name=='tt' and len(config.ttCalib)>0:
        det.setPars(config.ttCalib)
        if rank == 0:
            logger.info(f'Using user-defined tt parameters: {config.ttCalib}')
    userDataCfg[det.name] = det.params_as_dict()
for det in dets:
    try:
        userDataCfg[det._name] = det.params_as_dict()
    except:
        userDataCfg[det.name] = det.params_as_dict()

# is EOODet exists, save this later.
if EODet is None:
    Config={'UserDataCfg':userDataCfg}
    small_data.save(Config)

end_setup_dets = time.time()

if args.tiff: # this needs to be done for S3DF
    if onS3DF:
        dirname = S3DF_BASE / f"{exp[:3]}/{exp}/scratch/run{int(run)}"
    else:
        dirname = PSANA_BASE / '%s/%s/scratch/run%d'%(args.experiment[:3],args.experiment,int(args.run))
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

ds.break_after(args.nevents)
start_evt_loop = time.time()
for evt_num, evt in enumerate(ds.events()):
    if evt_num == 0 and EODet is not None:
        det_data = detOnceData(EODet, evt, args.noarch)
        if det_data[EODet.name] != {}:
            userDataCfg[EODet.name] = EODet.params_as_dict()
            Config={'UserDataCfg':userDataCfg}
            small_data.save(Config)
            small_data.save(det_data)
        else:
            Config={'UserDataCfg':userDataCfg}
            small_data.save(Config)

    def_data = detData(default_dets, evt)
    small_data.event(def_data)

    #detector data using DetObject 
    userDict = {}
    for det in dets:
        try:
            #this should be a plain dict. Really.
            det.getData(evt)
            det.processFuncs()
            userDict[det._name]=getUserData(det)
            try:
                envData=getUserEnvData(det)
                if len(envData.keys())>0:
                    userDict[det._name+'_env']=envData
            except:
                pass
            det.processSums(dropped=isDropped(def_data))
#             print(userDict[det._name])
        except:
            # handle when sum is bad for all shots on a rank (rare, but happens)
            for key in det._storeSum.keys():
                if det._storeSum[key] is None:
                    det._storeSum[key] = 0
                else:
                    det._storeSum[key] += 0

    small_data.event(userDict)

    if args.tiff:
        for key in userDict:
            for skey in userDict[key]:
                if skey.find('area')>=0 or skey.find('img')>=0:
                   if len(userDict[key][skey].shape)==2:
                       image = userDict[key][skey]
                       try:
                           mask_key = [ key for key in userDataCfg[key].keys() if key.find('mask_img')>=0 ]
                           if len(mask_key)>0:
                               maskImg = userDataCfg[key][mask_key[0]]
                               imageMasked = np.ma.array(image, mask=maskImg)
                               image = imageMasked.filled(fill_value=0)
                       except:
                           pass
                       im = Image.fromarray(image)
                       tiff_file = dirname / f"Run_{int(run)}_evt_{evt_num+1}_{key}.tiff"
                       im.save(tiff_file)

    #here you can add any data you like: example is a product of the maximumof two area detectors.
    #try:
    #    jungfrau1MMax = jungfrau1M.evt.dat.max()
    #    epix_vonHamosMax = epix_vonHamos.evt.dat.max()
    #    combDict = {'userValue': jungfrau1MMax*epix_vonHamosMax}
    #    small_data.event(combDict)
    #except:
    #    pass


    #the ARP will pass run & exp via the enviroment, if I see that info, the post updates
    if ( (evt_num<100 and evt_num%10==0) or (evt_num<1000 and evt_num%100==0) or (evt_num%1000==0)):
        if os.environ.get('ARP_JOB_ID', None) is not None:
            if ds.size == 1:
                requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Current Event</b>", "value": evt_num+1}])
            elif ds.rank == 0:
                requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Current Event / rank </b>", "value": evt_num+1}])
        else:
            if ds.size == 1:
                print('Current Event:', evt_num+1)
            elif ds.rank == 0:
                print('Current Event / rank :', evt_num+1)

end_evt_loop = time.time()

sumDict={'Sums': {}}
for det in dets:
    for key in det.storeSum().keys():
        if "max" in key:
            sumData = small_data.max(det.storeSum()[key])
        else:
            sumData = small_data.sum(det.storeSum()[key])
        sumDict['Sums']['%s_%s'%(det._name, key)]=sumData
if len(sumDict['Sums'].keys())>0:
#     print(sumDict)
    small_data.save(sumDict)

small_data.save()
small_data.close()
COMM.Barrier()

# Print duration summary
dets_time_start = (start_setup_dets-start_job)/60
dets_time_end = (end_setup_dets-start_job)/60
evt_time_start = (start_evt_loop-start_job)/60
evt_time_end = (end_evt_loop-start_job)/60
logger.debug(f"##### Timing benchmarks core {ds.rank}: ##### """)
logger.debug(f'Setup dets: \n\tStart: {dets_time_start:.2f} min\n\tEnd: {dets_time_end:.2f} min')
logger.debug(f'\tDuration:{dets_time_end-dets_time_start:.2f}')
logger.debug(f'Event loop: \n\tStart: {evt_time_start:.2f} min\n\tEnd: {evt_time_end:.2f} min')
logger.debug(f'\tDuration:{evt_time_end-evt_time_start:.2f}')
logger.debug('\n')

end_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
end_job = time.time()
if ds.rank==0:
    prod_time = (end_job-start_job)/60
    print('########## JOB TIME: {:03f} minutes ###########'.format(prod_time))

#finishing up here....
logger.debug('rank {0} on {1} is finished'.format(ds.rank, hostname))
if os.environ.get('ARP_JOB_ID', None) is not None:
    if ds.size > 1:
        if ds.rank == 0:
            requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": "~ %d cores * %d evts"%(ds.size,evt_num)},{"key": "<b>Duration</b>", "value": "%f min"%(prod_time)}])
    else:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": evt_num}])
logger.debug('Saved all small data')


# This is no broken. How to access the file under /cds/...?
# Should we put it under /sdf/ as well?
if args.postRuntable and ds.rank==0:
    print('Posting to the run tables.')
    locStr=''
    runtable_data = {"Prod%s_end"%locStr:end_prod_time,
                     "Prod%s_start"%locStr:begin_prod_time,
                     "Prod%s_jobstart"%locStr:begin_job_time,
                     "Prod%s_duration_mins"%locStr:prod_time,
                     "Prod%s_ncores"%locStr:ds.size}
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
    if os.environ.get("ARP_LOCATION", None) == "S3DF":
        with open('/sdf/group/lcls/ds/tools/forElogPost.txt') as reader:
            answer = reader.readline()
    else:
        with open('/cds/home/opr/%s/forElogPost.txt'%user,'r') as reader:
            answer = reader.readline()
    r = requests.post(ws_url, params={"run_num": args.run}, json=runtable_data, 
                      auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
    print(r)
    if det_presence!={}:
        rp = requests.post(ws_url, params={"run_num": args.run}, json=det_presence,
                           auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
        print(rp)
