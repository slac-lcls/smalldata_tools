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
from importlib import import_module
from mpi4py import MPI
COMM = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

logger = logging.getLogger(__name__)
# logging.basicConfig(level=logging.DEBUG)
log_level = 'INFO'
# log_level = DEBUG
log_format = '[ %(asctime)s | %(levelname)-8s | %(filename)s] %(message)s'
logging.basicConfig(format=log_format, level=log_level)


if rank == 0:
    logger.info(f"MPI size: {size}")
    logger.info('psana conda environment is {0}'.format(os.environ['CONDA_DEFAULT_ENV']))


# DEFINE DETECTOR AND ADD ANALYSIS FUNCTIONS
def define_dets(run, det_list):
    # Load DetObjectFunc parameters (if defined)
    # Assumes that the config file with function parameters definition
    # has been imported under "config"
    
    rois_args = [] # special implementation as a list to support multiple ROIs.
    fims_args = {}

    # Get the functions arguments from the production config
    if 'getROIs' in dir(config):
        rois_args = config.getROIs(run)
    if 'getFIMs' in dir(config):
        fims_args = config.getFIMs(run)

    dets = []

    for detname in det_list:
        havedet = (detname in thisrun.detnames)
        
        # Common mode (default: None)
        common_mode = None
        
        if not havedet:
            continue

        if detname.find('fim')>=0 or detname.find('w8')>=0: 
            # why different name for w8?
            det = DetObject(detname , thisrun, 
                            common_mode=common_mode, name=f'det_{detname}')
        else:
            det = DetObject(detname, thisrun, common_mode=common_mode)
        logger.debug(f'Instantiated det {detname}: {det}')

        # HSD need special treatment due to their data structure. 
        # TO REVIEW / REVISE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if detname.find('hsd')>=0:# and not args.nohsd:
            hsdsplit = hsdsplitFunc(writeHsd=False)
            if detname in rois_args:
                for sdetname in rois_args[detname]:
                    funcname='%s__%s'%(sdetname, 'ROI')
                    RF = hsdROIFunc(name='%s__%s'%(sdetname, 'ROI'), 
                                    writeArea=True, 
                                    ROI=ROIs[detname][sdetname])
                    hsdsplit.addFunc(RF)
            det.addFunc(hsdsplit)

        elif detname in rois_args:
            for iROI,ROI in enumerate(rois_args[detname]):
                try:
                    proj_dir = ROI.pop('proj_dir')
                except:
                    proj_dir = None

                thisROIFunc = ROIFunc(**ROI)
                if proj_dir is not None:
                    thsiROIFunc.addfunc(projectionFunc(axis=proj_dir))
                det.addFunc(thisROIFunc)


        if detname in fims_args:
            # Special treatment for the FIM mimicking a potential FEX config.
            fimFunc = fimSumFunc(**fims_args[detname])
            det.addFunc(fimFunc)
            
        det.storeSum(sumAlgo='calib')
        logger.debug(f'Rank {rank} Add det {detname}: {det}')
        dets.append(det)

    return dets



# General Workflow
# This is meant for arp which means we will always have an exp and run

fpath = os.path.dirname(os.path.abspath(__file__))
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


# Constants
HUTCHES = [
	'TMO',
	'RIX',
	'UED'
]

S3DF_BASE = Path('/sdf/data/lcls/ds/')
FFB_BASE = Path('/cds/data/drpsrcf/')
PSDM_BASE = Path(os.environ.get('SIT_PSDM_DATA', S3DF_BASE))
SD_EXT = Path('./hdf5/smalldata/')
logger.debug(f"PSDM_BASE={PSDM_BASE}")

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
                    default=0)
parser.add_argument('--directory',
                    help='directory for output files (def <exp>/hdf5/smalldata)')
parser.add_argument('--gather_interval',
                    help='gather interval',
                    type=int,
                    default=100)
parser.add_argument('--norecorder',
                    help='ignore recorder streams',
                    action='store_true',
                    default=False)
parser.add_argument('--url',
                    default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
parser.add_argument('--epicsAll', help='store all epics PVs',
                    action='store_true',
                    default=False)
parser.add_argument('--full',
                    help='store all data (please think before usig this)',
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
                    help='save all images also as single tiff (use with even more care)',
                    action='store_true',
                    default=False)
parser.add_argument('--postRuntable',
                    help="postTrigger for seconday jobs",
                    action='store_true',
                    default=False)
parser.add_argument('--wait',
                    help="wait for a file to appear",
                    action='store_true',
                    default=False)
parser.add_argument('--rawFim',
                    help="save raw Fim data",
                    action='store_true',
                    default=False)
parser.add_argument('--nohsd',
                    help="dont save HSD data",
                    action='store_true',
                    default=False)
parser.add_argument('--nosum',
                    help="dont save sums",
                    action='store_true',
                    default=False)
parser.add_argument('--noarch',
                    help="dont use archiver data",
                    action='store_true',
                    default=False)
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
        if onS3DF: # S3DF should now be the default
            write_dir = S3DF_BASE / hutch.lower() / exp / SD_EXT
        else:
            logger.error('On an unknown system, cannot figure where to save data.')
            logger.error('Please fix or pass a write_dir argument.')
            sys.exit()
    logger.debug(f'hdf5 directory: {write_dir}')

    write_dir = Path(write_dir)
    h5_f_name = write_dir / f'{exp}_Run{run.zfill(4)}.h5'
    if not write_dir.exists():
        logger.info(f'{write_dir} does not exist, creating directory now.')
        try:
            write_dir.mkdir(parents=True)
        except (PermissionError, FileNotFoundError) as e:
            logger.error(f'Unable to make directory {write_dir} for output' \
                        f'exiting on error: {e}')
            sys.exit()
    if rank==0:
        logger.info('Will write small data file to {0}'.format(h5_f_name))
    return h5_f_name

##### START SCRIPT ########
hostname = socket.gethostname()

# Parse hutch name from experiment and check it's a valid hutch
exp = args.experiment
run = args.run
station = args.stn
logger.debug('Analyzing data for EXP:{0} - RUN:{1}'.format(args.experiment, args.run))

begin_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')

hutch = exp[:3].upper()
if hutch not in HUTCHES:
	logger.error('Could not find {0} in list of available hutches'.format(hutch))
	sys.exit()

prod_cfg = f"prod_config_{hutch.lower()}"
if rank==0:
    logger.info(f"Producer cfg file: <{prod_cfg}>.")
config = import_module(prod_cfg)


# Figure out where we are running from and check that the 
# xtc files are where we expect them.
onS3DF = False
useFFB = False
xtc_files = []

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

else:
    logger.warning('On an unknow system, things may get weird.')

# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.directory, exp, hutch)

# Create data source.
datasource_args={
    'exp': exp,
    'run' : int(run)
}
if useFFB:
    datasource_args['live'] = True

if rank==0: print('Opening the data source:')
if args.nevents != 0:
    datasource_args['max_events'] = args.nevents

# Setup if integrating detectors are requested.
if len(config.integrating_detectors) > 0:
    datasource_args['intg_det'] = config.integrating_detectors[0]  # problem is we have more than 1 int det here?
    datasource_args['batch_size'] = 1
    os.environ['PS_SMD_N_EVENTS']='1'

ds = psana.DataSource(**datasource_args)

if ds.unique_user_rank():
    print("#### DATASOURCE AND PSANA ENV VAR INFO ####")
    print(f"Instantiated data source with arguments: {datasource_args}")
    print(f"MPI size: {size}")
    print(f"PS_EB_NODES={os.environ.get('PS_EB_NODES')}")
    print(f"PS_SRV_NODES={os.environ.get('PS_SRV_NODES')}")
    print(f"PS_SMD_N_EVENTS={os.environ.get('PS_SMD_N_EVENTS')}") # defaults to 1000
    print(f"DS batchsize: {ds.batch_size}")
    print("#### END DATASOURCE AND PSANA ENV VAR INFO ####\n")

thisrun = next(ds.runs())

# Generate smalldata object
if ds.unique_user_rank():
    print('Opening the h5file %s, gathering at %d'%(h5_f_name,args.gather_interval))
small_data = ds.smalldata(filename=h5_f_name, batch_size=args.gather_interval)
if ds.unique_user_rank():
    print(f"rank: {rank}")
    print('smalldata file has been successfully created.')



##########################################################
##
## Setting up the default detectors
##
##########################################################
if not ds.is_srv(): # srv nodes do not have access to detectors.
# if True:
    default_dets = defaultDetectors(hutch.lower(), thisrun)
    if ds.unique_user_rank():
        logger.info('Default detectors loaded.')

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
    elif len(config.epicsPV) > 0:
        default_dets.append(epicsDetector(PVlist=epicsPV, name='epicsUser', run=thisrun))

    if len(config.epicsOncePV) > 0:
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
        # print(f"This run: {thisrun}")
        dets = define_dets(int(args.run), config.detectors)
        int_dets = define_dets(int(args.run), config.integrating_detectors)
    logger.info(f'Rank {rank} detectors: {[det._name for det in dets]}')
    logger.info(f'Rank {rank} integrating detectors: {[det._name for det in int_dets]}')


    det_presence={}
    if args.full:
        try:
            aliases = [ dn for dn in thisrun.detnames ]
            vetoDets = ['epicsinfo'] #at least for run 339 of rixx43518

            for alias in aliases:
                det_presence[alias] = 1
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

    evt_num=-1 #set this to default until I have a useable rank for printing updates...
    if ds.unique_user_rank() == 0: print('And now the event loop....')

    normdict={}
    for det in intdets:
        normdict[det._name]={'count' : 0,
                            'timestamp_min' : 0,
                            'timestamp_max' : 0}

event_iter = thisrun.events()

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

    # Combine default data & user data into single dict.
    det_data.update(userDict)

    #timing - inhibit counts - collect in default data?
    if len(intdets) > 0:
        userDictInt = {}
        #userDict has keys we could sum & remove!
        #normdict['inhibit']+=det_data['timing']['inhibit'][2]
        for det in intdets:
            normdict[det._name]['count'] += 1
            # for now, sum up all default data....
            for k,v in det_data.items():
                if isinstance(v, dict):
                    for kk, vv in v.items():
                        sumkey = k+'_sum_'+kk
                        if k not in normdict[det._name]:
                            normdict[det._name][k] = {}
                            normdict[det._name][k][sumkey] = np.array(vv)
                        if sumkey in normdict[det._name][k].keys():
                            normdict[det._name][k][sumkey] += np.array(vv)
                        else:
                            normdict[det._name][k][sumkey] = np.array(vv)
                else:
                    sumkey = k+'_sum'
                    if sumkey in normdict[det._name]:
                        normdict[det._name][sumkey] += v
                    else:
                        normdict[det._name][sumkey] = v

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

    #store event-based data
    # if det_data is not None:
        # DO WE STILL WANT THAT???
        # #remove data fields from the save_def_dets list
        # if 'veto' in save_def_dets:
        #     for k in save_def_dets['veto']:
        #         v = det_data.pop(k, None)
        #     #for k,v in det_data.items():
        #     #    if k not in save_def_dets['veto']:
        #     #        save_det_data[k]=v
        # if 'save' in save_def_dets:
        #     save_det_data={}
        #     for k,v in det_data.items():
        #         if k in save_def_dets['save']:
        #             save_det_data[k]=v
        #     det_data = save_det_data
        # #save what was selected to be saved.
        # #print('SAVE ',det_data)
        
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

if not ds.is_srv():
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
            Config={'UserDataCfg' : userDataCfg}
    else:
        Config={'UserDataCfg' : userDataCfg}
    #if rank==0: print(Config)
    if small_data.summary:
        small_data.save_summary(Config) # this only works w/ 1 rank!

# Finishing up:
# The filesystem seems to make smalldata.done fail. Some dirty tricks
# are needed here.
# Hopefully this can be removed soon.
print(f"Smalldata type for rank {rank}: {small_data._type}")
import h5py

if small_data._type == 'server':
    print(f"Close smalldata server file on {rank}")
    # flush the data caches (in case did not hit cache_size yet)
    for dset, cache in small_data._server._cache.items():
        if cache.n_events > 0:
            small_data._server.write_to_file(dset, cache)
    time.sleep(1)
    small_data._server.file_handle.close()
    time.sleep(15)

print(f"smalldata.done() on rank {rank}")
small_data.done()


end_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
end_job = time.time()
prod_time = (end_job-start_job)/60
if ds.unique_user_rank():
    print('########## JOB TIME: {:03f} minutes ###########'.format(prod_time))
logger.debug('rank {0} on {1} is finished'.format(rank, hostname))



#if args.sort:
#    if ds.unique_user_rank():
#        import subprocess
#        cmd = ['timestamp_sort_h5', h5_f_name, h5_f_name]


if os.environ.get('ARP_JOB_ID', None) is not None:
    if size > 2 and rank == 2:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": "~ %d cores * %d evts"%(size,evt_num)}])
    else:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": evt_num}])

if args.postRuntable and ds.unique_user_rank():
    print('Posting to the run tables.')
    locStr=''
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
    logger.debug('URL: ', ws_url)
    user = (args.experiment[:3]+'opr').replace('dia','mcc')
    if os.environ.get("ARP_LOCATION", None) == "S3DF":
        with open('/sdf/group/lcls/ds/tools/forElogPost.txt') as reader: 
            answer = reader.readline()

    r = requests.post(ws_url, params={"run_num": args.run}, json=runtable_data, auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
    logger.debug(r)
    if det_presence != {}:
        rp = requests.post(ws_url, params={"run_num": args.run}, json=det_presence, auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
        logger.debug(rp)

