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
    ret_dict['tmoopal'] = []
    ret_dict['manta'] = []
    ret_dict['hsd'] = []
    ret_dict['andor'] = []
    if run>0:
        roi_dict = {}
        roi_dict['ROIs'] = [ [[157,487], [294,598]] ] # can define more than one ROI
        ROIs =  [ [[157,487], [294,598]] ] # can define more than one ROI
        ROI_area_dict={'thresADU': None}
        ROI_area_dict['writeArea']=True
        for iROI,ROI in enumerate(ROIs):
            ROI_area_dict['ROI']=ROI
            ROI_area_dict['name']='ROI_%d'%iROI
            ret_dict['tmoopal'].append(ROI_area_dict)
            ret_dict['manta'].append(ROI_area_dict)

        #this is currently saving the full traces for input detectors.
        ROI_area_dict2={'thresADU': None}
        ROI_area_dict2['writeArea']=True
        ROI_area_dict2['calcPars']=False
        ROI_area_dict2['ROI']=None
        ret_dict['hsd'].append(ROI_area_dict2)
        ret_dict['andor'].append(ROI_area_dict2)

    return ret_dict


# DEFINE DETECTOR AND ADD ANALYSIS FUNCTIONS
def define_dets(run):
    detnames = ['tmoopal', 'hsd', 'andor']
    #detnames = ['manta']
    #detnames = []
    dets = []
    # Load DetObjectFunc parameters (if defined)
    try:
        ROIs = getROIs(run)
    except:
        ROIs = []
    try:
        az = getAzIntParams(run)
    except:
        az = []
    try:
        phot = getPhotonParams(run)
    except:
        phot = []
    try:
        svd = getSvdParams(run)
    except:
        svd = []

    for detname in detnames:
        havedet = (detname in thisrun.detnames)
        # Common mode
        common_mode=0
        if havedet:
            print('have detname ',detname)
            det = DetObject_lcls2(detname , thisrun, common_mode=common_mode)
            
            # Analysis functions
            # ROIs:
            if detname in ROIs:
                for iROI,ROI in enumerate(ROIs[detname]):
                    det.addFunc(ROIFunc(**ROI))
            # Azimuthal binning
            if detname in az:
                det.addFunc(azimuthalBinning(**az[detname]))
            # Photon count
            if detname in phot:
                det.addFunc(photonFunc(**phot[detname]))
            # SVD waveform analysis
            if detname in svd:
                det.addFunc(svdFit(**svd[detname]))

            det.storeSum(sumAlgo='calib')
            #det.storeSum(sumAlgo='calib_img')
            #det.storeSum(sumAlgo='square_img')
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
print(fpathup)

from smalldata_tools.utilities import printMsg, checkDet
from smalldata_tools.SmallDataUtils import setParameter, getUserData, getUserEnvData
from smalldata_tools.SmallDataUtils import defaultDetectors, detData
from smalldata_tools.SmallDataDefaultDetector import lcls2_epicsDetector, genlcls2Detector
from smalldata_tools.DetObject import DetObject_lcls2
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, sparsifyFunc, imageFunc
from smalldata_tools.ana_funcs.waveformFunc import getCMPeakFunc, templateFitFunc
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
            write_dir = ''.join([FFB_BASE, '/', hutch, '/', exp, '/scratch', SD_EXT])
        else:
            write_dir = ''.join([PSDM_BASE, '/', hutch, '/', exp, SD_EXT])
    if args.default:
        if useFFB:
            write_dir = write_dir.replace('hdf5','hdf5_def')
        else:
            write_dir = write_dir.replace('hdf5','scratch')
    h5_f_name = ''.join([write_dir, '/', exp, '_Run', run.zfill(4), '.h5'])
    if not os.path.isdir(write_dir):
        logger.debug('{0} does not exist, creating directory'.format(write_dir))
        try:
            os.mkdir(write_dir)
        except OSError as e:
            logger.debug('Unable to make directory {0} for output, exiting: {1}'.format(write_dir, e))
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
    print('Files appeared after %s seconds'%(str(waitFilesStart-waitFilesEnd)))
    useFFB = True

# If not a current experiment or files in ffb, look in psdm
else:
    logger.debug('Not on FFB, use offline system')
    xtc_files = get_xtc_files(PSDM_BASE, hutch, run)
    if len(xtc_files)==0:
        print('We have no xtc files for run %s in %s in the offline system'%(run,exp))
        sys.exit()

# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.directory, exp, hutch)
#if args.default:
#    if useFFB:
#        h5_f_name = h5_f_name.replace('hdf5','hdf5_def')
#    else:
#        h5_f_name = h5_f_name.replace('hdf5','scratch')

# Define data source name and generate data source object, don't understand all conditions yet
ds_name = ''.join(['exp=', exp, ':run=', run, ':smd:live'])
#ds_name = ''.join(['exp=', exp, ':run=', run, ':smd'])
try:
    if useFFB:
        xtcdir = '/cds/data/drpsrcf/%s/%s/xtc'%(exp[0:3],exp)
        ds = psana.DataSource(exp=exp, run=int(run), dir=xtcdir)
    else:
        ds = psana.DataSource(exp=exp, run=int(run))
except Exception as e:
    logger.debug('Could not instantiate MPIDataSource with {0}: {1}'.format(ds_name, e))
    sys.exit()
#LCLS-2: need to get run to get detectors....
thisrun = next(ds.runs())

# Generate smalldata object
print('Opening the h5file %s'%h5_f_name)
small_data = ds.smalldata(filename=h5_f_name, batch_size=args.gather_interval)

# Not sure why, but here
if rank is 0:
    logger.debug('psana conda environment is {0}'.format(os.environ['CONDA_DEFAULT_ENV']))
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
    logger.debug('epicsStore names', epicsPV)
    if len(epicsPV)>0:
        logger.debug('adding all epicsPVs....')
        default_dets.append(lcls2_epicsDetector(PVlist=epicsPV, name='epicsAll', run=thisrun))
elif len(epicsPV)>0:
    default_dets.append(lcls2_epicsDetector(PVlist=epicsPV, name='epicsUser', run=thisrun))

default_det_aliases = [det.name for det in default_dets]
print(default_det_aliases)
if args.rawFim:
    for fim in default_det_aliases:
        if fim.find('fim')>=0: #are you really a FIM?
            default_dets.append(genlcls2Detector(fim, run=thisrun, h5name='%s_raw'%fim))

if not args.default:
    dets = define_dets(args.run)
else:
    dets = []

det_presence={}
if args.full:
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
        if alias=='MEC-XT2-BMMON-02' and args.experiment.find('mec')>=0:  continue
        if alias=='MEC-XT2-BMMON-03' and args.experiment.find('mec')>=0:  continue

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
            thisDet = DetObject(alias, ds.env(), int(run), name=alias)
            hasGeom=False
            for keyword in ['cs','Cs','epix','Epix','jungfrau','Jungfrau']:
                if alias.find(keyword)>=0 and args.image: hasGeom=True
            if hasGeom:
                fullROI=ROIFunc()
                fullROI.addFunc(imageFunc(coords=['x','y']))
            else:    
                fullROI = ROIFunc(writeArea=True)
            thisDet.addFunc(fullROI)
            dets.append(thisDet)
        except:
            pass
print('postfull')

if args.tiff:
    dirname = '/reg/d/psdm/%s/%s/scratch/run%d'%(args.experiment[:3],args.experiment,int(args.run))
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

#FIX: MPI size for LCLS-2
max_iter = args.nevents / size
event_iter = thisrun.events()

for evt_num, evt in enumerate(event_iter):
    print(evt_num)
    if evt_num > max_iter:
        break

    det_data = detData(default_dets, evt)
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
            print(userDict[det._name])
        except:
            pass
    try:
        small_data.event(userDict)
    except:
        small_data.event(evt,userDict)

    if args.tiff:
        for key in userDict:
            for skey in userDict[key]:
                if skey.find('area')>=0 or skey.find('img')>=0:
                   if len(userDict[key][skey].shape)==2:
                       im = Image.fromarray(userDict[key][skey])
                       tiff_file = dirname+'/Run_%d_evt_%d_%s.tiff'%(int(args.run), evt_num+1, key)
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
    print('counter update if ARP')
    if (int(os.environ.get('RUN_NUM', '-1')) > 0) and ((evt_num<100&evt_num%10==0) or (evt_num<1000&evt_num%100==0) or (evt_num%1000==0)):
        try:
            if size == 1:
                requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Current Event</b>", "value": evt_num+1}])
            elif rank == 0:
                requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Current Event / rank </b>", "value": evt_num+1}])
        except:
                requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Current Event</b>", "value": evt_num+1}])

print('sums')
sumDict={'Sums': {}}
for det in dets:
    for key in det.storeSum().keys():
        sumData=small_data.sum(det.storeSum()[key])
        sumDict['Sums']['%s_%s'%(det._name, key)]=sumData
if len(sumDict['Sums'].keys())>0:
#     print(sumDict)
    try:
        small_data.save(sumDict)
    except:
        small_data.save_summary(sumDict)

userDataCfg={}
for det in default_dets:
    #make a list of configs not to be saved as lists of strings don't work in ps-4.2.5
    noConfigSave = ['epics', 'scan']
    if det.name not in noConfigSave:
        userDataCfg[det.name] = det.params_as_dict()
for det in dets:
    try:
        userDataCfg[det._name] = det.params_as_dict()
    except:
        userDataCfg[det.name] = det.params_as_dict()
Config={'UserDataCfg':userDataCfg}
print('preconfig', Config)
try:
    small_data.save(Config)
except:
    small_data.save_summary(Config)
print('postconfig')

end_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
end_job = time.time()
prod_time = (end_job-start_job)/60
if rank==0:
    print('########## JOB TIME: {:03f} minutes ###########'.format(prod_time))
logger.debug('rank {0} on {1} is finished'.format(rank, hostname))

#finishing up here....
try:
    small_data.save()
except:
    print('done')
    small_data.done()
    pass

if (int(os.environ.get('RUN_NUM', '-1')) > 0):
    if size > 1 and rank == 0:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": "~ %d cores * %d evts"%(size,evt_num)}])
    else:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": evt_num}])
logger.debug('Saved all small data')

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
