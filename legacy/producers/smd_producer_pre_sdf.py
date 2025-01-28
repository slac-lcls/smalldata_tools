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
    if run>0:
        roi_dict = {}
        roi_dict['ROIs'] = [ [[1,2], [157,487], [294,598]] ] # can define more than one ROI
        roi_dict['writeArea'] = True
        roi_dict['thresADU'] = None
        ret_dict['jungfrau1M'] = roi_dict
    return ret_dict

def isDropped(def_data):
    if def_data['lightStatus']['xray'] == 0: 
        return True
    return False

##########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
# These lists are either PV names, aliases, or tuples with both.
#epicsPV = ['gon_h'] 
#epicsOncePV = [("XPP:GON:MMS:01.RBV", 'MyAlias'), 'gon_v', "XPP:GON:MMS:03.RBV",
#               "FOO:BAR:BAZ", ("X:Y:Z", "MCBTest"), 'A:B:C']
#epicsOncePV = [('GDET:FEE1:241:ENRC', "MyTest"), 'GDET:FEE1:242:ENRC', "FOO:BAR:BAZ"]
epicsOncePV = []
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


# DEFINE DETECTOR AND ADD ANALYSIS FUNCTIONS
def define_dets(run):
    detnames = ['jungfrau1M'] # add detector here
    dets = []
    
    # Load DetObjectFunc parameters (if defined)
    try:
        ROIs = getROIs(run)
    except Exception as e:
        print(f'Can\'t instantiate ROI args: {e}')
        ROIs = []
    try:
        az = getAzIntParams(run)
    except Exception as e:
        print(f'Can\'t instantiate azimuthalBinning args: {e}')
        az = []
    try:
        az_pyfai = getAzIntPyFAIParams(run)
    except Exception as e:
        print(f'Can\'t instantiate AzIntPyFAI args: {e}')
        az_pyfai = []
    try:
        phot = getPhotonParams(run)
    except Exception as e:
        print(f'Can\'t instantiate Photon args: {e}')
        phot = []
    try:
        drop = getDropletParams(run)
    except Exception as e:
        print(f'Can\'t instantiate Droplet args: {e}')
        drop = []
    try:
        drop2phot = getDroplet2Photons(run)
    except Exception as e:
        print(f'Can\'t instantiate Droplet2Photons args: {e}')
        drop2phot = []
    try:
        auto = getAutocorrParams(run)
    except Exception as e:
        print(f'Can\'t instantiate Autocorrelation args: {e}')
        auto = []
    try:
        svd = getSvdParams(run)
    except Exception as e:
        print(f'Can\'t instantiate SVD args: {e}')
        svd = []
        
    # Define detectors and their associated DetObjectFuncs
    for detname in detnames:
        havedet = checkDet(ds.env(), detname)
        # Common mode
        if havedet:
            if detname=='': 
                # change here to specify common mode for detname if desired. Else default is used
                common_mode=0
            else:
                common_mode=None
            det = DetObject(detname ,ds.env(), int(run), common_mode=common_mode)
            
            # Analysis functions
            # ROIs:
            if detname in ROIs:
                for iROI,ROI in enumerate(ROIs[detname]['ROIs']):
                    det.addFunc(ROIFunc(name='ROI_%d'%iROI,
                                        ROI=ROI,
                                        writeArea=ROIs[detname]['writeArea'],
                                        thresADU=ROIs[detname]['thresADU']))
            # Azimuthal binning
            if detname in az:
                det.addFunc(azimuthalBinning(**az[detname]))
            if detname in az_pyfai:
                det.addFunc(azav_pyfai(**az_pyfai[detname]))
            # Photon count
            if detname in phot:
                det.addFunc(photonFunc(**phot[detname]))
            # Droplet algo
            if detname in drop:
                if nData in drop:
                    nData = drop.pop('nData')
                else:
                    nData = None
                func = dropletFunc(**drop[detname])
                func.addFunc(roi.sparsifyFunc(nData=nData))
                det.addFunc(func)
            # Droplet to photons
            if detname in drop2phot:
                if 'nData' in drop2phot[detname]:
                    nData = drop2phot[detname].pop('nData')
                else:
                    nData = None
                # getp droplet dict
                droplet_dict = drop2phot[detname]['droplet']
                #get droplet2Photon dict
                d2p_dict = drop2phot[detname]['d2p']
                dropfunc = dropletFunc(**droplet_dict)
                drop2phot = droplet2Photons(**d2p_dict)
                sparsify = sparsifyFunc(nData=nData)
                drop2phot.addFunc(sparsify)
                dropfunc.addFunc(drop2phot)
                det.addFunc(dropfunc)
            # Autocorrelation
            if detname in auto:
                det.addFunc(Autocorrelation(**auto[detname]))
            # SVD waveform analysis
            if detname in svd:
                det.addFunc(svdFit(**svd[detname]))

            det.storeSum(sumAlgo='calib')
            det.storeSum(sumAlgo='calib_dropped')
            det.storeSum(sumAlgo='calib_dropped_square')
            #det.storeSum(sumAlgo='calib_img')
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
print(fpathup)

from smalldata_tools.utilities import printMsg, checkDet
from smalldata_tools.SmallDataUtils import setParameter, getUserData, getUserEnvData
from smalldata_tools.SmallDataUtils import defaultDetectors, detData, detOnceData
from smalldata_tools.SmallDataDefaultDetector import epicsDetector, eorbitsDetector
from smalldata_tools.SmallDataDefaultDetector import bmmonDetector, ipmDetector
from smalldata_tools.SmallDataDefaultDetector import encoderDetector, adcDetector
from smalldata_tools.SmallDataDefaultDetector import xtcavDetector
from smalldata_tools.DetObject import DetObject
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, imageFunc
from smalldata_tools.ana_funcs.sparsifyFunc import sparsifyFunc
from smalldata_tools.ana_funcs.waveformFunc import getCMPeakFunc, templateFitFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.droplet2Photons import droplet2Photons
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning
from smalldata_tools.ana_funcs.azav_pyfai import azav_pyfai
from smalldata_tools.ana_funcs.smd_svd import svdFit
from smalldata_tools.ana_funcs.correlations.smd_autocorr import Autocorrelation

logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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

S3DF_BASE = '/sdf/data/lcls/ds'
FFB_BASE = '/cds/data/drpsrcf'
PSDM_BASE = os.environ.get('SIT_PSDM_DATA', '/reg/d/psdm')
SD_EXT = '/hdf5/smalldata'
logger.debug(f"PSDM_BASE={PSDM_BASE}")

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run', type=str, default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--stn', help='hutch station', type=int, default=0)
parser.add_argument('--nevents', help='number of events', type=int, default=1e9)
parser.add_argument('--directory', help='directory for output files (def <exp>/hdf5/smalldata)')
parser.add_argument('--gather_interval', help='gather interval', type=int, default=25)
parser.add_argument('--norecorder', help='ignore recorder streams', action='store_true', default=False)
parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
parser.add_argument('--epicsAll', help='store all epics PVs', action='store_true', default=False)
parser.add_argument('--full', help='store all data (please think before using this)', action='store_true', default=False)
parser.add_argument('--fullSum', help='store sums for all area detectors', action='store_true', default=False)
parser.add_argument('--default', help='store only minimal data', action='store_true', default=False)
parser.add_argument('--image', help='save everything as image (use with care)', action='store_true', default=False)
parser.add_argument('--tiff', help='save all images also as single tiff (use with even more care)', action='store_true', default=False)
parser.add_argument('--centerpix', help='do not mask center pixels for epix10k detectors.', action='store_true', default=False)
parser.add_argument("--postRuntable", help="postTrigger for seconday jobs", action='store_true', default=False)
parser.add_argument("--wait", help="wait for a file to appear", action='store_true', default=False)
parser.add_argument("--xtcav", help="add xtcav processing", action='store_true', default=False)
parser.add_argument("--noarch", help="dont use archiver data", action='store_true', default=False)
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
if hostname.find('sdf')>=0:
    logger.debug('On S3DF')


elif hostname.find('drp')>=0:
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

# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.directory, exp, hutch)
#if args.default:
#    if useFFB:
#        h5_f_name = h5_f_name.replace('hdf5','hdf5_def')
#    else:
#        h5_f_name = h5_f_name.replace('hdf5','scratch')

# Define data source name and generate data source object, don't understand all conditions yet
ds_name = ''.join(['exp=', exp, ':run=', run, ':smd'])
if args.norecorder:
        ds_name += ':stream=0-79'
if useFFB:
        ds_name += ':live:dir=/cds/data/drpsrcf/%s/%s/xtc'%(exp[0:3],exp)
# try this: live & all streams (once I fixed the recording issue)
#ds_name = ''.join(['exp=', exp, ':run=', run, ':smd:live'])
try:
    ds = psana.MPIDataSource(ds_name)
except Exception as e:
    logger.debug('Could not instantiate MPIDataSource with {0}: {1}'.format(ds_name, e))
    sys.exit()

# Generate smalldata object
small_data = ds.small_data(h5_f_name, gather_interval=args.gather_interval)

# Not sure why, but here
if ds.rank == 0:
	logger.debug('psana conda environment is {0}'.format(os.environ['CONDA_DEFAULT_ENV']))

########################################################## 
##
## Setting up the default detectors
##
########################################################## 

start_setup_dets = time.time()

default_dets = defaultDetectors(hutch.lower())
if args.xtcav and not args.norecorder:
    #default_dets.append(xtcavDetector('xtcav','xtcav',method='COM'))
    default_dets.append(xtcavDetector('xtcav','xtcav'))

#
# add stuff here to save all EPICS PVs.
#
if args.full or args.epicsAll:
    logger.debug('epicsStore names', ds.env().epicsStore().pvNames())
    if args.experiment.find('dia')>=0:
        epicsPV=ds.env().epicsStore().pvNames()
    else:
        epicsPV=ds.env().epicsStore().aliases()
    if len(epicsPV)>0:
        logger.debug('adding all epicsPVs....')
        default_dets.append(epicsDetector(PVlist=epicsPV, name='epicsAll'))
elif len(epicsPV)>0:
    default_dets.append(epicsDetector(PVlist=epicsPV, name='epicsUser'))

if len(epicsOncePV)>0:
    EODet = epicsDetector(PVlist=epicsOncePV, name='epicsOnce')
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
                thisDet.storeSum(sumAlgo='calib_dropped_square')

            dets.append(thisDet)
        except:
           pass


# save detector config data
userDataCfg={}
for det in default_dets:
    if det.name=='tt' and len(ttCalib)>0:
        det.setPars(ttCalib)
        logger.info(f'Using user-defined tt parameters: {ttCalib}')
    userDataCfg[det.name] = det.params_as_dict()
for det in dets:
    try:
        userDataCfg[det._name] = det.params_as_dict()
    except:
        userDataCfg[det.name] = det.params_as_dict()

#is EOODet exists, save this later.
if EODet is None:
    Config={'UserDataCfg':userDataCfg}
    small_data.save(Config)

end_setup_dets = time.time()

if args.tiff:
    dirname = '/cds/data/psdm/%s/%s/scratch/run%d'%(args.experiment[:3],args.experiment,int(args.run))
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
    if ( (evt_num<100 and evt_num%10==0) or (evt_num<1000 and evt_num%100==0) or (evt_num%1000==0)):
        if (int(os.environ.get('RUN_NUM', '-1')) > 0):
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
        sumData=small_data.sum(det.storeSum()[key])
        sumDict['Sums']['%s_%s'%(det._name, key)]=sumData
if len(sumDict['Sums'].keys())>0:
#     print(sumDict)
    small_data.save(sumDict)

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
small_data.save()
if (int(os.environ.get('RUN_NUM', '-1')) > 0):
    if ds.size > 1:
        if ds.rank == 0:
            requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": "~ %d cores * %d evts"%(ds.size,evt_num)},{"key": "<b>Duration</b>", "value": "%f min"%(prod_time)}])
    else:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": evt_num}])
logger.debug('Saved all small data')

if args.postRuntable and ds.rank==0:
    print('Posting to the run tables.')
    locStr=''
    if useFFB:
        locStr='_ffb'
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
