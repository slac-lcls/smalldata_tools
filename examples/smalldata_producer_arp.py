#!/usr/bin/env python

import numpy as np
import psana
import time
import argparse
import socket
import os
import logging 
import requests
import sys
from glob import glob

# General Workflow
# This is meant for arp which means we will always have an exp and run
# Check if this is a current experiment
# If it is current, check in ffb for xtc data, if not there, default to psdm

# TODO: Fix this
#fpath=os.path.dirname(os.path.abspath(__file__))
#fpathup = '/'.join(fpath.split('/')[:-1])
#sys.path.append(fpathup)
#print(fpathup)
#this did not work for me, but test again...
#sys.path.append('../smalldata_tools')
sys.path.append('/reg/g/psdm/sw/tools/smalldata_tools')
from smalldata_tools.utilities import printMsg
from smalldata_tools.SmallDataUtils import setParameter, defaultDetectors, detData
from smalldata_tools.SmallDataUtils import getUserData, getUserEnvData
from smalldata_tools.SmallDataDefaultDetector import ttRawDetector, wave8Detector, epicsDetector

from smalldata_tools.SmallDataDefaultDetector import epicsDetector, eorbitsDetector
from smalldata_tools.SmallDataDefaultDetector import bmmonDetector, ipmDetector
from smalldata_tools.SmallDataDefaultDetector import encoderDetector
from smalldata_tools.roi_rebin import ROIFunc, imageFunc
from smalldata_tools.DetObject import DetObject

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.debug('Directory is {0}'.format(os.getcwd()))

# Constants
HUTCHES = [
	'AMO',
	'SXR',
	'XPP',
	'XCS',
	'MFX',
	'CXI',
	'MEC'
]

WS_URL = 'https://pswww.slac.stanford.edu/ws/lgbk'
ACTIVE_EXP_EXT = '/lgbk/ws/activeexperiment_for_instrument_station'
FFB_BASE = '/reg/d/ffb'
PSDM_BASE = '/reg/d/psdm'
SD_EXT = '/hdf5/smalldata'
WS_CUR_RUN = '/ws/current_run'

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run', type=str, default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--exp', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--stn', help='hutch station', type=int, default=0)
parser.add_argument('--nevents', help='number of events', type=int, default=1e9)
parser.add_argument('--dir', help='directory for output files (def <exp>/hdf5/smalldata)')
parser.add_argument('--offline', help='run offline (def for current exp from ffb)', action='store_true')
parser.add_argument('--gather', help='gather interval', type=int, default=100)
parser.add_argument('--epicsAll', help='save all EPICS PVs', action='store_true')
parser.add_argument('--full', help='save everything (use with care)', action='store_true')
parser.add_argument("--norecorder", help="ignore recorder streams", action='store_true')
parser.add_argument('--image', help='save everything as image (use with care)', action='store_true')
args = parser.parse_args()
logger.debug('Args to be used for small data run: {0}'.format(args))

###### Helper Functions ##########

def get_cur_exp(hutch, station):
	"""Get the active experiment for the given hutch, returns '' if no active
	experimets for given hutch
	"""
	endpoint = ''.join([WS_URL, ACTIVE_EXP_EXT])
	args = {'instrument_name': hutch, 'station': station}
	r = requests.get(endpoint, args)
	active_exp = str(r.json().get('value', {'name':''}).get('name'))
	
	return active_exp

def get_cur_run(exp):
	"""Get current run for current experiment"""
	endpoint = ''.join([WS_URL, '/lgbk/', exp, WS_CUR_RUN])
	run_info = requests.get(endpoint).json()['value']
	
	return str(int(run_info['num']))

def get_xtc_files(base, hutch, run):
	"""File all xtc files for given experiment and run"""
	run_format = ''.join(['r', run.zfill(4)])
	data_dir = ''.join([base, '/', hutch.lower(), '/', exp, '/xtc'])
	xtc_files = glob(''.join([data_dir, '/', '*', '-', run_format, '*']))

	return xtc_files

def get_sd_file(write_dir, exp, hutch):
	"""Generate directory to write to, create file name"""
	if write_dir is None:
		write_dir = ''.join([PSDM_BASE, '/', hutch, '/', exp, SD_EXT])
	h5_f_name = ''.join([write_dir, '/', exp, '_Run', run.zfill(4), '.h5'])
	logger.debug('Will write small data file to {0}'.format(h5_f_name))

	if not os.path.isdir(write_dir):
		logger.debug('{0} does not exist, creating directory'.format(write_dir))
		try:
			os.mkdir(write_dir)
		except OSError as e:
			logger.debug('Unable to make directory {0} for output, exiting: {1}'.format(write_dir, e))
			sys.exit()

	return h5_f_name

##### START SCRIPT ########

# Define hostname
hostname = socket.gethostname()

# Parse hutch name from experiment and check it's a valid hutch
exp = args.exp
run = args.run
station = args.stn
logger.debug('Analyzing data for EXP:{0} - RUN:{1}'.format(args.exp, args.run))

hutch = exp[:3].upper()
if hutch not in HUTCHES:
	logger.debug('Could not find {0} in list of available hutches'.format(hutch))
	sys.exit()	

# Get current exp and run
cur_exp = get_cur_exp(hutch, station)

xtc_files = []
# If experiment matches, check for files in ffb
ds_name = ''.join(['exp=', exp, ':run=', run, ':smd'])
if cur_exp == exp:
        xtc_files = get_xtc_files(FFB_BASE, hutch, run)

if xtc_files and not args.offline:
        ds_name = ds_name.replace(':smd',':dir=/reg/d/ffb/%s/%s/xtc:live:smd'%(exp[:3],exp))
# If not a current experiment or files in ffb, look in psdm
else:
	logger.debug('Either not a current exp or files not in ffb, looking in psdm')
	xtc_files = get_xtc_files(PSDM_BASE, hutch, run)
        ds_name = ''.join(['exp=', exp, ':run=', run, ':smd', ':live'])

if not xtc_files:
	logger.debug('XTC files not available')
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>XTC files:</b>", "value": 'not available'}])
        sys.exit()

# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.dir, exp, hutch)

# Define data source name and generate data source object, don't understand all conditions yet
if args.norecorder:
        ds_name += 'stream=0-79'

try:
	ds = psana.MPIDataSource(ds_name)
except Exception as e:
	logger.debug('Could not instantiate MPIDataSource with {0}: {1}'.format(ds_name, e))
	sys.exit()

# Generate smalldata object
small_data = ds.small_data(h5_f_name, gather_interval=args.gather)

# Not sure why, but here
if ds.rank is 0:
	logger.debug('psana conda environment is {0}'.format(os.environ['CONDA_DEFAULT_ENV']))

# gather default dets and add to data
default_dets = defaultDetectors(hutch.lower())
config_data = {det.name: det.params_as_dict() for det in default_dets}
small_data.save({'UserDataCfg': config_data})

#
# add stuff here to save all EPICS PVs.
#
logger.debug('epicsStore names', ds.env().epicsStore().pvNames())
epicsPVall=[]
if args.exp.find('dia')>=0:
    epicsPVall=ds.env().epicsStore().pvNames()
    default_dets.append(eorbitsDetector())
elif args.epicsAll:
    epicsPVall=ds.env().epicsStore().aliases()
elif args.full:
    epicsPVall=ds.env().epicsStore().aliases()
if len(epicsPVall)>0:
    logger.debug('adding all epicsPVs....total of %d PVs'%(len(epicsPVall)))
    default_dets.append(epicsDetector(PVlist=epicsPVall,name='epicsAll'))
    logger.debug('default dets: ',[det.name for det in default_dets])

dets=[]
if args.full:
        default_det_aliases = [det.name for det in default_dets]

        aliases = []
        for dn in psana.DetNames():
            if dn[1]!='':
                aliases.append(dn[1])
            else:
                aliases.append(dn[0])

        for alias in aliases:
            if alias in default_det_aliases: continue
            if alias=='FEEGasDetEnergy': continue #done by mpidatasource
            if alias=='PhaseCavity':     continue #done by mpidatasource
            if alias.find('evr')>=0:     continue #done by mpidatasource
            if alias=='ControlData':     continue #done by my code
            if alias.find('BMMON')>=0:
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
            try:
                thisDet = DetObject(alias, ds.env(), int(run), name=alias)
                hasGeom=False
                for keyword in ['cs','Cs','epix','Epix','jungfrau','Jungfrau']:
                        if alias.find(keyword)>=0 and args.image: hasGeom=True
                if hasGeom:
                        fullROI = ROIFunc()
                        fullROI.addFunc(imageFunc(coords=['x','y']))
                else:
                        fullROI = ROIFunc(writeArea=True)
                thisDet.addFunc(fullROI)
                dets.append(thisDet)
            except:
                pass


max_iter = args.nevents / ds.size
evt_num = -1
for evt_num, evt in enumerate(ds.events()):
	if evt_num >= max_iter:
		break

	det_data = detData(default_dets, evt)
	small_data.event(det_data)
        if (evt_num==0): 
            print('default det data: ',[k for k in det_data])
        if ((evt_num<100&evt_num%10==0) or (evt_num<1000&evt_num%100==0) or (evt_num%1000==0)):
            try:
                requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Current Event</b>", "value": evt_num}])
            except:
                pass

        #detector data using DetObject 
        if len(dets)>0:
            userDict = {}
            for det in dets:
                try:
                    det.getData(evt)
                    det.processFuncs()
                    userDict[det._name]=0
                    userDict[det._name]=getUserData(det)
                    #print(userDict[det._name].keys())
                    try:
                        envData=getUserEnvData(det)
                        if len(envData.keys())>0:
                            userDict[det._name+'_env']=envData
                    except:
                        pass
                except:
                    pass
            small_data.event(userDict)
        if ((evt_num<100&evt_num%10==0) or (evt_num<1000&evt_num%100==0) or (evt_num%1000==0)):
            requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Current Event</b>", "value": evt_num}])

logger.debug('rank {0} on {1} is finished'.format(ds.rank, hostname))
small_data.save()
try:
    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Last Event</b>", "value": evt_num}])
except:
    pass
logger.debug('Saved all small data')

