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
sys.path.append('/reg/g/psdm/sw/tools/smalldata_tools')
from smalldata_tools.utilities import printMsg
from smalldata_tools.SmallDataUtils import setParameter, defaultDetectors, detData
from smalldata_tools.SmallDataDefaultDetector import ttRawDetector, wave8Detector, epicsDetector

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

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
parser.add_argument('--run', help='run', type=str, default=os.environ['RUN_NUM'])
parser.add_argument('--exp', help='experiment name', type=str, default=os.environ['EXPERIMENT'])
parser.add_argument('--stn', help='hutch station', type=int, default=0)
parser.add_argument('--nevt', help='number of events', type=int, default=1e9)
parser.add_argument('--dir', help='directory for output files (def <exp>/hdf5/smalldata)')
parser.add_argument('--offline', help='run offline (def for current exp from ffb)')
parser.add_argument('--gather', help='gather interval', type=int, default=100)
parser.add_argument("--norecorder", help="ignore recorder streams", action='store_true')
args = parser.parse_args()

###### Helper Functions ##########

def get_cur_exp(hutch, station):
	"""Get the active experiment for the given hutch"""
	endpoint = ''.join([WS_URL, ACTIVE_EXP_EXT])
	args = {'instrument_name': hutch, 'station': station}
	r = requests.get(endpoint, args)
	active_exp = str(r.json().get('value', {}).get('name'))
	
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
logger.debug('Analayzing data for EXP:{0} - RUN:{1}'.format(args.exp, args.run))

hutch = exp[:3].upper()
if hutch not in HUTCHES:
    logger.debug('Could not find {0} in list of available hutches'.format(hutch))

# Get current exp and run
cur_exp = get_cur_exp(hutch, station)
cur_run = get_cur_run(cur_exp)  # This might be unecessary

xtc_files = []
# If experiment matches, check for files in ffb
if cur_exp == exp:
	xtc_files = get_xtc_files(FFB_BASE, hutch, cur_run)

# If not a current experiment or files in ffb, look in psdm
if not xtc_files:
	logger.debug('Either not a current exp or files not in ffb, looking in psdm')
	xtc_files = get_xtc_files(PSDM_BASE, hutch, cur_run)

# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.dir, exp, hutch)

# Define data source name and generate data source object, don't understand all conditions yet
ds_name = ''.join(['exp=', exp, ':run=', run, ':smd', ':stream=0-79'])
try:
	ds = psana.MPIDataSource(ds_name)
except Exception as e:
	logger.debug('Could not instantiate MPIDataSource with {0): {1}'.format(ds_name, e))
	sys.exit()

# Generate smalldata object
small_data = ds.small_data(h5_f_name, gather_interval=args.gather)

# Not sure why, but here
if ds.rank is 0:
	logger.debug('psana conda environemnt is {0}'.format(os.environ['CONDA_DEFAULT_ENV']))

# gather default dets and add to data
default_dets = defaultDetectors(hutch.lower())
config_data = {det.name: det.params_as_dict() for det in default_dets}
small_data.save({'UserDataCfg': config_data})

max_iter = args.nevt / ds.size
for evt_num, evt in enumerate(ds.events()):
	if evt_num > max_iter:
		break

	det_data = detData(default_dets, evt)
	small_data.event(det_data)

logger.debug('rank {0} on {1} is finished'.format(ds.rank, hostname))
small_data.save()
logger.debug('Saved all small data')

