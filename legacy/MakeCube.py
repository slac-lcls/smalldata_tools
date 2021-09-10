#!/usr/bin/env python

###
# run this like:
# ./examples/MakeCube.py --run 302 --exp xppo5616
#
###
import sys
import time
from datetime import datetime
import numpy as np
import argparse
import socket
import logging 
import os
import requests
import tables
from PIL import Image
from requests.auth import HTTPBasicAuth
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


fpath=os.path.dirname(os.path.abspath(__file__))
fpathup = '/'.join(fpath.split('/')[:-1])
sys.path.append(fpathup)
print(fpathup)

from smalldata_tools.SmallDataAna import SmallDataAna
from smalldata_tools.SmallDataAna_psana import SmallDataAna_psana

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
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

FFB_BASE = '/cds/data/drpsrcf'
PSDM_BASE = '/reg/d/psdm'
SD_EXT = '/hdf5/smalldata'

begin_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')

parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run', type=str, default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument("--indirectory", help="directory w/ smallData file if not default")
parser.add_argument("--outdirectory", help="directory w/ smallDatafor cube if not same as smallData", default='')
parser.add_argument("--nevents", help="number of events/bin")
parser.add_argument("--postRuntable", help="postTrigger for seconday jobs", action='store_true')
parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
args = parser.parse_args()

##### START SCRIPT ########

# Define hostname
hostname = socket.gethostname()

# Parse hutch name from experiment and check it's a valid hutch
exp = args.experiment
run = args.run
logger.debug('Analyzing data for EXP:{0} - RUN:{1}'.format(args.experiment, args.run))

hutch = exp[:3].upper()
if hutch not in HUTCHES:
	logger.debug('Could not find {0} in list of available hutches'.format(hutch))
	sys.exit()	


dirname=''
if args.indirectory:
    dirname = args.indirectory
    if dirname.find('/')<=0:
        dirname+='/'

ana = None
anaps = SmallDataAna_psana(exp,run,dirname)
if anaps and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
    print('create ana module from anaps')
    ana = anaps.sda
else:
    print('we will now try to open the littleData file directly')
    ana = SmallDataAna(exp,run, dirname, fname)
    if 'fh5' not in ana.__dict__.keys():
        ana = None

if ana is not None:
    ana.printRunInfo()

    filterName='PIPS'
    ####
    #set up different event filters.    
    ####
    PIPS = ana.getVar('pips-diode-air/channels')
    ana.addVar('PIPS',PIPS[:,0])
    ana.addCut('PIPS',1.5,2.0,filterName)
    nSel = ana.getFilter(filterName).sum()

    ####
    #(if you want a full detector or many)
    ####
    #save image data
    det1Dict = {'source':'Epix10kaQuad2','full':1, 'image':1, 'thresADU':1.}
    det2Dict = {'source':'epix100a_1_IXS','full':1, 'image':1, 'thresADU':15.}
    det3Dict = {'source':'epix100a_2_XRTS','full':1, 'image':1, 'thresADU':15.}
    #save photon images -- not used recently.
    #detDict = {'source':'jungfrau512k','photon_0p85':3.0, 'image':1}

    ####
    # list of variables to bin
    ####
    varList = ['PIPS', det1Dict, det2Dict, det3Dict]


    cubeName='cube' #initial name
    binName='evr/code_140'
    binSteps=[0.,1.1]

    if args.nevents:
        cubeName+='_%sEvents'%args.nevents

    ana.addCube(cubeName,binName,binSteps,filterName)
    ana.addToCube(cubeName,varList)


    if (int(os.environ.get('RUN_NUM', '-1')) > 0) and rank==0:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Cube </b>", "value": "Prepared"}])

    if args.nevents:
        print('make cube with fewer events/bin')
        if nSel>0:
            anaps.makeCubeData(cubeName, nEvtsPerBin=int(args.nevents),dirname=args.outdirectory)
    else:
        if nSel>0:
            anaps.makeCubeData(cubeName, dirname=args.outdirectory)
   
end_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')

#read hdf5 file and make tiff files....
if rank==0:
    tiffdir='/cds/data/drpsrcf/mec/meclu9418/scratch/run%s'%args.run
    cubeFName = 'Cube_%s_Run%04d_%s'%(args.experiment, int(run), cubeName )
    cubedirectory= args.outdirectory
    if cubedirectory=='': cubedirectory = '/cds/data/drpsrcf/mec/meclu9418/scratch/hdf5/smalldata'
    dat = tables.open_file('%s/%s.h5'%(cubedirectory, cubeFName )).root
    detnames = ['Epix10kaQuad2','epix100a_1_IXS','epix100a_2_XRTS']
    for detname in detnames:
        cubedata = np.squeeze(getattr(dat, detname))
        if cubedata.ndim==2:
            im = Image.fromarray(cubedata)
            tiff_file = '%s/Run_%s_%s_filter_%s.tiff'%(tiffdir, run, detname, filterName)
            im.save(tiff_file)

if (int(os.environ.get('RUN_NUM', '-1')) > 0) and rank==0:
    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Cube </b>", "value": "Done"}])

if args.postRuntable and rank==0:
    print('posting to the run tables.')
    runtable_data = {"Prod_cube_end":end_prod_time,
                     "Prod_cube_start":begin_prod_time,
                     "Prod_cube_ncores":ds.size}
    time.sleep(5)
    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(args.experiment)
    #print('URL:',ws_url)

    user=args.experiment[:3]+'opr'
    with open('/cds/home/opr/%s/forElogPost.txt'%user,'r') as reader:
        answer = reader.readline()
    r = requests.post(ws_url, params={"run_num": args.run}, json=runtable_data, auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
    print(r)
