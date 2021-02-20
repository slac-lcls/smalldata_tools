#!/usr/bin/env python

###
# run this like:
# ./examples/MakeCube.py --run 302 --exp xppo5616
#
# submit to queue like:
# before submitting to the queue, you'll need to unset the DISPLAY variable. 
# bsub -q psanaq -n <njobs> -o <path_to_logfile> python ./examples/MakeCube.py --run 302 --exp xppo5616
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
from requests.auth import HTTPBasicAuth

########################################################## 
##
## User Input start --> 
##
########################################################## 
##########################################################
# functions for run dependant parameters
##########################################################
def binBoundaries(run):
    if isinstance(run,basestring):
        run=int(run)
    if run == 6:
        return np.arange(-5.,5.,0.1)
    return None
########################################################## 
##
## <-- User Input end
##
########################################################## 



fpath=os.path.dirname(os.path.abspath(__file__))
fpathup = '/'.join(fpath.split('/')[:-1])
sys.path.append(fpathup)
print(fpathup)

from smalldata_tools.SmallDataAna import SmallDataAna
from smalldata_tools.SmallDataAna_psana import SmallDataAna_psana

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

    #CHANGE ME!
    ####
    #set up different event filters.    
    ####
    ana.addCut('lightStatus/xray',0.5,1.5,'filter1')

    ana.addCut('lightStatus/xray',0.5,1.5,'filter2')
    ana.addCut('evr/code_41',0.5,1.5,'filter2')


    #CHANGE ME 
    ####
    #(if you want a full detector or many)
    ####
    #save image data
    #detDict = {'source':'jungfrau1M','full':1, 'image':1, 'common_mode':7}
    detDict = {'source':'jungfrau1M','full':1, 'image':1}
    #save photon images -- not used recently.
    #detDict = {'source':'jungfrau512k','photon_0p85':3.0, 'image':1}
    #zylaDict = {'source':'zyla','full':1, 'common_mode':0}

    #CHANGE ME!
    ####
    # list of variables to bin
    ####
    #varList = ['ipm5/sum','snd_dio/t4d','snd_dio/dco','scan/_delay', zylaDict, detDict]
    varList = ['ipm2/sum','ipm3/sum','diodeU/channels', detDict]
    #varList = ['ipm2/sum','ipm3/sum','diodeU/channels']

    #CHANGE ME
    ####
    # if you are interested in laser-xray delay, please select the delay of choice here!
    ####
    ana.setDelay(use_ttCorr=True, addEnc=False, addLxt=False, reset=True)

    cubeName='cube' #initial name
    scanName, scanValues = ana.getScanValues()
    binSteps=[]
    binName=''
    filterName='filter1'
    if scanName!='':
        if scanName.find('delay')<0:
            scanSteps = np.unique(scanValues)
            scanSteps = np.append(scanSteps, scanSteps[-1]+abs(scanSteps[1]-scanSteps[0]))#catch value at right edge?
            scanSteps = np.append(scanSteps[0]-abs(scanSteps[1]-scanSteps[0]),scanSteps) #catch values at left edge
            binSteps = scanSteps
            cubeName = scanName
            if scanName == 'lxt':    
                print('bin data using ',scanName,' and bins: ',scanSteps)
                binName='delay'
            else:
                print('bin data using ',scanName,' and bins: ',scanSteps)
                binName='scan/%s'%scanName
        else:
            binSteps = binBoundaries(run)
            if binSteps is None:
                #assume a fast delay scan here.
                enc = ana.getVar('enc/lasDelay')
                binSteps = np.arange( (int(np.nanmin(enc*10)))/10., 
                                      (int(np.nanmax(enc*10)))/10., 0.1)
            binName = 'delay'
    else:
        cubeName='randomTest'
        binName='ipm2/sum'
        binVar=ana.getVar(binName)
        binSteps=np.percentile(binVar,[0,25,50,75,100])

    if args.nevents:
        cubeName+='_%sEvents'%args.nevents

    ana.addCube(cubeName,binName,binSteps,filterName)
    ana.addToCube(cubeName,varList)

    ####
    #there are two ways to get multiple cubes: if you want a standard laser on/off set, use the onoff flag
    #if you e.g. have data with an external magnet switched on/off, use the 'Sel2' approach.
    #both are not necessary.
    ####
    ana.addCube(cubeName+'Sel2',binName,binSteps,'filter2')
    ana.addToCube(cubeName+'Sel2',varList)
    nSel = ana.getFilter(filterName).sum()
    nSel2 = ana.getFilter('filter2').sum()

    if (int(os.environ.get('RUN_NUM', '-1')) > 0):
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Cube: </b>", "value": "Prepared"}])

    if args.nevents:
        print('make cube with fewer events/bin')
        if nSel>0:
            anaps.makeCubeData(cubeName, nEvtsPerBin=int(args.nevents),dirname=args.outdirectory)
        #if nSel2>0:
       #     anaps.makeCubeData(cubeName+'Sel2',  nEvtsPerBin=int(args.nevents))
    else:
        if nSel>0:
            anaps.makeCubeData(cubeName, dirname=args.outdirectory)
    #        #anaps.makeCubeData(cubeName, onoff=1) #request 'on' events base on input filter (add optical laser filter)
    #        #anaps.makeCubeData(cubeName, onoff=0) #request 'off' events base on input filter (switch optical laser filter, drop tt requirements)
    #    if nSel2>0:
    #        anaps.makeCubeData(cubeName+'Sel2')

end_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')

if (int(os.environ.get('RUN_NUM', '-1')) > 0):
    requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Cube: </b>", "value": "Done"}])

if args.postRuntable:
    print('posting to the run tables.')
    runtable_data = {"Prod_cube_end":end_prod_time,
                     "Prod_cube_start":begin_prod_time,
                     "Prod_cube_ncores":ds.size}
    time.sleep(5)
    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(args.experiment)
    print('URL:',ws_url)
    #krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    #r = requests.post(ws_url, headers=krbheaders, params={"run_num": args.run}, json=runtable_data)
    user=args.experiment[:3]+'opr'
    with open('/cds/home/opr/%s/forElogPost.txt'%user,'r') as reader:
        answer = reader.readline()
    r = requests.post(ws_url, params={"run_num": args.run}, json=runtable_data, auth=HTTPBasicAuth(args.experiment[:3]+'opr', answer[:-1]))
    print(r)
