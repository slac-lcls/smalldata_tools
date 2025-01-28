import sys
import argparse
import socket
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run")
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--cube", help="cube name", default='sort')
args = parser.parse_args()

if not args.run:
    run=raw_input('provide a run number of experiment %s:'%expname)
else:
    run=args.run

runList=[]
runs = run.split(',')
for trun in runs:
    if trun.find('-')==-1:
        runList.append(int(trun))
    else:
        for ttrun in range(int(trun.split('-')[0]), int(trun.split('-')[1])+1):
            runList.append(ttrun)

hutches=['amo','sxr','xpp','xcs','mfx','cxi','mec']
if not args.exp:
    hostname=socket.gethostname()
    foundHutch=False
    for ihutch in hutches:
        if hostname.find(ihutch)>=0:
            hutch=ihutch.upper()
            foundHutch=True
            break
    if not foundHutch and hostname.find('psusr')>=0:
        if hostname.find('psusr11')>=0:
            hutch='AMO'
        elif hostname.find('psusr12')>=0:
            hutch='SXR'
        elif hostname.find('psusr13')>=0:
            hutch='XPP'
        elif hostname.find('psusr21')>=0:
            hutch='XCS'
        elif hostname.find('psusr22')>=0:
            hutch='CXI'
        elif hostname.find('psusr23')>=0:
            hutch='MEC'
        elif hostname.find('psusr24')>=0:
            hutch='MFX'
        if hutch!='':
            foundHutch=True
    else:
        #then check current path
        path=os.getcwd()
        for ihutch in hutches:
            if path.find(ihutch)>=0:
                hutch=ihutch.upper()
                foundHutch=True
                break
        if not foundHutch:
            print('I cannot figure out which hutch we are in, so cannot determine the current experiment')
            sys.exit()

    try:
        import logging
        import requests
        ws_url = "https://pswww.slac.stanford.edu/ws/lgbk"
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        if hutch == 'CXI':
            print('Will assume the first CXI station, if this is wrong, please  -e <expname> on commandline')
        resp = requests.get(ws_url + "/lgbk/ws/activeexperiment_for_instrument_station", {"instrument_name": hutch, "station": 0})
        expname = resp.json().get("value", {}).get("name")
    except:
        print('could not determine experiment name, will quit')
        sys.exit()
else:
    expname=args.exp

#from smalldata_tools import SmallDataAna
sys.path.append('./smalldata_tools')
from CubeAna import CubeAna
t0 = time.time()

cube = CubeAna(expname, runList[0], cubeName=args.cube, dirname='./output', debug=True)
#cube.addAzInt(detname='cspad', center=[188752.179175, 104154.2268549],dis_to_sam=146., eBeam=9.5)
for thisrun in runList:
    if thisrun!=runList[0]:
        cube.addDataFromRun(thisrun)
#cube.applyAzav()
cube.makeImg()
cube.cubeSumToHdf5()
    
t1 = time.time()
print('this took %g seconds '%(t1-t0))
