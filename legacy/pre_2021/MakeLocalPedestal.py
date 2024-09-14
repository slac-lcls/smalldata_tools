import sys
import numpy as np
import argparse
from matplotlib import pyplot as plt
import socket
import os
#from smalldata_tools.utilties import hist2d,getUserData,rebin,addToHdf5
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run",type=int)
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--dir", help="output directory")
parser.add_argument("--file", help="file name")
args = parser.parse_args()

if not args.run:
    run=int(raw_input('provide a run number of experiment %s:'%expname))
else:
    run=args.run

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

dirname=''
if args.dir:
    dirname = args.dir
    if dirname.find('/')<=0:
        dirname+='/'
fname=''
if args.file:
    fname = args.file

#from smalldata_tools import SmallDataAna
sys.path.append('./smalldata_tools')
from SmallDataAna import SmallDataAna
ana = None
anaps = None
try:
    from SmallDataAna_psana import SmallDataAna_psana
    #from smalldata_tools import SmallDataAna_psana
    anaps = SmallDataAna_psana(expname,run,dirname,fname)
except:
    print('failed to import & create SmallDataAna_psana object')
    pass

if anaps is not None:
    #anaps.makePedestal('epix10ka2m',i0Check=None,numEvts=1000)
    ana.addCut('lightStatus/xray',-0.5,0.5,'xoff')
    ana.addCut('gas_detector/f_22_ENRC',-0.5,0.15,'xoff_gdet')
    anaps.makePedestal('epix',i0Check='ipm',dirname='/reg/d/psdm/xpp/xpplp7515/results/detector_monitoring_new/', useFilter='xoff_gdet', numEvts=1000)
    
