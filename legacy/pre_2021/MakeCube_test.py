###
# run this like:
# python ./examples/MakeCube.py --run 302 --exp xppo5616
#
# submit to queue like:
# before submitting to the queue, you'll need to unset the DISPLAY variable. 
# bsub -q psanaq -n <njobs> -o <path_to_logfile> python ./examples/MakeCube.py --run 302 --exp xppo5616
###
import sys
import numpy as np
import argparse
import socket
import os

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run",type=int)
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--dir", help="output directory")
parser.add_argument("--nev", help="number of events/bin")
parser.add_argument("--file", help="file name")
args = parser.parse_args()

if not args.run:
    run=int(raw_input('provide a run number of experiment %s:'%expname))
else:
    run=args.run

if not args.exp:
    hostname=socket.gethostname()
    if hostname.find('xpp')>=0:
        hutch='XPP'
    elif hostname.find('xcs')>=0:
        hutch='XCS'
    elif hostname.find('mfx')>=0:
        hutch='MFX'
    else:
        #then check current path
        path=os.getcwd()
        if path.find('xcs')>=0:
            hutch='XCS'
        elif path.find('xpp')>=0:
            hutch='XPP'
        elif path.find('mfx')>=0:
            hutch='MFX'
        else:
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

sys.path.append('./smalldata_tools')
print('try to import smalldata')
from SmallDataAna import SmallDataAna
print('imported smalldata')
from SmallDataAna_psana import SmallDataAna_psana
print('imported smalldata_psana')

ana = None
anaps = SmallDataAna_psana(expname,run,dirname,fname)
if anaps and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
    print('create ana module from anaps')
    ana = anaps.sda
else:
    print('we will now try to open the littleData file directly')
    ana = SmallDataAna(expname,run, dirname, fname)
    if 'fh5' not in ana.__dict__.keys():
        ana = None

if ana is not None:
    ana.printRunInfo()

    #somewhat non-obvious: if lower==upper; REJECT this point (e.g. FLTPOS==0)
    ana.addCut('lightStatus/xray',0.5,1.5,'xon')
    ana.addCut('lightStatus/xray',0.5,1.5,'on')
    ana.addCut('lightStatus/laser',0.5,1.5,'on')
    ana.addCut('lightStatus/xray',0.5,1.5,'off')
    ana.addCut('lightStatus/laser',-0.5,0.5,'off')
    
    ana.addCut('lightStatus/xray',0.5,1.5,'onSel')
    ana.addCut('lightStatus/laser',0.5,1.5,'onSel')
    ana.addCut('ipm3/sum',0.25,10.,'onSel')
    ana.addCut('tt/AMPL',0.025,10.,'onSel')
    ana.addCut('tt/FLTPOSFWHM',80.,130.,'onSel')

    #save raw data
    #cspadDict = {'source':'cspad','full':1}
    #save image data
    cspadDict = {'source':'cspad','full':1,'image':1}
    ana.addCube('cubeDelay','delay',np.linspace(-12.8,-10.8,24),'on')
    #ana.addCube('cubeDelay','delay',np.linspace(-12.8,-10.8,4),'on')
    ana.addToCube('cubeDelay',['ipm2/sum','ipm3/sum','diodeU/channels', cspadDict])

    #ana.addCube('cube','scan/varStep',np.arange(-0.5,ana.xrData.scan__varStep.max()+0.5,1),'onSel')
    #ana.addToCube('cube',['ipm2/sum','ipm3/sum', cspadDict])

    #anaps.makeCubeData('cube',  nEvtsPerBin=3, dirname='./output')
    if args.nev:
        anaps.makeCubeData('cubeDelay',  nEvtsPerBin=args.nev)
    else:
        anaps.makeCubeData('cubeDelay')
