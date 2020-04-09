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

################
def get_i0_filter(run):
    if int(run)<10:
        return -10.,100
    elif int(run)<100:
        return 10.,1e6
    else:
        return None,None
##################

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run",type=int)
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--dir", help="output directory")
parser.add_argument("--nev", help="number of events/bin", type=int)
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
        import RegDB.experiment_info
        expname=RegDB.experiment_info.active_experiment(hutch)[1]
    except:
        print 'could not determine experiment name, will quit'
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
print 'try to import smalldata'
from SmallDataAna import SmallDataAna
print 'imported smalldata'
from SmallDataAna_psana import SmallDataAna_psana
print 'imported smalldata_psana'

ana = None
anaps = SmallDataAna_psana(expname,run,dirname,fname)
if anaps and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
    print 'create ana module from anaps'
    ana = anaps.sda
else:
    print 'we will now try to open the littleData file directly'
    ana = SmallDataAna(expname,run, dirname, fname)
    if 'fh5' not in ana.__dict__.keys():
        ana = None

if ana is not None:
    ana.printRunInfo()

    #somewhat non-obvious: if lower==upper; REJECT this point (e.g. FLTPOS==0)
    ana.addCut('lightStatus/xray',0.5,1.5,'Sel')
    #--> i0 threshold here
    cutLow,cutHigh =  get_i0_filter(run)
    if cutLow is not None:
        ana.addCut('ipm3/sum',cutLow,cutHigh,'Sel')
    #<-- i0 threshold here

    #save raw data
    #cspadDict = {'source':'cspad','full':1, 'thresAdu':10.}
    #cspadDict = {'source':'cspad','full':1, 'thresAdu':10., 'image':1}
    #cspadDict = {'source':'cspad','full':1}
    #see if this actually works.....
    cspadDict = {'source':'cspad','full':1, 'photon_0p9':42.5}

    ls = ana.getVar('lightStatus/laser')
    postLaser = ana.getVar('evr/code_92')
    ana.addVar('eventBucket',ls+postLaser*2)
    ana.addCube('snelsonTest','eventBucket',[-0.5,2.5,1.],'Sel')
    ##if we want to bin each location separately.
    #scanStepVar = ana.getVar('scan/varStep')
    #scanSteps = np.unique(scanStepVar)
    #scanSteps = np.append(scanSteps, max(scanSteps+1))
    #ana.addCube('snelsonTest','eventBucket',[-0.5,2.5,1.],'Sel', addBinVars={'scan/varStep':scanSteps})
    ana.addToCube('snelsonTest',['ipm2/sum','ipm3/sum','diodeU/channels', 'cspad/azav', 'cspad/ROI_0_sum',  'cspad/ROI_0_com', 'UserDataCfg/cspad/azav_center',cspadDict])
    #ana.addToCube('snelsonTest',['ipm2/sum','ipm3/sum','diodeU/channels', 'cspad/azav', 'cspad/ROI_0_sum'])

    if args.nev:
        anaps.makeCubeData('snelsonTest',  nEvtsPerBin=args.nev)
    else:
        anaps.makeCubeData('snelsonTest')
