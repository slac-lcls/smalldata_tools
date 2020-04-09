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
parser.add_argument("--dir", help="directory w/ smallData file if not default")
parser.add_argument("--nev", help="number of events/bin")
parser.add_argument("--file", help="file name")
parser.add_argument("--cube", help="cube Name")
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

#xcslr6316, run38: t1d from 4000 to 16000
#xcslr6316, run38: dcc from 1000 to 25000
if ana is not None:
    ana.printRunInfo()

    #set up different event selection.
    ana.addCut('lightStatus/xray',0.5,1.5,'on')
    ana.addCut('lightStatus/laser',0.5,1.5,'on')
    ana.addCut('lightStatus/xray',0.5,1.5,'off')
    ana.addCut('lightStatus/laser',-0.5,0.5,'off')
    
    ana.addCut('lightStatus/xray',0.5,1.5,'good')
    ana.addCut('lightStatus/laser',0.5,1.5,'good')
    ana.addCut('ipm5/sum',7.5, 10.,'good')
    #ana.addCut('ipm5/sum',0.0, 10.,'good')

    #save raw data
    #detDict = {'source':'jungfrau512k','full':1}
    #save image data
    detDict = {'source':'jungfrau512k','full':1,'image':1}
    #save photon images
    #detDict = {'source':'jungfrau512k','photon_0p85':3.0, 'image':1}

    zylaDict = {'source':'zyla_ladm','full':1}
    #epixDict = {'source':'epix_ladm_1','full':1}
    #epixDict = {'source':'epix_ladm_1','full':1, 'image':1}
    #snelson: this does not add the second plot I was hoping it would (and though I tested it did..
    epixDict = {'source':'epix_ladm_1','full':1,'photon_0p85':155.0, 'image':1
    #yes, this works.
    #epixDict = {'source':'epix_ladm_1','photon_0p85':155.0,'image':1}

    varList = ['ipm5/sum','snd_dio/t4d','snd_dio/dco','scan/_delay', epixDict]

    cubeName=''
    scanName, scanValues = ana.getScanValues()
    if scanName!='':# and scanName.find('delay')<0:
        scanSteps = np.unique(scanValues)
        scanSteps = np.append(scanSteps, scanSteps[-1]+(scanSteps[1]-scanSteps[0]))#catch value at right edge?
        scanSteps = np.append(scanSteps[0]-(scanSteps[1]-scanSteps[0]),scanSteps) #catch values at left edge
        cubeName=scanName
        ana.addCube(cubeName,'scan/%s'%scanName,scanSteps,'good')
        ana.addToCube(cubeName,varList)

    ##artificial cube to speed data reading up.
    ana.addCube('random','user/random',np.linspace(0.,1.,100),'good') #set event selection
    ana.addToCube('random',varList) #list variables to be binned.

    ana.addCube('DelayLineScan','scan/_delay',np.linspace(0,0.1,11),'good')
    ana.addToCube('DelayLineScan',varList)

    #snelson: this does not work yet.
    ana.addCube('DelayLineScan_2d','scan/_delay',np.linspace(0,0.1,11),'good',addBinVars={'user/random': np.linspace(0.,1.,5)})
    ana.addToCube('DelayLineScan_2d',varList)

    #scanVal = ana.getVar('scan/delay')
    #ana.addCube('delay','scan/delay',np.unique(scanVal),'good')
    #ana.addToCube('delay',['ipm5/sum','zyla_ladm/ROI0',zylaDict])

    cubeName = 'DelayLineScan'
    if args.cube:
        cubeName = args.cube
    if args.nev:
        print 'make cube with fewer events/bin'
        #anaps.makeCubeData('cubeDelay',  nEvtsPerBin=int(args.nev))
        anaps.makeCubeData(cubeName,  nEvtsPerBin=int(args.nev))
    else:
        #anaps.makeCubeData('cubeDelay')
        anaps.makeCubeData(cubeName)
