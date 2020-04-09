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
parser.add_argument("--cube", help="cube name")
parser.add_argument("--onoff", help="cube name", action='store_true')
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
    ana.addCut('lightStatus/xray',0.5,1.5,'xon')
    ana.addCut('lightStatus/xray',0.5,1.5,'on')
    ana.addCut('lightStatus/laser',0.5,1.5,'on')
    ana.addCut('lightStatus/xray',0.5,1.5,'off')
    ana.addCut('lightStatus/laser',-0.5,0.5,'off')
    
    ana.addCut('lightStatus/xray',0.5,1.5,'onSel')
    ana.addCut('lightStatus/laser',0.5,1.5,'onSel')
    #ana.addCut('ipm3/sum',100,1.e6,'onSel')
    ana.addCut('ipm3/sum',0,1.e6,'onSel')
    #ana.addCut('tt/AMPL',0.01,10.,'onSel')
    #ana.addCut('tt/FLTPOSFWHM',80.,130.,'onSel')

    ana.setDelay(use_ttCorr=True, addEnc=False, addLxt=False, reset=True)

    #save raw data
    detDict = {'source':'jungfrau512k','full':1, 'image':1, 'common_mode':0}#maybe 7?

    #varList = ['ipm2/sum','ipm3/sum','diodeU/channels']
    varList = ['ipm2/sum','ipm3/sum','diodeU/channels', 'jungfrau512k/ROIsum_0_proj_sum','jungfrau512k/ROIsum_1_proj_sum','jungfrau512k/ROI_0_projax0_sum','jungfrau512k/ROI_0_projax1_sum', detDict]

    cubeName='cube' #initial name
    scanName, scanValues = ana.getScanValues()
    binSteps=[]
    binName=''
    filterName='onSel'
    stepSize=0.1
    if scanName!='':# and scanName.find('delay')<0:
        if scanName.find('lxt')>=0:
            binName='delay'
            cubeName='delay'
            binSteps = np.arange(scanValues.min(), scanValues.max(), stepSize)
            print 'bin data using ',scanName,' and bins: ',scanSteps
        elif scanName.find('delay')>=0:
            lasDelay = ana.getVar('enc/lasDelay')
            binName='delay'
            cubeName='delay'
            binSteps = np.arange(int(np.nanmin(lasDelay)*10.)/10., 
                                 int(np.nanmax(lasDelay)*10.)/10., stepSize)
            #print 'delay...'
            #sys.exit()
        else:
            scanSteps = np.unique(scanValues)
            scanSteps = np.append(scanSteps, scanSteps[-1]+abs(scanSteps[1]-scanSteps[0]))#catch value at right edge?
            scanSteps = np.append(scanSteps[0]-abs(scanSteps[1]-scanSteps[0]),scanSteps) #catch values at left edge
            binSteps = scanSteps
            cubeName = scanName
            print 'bin data using ',scanName,' and bins: ',scanSteps
            binName='scan/%s'%scanName

    else:
        print 'this run is not a scan, skip making a cube!'
        sys.exit()

    if args.nev:
        cubeName+='_%sEvents'%args.nev
    ana.addCube(cubeName,binName,binSteps,filterName)
    ana.addToCube(cubeName,varList)

    if args.cube:
        cubeName = args.cube

    #anaps.makeCubeData('cube',  nEvtsPerBin=3, dirname='./output')
    if args.nev:
        if args.onoff:
            anaps.makeCubeData(cubeName,  nEvtsPerBin=args.nev, onoff=0)
            anaps.makeCubeData(cubeName,  nEvtsPerBin=args.nev, onoff=1)
        else:
            anaps.makeCubeData(cubeName,  nEvtsPerBin=args.nev)
    else:
        if args.onoff:
            anaps.makeCubeData(cubeName,onoff=0)
            anaps.makeCubeData(cubeName,onoff=1)
        else:
            anaps.makeCubeData(cubeName)
