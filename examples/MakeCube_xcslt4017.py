###
# run this like:
# python ./examples/MakeCube.py --run 302 --exp xppo5616
#
# submit to queue like:
# before submitting to the queue, you'll need to unset the DISPLAY variable. 
# bsub -q psanaq -n <njobs> -o <path_to_logfile> python ./examples/MakeCube.py --run 302 --exp xppo5616
#
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

if ana is not None:
    ana.printRunInfo()

    #set up different event selection.    
    ana.addCut('lightStatus/xray',0.5,1.5,'xon')

    ana.addCut('lightStatus/xray',0.5,1.5,'good')
    ana.addCut('lightStatus/laser',0.5,1.5,'good')
    ana.addCut('ipm5/sum',1000.0, 1e6,'good')

    ana.addCut('lightStatus/xray',0.5,1.5,'off')
    ana.addCut('lightStatus/laser',-0.5,0.5,'off')
    ana.addCut('ipm5/sum',1000.0, 1e6,'off')

    ana.setDelay(use_ttCorr=True, addEnc=False, addLxt=False, reset=True)

    #save raw data
    detDict = {'source':'epix10ka2m','full':1,'common_mode':80}
    #save image data
    #detDict = {'source':'jungfrau512k','full':1, 'image':1, 'common_mode':7}
    #save photon images
    #detDict = {'source':'jungfrau512k','photon_0p85':3.0, 'image':1}

    varList = ['ipm5/sum','epix10ka2m/azav', detDict]

    cubeName='scan'
    scanName, scanValues = ana.getScanValues()
    stepSize=0.5
    if scanName!='':# and scanName.find('delay')<0:
        if scanName == 'lxt':    
            cubeName = 'lxt'
            scanSteps = np.unique(scanValues)
            scanSteps = np.append(scanSteps, scanSteps[-1]+abs(scanSteps[1]-scanSteps[0]))#catch value at right edge?
            scanSteps = np.append(scanSteps[0]-abs(scanSteps[1]-scanSteps[0]),scanSteps) #catch values at left edge
            #cubeName=scanName
            ana.setDelay(use_ttCorr=False, addEnc=False, addLxt=False)
            print 'bin data using ',scanName,' and bins: ',scanSteps
            ana.addCube(cubeName,'delay',scanSteps,'good')
            ana.addToCube(cubeName,varList)

            ana.addCube(cubeName+'Off','delay',scanSteps,'off')
            ana.addToCube(cubeName+'Off',varList)
            
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
            #cubeName=scanName
            print 'bin data using ',scanName,' and bins: ',scanSteps
            ana.addCube(cubeName,'scan/%s'%scanName,scanSteps,'good')
            ana.addToCube(cubeName,varList)

            ana.addCube(cubeName+'Off','scan/%s'%scanName,scanSteps,'off')
            ana.addToCube(cubeName+'Off',varList)

    else:
        cubeName='delay'
        
        ##should be fine to use all events as pre-split&delay...
        #binVarName='gas_detector/f_22_ENRC' #no matter which variable I use, this fails for some runs...
        #binBounds=np.percentile(binVar,[0,25,50,75,100])
        binVarName='delay' #always there. Need to make two cubes. Single bin
        binBounds=[-0.5,1.5]
        binVar=ana.getVar(binVarName)

        ana.addCube(cubeName,binVarName,binBounds,'good') #set event selection
        ana.addToCube(cubeName,varList) #list variables to be binned.
        
        ana.addCube(cubeName+'Off',binVarName,binBounds,'good')#set event selection
        ana.addToCube(cubeName+'Off',varList) #list variables to be binned.
        #ana.printCubes('random')

    if args.cube:
        cubeName = args.cube

    nGood = ana.getFilter('good').sum()
    nOff = ana.getFilter('off').sum()
    if run >= 146:# and run <=154:
        if args.nev:
            anaps.makeCubeData(cubeName,  nEvtsPerBin=int(args.nev))
        else:
            anaps.makeCubeData(cubeName)
        import sys
        sys.exit()
    if args.nev:
        print 'make cube with fewer events/bin'
        if nGood>0:
            anaps.makeCubeData(cubeName,  nEvtsPerBin=int(args.nev))
        #if nOff>0:
        #    anaps.makeCubeData(cubeName+'Off',  nEvtsPerBin=int(args.nev))
    else:
        if nGood>0:
            anaps.makeCubeData(cubeName)
        #if nOff>0:
        #    anaps.makeCubeData(cubeName+'Off')
