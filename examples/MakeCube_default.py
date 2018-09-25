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
    hutches=['amo','sxr','xpp','xcs','mfx','cxi','mec']
    hutch=None
    for thisHutch in hutches:
        if hostname.find(thisHutch)>=0:
            hutch=thisHutch.upper()
    if hutch is None:
        #then check current path
        path=os.getcwd()
        for thisHutch in hutches:
            if path.find(thisHutch)>=0:
                hutch=thisHutch.upper()
    if hutch is None:
        print 'cannot figure out which experiment to use, please specify -e <expname> on commandline'
        import sys
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
    detDict = {'source':'jungfrau512k','full':1, 'image':1, 'common_mode':7}
    #save photon images -- not used recently.
    #detDict = {'source':'jungfrau512k','photon_0p85':3.0, 'image':1}
    zylaDict = {'source':'zyla','full':1, 'common_mode':0}

    #CHANGE ME!
    ####
    # list of variables to bin
    ####
    #varList = ['ipm5/sum','snd_dio/t4d','snd_dio/dco','scan/_delay', zylaDict, detDict]
    varList = ['ipm5/sum','snd_dio/t4d','snd_dio/dco']

    #CHANGE ME
    ####
    # if you are interested in laser-xray delay, please select the delay of choice here!
    ####
    ana.setDelay(use_ttCorr=False, addEnc=False, addLxt=False)


    cubeName='cube' #initial name
    scanName, scanValues = ana.getScanValues()
    binSteps=[]
    binName=''
    filterName='filter1'
    if scanName!='':# and scanName.find('delay')<0:
        scanSteps = np.unique(scanValues)
        scanSteps = np.append(scanSteps, scanSteps[-1]+abs(scanSteps[1]-scanSteps[0]))#catch value at right edge?
        scanSteps = np.append(scanSteps[0]-abs(scanSteps[1]-scanSteps[0]),scanSteps) #catch values at left edge
        binSteps = scanSteps
        cubeName = scanName
        if scanName == 'lxt':    
            print 'bin data using ',scanName,' and bins: ',scanSteps
            binName='delay'
        else:
            print 'bin data using ',scanName,' and bins: ',scanSteps
            binName='scan/%s'%scanName
    else:
        cubeName='randomTest'
        binName='ipm4/sum'
        binVar=ana.getVar(binVarName)
        binSteps=np.percentile(binVar,[0,25,50,75,100])

    if args.nev:
        cubeName+='_%sEvents'%args.nev
    ana.addCube(cubeName,binName,binSteps,filterName)
    ana.addToCube(cubeName,varList)

    if args.cube:
        cubeName = args.cube

    ####
    #there are two ways to get multiple cubes: if you want a standard laser on/off set, use the onoff flag
    #if you e.g. have data with an external magnet switched on/off, use the 'Sel2' approach.
    #both are not necessary.
    ####
    ana.addCube(cubeName+'Sel2',binName,binSteps,'filter2')
    ana.addToCube(cubeName+'Sel2',varList)
    nSel = ana.getFilter(filterName).sum()
    nSel2 = ana.getFilter('filter2').sum()
    if args.nev:
        print 'make cube with fewer events/bin'
        if nSel>0:
            anaps.makeCubeData(cubeName,  nEvtsPerBin=int(args.nev))
        if nSel2>0:
            anaps.makeCubeData(cubeName+'Sel2',  nEvtsPerBin=int(args.nev))
    else:
        if nSel>0:
            anaps.makeCubeData(cubeName)
            #anaps.makeCubeData(cubeName, onoff=1) #request 'on' events base on input filter (add optical laser filter)
            #anaps.makeCubeData(cubeName, onoff=0) #request 'off' events base on input filter (switch optical laser filter, drop tt requirements)
        if nSel2>0:
            anaps.makeCubeData(cubeName+'Sel2')
