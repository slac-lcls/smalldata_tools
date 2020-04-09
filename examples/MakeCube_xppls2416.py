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
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run",type=int)
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--dir", help="directory w/ smallData file if not default")
parser.add_argument("--nev", help="number of events/bin")
parser.add_argument("--cube", help="cube name")
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
if rank==0: print 'try to import smalldata'
from SmallDataAna import SmallDataAna
if rank==0: print 'imported smalldata'
from SmallDataAna_psana import SmallDataAna_psana
if rank==0: print 'imported smalldata_psana'

ana = None
anaps = SmallDataAna_psana(expname,run,dirname,fname)
if anaps and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
    if rank==0: print 'create ana module from anaps'
    ana = anaps.sda
else:
    if rank==0: print 'we will now try to open the littleData file directly'
    ana = SmallDataAna(expname,run, dirname, fname)
    if 'fh5' not in ana.__dict__.keys():
        ana = None

if ana is not None:
    if rank==0: ana.printRunInfo()

    #somewhat non-obvious: if lower==upper; REJECT this point (e.g. FLTPOS==0)
    ana.addCut('lightStatus/xray',0.5,1.5,'good')
    ana.addCut('lightStatus/laser',0.5,1.5,'good')
    ana.addCut('ipm2/sum',0.05,1.,'good')
    #ana.addCut('ipm3/sum',0.25,10.,'onSel')

    if rank==0: ana.printSelections('good')
    #save raw data
    #detDict = {'source':'jungfrau1M','full':1}
    #save image data
    #detDict = {'source':'jungfrau1M','full':1,'image':1}
    #detDict = {'source':'jungfrau1M', 'full':1, 'image':1} #add all images in bin
    #detDict = {'source':'jungfrau1M', 'full':1, 'image':1, 'common_mode': 72} #use common_mode
    detDict = {'source':'jungfrau1M', 'full':1, 'image':1, 'common_mode': 72, 'thresADU': 2.} #threshold on single images before adding
    #detDict = {'source':'jungfrau1M', 'photon_09':9.5, 'image':1, 'common_mode': 72} #photonize image.

    #detDict = {'source':'jungfrau1M', 'full':1, 'image':1, 'common_mode': 7, 'thresADU': 2.} #threshold on single images before adding
    #detDict = {'source':'jungfrau1M', 'photon_09':9.5, 'image':1, 'common_mode': 7} #photonize image.

    #example for test runs w/ noise.
    #ana.addCube('cubeDelay','delay',np.linspace(-5, 155, 33),'onSel')

    cubeName=''
    scanName, scanValues = ana.getScanValues()
    if scanName!='' and scanName.find('delay')<0:
        scanSteps = np.unique(scanValues)
        scanSteps = np.append(scanSteps[0]-(scanSteps[1]-scanSteps[0]),scanSteps)
        scanSteps = np.append(scanSteps, scanSteps[-1]+(scanSteps[1]-scanSteps[0]))
        cubeName=scanName
        ana.addCube(cubeName,'scan/%s'%scanName,scanSteps,'good')
    else:
        #here we care about the delay, so we ensure the timetool looks good.
        ana.addCut('tt/AMPL',0.0025,10.,'good')
        ana.addCut('tt/FLTPOSFWHM',80.,130.,'good')
        cubeName = 'delay'
        encDelay = ana.getVar('enc/lasDelay')
        lxt = ana.getVar('epics/lxt_ttc')*1e12
        encDelay = encDelay.copy()+lxt
        minEnc = (int(min(encDelay)*10.))
        if minEnc<0:
            minEnc+=-1
        minEnc /= 10.
        maxEnc = (int(max(encDelay)*10.)+1)/10.
        #binWidth=0.05
        binWidth=0.1
        #binWidth=5.
        #ana.addCube(cubeName,'delay',np.linspace(-2, 50, 26),'good')
        ana.addCube(cubeName,'delay',np.arange(minEnc, maxEnc, binWidth),'good')

    if args.cube:
        cubeName=args.cube
    ana.addToCube(cubeName,['ipm2/sum','ipm3/sum','diodeU/channels', detDict])
    #ana.addToCube(cubeName,['ipm2/sum','ipm3/sum','diodeU/channels'])

    if cubeName=='':
        print 'no cube defined, will exit'
        sys.exit()

    if args.nev:
        print 'make cube with fewer events/bin'
        #anaps.makeCubeData(cubeName,  nEvtsPerBin=int(args.nev))
        if cubeName.find('delay')>=0:
            anaps.makeCubeData(cubeName,  nEvtsPerBin=int(args.nev), onoff=0)
            #anaps.makeCubeData(cubeName,  nEvtsPerBin=int(args.nev), onoff=1)
        else:
            anaps.makeCubeData(cubeName,  nEvtsPerBin=int(args.nev), onoff=0)
            #anaps.makeCubeData(cubeName,  nEvtsPerBin=int(args.nev), onoff=1)
    else:
        if cubeName.find('delay')>=0:
            anaps.makeCubeData(cubeName, onoff=0)
            anaps.makeCubeData(cubeName, onoff=1)
        else:
            anaps.makeCubeData(cubeName, onoff=0)
            anaps.makeCubeData(cubeName, onoff=1)
