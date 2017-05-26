import sys
import numpy as np
import argparse
from matplotlib import pyplot as plt
import socket
import os
from smalldata_tools import droplets,hist2d,getUserData,rebin,addToHdf5
from smalldata_tools import dropObject
from IPython.terminal.prompts import Prompts,Token
import itertools

class myPrompt(Prompts):
    def in_prompt_tokens(self, cli=None):
        return [(Token.Prompt, 'SDAna In ['), (Token.PromptNum, str(self.shell.execution_count)), (Token.Prompt, ']: ' ),]
    def out_prompt_tokens(self):
        return [(Token.Prompt, 'SDAna Out ['), (Token.PromptNum, str(self.shell.execution_count)), (Token.Prompt, ']: ' ),]
        #return [(Token.Prompt, 'LDAna Out:' ),]

ip=get_ipython()
ip.prompts = myPrompt(ip)

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

from smalldata_tools import SmallDataAna
ana = None
anaps = None
try:
    from smalldata_tools import SmallDataAna_psana
    anaps = SmallDataAna_psana(expname,run,dirname,fname)
except:
    pass

if anaps is not None and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
    print 'create ana module from anaps'
    ana = anaps.sda
else:
    print 'we will now try to open the littleData file directly'
    ana = SmallDataAna(expname,run, dirname, fname)
    if 'fh5' not in ana.__dict__.keys():
        ana = None

if ana is not None:
    ana.printRunInfo()

    ana.addCut('lightStatus/xray',0.5,1.5,'xon')

    ana.addCut('lightStatus/xray',0.5,1.5,'on')
    ana.addCut('lightStatus/laser',0.5,1.5,'on')

    ana.addCut('lightStatus/xray',0.5,1.5,'off')
    ana.addCut('lightStatus/laser',-0.5,0.5,'off')
    
    #delay=ana.getDelay()
    ana.addCube('cube','scan/varStep',np.arange(-0.5,ana.xrData.scan__varStep.max()+0.5,1),'on')

    #make an "opt" cube using the best image in a set of event indices rather than the sum
    addVar='userValues/mycalc'
    ana.addToCube('cube',['scan/nano_y','scan/nano_z',addVar])

    cubeData, eventIdxDict = ana.makeCubeData('cube', returnIdx=True, addIdxVar=addVar)

    nano_y=cubeData['scan__nano_y'].values/cubeData['nEntries'].values
    nano_z=cubeData['scan__nano_z'].values/cubeData['nEntries'].values

    import psana
    from smalldata_tools import dropObject

    #selecting the highest userVal event here. Can do whatever though.
    maxVar_fiducials=[]
    maxVar_evt_times=[]
    targetVal = np.nanpercentile(ana.getVar(addVar),50)
    for fids,times,addVars in zip(eventIdxDict['fiducial'],eventIdxDict['evttime'],eventIdxDict[addVar]):
        #maxUserVal = np.argmax(addVars).values.tolist()        
        maxUserVal = np.argmax(np.abs(addVars-targetVal)).values.tolist()
        maxVar_fiducials.append(fids[maxUserVal].values.tolist())
        maxVar_evt_times.append(times[maxUserVal].values.tolist())
    maxVar_fiducials=np.array(maxVar_fiducials)
    maxVar_evt_times=np.array(maxVar_evt_times)



    runIdx = anaps.dsIdxRun
    cs140_dict = {'source': 'cs140_dump', 'full': 1, 'common_mode':1}
    anaps.addDetInfo(cs140_dict)
    det = anaps.cs140_dump
    det.evt = dropObject()
    cs140_data = []
    cs140_img = []
    for evtfid,evttime in itertools.izip(maxVar_fiducials, maxVar_evt_times):
        evtt = psana.EventTime(evttime,evtfid)
        evt = runIdx.event(evtt)
        #now loop over detectors in this event
        det.getData(evt)
        det.processDetector()
        cs140_data.append(det.evt.write_full)
        cs140_img.append(det.det.image(run,det.evt.write_full))

    #now write hdf5 file & add detector info to it.
    #note that cs140_data is a masked array.....
    cs140_data = np.array(cs140_data)
    cs140_img = np.array(cs140_img)
    import h5py
    fh5 = h5py.File('test_out.h5','w')
    fh5['cs140_raw']=cs140_data.astype(float)
    fh5['nano_y']=nano_y.astype(float)
    fh5['nano_z']=nano_z.astype(float)
    fh5.close()
