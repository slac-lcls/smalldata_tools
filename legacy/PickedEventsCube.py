import sys
import numpy as np
import argparse
from matplotlib import pyplot as plt
import socket
import os
import itertools    
import psana
import h5py
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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

from smalldata_tools.SmallDataAna import SmallDataAna
ana = None
anaps = None
try:
    from smalldata_tools.SmallDataAna_psana import SmallDataAna_psana
    anaps = SmallDataAna_psana(expname,run,dirname,fname)
except:
    pass

if anaps and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
    print('create ana module from anaps')
    ana = anaps.sda
else:
    print('we will now try to open the littleData file directly')
    ana = SmallDataAna(expname,run, dirname, fname)
    if 'fh5' not in ana.__dict__.keys():
        ana = None

if ana is not None:
    #build up your filter conditions. You can also use variables that can be build up from variables in the smallData file
    #see xpp wiki for description.
    ana.addCut('lightStatus/xray',0.5,1.5,'on')
    #add other selection cuts like this:
    #ana.addCut('ipm3/sum',0.2,10.,'on')
    
    #bin data in scan steps using the selection defined above.
    ana.addCube('cube','scan/varStep',np.arange(-0.5,ana.xrData.scan__varStep.max()+0.5,1),'on')

    #define the binning variable here as it'll appear multiple places
    selectionVar='userValues/mycalc'
    #make an "opt" cube using the best image in a set of event indices rather than the sum
    ana.addToCube('cube',['scan/nano_y','scan/nano_z',selectionVar])

    #cubeData will contain the "standard" cube, which means it'll have sums over all events in the bins (steps)
    #eventIdxDict will contain a dict of lists for each scan step that contain the fiducial&event_time for
    #    picking the events for the "single event cube". You add other variables in addIdxVar. This is used to select
    #    which event to pick
    cubeData, eventIdxDict = ana.makeCubeData('cube', returnIdx=True, addIdxVar=selectionVar)

    #setup to go back and get data from the "big" xtc files.
    runIdx = anaps.dsIdxRun
    cs140_dict = {'source': 'cs140_dump', 'full': 1, 'common_mode':1}
    anaps.addDetInfo(cs140_dict)
    det = anaps.cs140_dump

    #get the average motor positions at scan step to be saved later in hdf5 files w/ cs140 data
    nano_y=cubeData['scan__nano_y'].values/cubeData['nEntries'].values
    nano_z=cubeData['scan__nano_z'].values/cubeData['nEntries'].values
    #this should go in loop.

    #pickVals is supposed to be a list of the pickVals you want a single event image cubes for
    pickVals = [np.nanpercentile(ana.getVar(selectionVar),40),np.nanpercentile(ana.getVar(selectionVar),60)]
    for ipV, pickVal in enumerate(pickVals):
        #this code here should allow the code to run MPI'ed
        if size>= len(pickVals):
            if rank!=ipV:
                continue

        #make a list of timestamps for events to get big data from.
        selEvt_fiducials=[]
        selEvt_evt_times=[]
        selEvt_selVar=[]
        for fids,times,selVars in zip(eventIdxDict['fiducial'],eventIdxDict['evttime'],eventIdxDict[selectionVar]):
            #get the position of the "best" event at this step
            #selecting the highest userVal event here. Can do whatever though.
            selUserVal = np.argmin(np.abs(selVars-pickVal)).values.tolist()
            #and now get the event ID & selection variable at this point
            selEvt_fiducials.append(fids[selUserVal].values.tolist())
            selEvt_evt_times.append(times[selUserVal].values.tolist())
            selEvt_selVar.append(selVars[selUserVal].values.tolist())
        selEvt_fiducials=np.array(selEvt_fiducials)
        selEvt_evt_times=np.array(selEvt_evt_times)
        selEvt_selVar=np.array(selEvt_selVar)

        #now loop over these events to get data. Get both raw data & image here.
        cs140_data = []
        cs140_img = []
        for evtfid,evttime in itertools.izip(selEvt_fiducials, selEvt_evt_times):
            evtt = psana.EventTime(evttime,evtfid)
            evt = runIdx.event(evtt)
            #now loop over detectors in this event
            det.getData(evt)
            det.processFuncs()
            cs140_data.append(det.evt.write_full)
            cs140_img.append(det.det.image(run,det.evt.write_full))

        #now write hdf5 file & add detector info to it.
        cs140_data = np.array(cs140_data)
        cs140_img = np.array(cs140_img)
        fh5 = h5py.File('%s/%s_%d_%s_%d.h5'%(dirname,expname, run, selectionVar.replace('/','_'), ipV),'w')
        fh5['cs140_raw']=cs140_data.astype(float)
        fh5['cs140_img']=cs140_img.astype(float)
        fh5['nano_y']=nano_y.astype(float)
        fh5['nano_z']=nano_z.astype(float)
        fh5['selVar']=selEvt_selVar.astype(float)
        fh5['pickValues']=np.array(pickVals).astype(float)
        fh5.close()
