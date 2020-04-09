import numpy as np
import sys
import argparse
import socket
import os
import itertools
import psana
import h5py
from smalldata_tools import dropObject
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
fname=None
if args.file:
    fname = args.file

dirname = '/reg/d/psdm/%s/%s/hdf5/smalldata'%('xpp',expname)

if fname is None:
    fname = '%s_Run%03d.h5'%(expname,int(run))

sys.path.append('./smalldata_tools')
from SmallDataAna import SmallDataAna
ana = None
anaps = None
try:
    from SmallDataAna_psana import SmallDataAna_psana
    anaps = SmallDataAna_psana(expname,run,dirname,fname)
except:
    print "I have an error"
    pass

if anaps and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
    print 'create ana module from anaps'
    ana = anaps.sda
else:
    print 'we will now try to open the littleData file directly'
    ana = SmallDataAna(expname,run, dirname, fname)
    if 'fh5' not in ana.__dict__.keys():
        ana = None

if ana is not None:
    #build up your filter conditions. You can also use variables that can be build up from variables in the smallData file
    #see xpp wiki for description.
    ana.addCut('lightStatus/xray',0.5,1.5,'on')
    #add other selection cuts like this:
    #ana.addCut('ipm2/sum',0.2,10.,'on')

    #bin data in scan steps using the selection defined above.
    ana.addCube('cube','scan/varStep',np.arange(-0.5,ana.xrData.scan__varStep.max()+1+0.5,1),'on')

    #define the binning variable here as it'll appear multiple places
    #selectionVar='userValues/userCOM0'
    selectionVar='userValues/userError0'

    #try binning variable constructed of two userVars
    # selectionVarH = 'userValues/userHori'
    # selectionVarV = 'userValues/userVert'
    # selectionVar = 'pcDist'

    # userHori = ana.getVar(selectionVarH)
    # userVert = ana.getVar(selectionVarV)
    # ana.addVar(selectionVar, (userHori - 5)**2 + (userVert - 0)**2)

    #make an "opt" cube using the best image in a set of event indices rather than the sum
    ana.addToCube('cube',['scan/nano_x','scan/nano_y','ipm2/sum','ipm3/sum',selectionVar])
    #ana.addToCube('cube',['scan/nano_y','scan/nano_z',selectionVarH,selectionVarV,selectionVar])

    #cubeData will contain the "standard" cube, which means it'll have sums over all events in the bins (steps)
    #eventIdxDict will contain a dict of lists for each scan step that contain the fiducial&event_time for
    #    picking the events for the "single event cube". You add other variables in addIdxVar. This is used to select
    #    which event to pick
    cubeData, eventIdxDict = ana.makeCubeData('cube', returnIdx=True, addIdxVar=selectionVar)

    #setup to go back and get data from the "big" xtc files.
    runIdx = anaps.dsIdxRun
    cs140_dict = {'source': 'cs140_front', 'full': 1, 'common_mode':1}
    anaps.addDetInfo(cs140_dict)
    det = anaps.cs140_front
    det.evt = dropObject()

    # zyla_dict = {'source': 'zyla_1', 'full': 1, 'common_mode':-1}
    # anaps.addDetInfo(zyla_dict)
    # detzyla = anaps.zyla_1
    # detzyla.evt = dropObject()

    #get the average motor positions at scan step to be saved later in hdf5 files w/ cs140 data
    nano_x=cubeData['scan__nano_x'].values/cubeData['nEntries'].values
    nano_y=cubeData['scan__nano_y'].values/cubeData['nEntries'].values
    ipm2=cubeData['ipm2__sum'].values
    #this should go in loop.

    #pickVals is supposed to be a list of the pickVals you want a single event image cubes for
    #pickVals = [np.nanpercentile(ana.getVar(selectionVar),40), np.nanpercentile(ana.getVar(selectionVar),50), np.nanpercentile(ana.getVar(selectionVar),60)]
    pickVals = [ana.getVar(selectionVar).min()]
    print pickVals
    for ipV, pickVal in enumerate(pickVals):
        #this code here should allow the code to run MPI'ed
        if size>= len(pickVals):
            if rank!=ipV:
                continue
        print 'look at events close to ',pickVal

        #make a list of timestamps for events to get big data from.
        selEvt_fiducials=[]
        selEvt_evt_times=[]
        selEvt_selVar=[]
        for fids,times,selVars in zip(eventIdxDict['fiducial'],eventIdxDict['evttime'],eventIdxDict[selectionVar]):
            #get the position of the "best" event at this step
            #selecting the highest userVal event here. Can do whatever though.
            #shift selVar array as Zyla has information from previous shot.
            selVars = selVars[1:]
            fids = fids[:-1]
            times = times[:-1]
            ipm2s = ipm2s[:-1]
            #now remove bad ipm shots
            ipmGood = (ipm2s<10.).
            fids_filter = fids[ipm]
            #go on as before, event ID and zyla information are aligned now again.
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
        cs140_sum = []
        eventNr = []
        for evtfid,evttime in itertools.izip(selEvt_fiducials, selEvt_evt_times):
            evtt = psana.EventTime(evttime,evtfid)
            evt = runIdx.event(evtt)
            #now loop over detectors in this event
            det.getData(evt)
            det.processDetector()
            cs140_data.append(det.evt.write_full)
            cs140_img.append(np.rot90(det.det.image(run,det.evt.write_full), 2))
            cs140_sum.append(np.nansum(det.det.image(run,det.evt.write_full)))

        #now write hdf5 file & add detector info to it.
        cs140_data = np.array(cs140_data)
        cs140_img = np.array(cs140_img)
        hdf_filename = '%s/%s_%d_%s_%d.h5'%(dirname, expname, run, selectionVar.replace('/','_'), ipV)
        fh5 = h5py.File(hdf_filename,'w')
        print 'now save as file:  %s/%s_%d_%s_%d.h5'%(dirname,expname, run, selectionVar.replace('/','_'), ipV)
        fh5['cs140_raw']=cs140_data.astype(float)
        fh5['cs140_img']=np.maximum( 0, cs140_img.astype(float) )
        data_group=fh5.create_group('/entry/data')
        data_group['data']=fh5['cs140_img']
        fh5['nano_x']=nano_x.astype(float)
        fh5['nano_y']=nano_y.astype(float)
        fh5['cs140_sum']=ipm2.astype(float)
        fh5['selVar']=selEvt_selVar.astype(float) #<--
        fh5['pickValues']=np.array(pickVals).astype(float)
        fh5.close()

        ptychopos = []
        ptychopos_filename = '%s/positions_%s_%d_%s_%d.txt'%(dirname, expname, run, selectionVar.replace('/','_'), ipV)
        posX = nano_x.astype(float) / 0.2 / 1000.
        posY = nano_y.astype(float) / 0.2 / 1000.
        ints = [float(x) for x in cs140_sum]
        for i, (Int,(x,y)) in enumerate(zip(ints,zip(posX, posY))):
            line = "%d  %f  %f  %s %f\n"%(i, x, y, hdf_filename, Int)
            ptychopos.append(line)

        with open(ptychopos_filename, 'w') as G:
            G.writelines(ptychopos)


# TO DO write ptycho config file
# TO DO adapt path to XCS

