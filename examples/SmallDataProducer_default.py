# importing generic python modules
import numpy as np
import h5py
import psana
import time
import argparse
import socket
import os
import RegDB.experiment_info
from smalldata_tools import defaultDetectors,epicsDetector,printMsg,detData
from smalldata_tools import ttRawDetector,wave8Detector
########################################################## 
##
## User Input start --> 
##
########################################################## 
##########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw'] 
epicsPV = [] 
#fix timetool calibration if necessary
#ttCalib=[0.,2.,0.]
ttCalib=[]
#decide which analog input to save & give them nice names
#aioParams=[[1],['laser']]
aioParams=[]
########################################################## 
##
## <-- User Input end
##
########################################################## 


##########################################################
#command line input parameter: definitions & reading
##########################################################
maxNevt=1e9
gatherInterval=100
dirname = None
parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run")
parser.add_argument("--exp", help="expeariment name")
parser.add_argument("--nevt", help="number of events", type=int)
parser.add_argument("--dir", help="directory for output files (def <exp>/hdf5/smalldata)")
parser.add_argument("--offline", help="run offline (def for current exp from ffb)")
parser.add_argument("--gather", help="gather interval (def 100)", type=int)
parser.add_argument("--live", help="add data to redis database (quasi-live feedback)", action='store_true')
hostname=socket.gethostname()
args = parser.parse_args()
if not args.run:
    run=raw_input("Run Number:\n")
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
    expname=RegDB.experiment_info.active_experiment(hutch)[1]
    dsname='exp='+expname+':run='+run+':smd:dir=/reg/d/ffb/%s/%s/xtc:live'%(hutch.lower(),expname)
    #data gets removed from ffb faster now, please check if data is still available
    isLive = (RegDB.experiment_info.experiment_runs(hutch)[-1]['end_time_unix'] is None)
    if not isLive:
        dirname = '/reg/d/ffb/%s/%s/xtc'%(hutch.lower(),expname)
        xtcname=dirname+'/e*-r%04d-*'%int(run)
        import glob
        presentXtc=glob.glob('%s'%xtcname)
        if len(presentXtc)==0:
            dsname='exp='+expname+':run='+run+':smd'
else:
    expname=args.exp
    hutch=expname[0:3]
    dsname='exp='+expname+':run='+run+':smd'
if args.offline:
    dsname='exp='+expname+':run='+run+':smd'
if args.gather:
    gatherInterval=args.gather
if args.nevt:
    maxNevt=args.nevt
if args.dir:
    dirname=args.dir
    if dirname[-1]=='/':
        dirname=dirname[:-1]

debug = True
time_ev_sum = 0.
try:
    ds = psana.MPIDataSource(dsname)
    if ds.rank==0:
        print 'looking at data: ',dsname
    if dirname is None:
        dirname = '/reg/d/psdm/%s/%s/hdf5/smalldata'%(hutch.lower(),expname)
    smldataFile = '%s/%s_Run%03d.h5'%(dirname,expname,int(run))
    if ds.rank==0:
        print 'saving data in file: %s with gatherINterval %f'%(smldataFile,gatherInterval)
    smldata = ds.small_data(smldataFile,gather_interval=gatherInterval)
    if args.live:
        smldata.connect_redis()
except:
    print 'we seem to not have small data, you need to use the idx file based Ldat_standard code....',dsname
    import sys
    sys.exit()


defaultDets = defaultDetectors(hutch)
if len(epicsPV)>0:
    defaultDets.append(epicsDetector(PVlist=epicsPV, name='epicsUser'))
if len(ttCalib)>0:
    setParameter(defaultDets, ttCalib)
if len(aioParams)>0:
    setParameter(defaultDets, aioParams, 'ai')

lasttime=time.time()
for eventNr,evt in enumerate(ds.events()):
    #printMsg(eventNr, evt.run(), ds.rank, ds.size)
    printMsg(eventNr, evt.run(), ds.rank)
    if (eventNr%int(1000/ds.size)==0):
        nowtime=time.time()
        if ds.rank==0:
            print 'in run %g needed %f for 1k events at %d events'%(evt.run(),(nowtime-lasttime)*1000./((1000/ds.size)*float(ds.size)),ds.size*eventNr)
        lasttime=nowtime
    #time_ev_start = MPI.Wtime()

    if eventNr >= maxNevt:
        break

    #add default data
    defData = detData(defaultDets, evt)
    smldata.event(defData)

print 'rank %d on %s is finished'%(ds.rank, hostname)
smldata.save()
