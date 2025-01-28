# importing genereric python modules
import numpy as np
import h5py
import psana
import time
import argparse
import socket
import os

from smalldata_tools import defaultDetectors,epicsDetector,printMsg,detData,checkDet,getCfgOutput,getUserData,getUserEnvData
from smalldata_tools.DetObject import DetObject
from smalldata_tools.roi_rebin import ROIFunc
########################################################## 
##
## User Input start --> 
##
########################################################## 
##########################################################
# functions for run dependant parameters
##########################################################

def getROIs(run):
    if run>=210 and run<220:
	sigROI = [] #no signal apparent for run 210
        #np.append(sigROI,[[0,1], [150,200], [150,200]],axis=0)
    elif run>=220 and run<230:
	sigROI = [[1,2], [87,146], [6,384]]
    else:
	sigROI = []

    if len(sigROI)>1:
        return sigROI, 
    else:
        return sigROI

##########################################################
# run independent parameters 
##########################################################
#event codes which signify no xray/laser
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw'] 
epicsPV = [] 
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
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--nevt", help="number of events", type=int)
parser.add_argument("--dir", help="directory for output files (def <exp>/hdf5/smalldata)")
parser.add_argument("--offline", help="run offline (def for current exp from ffb)")
parser.add_argument("--gather", help="gather interval (def 100)", type=int)
parser.add_argument("--live", help="add data to redis database (quasi-live feedback)", action='store_true')
args = parser.parse_args()
hostname=socket.gethostname()
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
        sys.exit()
    expname=RegDB.experiment_info.active_experiment(hutch)[1]
    dsname='exp='+expname+':run='+run+':smd:dir=/reg/d/ffb/%s/%s/xtc:live'%(hutch.lower(),expname)
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
except:
    print 'failed to make MPIDataSource for ',dsname
    import sys
    sys.exit()

try:    
    if dirname is None:
        dirname = '/reg/d/psdm/%s/%s/hdf5/smalldata'%(hutch.lower(),expname)
    smldataFile = '%s/%s_Run%03d.h5'%(dirname,expname,int(run))

    smldata = ds.small_data(smldataFile,gather_interval=gatherInterval)
    if args.live:
        smldata.connect_redis()
except:
    print 'failed making the output file ',smldataFile
    import sys
    sys.exit()

########################################################## 
##
## User Input start --> 
##
########################################################## 
dets=[]
ROIs = getROIs(int(run))
haveCspad = checkDet(ds.env(), 'cs140_0')
if haveCspad:
    cspad = DetObject('cs140_0' ,ds.env(), int(run), name='cs140_0')
    for iROI,roi in enumerate(ROIs):
        print 'adding func'
        cspad.addFunc(ROIFunc(name='ROI_%d'%iROI, ROI=roi, writeArea=True))
    dets.append(cspad)

########################################################## 
##
## <-- User Input end
##
########################################################## 
#dets = [ det for det in dets if checkDet(ds.env(), det._srcName)]
dets = [ det for det in dets if checkDet(ds.env(), det.det.alias)]
#for now require all area detectors in run to also be present in event.

defaultDets = defaultDetectors(hutch)
#ttCalib=[0.,2.,0.]
#setParameter(defaultDets, ttCalib)
#aioParams=[[1],['laser']]
#setParameter(defaultDets, aioParams, 'ai')
if len(epicsPV)>0:
    defaultDets.append(epicsDetector(PVlist=epicsPV, name='epicsUser'))

#add config data here
userDataCfg={}
for det in dets:
    userDataCfg[det._name] = det.params_as_dict()
    #print(userDataCfg[det._name].keys())
Config={'UserDataCfg':userDataCfg}
smldata.save(Config)

for eventNr,evt in enumerate(ds.events()):
    printMsg(eventNr, evt.run(), ds.rank)

    if eventNr >= maxNevt/ds.size:
        break

    #add default data
    defData = detData(defaultDets, evt)
    #for key in defData.keys():
    #    print eventNr, key, defData[key]
    smldata.event(defData)

    #detector data using DetObject 
    userDict = {}
    for det in dets:
        try:
            det.getData(evt)
            det.processFuncs()
            userDict[det._name]=getUserData(det)
            try:
                userDict[det._name+'_env']=getUserEnvData(det) #need epix data to try
            except:
                pass
            #print userDict[det._name]
        except:
            pass
    smldata.event(userDict)

#    if args.live:
#        import time
#        time.sleep(0.1)

#gather whatever event did not make it to the gather interval
print 'rank %d on %s is finished'%(ds.rank, hostname)
smldata.save()
