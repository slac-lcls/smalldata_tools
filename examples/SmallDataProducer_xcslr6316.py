# importing generic python modules
import numpy as np
import psana
import time
import argparse
import socket
import os
import RegDB.experiment_info
from smalldata_tools import defaultDetectors,epicsDetector,printMsg,detData,DetObject
from smalldata_tools import checkDet,getCfgOutput,getUserData,getUserEnvData,dropObject
from smalldata_tools import ttRawDetector,wave8Detector,defaultRedisVars,setParameter
########################################################## 
##
## User Input start --> 
##
########################################################## 
##########################################################
# functions for run dependant parameters
##########################################################d
def getROI_jungfrau(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <= 11:
        return  [ [[271,419], [645,748]] ]
    elif run <= 21:
        return  [ [[0,512], [0,1024]] ]
    else:
        return  [ [[224,281], [231,341]] ]

def getROI_zyla(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <= 11:
        return  [ [[224,281], [231,341]] ]
    elif run <= 21:
        return  [ [[0,5000], [0,5000]] ]  #this will use the whole detector
    else:
        return  [ [[224,281], [231,341]] ]

def getROI_zyla_ladm(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <= 11:
        return  [ [[224,281], [231,341]] ]
    elif run <= 21:
        return  [ [[0,5000], [0,5000]] ]  #this will use the whole detector
    else:
        return  [ [[224,281], [231,341]] ]


##########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw'] 
epicsPV = ['bsx','bsy','bsth','bmx','bmth','bmchi','samx','samy','samth','sam2th','samh','samn','ath','my ','diode_x','jjvg','jjvo','jjhg','jjho','dety','detz','lens1x','lens1y','lens1z','lens1thx','lens1thy','lens1f','lens2x','lens2y','lens2z','lens2thx','lens2thy','lens2f','snd_t1_l']
#fix timetool calibration if necessary
ttCalib=[]
#ttCalib=[1.860828, -0.002950]
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
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--nevt", help="number of events", type=int)
parser.add_argument("--dir", help="directory for output files (def <exp>/hdf5/smalldata)")
parser.add_argument("--offline", help="run offline (def for current exp from ffb)")
parser.add_argument("--gather", help="gather interval (def 100)", type=int)
parser.add_argument("--live", help="add data to redis database (quasi-live feedback)", action='store_true')
parser.add_argument("--liveFast", help="add data to redis database (quasi-live feedback)", action='store_true')
parser.add_argument("--norecorder", help="ignore recorder streams", action='store_true')
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
        import sys
        sys.exit()
    expname=RegDB.experiment_info.active_experiment(hutch)[1]
    dsname='exp='+expname+':run='+run+':smd:dir=/reg/d/ffb/%s/%s/xtc:live'%(hutch.lower(),expname)
    #data gets removed from ffb faster now, please check if data is still available
    isLive = (RegDB.experiment_info.experiment_runs(hutch)[-1]['end_time_unix'] is None)
    if not isLive:
        xtcdirname = '/reg/d/ffb/%s/%s/xtc'%(hutch.lower(),expname)
        xtcname=xtcdirname+'/e*-r%04d-*'%int(run)
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

if args.norecorder:
    dsname=dsname+':stream=0-79'

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

except:
    print 'failed making the output file ',smldataFile
    import sys
    sys.exit()

if ds.rank==0:
    version='unable to detect psana version'
    for dirn in psana.__file__:
        if dirn.find('ana-')>=0:
            version=dirn
    print 'Using psana version ',version
########################################################## 
##
## User Input start --> 
##
########################################################## 
dets=[]

zylanames=['zyla_ladm','zyla']
zROIs = [ getROI_zyla_ladm(run),  getROI_zyla(run)]
for zylaname,zROI in zip(zylanames,zROIs):
    have_zyla = checkDet(ds.env(), zylaname)
    if have_zyla:
#    #use anaps.AvImage('zyla'); anaps.SelectRegion('zyla') to get the tuple for the ROI selection
        zyla = DetObject(zylaname ,ds.env(), int(run)) #subtract pedestal if possible.
        for iROI,ROI in enumerate(zROI):
            zyla.addROI('ROI%d'%iROI, ROI)
            #zyla.addROI('ROI%d'%iROI, ROI, writeArea=True)
        dets.append(zyla)

jROI = getROI_jungfrau(run)
jungfrauname='jungfrau512k'
have_jungfrau = checkDet(ds.env(), jungfrauname)
if have_jungfrau:
    #create detector object. needs run for calibration data
    #optinally try common mode method 7
    jungfrau = DetObject(jungfrauname ,ds.env(), int(run), common_mode=0)
    for iROI,ROI in enumerate(zROI):
        jungfrau.addROI('ROI_%d'%iROI, ROI)
        jungfrau.__dict__['ROI_%d'%iROI].addProj(pjName='thresAdu10', axis=-1, thresADU=1., mean=True)
        jungfrau.__dict__['ROI_%d'%iROI].addProj(pjName='thresRms5', axis=-1, thresRms=5, mean=True)
        jungfrau.__dict__['ROI_%d'%iROI].addProj(pjName='sum', axis=-1, mean=True)
    jungfrauImg = None
    dets.append(jungfrau)

##adding raw timetool traces:
#defaultDets.append(ttRawDetector(env=ds.env()))
##adding wave8 traces:
#defaultDets.append(wave8Detector('Wave8WF'))
########################################################## 
##
## <-- User Input end
##
########################################################## 
dets = [ det for det in dets if checkDet(ds.env(), det._srcName)]
#for now require all area detectors in run to also be present in event.

defaultDets = defaultDetectors(hutch)
if len(ttCalib)>0:
    setParameter(defaultDets, ttCalib)
if len(aioParams)>0:
    setParameter(defaultDets, aioParams, 'ai')
if len(epicsPV)>0:
    defaultDets.append(epicsDetector(PVlist=epicsPV, name='epicsUser'))

#add config data here
userDataCfg={}
for det in dets:
    userDataCfg[det._name]=getCfgOutput(det)
Config={'UserDataCfg':userDataCfg}
smldata.save(Config)

for eventNr,evt in enumerate(ds.events()):
    printMsg(eventNr, evt.run(), ds.rank, ds.size)

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
            #this should be a plain dict. Really.
            det.evt = dropObject()
            det.getData(evt)
            det.processDetector()
            userDict[det._name]=getUserData(det)
            try:
                envData=getUserEnvData(det)
                if len(envData.keys())>0:
                    userDict[det._name+'_env']=envData
            except:
                pass
            #print userDict[det._name]
        except:
            pass
    smldata.event(userDict)

    try:
        if jungfrauImg is None:
            jungfrauImg = jungfrau.evt.dat.squeeze()
        else:
            jungfrauImg += jungfrau.evt.dat.squeeze()
    except:
        pass
            
runsum=smldata.sum(jungfrauImg)
print 'rank %d on %s is finished'%(ds.rank, hostname)
smldata.save(jungfrauImg=runsum)
