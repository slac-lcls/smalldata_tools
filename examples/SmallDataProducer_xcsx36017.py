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
##########################################################

def getROI_jungfrau(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <= 11:
        return  [ [[0,1], [271,419], [645,748]] ]
    else:
        return  [ [[0,1], [0,512], [0,1024]] ]

def getNmaxPhoton(run):
    if isinstance(run,basestring):
        run=int(run)

    if run <= 10:
        return 5000
    elif run >=145 and run <= 162:
        return 10000
    else:
        return 4000    

##########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw'] 
epicsPV = []
#fix timetool calibration if necessary
#ttCalib=[0.,2.,0.]
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

#epixnames = ['epix_ladm_1','epix_ladm_2']
#epixnames = ['epix_ladm_1']
epixnames = []
nPhoton = getNmaxPhoton(int(run)) #could use this, values hardcoded below for now.
for iepix,epixname in enumerate(epixnames):
    have_epix = checkDet(ds.env(), epixname)
    if have_epix:
        #create detector object. needs run for calibration data
        #common mode: 46 is the "original" to psana method 6(?)
        #row & column correction on data w/ photons&neighbors removed.
        #9sec/800events CM6
        #epix = DetObject(epixname ,ds.env(), int(run), common_mode=46)
        #7sec/800events CM6
        epix = DetObject(epixname ,ds.env(), int(run), common_mode=6)
        #4sec/800events no CM1
        #epix = DetObject(epixname ,ds.env(), int(run), common_mode=0)

        ## this will saw the epix data in "raw" format w/o a pixel mask(!), but with calibration
        #epix.saveFull()

        #now add photon algorithms. Only works for single energy photon data
        # ADU_per_photon: expected ADU for photon of expected energy
        # thresADU: fraction of photon energy in photon candidate 
        #(2 neighboring pixel)
        #retImg: 0 (def): only histogram of 0,1,2,3,...,24 photons/pixel is returned
        #        1 return Nphot, x, y arrays
        #        2 store image using photons /event
        #        -1 store "trunkated" image after noise supression.

        #if epixname=='epix_3':
        epix.addPhotons(ADU_per_photon=165, thresADU=0.85, retImg=1, nphotMax=100, nphotRet=nPhoton)
        #else:
        #    epix.addPhotons(ADU_per_photon=160, thresADU=0.9, retImg=1, nphotMax=100, nphotRet=int(nPhoton/5))

        #if iepix==0 and run < 200:
        #    epix.addPhotons(ADU_per_photon=160, thresADU=0.9, retImg=2, nphotMax=100)
        #else:
        #    epix.addPhotons(ADU_per_photon=160, thresADU=0.9, retImg=1, nphotMax=100, nphotRet=5000)
        #make tiny small data files w/o epix for now.
        if (int(run) == 162 or int(run)== 146):
            epix.addDroplet(threshold=10., thresholdLow=3., thresADU=0.,name='droplet')
            epix['droplet'].addDropletSave(maxDroplets=10000)
        dets.append(epix)

if int(run) == 3  and int(run)<=232:
    dets=[]
    have_zyla = checkDet(ds.env(), 'zyla')
    if have_zyla:
        zyla = DetObject('zyla' ,ds.env(), int(run))
        zyla.saveFull()
        dets.append(zyla)

jROI = getROI_jungfrau(run)
jungfrauname='jungfrau1M'
#jungfrauname='jungfrau512k'
have_jungfrau = checkDet(ds.env(), jungfrauname)
if have_jungfrau and int(run)>3:
    #create detector object. needs run for calibration data
    #optinally try common mode method 0 (only pedestal subtraction)
    #jungfrau = DetObject(jungfrauname ,ds.env(), int(run), common_mode=7)
    jungfrau = DetObject(jungfrauname ,ds.env(), int(run), common_mode=0)
    #for iROI,ROI in enumerate(jROI):
    #    print 'adding an ROI: ',jROI
    #    jungfrau.addROI('ROI_%d'%iROI, ROI)
    #    jungfrau.__dict__['ROI_%d'%iROI].addProj(pjName='thresAdu10', axis=-1, thresADU=1., mean=True)
    #    jungfrau.__dict__['ROI_%d'%iROI].addProj(pjName='thresRms5', axis=-1, thresRms=5, mean=True)
    #    jungfrau.__dict__['ROI_%d'%iROI].addProj(pjName='sum', axis=-1, mean=True)

    ##save full image
    #jungfrau.saveFull()

    ##save photonized jungfrau - smaller as integers
    jungfrau.addPhotons(ADU_per_photon=9.5, thresADU=0.9, retImg=1, nphotMax=5000)
    #jungfrau.storeSum(sumAlgo='calib')
    jungfrau.storeSum(sumAlgo='img')
    #jungfrau.storeSum(sumAlgo='thresADU9')
    #jungfrau.storeSum(sumAlgo='nhits_thresADU4')
    #jungfrau.storeSum(sumAlgo='square')
    
    dets.append(jungfrau)

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
##adding raw timetool traces:
#defaultDets.append(ttRawDetector(env=ds.env()))
##adding wave8 traces:
#defaultDets.append(wave8Detector('Wave8WF'))

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
            #print userDict[det._name].keys()
        except:
            pass

    smldata.event(userDict)


    #here you can add any data you like: example is a product of the maximumof two area detectors.
    #try:
    #    cspadMax = cspad.evt.dat.max()
    #    epix_vonHamosMax = epix_vonHamos.evt.dat.max()
    #    combDict = {'userValue': cspadMax*epix_vonHamosMax}
    #    smldata.event(combDict)
    #except:
    #    pass

sumDict={'Sums': {}}
for det in dets:
    for key in det.storeSum().keys():
        sumData=smldata.sum(det.storeSum()[key])
        sumDict['Sums']['%s_%s'%(det._name, key)]=sumData
if len(sumDict['Sums'].keys())>0:
    smldata.save(sumDict)

print 'rank %d on %s is finished'%(ds.rank, hostname)
smldata.save()
