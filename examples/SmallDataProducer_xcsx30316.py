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
def getROI_zyla(run):
    if isinstance(run,basestring):
        run=int(run)

    if run <= 6:
        return [[[0,512], [25, 275]],
                [[0,512], [25, 275]]]
    else:
        return [[[0,512], [0, 512]],
                [[0,512], [0, 512]]]

##########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw'] 
epicsPV = ['snd_t2_th']
#fix timetool calibration if necessary
#ttCalib=[0.,2.,0.]
ttCalib=[]#1.860828, -0.002950]
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

########################################################## 
##
## User Input start --> 
##
########################################################## 
dets=[]

#ROI_zyla = getROI_zyla(int(run))
#if checkDet(ds.env(), 'zyla_ladm'):
#    zyla = DetObject('zyla_ladm' ,ds.env(), int(run))
#    for iROI,ROI in enumerate(ROI_zyla):
#        zyla.addROI('ROI_%d'%iROI, ROI)
#    dets.append(zyla)

#epixnames=['epix_ladm_1','epix_ladm_2']
epixnames=['epix_ladm_1']
nPhot = 200
if int(run)==44:
    nDrop = [16000,300]
elif int(run)==46:
    nDrop = [6000,50]
elif int(run)==43:
    nDrop = [6000,4000]

if int(run)==43:
    ROI_low = [[260,320], [10,60]]
    ROI_med = [[260,320], [150,200]]
    ROI_high = [[260,320], [250,300]]
else:
    ROI_low = [[175,225], [290,340]]
    ROI_high = [[175,225], [380,430]]
for iepix,epixname in enumerate(epixnames):
    have_epix = checkDet(ds.env(), epixname)
    if have_epix:
        #create detector object. needs run for calibration data
        #common mode: 46 is the "original" to psana method 7(?)
        #row & column correction on data w/ photons&neighbors removed.
        epix = DetObject(epixname ,ds.env(), int(run), common_mode=46)

        #now add photon algorithms. Only works for single energy photon data
        # ADU_per_photon: expected ADU for photon of expected energy
        # thresADU: fraction of photon energy in photon candidate 
        #(2 neighboring pixel)
        #retImg: 0 (def): only histogram of 0,1,2,3,...,24 photons/pixel is returned
        #        1 return Nphot, x, y arrays
        #        2 store image using photons /event
        # addPhotons: use official algorithm, photons2&3 are own. 
        #             photon3 w/maxMethod=1 is very similar to official,
        #             official used to be slower. Need to verify performance &
        #                   provide example images
        epix.addPhotons(ADU_per_photon=160, thresADU=0.85, nphotRet=nPhot, retImg=2, ROI=ROI_low, name='ROI_low')
        epix.addPhotons(ADU_per_photon=160, thresADU=0.85, nphotRet=nPhot, retImg=2, ROI=ROI_med, name='ROI_med')
        epix.addPhotons(ADU_per_photon=160, thresADU=0.85, nphotRet=nPhot, retImg=2, ROI=ROI_high, name='ROI_high')

        #two threshold droplet finding. Tends to add photons together into single droplet if occupancy 
        #is not low, might need photonizing step to get single photon positions
        if (int(run)==43 or int(run)==38 or int(run)==39):
            epix.addDroplet(threshold=10., thresholdLow=3., thresADU=0.,name='droplet')
            epix['droplet'].addDropletSave(maxDroplets=300, ROI=ROI_low, name='ROI_low')
            epix['droplet'].addDropletSave(maxDroplets=300, ROI=ROI_med, name='ROI_med')
            epix['droplet'].addDropletSave(maxDroplets=300, ROI=ROI_high, name='ROI_high')
            epix.addROI('ROI_low',ROI_low, writeArea=True)
            epix.addROI('ROI_med',ROI_med, writeArea=True)
            epix.addROI('ROI_high',ROI_high, writeArea=True)

        dets.append(epix)


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

    #first event.
    if ds.rank==0 and eventNr==0 and (args.live or args.liveFast):
        if not args.liveFast:
            #this saves all fields
            smldata.connect_redis()
        else:
            redisKeys = defaultRedisVars(hutch)
            redisList=['fiducials','event_time']
            for key in redisKeys:
                if key.find('/')>=0 and key in smldata._dlist.keys():
                    redisList.append(key)
                else:
                    for sdkey in smldata._dlist.keys():
                        if sdkey.find(key)>=0:
                            redisList.append(sdkey)
            print 'Saving in REDIS: ',redisList
            smldata.connect_redis(redisList)

print 'rank %d on %s is finished'%(ds.rank, hostname)
smldata.save()
