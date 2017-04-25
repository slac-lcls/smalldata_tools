# importing generic python modules
import numpy as np
import h5py
import psana
import time
import argparse
import socket
import os

from smalldata_tools import defaultDetectors,epicsDetector,printMsg,detData,DetObject,checkDet,getCfgOutput,getUserData,getUserEnvData
#from SmallDataUtils import *
#from utilities_wSmallData import *
#from DetObject import DetObject
########################################################## 
##
## User Input start --> 
##
########################################################## 
##########################################################
# functions for run dependant parameters
##########################################################
# none for now, start w/ full image saving to see if 
# start works with the new smallData
#
def getAzIntParams(run):
    ret_dict = {'eBeam': 8.015}
    ret_dict['cspad_center'] = [87697.946016760892, 94865.383526655729]
    ret_dict['cspad_dis_to_sam'] = 110.
    return ret_dict

def getROIs(run):
    return [[10, 20], [10, 20]]

def getNmaxDrop(run):
    if run >= 10:
        return 2000, 100
    else:
        return 400,400
##########################################################
# run independent parameters 
##########################################################
#event codes which signify no xray/laser
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw'] 
epicsPV = [] 
acqROI = [0,8000]
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
args = parser.parse_args()
if not args.run:
    run=raw_input("Run Number:\n")
else:
    run=args.run
if not args.exp:
    hutches=['amo','sxr','xpp','xcs','mfx','cxi','mec']
    hostname=socket.gethostname()
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
except:
    print 'failed making the output file ',smldataFile
    import sys
    sys.exit()

########################################################## 
##
## User Input start --> 
##
########################################################## 
nDrop = getNmaxDrop(int(run))
azIntParams = getAzIntParams(run)
epixnames = ['epix_vonHamos']
dets=[]
for iepix,epixname in enumerate(epixnames):
    have_epix = checkDet(ds.env(), epixname)
    if have_epix:
        print 'creating epix detector object  for epix ',epixname
        epix = DetObject(epixname ,ds.env(), int(run), name=epixname,common_mode=46)

        epix.addDroplet(threshold=10., thresholdLow=3., thresADU=0.,name='droplet')
        epix['droplet'].addAduHist([0.,1500.])
        epix['droplet'].addDropletSave(maxDroplets=nDrop[iepix])

        dets.append(epix)

ROIs = getROIs(int(run))
#if checkDet(env, 'epix_Rowland']:
#    epix = DetObject(epixname ,env, int(run), name='epix_Rowland', common_mode=46)
#    epix.addROI('ROI', ROIs, rm=30.0)
#    dets.append(epix)

haveCspad = checkDet(ds.env(), 'cspad')
if haveCspad:
    cspad = DetObject('cspad' ,ds.env(), int(run), name='cspad')
    for iROI,ROI in enumerate(ROIs):
        cspad.addROI('ROI_%d'%iROI, ROI)

    cspad.azav_eBeam=azIntParams['eBeam']
    if azIntParams.has_key('cspad_center'):
        cspad.azav_center=azIntParams['cspad_center']
        cspad.azav_dis_to_sam=azIntParams['cspad_dis_to_sam']
        try:
            cspad.addAzAv(phiBins=7)
        except:
            pass
    dets.append(cspad)


########################################################## 
##
## <-- User Input end
##
########################################################## 
dets = [ det for det in dets if checkDet(ds.env(), det._srcName)]
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
    userDataCfg[det._name]=getCfgOutput(det)
Config={'UserDataCfg':userDataCfg}
smldata.save(Config)

reqDet=[]
if expname=='xcsl2516':
    reqDet = ['acq01','epix_1','epix_2','epix_3','epix_4']
elif expname=='xcsm9816' or expname=='xcs06016':
    reqDet = ['epix_1','epix_2','epix_3','epix_4']
elif expname=='xcs12316':
    reqDet = ['epix_1','epix_2','epix_ladm']
elif expname=='xppn4516':
    reqDet = ['epix','epix_alcove']
for eventNr,evt in enumerate(ds.events()):
    printMsg(eventNr, evt.run(), ds.rank)

    if eventNr >= maxNevt/ds.size:
        break

    allDets = True
    aliases = [ k.alias() for k in evt.keys() ]
    for detn in reqDet:
        if not detn in aliases:
            allDets = False
    if not allDets:
        continue

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
                userDict[det._name+'_env']=getUserEnvData(det)
            except:
                pass
            #print userDict[det._name]
        except:
            pass
    smldata.event(userDict)

#gather whatever event did not make it to the gather interval
print 'I am finished, call save to be safe...'    
smldata.save()
