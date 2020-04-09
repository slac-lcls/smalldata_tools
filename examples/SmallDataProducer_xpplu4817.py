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
########################################################## 
##
## User Input start --> 
##
########################################################## 
##########################################################
# functions for run dependant parameters
##########################################################
def getAzIntParams(run):
    if isinstance(run,basestring):
        run=int(run)
        
    ret_dict = {'eBeam': 8.35}
    ret_dict['cspad_center'] = [90201.147384, 93557.517461]
    ret_dict['cspad_dis_to_sam'] = 70.
    return ret_dict

def getCspadROIs(run):
    if isinstance(run,basestring):
        run=int(run)

    return [  [[1,2], [5,174], [6,300]], 
             [[9,10], [0,175], [3,348]] ]

def getNmaxDrop(run):
    if isinstance(run,basestring):
        run=int(run)

    if run <= 2:
        return 10
    elif run <= 5:
        return 1500
    elif run == 80:
        return 35000
    elif run == 81:
        return 30000
    else:
        return 40
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
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--nevt", help="number of events", type=int)
parser.add_argument("--dir", help="directory for output files (def <exp>/hdf5/smalldata)")
parser.add_argument("--offline", help="run offline (def for current exp from ffb)")
parser.add_argument("--gather", help="gather interval (def 100)", type=int)
parser.add_argument("--live", help="add data to redis database (quasi-live feedback)", action='store_true')
parser.add_argument("--skipDrop", help="save droplet information", action='store_true')

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
    dsname='exp='+expname+':run='+run+':smd:live'
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
#ttCalib=[0.,2.,0.]
#setParameter(defaultDets, ttCalib, 'tt')
##this gives the analog input channels friendlier names
#aioParams=[[1],['laser']]
#setParameter(defaultDets, aioParams, 'ai')


nDrop = getNmaxDrop(run)
epixnames = ['epix_1']
dets=[]
#epixROI=[ [[250,580],[250,580]], 
#          [[330,500],[330,500]] ]
if run <= 40: # filt and nofilt
	epixROI=[ [[68,580],[38,262]], 
                  [[117,624],[318,605]] ]
elif run<=79: # based on 83 full filtered epix
        epixROI=[[[94,616],[190,665]]]
else: 
        epixROI=[[[94,616],[190,665]]]

#epixROI_f=[[68,580],[38,262]]
#epixROI_nf=[[117,624],[318,605]]



for iepix,epixname in enumerate(epixnames):
    have_epix = checkDet(ds.env(), epixname)
    if have_epix:
        #print 'creating epix detector object  for epix ',epixname
        epix = DetObject(epixname ,ds.env(), int(run), name=epixname,common_mode=6)
        #epix.addDroplet(threshold=10., thresholdLow=3., thresADU=0.,name='droplet')
        ##save a histogram of ADU values in detector/event
        #if dirname.find('output')>=0:
        #    epix['droplet'].add_aduHist([50.,800.,1.])
        #else:
        #    epix['droplet'].add_aduHist([50.,600.,10.])
        ##same, but with ROI
        #epix['droplet'].add_aduHist([50.,600., 2],ROI=[59,538,250,675],name='ROI_')
        #epix['droplet'].add_aduHist([50.,600., 2],name='noROI_')
        #if dirname.find('output')>=0:
        #    epix['droplet'].addDropletSave(maxDroplets=35000)
        #epix['droplet'].checkADUHist()
        
        for iROI,ROI in enumerate(epixROI):
            epix.addROI('ROI_%d'%iROI, ROI)
            epix['ROI_%d'%iROI].addProj(axis=-1, thresADU=50.)            
        dets.append(epix)

        epix4 = DetObject(epixname ,ds.env(), int(run), name=epixname+'_cm4',common_mode=4)
        for iROI,ROI in enumerate(epixROI):
            epix4.addROI('ROI_%d'%iROI, ROI)
            epix4['ROI_%d'%iROI].addProj(axis=-1, thresADU=50.)
        dets.append(epix4)

        epix34 = DetObject(epixname ,ds.env(), int(run), name=epixname+'_cm34',common_mode=34)
        for iROI,ROI in enumerate(epixROI):
            epix34.addROI('ROI_%d'%iROI, ROI)
            epix34['ROI_%d'%iROI].addProj(axis=-1, thresADU=50.)
        dets.append(epix34)

        epix6 = DetObject(epixname ,ds.env(), int(run), name=epixname+'_cm6',common_mode=6)
        for iROI,ROI in enumerate(epixROI):
            epix6.addROI('ROI_%d'%iROI, ROI)
            epix6['ROI_%d'%iROI].addProj(axis=-1, thresADU=50.)
        dets.append(epix6)

        epixps = DetObject(epixname ,ds.env(), int(run), name=epixname+'_pedSub',common_mode=0)
        for iROI,ROI in enumerate(epixROI):
            epixps.addROI('ROI_%d'%iROI, ROI)
            epixps['ROI_%d'%iROI].addProj(axis=-1, thresADU=50.)
        dets.append(epixps)

        epixc = DetObject(epixname ,ds.env(), int(run), name=epixname+'_calib',common_mode=30)
        for iROI,ROI in enumerate(epixROI):
            epixc.addROI('ROI_%d'%iROI, ROI)
            epixc['ROI_%d'%iROI].addProj(axis=-1, thresADU=50.)
        dets.append(epixc)


azIntParams = getAzIntParams(run)
ROI_cspad = getCspadROIs(int(run))
haveCspad = checkDet(ds.env(), 'cspad')
if haveCspad:
    cspad = DetObject('cspad' ,ds.env(), int(run), name='cspad')
    for iROI,ROI in enumerate(ROI_cspad):
        cspad.addROI('ROI_%d'%iROI, ROI)
        #cspad.addROI('ROI_%d'%iROI, ROI, writeArea=True)

    #cspad.azav_eBeam=azIntParams['eBeam']
    #if azIntParams.has_key('cspad_center'):
    #    cspad.azav_center=azIntParams['cspad_center']
    #    cspad.azav_dis_to_sam=azIntParams['cspad_dis_to_sam']
    #    try:
    #        cspad.addAzAv(phiBins=11, Pplane=0)
    #    except:
    #        pass
    dets.append(cspad)

    cspad0 = DetObject('cspad' ,ds.env(), int(run), name='cspad_0',common_mode=0)
    for iROI,ROI in enumerate(ROI_cspad):
        cspad0.addROI('ROI_%d'%iROI, ROI)
    dets.append(cspad0)

    cspad1 = DetObject('cspad' ,ds.env(), int(run), name='cspad_1',common_mode=1)
    for iROI,ROI in enumerate(ROI_cspad):
        cspad1.addROI('ROI_%d'%iROI, ROI)
    dets.append(cspad1)

    cspad5 = DetObject('cspad' ,ds.env(), int(run), name='cspad_5',common_mode=5)
    for iROI,ROI in enumerate(ROI_cspad):
        cspad5.addROI('ROI_%d'%iROI, ROI)
    dets.append(cspad5)
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

    #here you can add any data you like: example is a product of the maximumof two area detectors.
    try:
        cspadMax = cspad.evt.dat.max()
        epix_vonHamosMax = epix_vonHamos.evt.dat.max()
        combDict = {'userValue': cspadMax*epix_vonHamosMax}
        smldata.event(combDict)
    except:
        pass

print 'rank %d on %s is finished'%(ds.rank, hostname)
smldata.save()
