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
def getAzIntParams(run):
    if isinstance(run,basestring):
        run=int(run)
        
    ret_dict = {'eBeam': 9.5}
    ret_dict['epix10k_center'] = [82.20, -65.675]
    ret_dict['epix10k_dis_to_sam'] = 56.
    ret_dict['Pplane'] = 1
    return ret_dict

def getROI_epix10k(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <=6:
        return [ [[0,1], [1,74], [312,381]] ]
    else:
        return [ [[0,1], [1,74], [312,381]] ]

def getPed_epix10k(run):
    ped = None
    maxValid=-1
    fileValid=''
    for f in os.listdir('./pedestals/'):
        validRange=f.replace('.data','').split('-')
        if (validRange[1]=='end' or int(valueRange[1])>int(run)) and int(validRange[0]) > maxValid:
            maxValid = int(validRange[0])
            fileValid = f
    print('I will use pedestal /reg/d/psdm/xcs/xcslt5117/results/smalldata_tools_wEpix10k/pedestals/%s for run %d'%(fileValid,int(run)))
    ped = np.loadtxt('./pedestals/%s'%fileValid).reshape(16,352,384)
    return ped

def get_radav_mask(run):
    if isinstance(run,basestring):
        run=int(run)
    other_mask = np.loadtxt('Mask_AvImg_epix10ka2m_xcslt5117_Run065.data').reshape(16,352,384)
    try:
        cen_edge_mask = np.loadtxt('Mask_epix10ka2m_edge_center_2.data').reshape(16,352,384)
        #cen_edge_mask = np.loadtxt('Mask_epix10ka2m_edge_center_5.data').reshape(16,352,384)
        #cen_edge_mask = np.loadtxt('Mask_epix10ka2m_edge_center_x.data').reshape(16,352,384)
    except:
        cen_edge_mask = np.ones((16,352,384))
    if run <= 30:
        fname='Mask_AvImg_epix10ka2m_xcslt5117_Run016.data'
    else:
        fname='Mask_AvImg_pedSub_epix10ka2m_xcslt5117_Run064.data'
    radav_mask=np.loadtxt(fname).reshape(16,352,384)

    return radav_mask*cen_edge_mask#*other_mask

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
    for dirn in psana.__file__.split('/'):
        if dirn.find('ana-')>=0:
            version=dirn
    print 'Using psana version ',version

########################################################## 
##
## User Input start --> 
##
########################################################## 
dets=[]

azIntParams = getAzIntParams(run)

ROI_epix10k = getROI_epix10k(int(run))
#Ped_epix10k = getPed_epix10k(int(run))
radav_mask = get_radav_mask(run)
haveEpix10k = checkDet(ds.env(), 'epix10ka2m')
if haveEpix10k:
    #common_mode=-1 # raw
    #common_mode=0 # pedSub 
    #common_mode=81 # pedSub w/ gain
    common_mode=80 # somewhat more official version of calib. Similar to 81.
    epix10k = DetObject('epix10ka2m' ,ds.env(), int(run), name='epix10ka2m', common_mode=common_mode)
    #if Ped_epix10k is not None and (common_mode==0 or common_mode==81):
    #    epix10k.setPed(Ped_epix10k)
    for iROI,ROI in enumerate(ROI_epix10k):
        epix10k.addROI('ROI_%d'%iROI, ROI)

    ##saves the full detector in raw data shape
    #epix10k.saveFull()

    epix10k.azav_eBeam=azIntParams['eBeam']
    if azIntParams.has_key('epix10k_center'):
        epix10k.azav_center=azIntParams['epix10k_center']
        epix10k.azav_dis_to_sam=azIntParams['epix10k_dis_to_sam']
        epix10k.addAzAv(phiBins=11, Pplane=azIntParams['Pplane'])

        #epix10k.radav_eBeam=azIntParams['eBeam']
        #epix10k.radav_center=azIntParams['epix10k_center']
        #epix10k.radav_dis_to_sam=azIntParams['epix10k_dis_to_sam']
        #epix10k.addAzAv(phiBins=720, Pplane=azIntParams['Pplane'], userMask=radav_mask, azavName='radav', qBin=0.2)

        #other_mask = np.loadtxt('Mask_AvImg_epix10ka2m_xcslt5117_Run065.data').reshape(16,352,384)
        #cen_edge_mask = np.loadtxt('Mask_epix10ka2m_edge_center_2.data').reshape(16,352,384)
        #epix10k.radav2_eBeam=azIntParams['eBeam']
        #epix10k.radav2_center=azIntParams['epix10k_center']
        #epix10k.radav2_dis_to_sam=azIntParams['epix10k_dis_to_sam']
        #epix10k.addAzAv(phiBins=720, Pplane=azIntParams['Pplane'], userMask=other_mask*cen_edge_mask,azavName='radav2', qBin=0.2)

        #cen_edge_mask = np.loadtxt('Mask_epix10ka2m_edge5_center1.data').reshape(16,352,384)
        cen_edge_mask = np.loadtxt('Mask_epix10ka2m_edge3_center1.data').reshape(16,352,384)
        epix10k.radav3_eBeam=azIntParams['eBeam']
        epix10k.radav3_center=azIntParams['epix10k_center']
        epix10k.radav3_dis_to_sam=azIntParams['epix10k_dis_to_sam']
        epix10k.addAzAv(phiBins=720, Pplane=azIntParams['Pplane'], userMask=cen_edge_mask,azavName='radav3', qBin=0.2)

        epix10k.radav4_eBeam=azIntParams['eBeam']
        epix10k.radav4_center=azIntParams['epix10k_center']
        epix10k.radav4_dis_to_sam=azIntParams['epix10k_dis_to_sam']
        epix10k.addAzAv(phiBins=720, Pplane=azIntParams['Pplane'], azavName='radav4', qBin=0.2)

        ##pPlane IS 1 (as used above)
        #if int(run)==63:
        #    epix10k.radav4_eBeam=azIntParams['eBeam']
        #    epix10k.radav4_center=azIntParams['epix10k_center']
        #    epix10k.radav4_dis_to_sam=azIntParams['epix10k_dis_to_sam']
        #    if azIntParams['Pplane'] == 0:            
        #        epix10k.addAzAv(phiBins=720, Pplane=1, userMask=cen_edge_mask,azavName='radav4', qBin=0.2)
        #    else:
        #        epix10k.addAzAv(phiBins=720, Pplane=0, userMask=cen_edge_mask,azavName='radav4', qBin=0.2)
    epix10k.storeSum(sumAlgo='calib')
    #epix10k.storeSum(sumAlgo='square')
    dets.append(epix10k)

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

##adding wave8 traces:
#defaultDets.append(wave8Detector('Wave8WF'))
##adding raw timetool traces:
#try:
#    defaultDets.append(ttRawDetector(env=ds.env()))
#except:
#    pass

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

    #here you can add any data processing you like, just add a dictionay with the relevant information at the end
    try:
        mydata = epix10k.evt.dat
        combDict = {'userValue': np.nanmax(mydata)}
        smldata.event(combDict)
    except:
        pass

sumDict={'Sums': {}}
for det in dets:
    for key in det.storeSum().keys():
        sumData=smldata.sum(det.storeSum()[key])
        sumDict['Sums']['%s_%s'%(det._name, key)]=sumData
if len(sumDict['Sums'].keys())>0:
    smldata.save(sumDict)

print 'rank %d on %s is finished'%(ds.rank, hostname)
smldata.save()
