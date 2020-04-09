# importing generic python modules
import numpy as np
import psana
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
    ret_dict['cspad_center'] = [87526.79161840, 92773.3296889500]
    ret_dict['cspad_dis_to_sam'] = 80.
    ret_dict['jungfrau_center'] = [87526.79161840, 92773.3296889500]
    ret_dict['jungfrau_dis_to_sam'] = 80.
    ret_dict['epix10k_center'] = [87526.79161840, 92773.3296889500]
    ret_dict['epix10k_dis_to_sam'] = 80.
    return ret_dict

def getROI_cspad(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <=6:
        return [ [[0,1], [1,74], [312,381]],
                 [[8,9], [8,89], [218,303]] ]
    else:
        return [ [[0,1], [1,74], [312,381]] ]

def getROI_jungfrau(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <=6:
        return [ [[0,1], [1,74], [312,381]] ]
    else:
        return [ [[0,1], [1,74], [312,381]] ]

def getROI_epix10k(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <=6:
        return [ [[0,1], [1,74], [312,381]] ]
    else:
        return [ [[0,1], [1,74], [312,381]] ]

def getPed_epix10k(run):
    ped = None
    if int(run) <= 60:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/14-end.data').reshape(16,352,384)
    elif int(run) <= 63:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/61-end.data').reshape(16,352,384)
    elif int(run) <= 72:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/64-end.data').reshape(16,352,384)
    elif int(run) <= 76:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/73-end.data').reshape(16,352,384)
    elif int(run) <= 95:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/77-end.data').reshape(16,352,384)
    elif int(run) <= 98:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/96-end.data').reshape(16,352,384)
    elif int(run) <= 106:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/99-end.data').reshape(16,352,384)
    elif int(run) <= 107:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/107-end.data').reshape(16,352,384)
    elif int(run) <= 110:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/108-end.data').reshape(16,352,384)
    elif int(run) <= 111:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/111-end.data').reshape(16,352,384)
    elif int(run) <= 233:
        ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/112-end.data').reshape(16,352,384)
    #else:# int(run) < 61:
    #    ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/254-end.data').reshape(16,352,384)
    #    #now recode Mikhail lookup....
    else:
        maxValid=-1
        for f in os.listdir('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/'):
            validRange=f.replace('.data','').split('-')
            if (validRange[1]=='end' or int(valueRange[1])>int(run)) and int(validRange[0]) > maxValid:
                maxValid = int(validRange[0])
                print('I will use pedestal /reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/%s for run %d'%(f,int(run)))
                ped = np.loadtxt('/reg/d/psdm/xcs/xcsx35617/results/smalldata_tools_wEpix10k/pedestals/%s'%f).reshape(16,352,384)
    return ped

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

ROI_cspad = getROI_cspad(int(run))
haveCspad = checkDet(ds.env(), 'cspad')
if haveCspad:
    cspad = DetObject('cspad' ,ds.env(), int(run), name='cspad')
    for iROI,ROI in enumerate(ROI_cspad):
        cspad.addROI('ROI_%d'%iROI, ROI)

    #saves the full detector in raw data shape
    #cspad.saveFull()

    cspad.azav_eBeam=azIntParams['eBeam']
    if azIntParams.has_key('cspad_center'):
        cspad.azav_center=azIntParams['cspad_center']
        cspad.azav_dis_to_sam=azIntParams['cspad_dis_to_sam']
        try:
            cspad.addAzAv(phiBins=11, Pplane=0)
        except:
            pass
    cspad.storeSum(sumAlgo='calib')
    cspad.storeSum(sumAlgo='square')
    dets.append(cspad)


ROI_jungfrau = getROI_jungfrau(int(run))
haveJungfrau = checkDet(ds.env(), 'jungfrau1M')
if int(run)==3:#debug cspad on this run
    haveJungfrau = False
if haveJungfrau:
    jungfrau = DetObject('jungfrau1M' ,ds.env(), int(run), name='jungfrau1M')
    for iROI,ROI in enumerate(ROI_jungfrau):
        jungfrau.addROI('ROI_%d'%iROI, ROI)

    jungfrau.azav_eBeam=azIntParams['eBeam']
    if azIntParams.has_key('jungfrau_center'):
        jungfrau.azav_center=azIntParams['jungfrau_center']
        jungfrau.azav_dis_to_sam=azIntParams['jungfrau_dis_to_sam']
        try:
            jungfrau.addAzAv(phiBins=11, Pplane=0)
        except:
            pass
    jungfrau.storeSum(sumAlgo='calib')
    jungfrau.storeSum(sumAlgo='square')
    dets.append(jungfrau)

Ped_epix10k = getPed_epix10k(int(run))
ROI_epix10k = getROI_epix10k(int(run))
haveEpix10k = checkDet(ds.env(), 'epix10ka2m')
###
# pedestal: we do not update pedestals yet. I might hack this later.
###
if haveEpix10k:
    #epix10k = DetObject('epix10ka2m' ,ds.env(), int(run), name='epix10ka2m', common_mode=-1)
    ##for iROI,ROI in enumerate(ROI_epix10k):
    ##    epix10k.addROI('ROI_%d'%iROI, ROI)

    #saves the full detector in raw data shape
    #epix10k.saveFull()

    ##epix10k.azav_eBeam=azIntParams['eBeam']
    ##if azIntParams.has_key('epix10k_center'):
    ##    epix10k.azav_center=azIntParams['epix10k_center']
    ##    epix10k.azav_dis_to_sam=azIntParams['epix10k_dis_to_sam']
    ##    try:
    ##        epix10k.addAzAv(phiBins=11, Pplane=0)
    ##    except:
    ##        pass

    #epix10k.storeSum(sumAlgo='calib')
    #epix10k.storeSum(sumAlgo='square')
    #dets.append(epix10k)

    ##ped subtracted.
    #epix10k_ps = DetObject('epix10ka2m' ,ds.env(), int(run), name='epix10ka2m_pedSub', common_mode=0)
    #if Ped_epix10k is not None:
    #    epix10k_ps.setPed(Ped_epix10k)
    #epix10k_ps.saveFull()
    #dets.append(epix10k_ps)

    #calib
    epix10k_calib = DetObject('epix10ka2m' ,ds.env(), int(run), name='epix10ka2m_calib', common_mode=80)
    #epix10k_calib.saveFull()
    dets.append(epix10k_calib)

    #ped subtracted, fixed gain "corrected"
    #epix10k_fg = DetObject('epix10ka2m' ,ds.env(), int(run), name='epix10ka2m_fixedGain', common_mode=81)
    #if Ped_epix10k is not None:
    #    epix10k_fg.setPed(Ped_epix10k)
    #epix10k_fg.saveFull()
    #dets.append(epix10k_fg)

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
    for k in userDict.keys():
        print k
        if isinstance(userDict[k], dict):
            for kk in userDict[k].keys():
                print kk, userDict[k][kk]
        else:
            print userDict[k].shape
    smldata.event(userDict)

    print 'test: ',epix10k_calib.evt.dat.shape
    smldata.event({'epd':epix10k_calib.evt.dat})
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
