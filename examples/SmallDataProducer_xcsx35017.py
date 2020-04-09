# importing generic python modules
import numpy as np
import psana
import time
import argparse
import socket
import os
import RegDB.experiment_info
from smalldata_tools import defaultDetectors,epicsDetector,printMsg,detData,DetObject
from smalldata_tools import checkDet,getCfgOutput,getUserData,getUserEnvData
from smalldata_tools import ttRawDetector,wave8Detector,setParameter
from smalldata_tools.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, sparsifyFunc
from smalldata_tools.waveformFunc import getCMPeakFunc, templateFitFunc
from smalldata_tools.droplet import dropletFunc
from smalldata_tools.photons import photon
#for reading the template for the acqiris fit.
import tables
########################################################## 
##
## User Input start --> 
##
########################################################## 
##########################################################
# functions for run dependant parameters
##########################################################d
def get_eventCodesToSave(run):
    codesToSave=[88]
    if run <=12:
        codesToSave.append(40)
    return  codesToSave

def getAzIntParams_jungfrau(run):
    #x,y:  414.709216692 -78256.9000939  R  71966.9906644 (this one!)
    #x,y:  3804.84897642 -77171.8796598  R  71387.0981453 (this one looks wierd as azav...)
    if isinstance(run,basestring):
        run=int(run)
        
    #jungfrau512k = by hand
    #x,y:  -44677.7513913 70700.0689637  R  55313.9486112
#x,y:  -43751.7088403 67815.1958565  R  67804.1963359

    if run <= 23: 
        ret_dict = {'eBeam': 9.5}
        ret_dict['cspad_center'] = [414.709216692, -78256.90]
        ret_dict['cspad_dis_to_sam'] = 600.
    else:
        ret_dict = {'eBeam': 9.5}
        ret_dict['cspad_center'] = [414.709216692, -78256.90]
        ret_dict['cspad_dis_to_sam'] = 600.

    return ret_dict

##########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw'] 
epicsPV = []
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
parser.add_argument("--peaks", help="number of peaks", type=int)
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
if args.peaks:
    peaks=args.peaks
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

#start by testing jungfrau stuff.
jungfrauName = 'jungfrau512kTest'
#jungfrauName = 'jungfrau512k'
haveJungfrau = checkDet(ds.env(), jungfrauName)
if haveJungfrau:
    jungfrau = DetObject.getDetObject(jungfrauName ,ds.env(), int(run), name=jungfrauName, common_mode=7)
    jungfrau.addFunc(ROIFunc(writeArea=True))
    dets.append(jungfrau)

#uxinames=['icarus_pink','icarus_yellow']
#uxinames=['icarus_pink']
uxinames=['icarus_pinkTest']
haveAnyUxi=False
for uxiname in uxinames:
    have_uxi = checkDet(ds.env(), uxiname)
    print 'have uxi: ',have_uxi
    if have_uxi:
        uxi = DetObject.getDetObject(uxiname ,ds.env(), int(run), common_mode=98) #subtract unb pixel, 98 just saves info
        #uxi = DetObject(uxiname ,ds.env(), int(run), common_mode=0) #subtract pedestal if possible.
        fullArea=ROIFunc(writeArea=True)
        specFunc_m100_200=spectrumFunc(name='spec_m100_200',bins=[-100,200,1.])
        specFunc_m100_200.addFunc(getCMPeakFunc(name='cmPk150',nPeakSel=4, minPeakNum=150))
        specFunc_m100_200.addFunc(getCMPeakFunc(name='cmPk1000',nPeakSel=4, minPeakNum=1000))
        specFunc_m100_200.addFunc(getCMPeakFunc(name='cmPk10000',nPeakSel=4, minPeakNum=10000))
        tileAreas=[ ROIFunc(ROI=[i,i+1,0,1e6], name='ROI_%d'%i) for i in range(uxi.ped.shape[0])]
        for tileROI in tileAreas:
            tileROI.addFunc(specFunc_m100_200)
            uxi.addFunc(tileROI)
        #fullArea.addFunc(specFunc_m100_200)
        uxi.addFunc(fullArea)
        
        dets.append(uxi)
        haveAnyUxi=True

acqirisName = 'acq01'
#acqirisName = 'acq01Test'
haveAcqiris = checkDet(ds.env(), acqirisName)
peakMain = None
fit_nPeaks=4
if haveAcqiris:
    acqiris = DetObject.getDetObject(acqirisName ,ds.env(), int(run), name=acqirisName)
    fullArea=ROIFunc(writeArea=True)

    #fit_nPeak=peaks #commandline argument. Should be argument in function def.
    peakTemplate = tables.open_file('/reg/d/psdm/xcs/xcsx31116/results/smalldata_tools/SingleAcqirisPeak.h5').root.acq_pk2_shift
    peakMain = peakTemplate[50:200].copy()
    peakTemplate = '/reg/d/psdm/xcs/xcsx35017/results/smalldata_tools/smalldata_tools/acqTemplate.npy'
    peakTemplate = np.load(peakTemplate)
    peakMain = peakTemplate[1460:1500].copy()
    peakMain = peakMain/peakMain.max()
    fullArea.addFunc(templateFitFunc(nPeaks=fit_nPeaks, fitMethod='sn_trf', name='sn_trf_%d'%fit_nPeaks, template=peakMain))
    fullArea.addFunc(templateFitFunc(nPeaks=fit_nPeaks, fitMethod='sn_lm', name='sn_lm_%d'%fit_nPeaks, template=peakMain))
    fullArea.addFunc(templateFitFunc(nPeaks=fit_nPeaks, fitMethod='sn_old', name='sn_old_%d'%fit_nPeaks, template=peakMain))
    fullArea.addFunc(templateFitFunc(nPeaks=fit_nPeaks, fitMethod='pah_trf', name='pah_trf_%d'%fit_nPeaks, template=peakMain))
    acqiris.addFunc(fullArea)
    dets.append(acqiris)

##adding raw timetoDetOol traces:
#defaultDets.append(ttRawDetector(env=ds.env()))
##adding wave8 traces:
#defaultDets.append(wave8Detector('Wave8WF'))
########################################################## 
##
## <-- User Input end
##
########################################################## 
dets = [ det for det in dets if checkDet(ds.env(), det.det.alias)]
#for now require all area detectors in run to also be present in event.

defaultDets = defaultDetectors(hutch)
if len(ttCalib)>0:
    setParameter(defaultDets, ttCalib)
if len(aioParams)>0:
    setParameter(defaultDets, aioParams, 'ai')
if len(epicsPV)>0:
    defaultDets.append(epicsDetector(PVlist=epicsPV, name='epicsUser'))

#OK, this needs debugging.
#add config data here
userDataCfg={}
for det in dets:
    userDataCfg[det._name] = det.params_as_dict()
Config={'UserDataCfg':userDataCfg}
smldata.save(Config)

#code to oinly select certain event codes for production
evr=psana.Detector('evr0')
codesToSave = get_eventCodesToSave(int(run))
if not haveAnyUxi:
    codesToSave.append(40)

for eventNr,evt in enumerate(ds.events()):
    printMsg(eventNr, evt.run(), ds.rank, ds.size)

    if eventNr >= maxNevt/ds.size:
        break

    skipEvent=True
    evtCodes = evr.eventCodes(evt)
    if evtCodes is not None:
        for code in codesToSave:
            if code in evtCodes:
                skipEvent=False
    if skipEvent and expname=='xcsx35017' and  int(run)!=67:
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
            det.getData(evt)
            det.processFuncs()
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

    ##here you can add any data you like: example is a product of the maximumof two area detectors.
    ##try:
    #trace = acqiris.evt.dat
    ##do stuff.
    ##    epix_vonHamosMax = epix_vonHamos.evt.dat.max()
    #try:
    #    acqDict = {'acqMax': trace.max()}
    #except:
    #    acqDict = {'acqMax': -666} ## missing trace in run 113 somewhere after event 900
    #smldata.event(acqDict)
    ##except:
    ##    pass

smldata.save()
print 'rank %d on %s is finished'%(ds.rank, hostname)
