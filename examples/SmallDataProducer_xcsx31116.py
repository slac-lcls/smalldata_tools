# importing generic python modules
import numpy as np
import psana
import time
import argparse
import socket
import os
import RegDB.experiment_info
from smalldata_tools.DetObject import DetObject
from smalldata_tools.utilities import checkDet, printMsg
from smalldata_tools.SmallDataUtils import setParameter, getUserData, getUserEnvData, detData, defaultDetectors
from smalldata_tools.SmallDataDefaultDetector import ttRawDetector, wave8Detector, epicsDetector
from smalldata_tools.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, sparsifyFunc
from smalldata_tools.waveformFunc import getCMPeakFunc, templateFitFunc
from smalldata_tools.droplet import dropletFunc
from smalldata_tools.photons import photon
from smalldata_tools.read_uxi import getUxiDict
from smalldata_tools.azimuthalBinning import azimuthalBinning
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

def getAzIntParams_uxi(run):
    if isinstance(run,basestring):
        run=int(run)
    if run >= 96 and run <= 101:
        return [4650, 123000]
    elif run >= 102: #center seems to be about the same.
        return [4650, 123000]
    else:
        return []

def cm_Mask(run):
    if isinstance(run,basestring):
        run=int(run)

    if run >= 25 and run <= 41:
        return [2,52,2,512-2]
    elif run >= 42 and run <=50: #small slit open.
        return [1024-52,1024-2,2,512-2]
    elif run >= 51 and run <=59: #larger slit.
        return [1024-52,1024-2,2,512-2]
    elif run >= 60 and run <=76: #could not see beam.
        return [1024-52,1024-2,2,512-2]
    elif run >= 85 and run <=95: # in 88 I see beam in top of camera, bottom still ok
        return [1024-52,1024-2,2,512-2]
    elif run >= 96 and run <= 101:
        return [2,52,2,512-2]
    elif run >= 102:
        return [1024-52,1024-2,512-102,512-2]
    else:
        return None

def acqROI(run):
    if isinstance(run,basestring):
        run=int(run)

    ##25 has no signal.
    #if run > 96:
    #    return [0,1,2750,4250]
    ##have not seen any signal elsewhere....    
    return [0,1,3100,3500]
    #return [0,1,2750,4250]

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
    #dsname='exp='+expname+':run='+run+':smd:stream=1' #for run 109, most xtc look broken.
    dsname='exp='+expname+':run='+run+':smd' #for run 109, most xtc look broken.
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
jungfrauName = 'jungfrau1M'
acqirisName = 'acq01'

#True:
#run over all event, do NOT store jungfrau & uxi data, only uxi presence in data
#False
#save only uxi events, full information
#runAllEvents = True
runAllEvents = False

if runAllEvents and int(run)!=74:
    jungfrauName = 'jungfrau1MTest'

#acqirisName = 'acq01Test'

haveJungfrau = checkDet(ds.env(), jungfrauName)
if haveJungfrau:
    #jungfrau = DetObject.getDetObject(jungfrauName ,ds.env(), int(run), name=jungfrauName, common_mode=7)
    jungfrau = DetObject(jungfrauName ,ds.env(), int(run), name=jungfrauName, common_mode=7)
    jungfrau.addFunc(ROIFunc(writeArea=True))
    print 'appended the jungfrau'
    dets.append(jungfrau)

haveAcqiris = checkDet(ds.env(), acqirisName)
peakMain = None
fit_nPeaks=2
if haveAcqiris:
    #acqiris = DetObject.getDetObject(acqirisName ,ds.env(), int(run), name=acqirisName)
    acqiris = DetObject(acqirisName ,ds.env(), int(run), name=acqirisName)
    ROI=acqROI(run)
    fullArea=ROIFunc(writeArea=True, ROI=ROI)

    peakTemplate = tables.open_file('/reg/d/psdm/xcs/xcsx31116/results/smalldata_tools/SingleAcqirisPeak.h5').root.acq_pk2_shift
    peakMain = peakTemplate[50:200].copy()
    #peakTemplate = '/reg/d/psdm/xcs/xcsx35017/results/smalldata_tools/smalldata_tools/acqTemplate.npy'
    #peakTemplate = np.load(peakTemplate)
    #peakMain = peakTemplate[1460:1500].copy()
    peakMain = peakMain/peakMain.max()
    invert=False
    if int(run) <= 84: #def neg. 87 def pos.
        invert=True
    fullArea.addFunc(templateFitFunc(nPeaks=fit_nPeaks, fitMethod='sn_trf', name='sn_trf_%d'%fit_nPeaks, template=peakMain, invert=invert, baseline=[300,400]))
    fullArea.addFunc(templateFitFunc(nPeaks=fit_nPeaks, fitMethod='sn_lm', name='sn_lm_%d'%fit_nPeaks, template=peakMain, invert=invert, baseline=[300,400]))
    fullArea.addFunc(templateFitFunc(nPeaks=fit_nPeaks, fitMethod='sn_old', name='sn_old_%d'%fit_nPeaks, template=peakMain, invert=invert, baseline=[300,400]))
    fullArea.addFunc(templateFitFunc(nPeaks=fit_nPeaks, fitMethod='pah_trf', name='pah_trf_%d'%fit_nPeaks, template=peakMain, invert=invert, baseline=[300,400]))
    acqiris.addFunc(fullArea)
    dets.append(acqiris)

uxiDict, uxiConfigDict = getUxiDict(int(run))

raredets=[]
raredetSave=False
if uxiDict != {}:
    #-1: raw data
    #0 : pedestal subtracted
    #1 : subtract masked area value
    #2 : zero photon peak for each row
    #3 : first masked area, then row-based.
    cm_maskROI =  cm_Mask(run)
    #uxi = DetObject.getDetObject('uxi_standalone',ds.env(), int(run), common_mode=81, uxiDict=uxiDict, uxiConfigDict=uxiConfigDict, cm_minFrac=0.05, cm_photonThres=20, cm_maskedROI=cm_maskROI, cm_maskNeighbors=1)
    uxi = DetObject('uxi_standalone',ds.env(), int(run), common_mode=81, uxiDict=uxiDict, uxiConfigDict=uxiConfigDict, cm_minFrac=0.05, cm_photonThres=20, cm_maskedROI=cm_maskROI, cm_maskNeighbors=1)
    fullArea=ROIFunc(writeArea=True)
    specFunc_m100_300=spectrumFunc(name='spec_m100_200',bins=[-100,300,1.])
#    #specFunc_m100_200.addFunc(getCMPeakFunc(name='cmPk150',nPeakSel=4, minPeakNum=150))
#    #specFunc_m100_200.addFunc(getCMPeakFunc(name='cmPk1000',nPeakSel=4, minPeakNum=1000))
#    #specFunc_m100_200.addFunc(getCMPeakFunc(name='cmPk10000',nPeakSel=4, minPeakNum=10000))
#    #for tileROI in tileAreas:
#    #    tileROI.addFunc(specFunc_m100_200)
#    #    uxi.addFunc(tileROI)
    if not runAllEvents:
        fullArea.addFunc(specFunc_m100_300)
        uxi.addFunc(fullArea)
        raredetSave=True
        tileAreas=[ ROIFunc(ROI=[i,i+1,0,1e6], name='ROI_%d'%i) for i in range(uxi.ped.shape[0])]
        projection = projectionFunc(axis=1)
        center=getAzIntParams_uxi(run)
        if center!=[]:
            azav    = azimuthalBinning(center=[center[1], center[0]], dis_to_sam=150.6, qbin=1e-3, phiBins=50, name='azav', eBeam=7.2)
            azav102 = azimuthalBinning(center=[center[1], center[0]], dis_to_sam=102., qbin=1e-3, phiBins=50, name='azav_102', eBeam=7.2)
            azav8955 = azimuthalBinning(center=[center[1], center[0]], dis_to_sam=89.55, qbin=1e-3, phiBins=50, name='azav_8955', eBeam=7.2)

            azav1d_102 = azimuthalBinning(center=[center[1], center[0]], dis_to_sam=102., qbin=1e-3, name='azav_1d_102', eBeam=7.2)
            azav1d_103 = azimuthalBinning(center=[center[1], center[0]], dis_to_sam=103., qbin=1e-3, name='azav_1d_103', eBeam=7.2)
            azav1d_104 = azimuthalBinning(center=[center[1], center[0]], dis_to_sam=104., qbin=1e-3, name='azav_1d_104', eBeam=7.2)
            azav1d_101 = azimuthalBinning(center=[center[1], center[0]], dis_to_sam=101., qbin=1e-3, name='azav_1d_101', eBeam=7.2)
            azav1d_1025 = azimuthalBinning(center=[center[1], center[0]], dis_to_sam=102.5, qbin=1e-3, name='azav_1d_1025', eBeam=7.2)
        for tileROI in tileAreas:
            tileROI.addFunc(projection)
            if center!=[]:
                tileROI.addFunc(azav)
                tileROI.addFunc(azav102)
                tileROI.addFunc(azav8955)

                tileROI.addFunc(azav1d_102)
                tileROI.addFunc(azav1d_103)
                tileROI.addFunc(azav1d_101)
                tileROI.addFunc(azav1d_1025)
                tileROI.addFunc(azav1d_104)
            uxi.addFunc(tileROI)
    raredets.append(uxi)

    if not runAllEvents:
        uxiPed = DetObject('uxi_standalone',ds.env(), int(run), common_mode=0, uxiDict=uxiDict, uxiConfigDict=uxiConfigDict, name='uxiPed')
        #uxiPed = DetObject.getDetObject('uxi_standalone',ds.env(), int(run), common_mode=0, uxiDict=uxiDict, uxiConfigDict=uxiConfigDict, name='uxiPed')
        uxiPed.addFunc(fullArea)
        raredets.append(uxiPed)

        #if int(run)==103 or int(run)==75:
            #uxiNoMask = DetObject.getDetObject('uxi_standalone',ds.env(), int(run), common_mode=81, uxiDict=uxiDict, uxiConfigDict=uxiConfigDict, name='uxiNoMask', cm_minFrac=0.05, cm_photonThres=20, cm_maskNeighbors=1)
            #uxiNoMask.addFunc(fullArea)
            #raredets.append(uxiNoMask)

            #uxiMask = DetObject.getDetObject('uxi_standalone',ds.env(), int(run), common_mode=80, uxiDict=uxiDict, uxiConfigDict=uxiConfigDict, name='uxiMask', cm_minFrac=0.05, cm_photonThres=20, cm_maskedROI=cm_maskROI)
            #uxiMask.addFunc(fullArea)
            #raredets.append(uxiMask)

            #uxi50 = DetObject.getDetObject('uxi_standalone',ds.env(), int(run), common_mode=1, uxiDict=uxiDict, uxiConfigDict=uxiConfigDict, name='uxi50', cm_minFrac=0.05, cm_photonThres=50, cm_maskedROI=cm_maskROI, cm_maskNeighbors=1)
            #uxi50.addFunc(fullArea)
            #raredets.append(uxi50)


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
for det in raredets:
    userDataCfg[det._name] = det.params_as_dict()
Config={'UserDataCfg':userDataCfg}
smldata.save(Config)

for eventNr,evt in enumerate(ds.events()):
    printMsg(eventNr, evt.run(), ds.rank, ds.size)

    if eventNr >= maxNevt/ds.size:
        break

    #detector data using DetObject 
    rarepresentDict = {'damage':{}}
    rareDict = {}
    for det in raredets:
        rarepresentDict['damage'][det._name]=0
        try:
            det.getData(evt)
            if det.evt.dat is None:
                continue
            rarepresentDict['damage'][det._name]=1
            det.processFuncs()
            rareDict[det._name]=getUserData(det)
            #print rareDict[det._name]
            try:
                envData=getUserEnvData(det)
                if len(envData.keys())>0:
                    rareDict[det._name+'_env']=envData
            except:
                pass
            #print rareDict[det._name]
        except:
            pass
        
    haveRareData=False
    for savedRareDet,rareDetData in rareDict.iteritems():
        haveRareData=True
        smldata.event({savedRareDet: rareDetData})
        
    #print 'w ',raredetSave, not(haveRareData)
    if raredetSave and not(haveRareData):
        continue
                      
    #print 'ww ', len(raredets)
    if len(raredets)>0:
        #print 'rarepresentDict ',rarepresentDict
        smldata.event(rarepresentDict)

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


smldata.save()
print 'rank %d on %s is finished'%(ds.rank, hostname)
