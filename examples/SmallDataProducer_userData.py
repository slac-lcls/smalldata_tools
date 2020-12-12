# importing generic python modules
import numpy as np
import psana
import time
import argparse
import socket
import os
import logging 
import requests

from smalldata_tools.DetObject import DetObject
from smalldata_tools.utilities import checkDet, printMsg
from smalldata_tools.SmallDataUtils import setParameter, getUserData, getUserEnvData, detData, defaultDetectors
from smalldata_tools.SmallDataDefaultDetector import ttRawDetector, wave8Detector, epicsDetector
from smalldata_tools.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, sparsifyFunc
from smalldata_tools.waveformFunc import getCMPeakFunc, templateFitFunc
from smalldata_tools.droplet import dropletFunc
from smalldata_tools.photons import photonFunc
from smalldata_tools.azimuthalBinning import azimuthalBinning

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
    ret_dict['epix10k2M_center'] = [87526.79161840, 92773.3296889500]
    ret_dict['epix10k2M_dis_to_sam'] = 80.
    return ret_dict

def getROI_epix10k2M(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <=6:
        return [ [[0,1], [1,74], [312,381]],
                 [[8,9], [8,89], [218,303]] ]
    else:
        return [ [[0,1], [1,74], [312,381]],
                 [[8,9], [8,89], [218,303]] ]

def getROI_rowland(run):
    if isinstance(run,basestring):
        run=int(run)

    if run <= 6:
        return [[[0,1], [25, 275], [516, 556]], 
                [[0,1], [25, 275], [460, 500]]]
    else:
        return [[[0,1], [25, 275], [516, 556]], 
                [[0,1], [25, 275], [460, 500]]]

def getNmaxDrop(run):
    if isinstance(run,basestring):
        run=int(run)

    if run >= 10:
        return 2000
    else:
        return 400

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
parser.add_argument("--norecorder", help="ignore recorder streams", action='store_true')
args = parser.parse_args()

ws_url = "https://pswww.slac.stanford.edu/ws/lgbk"
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
        print('cannot figure out which experiment to use, please specify -e <expname> on commandline')
        import sys
        sys.exit()
    if hutch == 'cxi':
        print('Will assume the first station, if this is wrong, please  -e <expname> on commandline')    
    resp = requests.get(ws_url + "/lgbk/ws/activeexperiment_for_instrument_station", {"instrument_name": hutch.upper(), "station": 0})
    expname = resp.json().get("value", {}).get("name")
    dsname='exp='+expname+':run='+run+':smd:dir=/reg/d/ffb/%s/%s/xtc:live'%(hutch.lower(),expname)
    #data gets removed from ffb faster now, please check if data is still available
    rundoc = requests.get(ws_url + "/lgbk/" + expname  + "/ws/current_run").json()["value"]
    if not rundoc:
        logger.error("Invalid response from server")
    isLive = rundoc.get('end_time', None)
    if not isLive:
        xtcdirname = '/reg/d/ffb/%s/%s/xtc'%(hutch.lower(),expname)
        xtcname=xtcdirname+'/*-r%04d-*'%int(run)
        import glob
        presentXtc=glob.glob('%s'%xtcname)
        if len(presentXtc)==0:
            dsname='exp='+expname+':run='+run+':smd'
else:
    expname=args.exp
    hutch=expname[0:3]
    dsname='exp='+expname+':run='+run+':smd'

    resp = requests.get(ws_url + "/lgbk/ws/activeexperiment_for_instrument_station", {"instrument_name": hutch.upper(), "station": 0})
    expnameC = resp.json().get("value", {}).get("name")
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
    ds = psana.MPIDataSource(str(dsname))
except:
    print('failed to make MPIDataSource for ',dsname)
    import sys
    sys.exit()

try:    
    if dirname is None:
        dirname = '/reg/d/psdm/%s/%s/hdf5/smalldata'%(hutch.lower(),expname)
    smldataFile = '%s/%s_Run%03d.h5'%(dirname,expname,int(run))

    smldata = ds.small_data(smldataFile,gather_interval=gatherInterval)

except:
    print('failed making the output file ',smldataFile)
    import sys
    sys.exit()

if ds.rank==0:
    version='unable to detect psana version'
    for dirn in psana.__file__.split('/'):
        if dirn.find('ana-')>=0:
            version=dirn
    print('Using psana version ',version)


########################################################## 
##
## Setting up the default detectors
## needs to be before the user detectors only for epix10k 
## data that needs ghost corrections and uses a psana 
## detector that does not take the event as inpu (EPICS PV, tt)
##
########################################################## 
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

########################################################## 
##
## User Input start --> 
##
########################################################## 
dets=[]

epixname='epix_diff'
nDrop = getNmaxDrop(int(run))
have_epix = checkDet(ds.env(), epixname)
if have_epix:
    #create detector object. needs run for calibration data
    #common mode: 46 is the "original" to psana method 7(?)
    #row & column correction on data w/ photons&neighbors removed.
    det = DetObject(epixname ,ds.env(), int(run), common_mode=46)

    #two threshold droplet finding.
    #for data w/ > 1 photon energy this is the only thing that will work.
    #Tends to add photons together into single droplet if occupancy
    #is not low, might need photonizing step to get single photon positions
    droplet = dropletFunc(threshold=10., thresholdLow=3., thresADU=0.,name='droplet')
    #spectrum of found droplets
    specFunc_300=spectrumFunc(name='spec_300',bins=[0,1000,2.5])
    droplet.addFunc(specFunc_300) 
    #save droplets as rectangular array.
    droplet.addFunc(sparsifyFunc(nData=10000))
    save droplets in image format - not working yet for epix10k: 
        might need geometry or hardcode stuff.
    droplet.addFunc(imageFunc())
    det.addFunc(droplet)

    #now add photon algorithms. Only works for single energy photon data
    # ADU_per_photon: expected ADU for photon of expected energy
    # thresADU: fraction of photon energy in photon candidate
                #(2 neighboring pixel)
    #retImg: 0 (def): only histogram of 0,1,2,3,...,24 photons/pixel is returned
    #        1 return Nphot, x, y arrays
    #        2 store image using photons /event
    photon = photonFunc(ADU_per_photon=10, thresADU=0.9, name='photon')
    #save number/photons per pixel
    nPhotPP=spectrumFunc(name='nPhotPP',bins=[0,30,1.])
    photon.addFunc(nPhotPP)
    #save sparsified array
    sparse = sparsifyFunc(nData=15000)
    photon.addFunc(sparse)
    #save as full area
    fullImg = ROIFunc(writeArea=True)
    photon.addFunc(fullImg)
    det.addFunc(photon)

    dets.append(det)

ROI_rowland = getROI_rowland(int(run))
if checkDet(ds.env(), 'cs140_diff'):
    cs140 = DetObject('cs140_diff' ,ds.env(), int(run))#, name='Rowland')
    for iROI,ROI in enumerate(ROI_rowland):
        cs140.addFunc(ROIFunc(ROI=ROI, name='ROI_%d'%iROI))
    dets.append(cs140)

azIntParams = getAzIntParams(run)
ROI_epix10k2M = getROI_epix10k2M(int(run))
haveEpix10k2M = checkDet(ds.env(), 'epix10k2M')
if haveEpix10k2M:
    epix10k2M = DetObject('epix10k2M' ,ds.env(), int(run), name='epix10k2M')
    for iROI,ROI in enumerate(ROI_epix10k2M):
        epix10k2M.addFunc(ROIFunc(ROI=ROI, name='ROI_%d'%iROI))

    epix10k2M.azav_eBeam=azIntParams['eBeam']
    if azIntParams.has_key('epix10k2M_center'):
        try:
            azav = azimuthalBinning(center=azIntParams['epix10k2M_center'], dis_to_sam=azIntParams['epix10k2M_dis_to_sam'], phiBins=11, Pplane=0)
            epix10k2M.addFunc(azav)
        except:
	        pass


    epix10k2M.storeSum(sumAlgo='calib')
    epix10k2M.storeSum(sumAlgo='square')
    dets.append(epix10k2M)

########################################################## 
##
## <-- User Input end
##
########################################################## 
#dets = [ det for det in dets if checkDet(ds.env(), det._src.__repr__())]
#for now require all area detectors in run to also be present in event.


#add config data here
userDataCfg={}
for det in defaultDets:
    userDataCfg[det.name] = det.params_as_dict()
for det in dets:
    userDataCfg[det._name] = det.params_as_dict()
#for det in raredets:
#    userDataCfg[det._name] = det.params_as_dict()
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
            det.getData(evt)
            det.processFuncs()
            userDict[det._name]=getUserData(det)
            try:
                envData=getUserEnvData(det)
                if len(envData.keys())>0:
                    userDict[det._name+'_env']=envData
            except:
                pass
            det.processSums()
            #print userDict[det._name]
        except:
            pass
    smldata.event(userDict)

    #here you can add any data you like: example is a product of the maximumof two area detectors.
    #try:
    #    epix10k2MMax = epix10k2M.evt.dat.max()
    #    epix_vonHamosMax = epix_vonHamos.evt.dat.max()
    #    combDict = {'userValue': epix10k2MMax*epix_vonHamosMax}
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

print('rank %d on %s is finished'%(ds.rank, hostname))
smldata.save()
