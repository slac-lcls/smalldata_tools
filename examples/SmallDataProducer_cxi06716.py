# importing generic python modules
import numpy as np
import h5py
import psana
import time
import argparse
import socket
import os
import sys
import RegDB.experiment_info
from smalldata_tools.DetObject import DetObject
from smalldata_tools.utilities import checkDet, printMsg
from smalldata_tools.SmallDataUtils import setParameter, getUserData, getUserEnvData, detData, defaultDetectors
from smalldata_tools.SmallDataDefaultDetector import epicsDetector
from smalldata_tools.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, sparsifyFunc, imageFunc
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
det_diffpat = 'Dg3CsPad2x2'

def error_function(image, master_file, background=1, max_noise=1, default_error=1E3):
    # Filter out to noisy images
    # got strange error "divide by zero", should not be possible
    noise_to_signal = np.nansum(image) *1.0 / (np.nansum( np.maximum(0,1.0*image-background) ) + 1)
    if (noise_to_signal > max_noise):
        return default_error #, noise_to_signal


    # subtract background and normalize
    cleaned_image = 1.0 * np.maximum(0.0, 1.0*image-background)
    cleaned_image /= np.nansum( cleaned_image )

    difference = np.nansum( (cleaned_image-master_file)**2 )
    return difference #, noise_to_signal

##########################################################
# run independent parameters
##########################################################
#event codes which signify no xray/laser
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw']
epicsPV = ['Sc1Sample_x','sc1diode_x']
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
dets=[]

have_diffpat = checkDet(ds.env(), det_diffpat)
if have_diffpat:
    #diffpat = DetObject(det_diffpat, ds.env(), int(run), name=det_diffpat, common_mode=-1)
    diffpat = DetObject(det_diffpat, ds.env(), int(run), name=det_diffpat, common_mode=0)
    fullROI_write = ROIFunc(writeArea=True)
    fullROI = ROIFunc()
    fullROI.addFunc(imageFunc(coords=['x','y']))
    #diffpat.addFunc(fullROI)
    diffpat.addFunc(fullROI_write)
    dets.append(diffpat)

##########################################################
##
## <-- User Input end
##
##########################################################
dets = [ det for det in dets if checkDet(ds.env(), det.det.alias)]

defaultDets = defaultDetectors(hutch)
if len(epicsPV)>0:
    defaultDets.append(epicsDetector(PVlist=epicsPV, name='epicsUser'))

#add config data here
userDataCfg={}
for det in defaultDets:
    userDataCfg[det.name] = det.params_as_dict()
for det in dets:
    userDataCfg[det._name] = det.params_as_dict()
Config={'UserDataCfg':userDataCfg}
smldata.save(Config)


noise_signal = 0
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
            #print userDict[det._name]
        except:
            pass
    smldata.event(userDict)

#gather whatever event did not make it to the gather interval
smldata.save()
