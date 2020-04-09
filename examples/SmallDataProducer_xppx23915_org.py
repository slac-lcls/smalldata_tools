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
from smalldata_tools import defaultDetectors,epicsDetector,printMsg,detData,DetObject
from smalldata_tools import checkDet,getCfgOutput,getUserData,getUserEnvData,dropObject
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
from scipy.ndimage.measurements import center_of_mass

from IPython.core.debugger import Pdb
pdb = Pdb()

#det_monitor = 'zyla'
#det_diffpat = 'cs140_ladm'
det_monitor = 'zyla_1'
det_diffpat = 'cs140_front'

def getROIs(run):
    if run > 7:
        return [[[0,512],[0,512]]]
    else:
        #return [[[0,10],[0,10]]]
        return [[[0,2160],[0,2560]]]
        #return []

##########################################################
# run independent parameters
##########################################################
#event codes which signify no xray/laser
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw']
epicsPV = []
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
ROIs = getROIs(int(run))
have_monitor = checkDet(ds.env(), det_monitor)
if have_monitor:
    #monitor = DetObject(det_monitor ,ds.env(), int(run), name=det_monitor, common_mode=-1)
    monitor = DetObject(det_monitor ,ds.env(), int(run), name=det_monitor, common_mode=0)
    for iROI,ROI in enumerate(ROIs):
        print 'debug setup: added ROI ',ROI
        monitor.addROI('ROI%d'%iROI,ROI,writeArea=True)
        #monitor.addROI('ROI%d'%iROI,ROI)

# have_diffpat = checkDet(ds.env(), det_diffpat)
# if have_diffpat:
#     diffpat = DetObject(det_diffpat, ds.env(), int(run), name=det_diffpat)

##########################################################
##
## <-- User Input end
##
##########################################################
#dets = [ det for det in dets if checkDet(ds.env(), det._srcName)]
#for now require all area detectors in run to also be present in event.

defaultDets = defaultDetectors(hutch)
if len(epicsPV)>0:
    defaultDets.append(epicsDetector(PVlist=epicsPV, name='epicsUser'))

#add config data here
userDataCfg={}
try:
    userDataCfg[det_monitor]=getCfgOutput(monitor)
    # userDataCfg[det_diffpat]=getCfgOutput(diffpat)
    Config={'UserDataCfg':userDataCfg}
except:
    pass
smldata.save(Config)

#load sample image and define error function
names = ['r81_1410.npy', 'r81_579.npy']
try:
    master_file = [np.load('/reg/d/psdm/xcs/xcslq2215/results/smalldata_tools/'+x) for x in names]
    back_ground = 105
    # check if it is already cleaned
    for i, mf in enumerate(master_file):
        if master_file[i].min()>0:
            master_file[i] = np.maximum(0, master_file[i]-back_ground)
            master_file[i] /= 1.0 * np.sum( master_file[i] )
    smldata.save({"ReferenceFiles":names})
except:
    master_file = None

def error_function(image, master_file, background=1, max_noise=1, default_error=1E3):

    # Filter out to noisy images
    # got strange error "divide by zero", should not be possible
    noise_to_signal = np.sum(image) *1.0 / (np.sum( np.maximum(0,1.0*image-background) ) + 1)
    if (noise_to_signal > max_noise):
        return default_error #, noise_to_signal


    # subtract background and normalize
    cleaned_image = 1.0 * np.maximum(0.0, 1.0*image-background)
    cleaned_image /= np.sum( cleaned_image )

    difference = np.sum( (cleaned_image-master_file)**2 )
    return difference #, noise_to_signal

noise_signal = 0
for eventNr,evt in enumerate(ds.events()):
    printMsg(eventNr, evt.run(), ds.rank)

    if eventNr >= maxNevt/ds.size:
        break

    #add default data
    defData = detData(defaultDets, evt)

    #for key in defData.keys():
    #    print eventNr, key, defData[key]
    smldata.event(defData)

    #detector data using DetObject
    userDict = {}
    try:
        #this should be a plain dict. Really.
        monitor.evt = dropObject()
        monitor.getData(evt)
        monitor.processDetector()
        userDict[monitor._name]=getUserData(monitor)

        #there is only user data if a ROI or something else is set
        # diffpat.evt = dropObject()
        # diffpat.getData(evt)
        # diffpat.processDetector()
        # userDict[diffpat._name] = getUserData(diffpat)
    except:
        pass
    smldata.event(userDict)

    #monitor.evt.dat is monitor data
    #this is where the calculation goes. Here get the standard deviation of opal camera.
    # userStd = monitor.ROI0.area.std()
    # userSum = monitor.ROI0.area.sum()
    # userCOM0 = center_of_mass(monitor.ROI0.area)[0]
    # userCOM1 = center_of_mass(monitor.ROI0.area)[1]

    userError = []
    try:
        for i, mf in enumerate(master_file):
            userErr = np.float64(error_function(monitor.ROI0.area, mf, background=back_ground, max_noise=20))
            userError.append(userErr)
            #helper.append(help_noise)
        #noise_signal = max(noise_signal, max(helper))
    except:
        pass

    if userError == []:
        userError = [0]

    #print np.array(userError)[0]

    # cspad140_max = diffpat.evt.dat.max()
    # cspad140_sum = diffpat.evt.dat.sum()

    try:
        combDict = {'userEvent': eventNr}
                    # 'userStd': userStd,
                    # 'userSum': userSum,
                    # 'userCOM0': userCOM0,
                    # 'userCOM1': userCOM1,
                    # 'userNoise': np.araray(helper),
                    # 'cspadMax': cspad140_max,
                    # 'cspadSum': cspad140_sum
        for i, ErrorList in enumerate(userError):
            combDict["userError%d"%(i)] = ErrorList
        userValDict={'userValues':combDict}
        smldata.event(userValDict)
    except:
        pass

    try:
        img_monitor = {'ImageMonitor': {monitor._name: np.array(monitor.ROI0.area)}}
        smldata.event(img_monitor)
    except:
        pass

#gather whatever event did not make it to the gather interval
smldata.save()
