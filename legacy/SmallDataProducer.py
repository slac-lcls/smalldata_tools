# importing generic python modules
import numpy as np
import psana
import time
import argparse
import socket
import os
import logging 
import requests

import sys
sys.path.append('/reg/g/xpp/xppcode/python/smalldata_tools/')

from smalldata_tools.utilities import printMsg
from smalldata_tools.SmallDataUtils import setParameter, defaultDetectors, detData
from smalldata_tools.SmallDataDefaultDetector import ttRawDetector, wave8Detector, epicsDetector


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

hostname=socket.gethostname()
if not args.run:
    run=raw_input("Run Number:\n")
else:
    run=args.run

ws_url = "https://pswww.slac.stanford.edu/ws/lgbk"
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
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
    lastRun = int(rundoc['num'])

    if (run < lastRun) or (run == lastRun and (not rundoc.get('end_time', None))):
        xtcdirname = '/reg/d/ffb/%s/%s/xtc'%(hutch.lower(),expname)
        xtcname=xtcdirname+'/*-r%04d-*'%int(run)
        import glob
        presentXtc=glob.glob('%s'%xtcname)
        if len(presentXtc)==0:
            dsname='exp='+expname+':run='+run+':smd'
    #this will fail when the data is not there after the timeout. On purpose to not fill up the queue.
else:
    expname=args.exp
    hutch=expname[0:3].upper()
    resp = requests.get(ws_url + "/lgbk/ws/activeexperiment_for_instrument_station", {"instrument_name": hutch, "station": 0})
    expnameCurr = resp.json().get("value", {}).get("name")
    dsname='exp='+expname+':run='+run+':smd'
    if expnameCurr == expname:
        dsname='exp='+expname+':run='+run+':smd'
        rundoc = requests.get(ws_url + "/lgbk/" + expname  + "/ws/current_run").json()["value"]
        if not rundoc:
            logger.error("Invalid response from server")
        lastRun = int(rundoc['num'])
        if (int(run) < int(lastRun)) or (int(run) == int(lastRun) and (not rundoc.get('end_time', None))):
            xtcdirname = '/reg/d/ffb/%s/%s/xtc'%(hutch.lower(),expname)
            xtcname=xtcdirname+'/e*-r%04d-*'%int(run)
            import glob
            presentXtc=glob.glob('%s'%xtcname)
            if len(presentXtc)>0:
                dsname='%s:dir=/reg/d/ffb/%s/%s/xtc'%(dsnam,hutch.lower(),expname)
                #check if we need live mode.
                for currXtc in presentXtc:
                    if currXtc.find('inprogress')>=0:
                        dsname=dsname+':live'
            else:
                xtcdirname = '/reg/d/psdm/%s/%s/xtc'%(hutch.lower(),expname)
                xtcname=xtcdirname+'/e*-r%04d-*'%int(run)
                files_there=False; nWait=0
                while not files_there:
                    presentXtc=glob.glob('%s'%xtcname)
                    if len(presentXtc)==0:
                        print('no files yet, wait')
                        sleep(30)
                        nWait=nWait+1
                    else:
                        files_there=True
                        for currXtc in presentXtc:
                            if currXtc.find('inprogress')>=0:
                                dsname=dsname+':live'
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

#for this, never wait for recorder streams....
#if args.norecorder:
#    dsname=dsname+':stream=0-79'
dsname=dsname+':stream=0-79'

debug = True
time_ev_sum = 0.
try:
    ds = psana.MPIDataSource(dsname)
except:
    import sys
    sys.exit()

try:    
    if dirname is None:
        dirname = '/reg/d/psdm/%s/%s/hdf5/smalldata'%(hutch.lower(),expname)
        #dirname = '/reg/d/psdm/%s/%s/results/arphdf5'%(hutch.lower(),expname)
    directory = os.path.dirname(dirname)
    #I think this is not actually working. Can't do it from script. Need to create first.
    if ds.rank==0 and not os.path.isdir(dirname):
        print('made directory for output files: %s'%dirname)
        os.mkdir(directory)

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

defaultDets = defaultDetectors(hutch)
epicsPV=[] #automatically read PVs from questionnaire/epicsArch file 
if len(epicsPV)>0:
    defaultDets.append(epicsDetector(PVlist=epicsPV, name='epicsUser'))
##adding raw timetool traces:
#try:
#    defaultDets.append(ttRawDetector(env=ds.env()))
#except:
#    pass

#add config data here
userDataCfg={}
#look for default data config?
for det in defaultDets:
    userDataCfg[det.name] = det.params_as_dict()
Config={'UserDataCfg':userDataCfg}
smldata.save(Config)

for eventNr,evt in enumerate(ds.events()):
    printMsg(eventNr, evt.run(), ds.rank, ds.size)

    if eventNr >= maxNevt/ds.size:
        break

    #add default data
    defData = detData(defaultDets, evt)
    #for key in defData.keys():
    #    print(eventNr, key, defData[key])
    smldata.event(defData)

print('rank %d on %s is finished'%(ds.rank, hostname))
smldata.save()
