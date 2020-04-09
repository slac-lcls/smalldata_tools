# importing generic python modules
import numpy as np
import psana
import time
import argparse
import socket
import os
import sys
import scipy
import tables
import RegDB.experiment_info

#sys.path.append('/reg/g/xpp/xppcode/python/smalldata_tools/')

from smalldata_tools import defaultDetectors,epicsDetector,printMsg,detData,DetObject
from smalldata_tools import checkDet,getCfgOutput,getUserData,getUserEnvData,dropObject
from read_uxi import read_uxi, get_uxi_timestamps, getDarks, getUnbonded, getZeroPeakFit
from smalldata_tools import findPeakSimple

def fillUxiDict(fname, xtcBeginTime, xtcEndTime, uxiDict, debug=False):
    if not os.path.isfile(fname): return
    dataDictTS = {}
    #read the uxi file, get dict of frame fids&evttimes & images & config.
    t_preUxiRead=time.time()
    dataDictTS['uxi'] = get_uxi_timestamps(fname)    
    
    print 'UXI file ',fname, dataDictTS['uxi']['lcls_ts_secs']
    minUxi = int(min(dataDictTS['uxi']['lcls_ts_secs']))
    maxUxi = int(max(dataDictTS['uxi']['lcls_ts_secs']))
    if debug: print 'min/max time: ',minUxi, maxUxi
    
    configDict=None
    #if uxi event time stamps lie in uxi range, get all uxi data
    if minUxi<=xtcEndTime and maxUxi>=xtcBeginTime:
        inXtc = np.array([(int(ts) >= xtcBeginTime and int(ts) <= xtcEndTime) for ts in dataDictTS['uxi']['lcls_ts_secs']])
        uxiDictAll,configDict = read_uxi(fname, returnConfig=True)    
        for key in uxiDictAll.keys():
            if key not in uxiDict.keys():
                uxiDict[key] = np.array((uxiDictAll[key]))[inXtc]
            else:
                uxiDict[key] = np.append(uxiDict[key], np.array((uxiDictAll[key]))[inXtc], axis=0)
        if debug: print 'check lengths: ',len(uxiDictAll['lcls_ts_secs']),' inXtc: ',inXtc.sum(),' out dict ',len(uxiDict['lcls_ts_secs'])
    return configDict

def findPars(trace, nPeak=2):
    difftrace = np.diff(trace)
    peakIdx = scipy.signal.find_peaks_cwt(difftrace, widths=[3])
    maxIdx = [ x for _,x in sorted(zip(difftrace[peakIdx],peakIdx), key=lambda pair: pair[0], reverse=True)]
    peakPos=[];peakVal=[]
    for iPk in range(nPeak):
       idx=maxIdx[iPk]
       while difftrace[idx]>0:
            idx+=1
       peakPos.append(idx)
       peakVal.append(trace[idx])
    #sort in peak position.
    pPos = [ pos for pos,_ in sorted(zip(peakPos,peakVal))]
    pPk = [ pk for _,pk in sorted(zip(peakPos,peakVal))]
    for pVal in pPk:
        pPos.append(pVal)
    return pPos
    #for pVal in peakVal:
    #    peakPos.append(pVal)
    #return peakPos

#add template for peak w/ peak shifted by 1 with fraction of peak position.
#now the optimization actually works!
def templateArray(args, templateShape, peakMain):
    if peakMain is None: return None
    template = peakMain#[10:110]
    templateMaxPos = np.argmax(template)
    pk1Pos = int(args[0])
    if (pk1Pos>templateMaxPos):
        templatePk1 = np.append(np.zeros(max(0,pk1Pos-templateMaxPos)), template)
    else:
        templatePk1 = template
    if (templateShape-templatePk1.shape[0])>0:
        templatePk1 = np.append(templatePk1, np.zeros(templateShape-templatePk1.shape[0]))
    elif (templateShape-templatePk1.shape[0])<0:
        templatePk1 = templatePk1[:templateShape]
    templatePkp1 = np.append(np.array([0]), templatePk1[:-1])

    pk2Pos = int(args[1])
    templatePk2 = np.append(np.zeros(max(0,pk2Pos-templateMaxPos)), template)
    if (templateShape-templatePk2.shape[0])>0:
        templatePk2 = np.append(templatePk2, np.zeros(templateShape-templatePk2.shape[0]))
    elif (templateShape-templatePk2.shape[0])<0:
        templatePk2 = templatePk2[:templateShape]
    templatePkp2 = np.append(np.array([0]), templatePk2[:-1])
    frac1 = args[0]-int(args[0])
    frac2 = args[1]-int(args[1])
    template1 = templatePk1*(1.-frac1)+templatePkp1*frac1
    template2 = templatePk2*(1.-frac2)+templatePkp2*frac2
    templateSum=template1*args[2]+template2*args[3]
    return templateSum

#make sure to exclude points > max from calc!
def fitTemplateLeastsq(trace, peakMain, saturationFrac=0.98, debug=False):
    args0 = findPars(trace)
    maxTrace = np.nanmax(trace)*saturationFrac
    if debug: print 'maxtrace: ',maxTrace,' n high pix ',(trace>maxTrace).sum() 
    if (trace>maxTrace).sum() > 2:
        maskSaturation=[trace<maxTrace]
    else:
        maskSaturation=np.ones_like(trace).astype(bool)
    errorfunction = lambda p: templateArray(p, trace.shape[0], peakMain)[maskSaturation]-trace[maskSaturation]
    p, success = scipy.optimize.leastsq(errorfunction, args0)
    return p

def templateSingleSum(fitPk, fitScale, templateShape):
    template = peakMain#[10:110]
    templateMaxPos = np.argmax(template)
    templatePk1 = np.append(np.zeros(fitPk-templateMaxPos), template)
    if (templateShape-templatePk1.shape[0])>0:
        templatePk1 = np.append(templatePk1, np.zeros(templateShape-templatePk1.shape[0]))
    templatePkp1 = np.append(np.array([0]), templatePk1[:-1])
    frac1 = fitPk-int(fitPk)
    template1 = templatePk1*(1.-frac1)+templatePkp1*frac1
    templateSum=(template1*fitScale).sum()
    return templateSum

def residuals(ar1, ar2, saturationFrac=0.98, error=9e-4):
    resid=0.
    maxTrace = np.nanmax(ar1)*saturationFrac
    if (ar1>maxTrace).sum() > 2:
        mask=np.array(ar1<maxTrace).astype(bool)
        #print 'have saturation! ', mask.shape
    else:
        mask=np.ones_like(ar1).astype(bool)
    ndof=0
    #print 'arresid: ', ar1.shape, ar2.shape, mask.shape
    for r1,r2,m in zip(ar1, ar2, mask):
        #print 'r1, r2. resid :', r1, r2, resid, ndof
        diffSq = (r1-r2)/error*(r1-r2)/error*m
        if m != 0:
            resid+=diffSq
            ndof+=1

    #print 'resid', resid, ndof    
    return resid/(ndof-1)


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
        sys.exit()
    expname=RegDB.experiment_info.active_experiment(hutch)[1]
    dsname='exp='+expname+':run='+run+':smd:dir=/reg/d/ffb/%s/%s/xtc:live'%(hutch.lower(),expname)
    #data gets removed from ffb faster now, please check if data is still available
    lastRun = RegDB.experiment_info.experiment_runs(hutch)[-1]['num'] 
    if (run < lastRun) or (run == lastRun and (RegDB.experiment_info.experiment_runs(hutch)[-1]['end_time_unix'] is not None)):
        xtcdirname = '/reg/d/ffb/%s/%s/xtc'%(hutch.lower(),expname)
        xtcname=xtcdirname+'/e*-r%04d-*'%int(run)
        import glob
        presentXtc=glob.glob('%s'%xtcname)
        if len(presentXtc)==0:
            dsname='exp='+expname+':run='+run+':smd'
    #this will fail when the data is not there after the timeout. On purpose to not fill up the queue.
else:
    expname=args.exp
    hutch=expname[0:3].upper()
    expnameCurr=RegDB.experiment_info.active_experiment(hutch)[1]
    dsname='exp='+expname+':run='+run+':smd'
    if expnameCurr == expname:
        dsname='exp='+expname+':run='+run+':smd'
        lastRun = RegDB.experiment_info.experiment_runs(hutch)[-1]['num'] 
        if (int(run) < int(lastRun)) or (int(run) == int(lastRun) and (RegDB.experiment_info.experiment_runs(hutch)[-1]['end_time_unix'] is not None)):
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
                        print 'no files yet, wait'
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
    sys.exit()

try:    
    if dirname is None:
        dirname = '/reg/d/psdm/%s/%s/hdf5/smalldata'%(hutch.lower(),expname)
        #dirname = '/reg/d/psdm/%s/%s/results/arphdf5'%(hutch.lower(),expname)
    directory = os.path.dirname(dirname)

    smldataFile = '%s/%s_Run%03d.h5'%(dirname,expname,int(run))
    smldata = ds.small_data(smldataFile,gather_interval=gatherInterval)

    #I think this is not actually working. Can't do it from script. Need to create first.
    if ds.rank==0 and not os.path.isdir(dirname):
        print 'made directory for output files: %s'%dirname
        os.mkdir(directory)

except:
    print 'failed making the output file ',smldataFile
    sys.exit()

if ds.rank==0:
    version='unable to detect psana version'
    for dirn in psana.__file__.split('/'):
        if dirn.find('ana-')>=0:
            version=dirn
    print 'Using psana version ',version

defaultDets = defaultDetectors(hutch)
epicsPV=[] #automatically read PVs from questionnaire/epicsArch file 
if len(epicsPV)>0:
    defaultDets.append(epicsDetector(PVlist=epicsPV, name='epicsUser'))

dets=[]
jungfrauName = 'jungfrau1M'
haveJungfrau = checkDet(ds.env(), jungfrauName)
if haveJungfrau:
    jungfrau = DetObject(jungfrauName ,ds.env(), int(run), name=jungfrauName, common_mode=7)
    jungfrau.addROI('ROI', [[0,2], [0,512], [0,1024]], writeArea=True) #fix ROI: full mode?!
    dets.append(jungfrau)#TURN BACK ON< OFF FOR DEBUG


acqirisName = 'acq01'
haveAcqiris = checkDet(ds.env(), acqirisName)
peakMain = None
if haveAcqiris:
    acqiris = DetObject(acqirisName ,ds.env(), int(run), name=acqirisName)
    acqiris.addROI('ROI', [0,1e6], writeArea=True)
    dets.append(acqiris)

    peakTemplate = tables.open_file('/reg/d/psdm/xcs/xcsx31116/results/smalldata_tools/SingleAcqirisPeak.h5').root.acq_pk2_shift
    peakMain = peakTemplate[50:200].copy()
    peakMain = peakMain/peakMain.max()

##
# get timestamp range for xtc file (???)
##
runsDict = RegDB.experiment_info.experiment_runs(hutch,expname)
xtcBeginTime=None
xtcEndTime=None
for thisRun in runsDict:
    if thisRun['num']==int(run):
        xtcBeginTime = int(thisRun['begin_time_unix'])
        xtcEndTime = int(thisRun['end_time_unix'])
##
# put some UXI stuff here.
##
uxiPath='/reg/d/psdm/xcs/xcsx31116/scratch/uxi/'
#start with actual run first.
fname='uxi_x311_%04d.dat'%int(run)
fnamePre='uxi_x311_%04d.dat'%(int(run)-1)
if int(run) == 19: #yes, this is special....
    fnamePre='uxi_x311_%04d.dat'%(int(run)-2)
fnamePost='uxi_x311_%04d.dat'%(int(run)+1)

print 'look for uxi: ',uxiPath+fname
minUxi=None; maxUxi=None
print 'xtc begin: ',xtcBeginTime, xtcEndTime

uxiDict={}
configDict=fillUxiDict(uxiPath+fname, xtcBeginTime, xtcEndTime, uxiDict)
#print configDict

#print 'DEBUG times: ',xtcBeginTime,minUxi
if minUxi is None:
    print 'here:',fnamePre
    configDict=fillUxiDict(uxiPath+fnamePre, xtcBeginTime, xtcEndTime, uxiDict)
if xtcBeginTime<minUxi and os.path.isfile(uxiPath+fnamePre):
    print 'cheking the run before....'
    fillUxiDict(uxiPath+fnamePre, xtcBeginTime, xtcEndTime, uxiDict)

if xtcEndTime>maxUxi and os.path.isfile(uxiPath+fnamePost):
    fillUxiDict(uxiPath+fnamePost, xtcBeginTime, xtcEndTime, uxiDict)

if 'lcls_ts_secs' not in uxiDict.keys():
    print 'did not find any icarus events, will quit'
    sys.exit()

print 'fuxiPost: ',len(uxiDict['lcls_ts_secs']),uxiDict.keys()

iDark, darkA, darkB =  getDarks(int(run))

if configDict is None:
    configDict={}
configDict['pedestal_A'] = darkA
configDict['pedestal_B'] = darkB
configDict['pedestal_run'] = [iDark]
uxiDict['frameA'] = np.array([fA-darkA for fA in uxiDict['frameA']])
uxiDict['frameB'] = np.array([fB-darkB for fB in uxiDict['frameB']])

fidmask = int('0x1ffff',16)
uxiSecs=[]; uxiNsecs=[]; uxiFids=[]
for tsec, tnsec, tfid in zip(uxiDict['lcls_ts_secs'], uxiDict['lcls_ts_necs'], uxiDict['lcls_ts_high']):
    uxiSecs.append(np.int64(tsec))
    uxiNsecs.append(np.int64(tnsec))
    uxiFids.append(int(tfid)&int(fidmask))

uxiDict['lcls_ts_secs'] = np.array(uxiSecs)
uxiDict['lcls_ts_necs'] = np.array(uxiNsecs)
uxiDict['lcls_ts_fids'] = np.array(uxiFids)


print 'DDD ',uxiDict['lcls_ts_secs']
evtNr_withUxi=0
evtFid_withUxi=0
for eventNr,evt in enumerate(ds.events()):
    printMsg(eventNr, evt.run(), ds.rank, ds.size)

    if eventNr >= maxNevt/ds.size:
        break

    #check if event is in uxi file.
    evttime = evt.get(psana.EventId).time()
    evtfid = evt.get(psana.EventId).fiducials()
    #find uxi pictures taken the same second & get nsec & fiducials
    uxiEvent=False; uxiIdx=-1
    if (len(np.argwhere(uxiDict['lcls_ts_secs']==evttime[0]))>0):
        secIdx=np.argwhere(uxiDict['lcls_ts_secs']==evttime[0])[0][0]
        #print 'DEBUG1 IDX: ',secIdx, evttime[0], uxiDict['lcls_ts_secs'][secIdx]
        #print 'DEBUG2 IDX: ',secIdx, evttime[1], uxiDict['lcls_ts_necs'][secIdx]
        #print 'DEBUG3 IDX: ',secIdx, evtfid, uxiDict['lcls_ts_fids'][secIdx]
        if len(np.argwhere(uxiDict['lcls_ts_necs']==evttime[1]))>0 and len(np.argwhere(uxiDict['lcls_ts_fids']==evtfid))>0:
            print eventNr,eventNr-evtNr_withUxi,' dFid:',evtfid-evtFid_withUxi,' IDX: ',secIdx, ' sec: ',evttime[0],' nsec: ',evttime[1],uxiDict['lcls_ts_necs'][secIdx], ' fid ', evtfid, uxiDict['lcls_ts_fids'][secIdx]
            evtNr_withUxi=eventNr
            evtFid_withUxi=evtfid
            if evttime[1]==uxiDict['lcls_ts_necs'][secIdx] and evtfid==uxiDict['lcls_ts_fids'][secIdx]:
                uxiEvent=True
                uxiIdx=secIdx
                #print 'exact'
            #else:
            #    print 'diff: ',secIdx==nsecIdx, '---', np.argwhere(uxiDict['lcls_ts_fids']==evtfid),' -- ', np.argwhere(uxiDict['lcls_ts_fids']==evtfid)[-1][0]

    if not uxiEvent:
        continue;

    #add default data
    defData = detData(defaultDets, evt)
    #for key in defData.keys():
    #    print eventNr, key, defData[key]
    smldata.event(defData)

    userData = {}
    for det in dets:
        try:
            #this should be a plain dict. Really.
            det.evt = dropObject()
            det.getData(evt)
            det.processDetector()
            userData[det._name]=getUserData(det)
            #print userData[det._name]
        except:
            pass

        ##make an image
        if det._name == 'jungfrau1M':
            #pass
            jungfrau_image = det.det.image(run, det.evt.dat)
            userData[det._name]['image'] = jungfrau_image 

        elif det._name == 'acq01':
            trace =  det.evt.dat.squeeze()
            retDict = findPeakSimple(trace, nMaxPeak=2, peakWid=20, peakOff=5)
            if retDict['peakIdx'][0]<retDict['peakIdx'][1]:
                userData[det._name]['peakIdx'] = np.array([retDict['peakIdx'][0], retDict['peakIdx'][1]])
                userData[det._name]['peakSum'] = np.array([retDict['peakSum'][0], retDict['peakSum'][1]])
                userData[det._name]['peakMax'] = np.array([retDict['peakMax'][0], retDict['peakMax'][1]])
            else:
                userData[det._name]['peakIdx'] = np.array([retDict['peakIdx'][1], retDict['peakIdx'][0]])
                userData[det._name]['peakSum'] = np.array([retDict['peakSum'][1], retDict['peakSum'][0]])
                userData[det._name]['peakMax'] = np.array([retDict['peakMax'][1], retDict['peakMax'][0]])
            #hardcode the positions.
            userData[det._name]['peakFixed'] = np.array([trace[3150:3170].sum(), trace[3185:3205].sum()])

            fitPars = fitTemplateLeastsq(trace, peakMain)
            fitRes = templateArray(fitPars, trace.shape[0], peakMain)
            userData[det._name]['resid'] = residuals(trace[3125:3300], fitRes[3125:3300])
            for ip,p in enumerate(fitPars):
                userData[det._name]['par%d'%ip] = p

    smldata.event(userData)

    #now add uxi data.
    evtUxiDict={'uxi':{}}
    for key in uxiDict.keys():
        if isinstance(uxiDict[key][uxiIdx],basestring):
            if uxiDict[key][uxiIdx].find('.')>=0:
                evtUxiDict['uxi'][key]=float(uxiDict[key][uxiIdx])
            else:
                evtUxiDict['uxi'][key]=int(uxiDict[key][uxiIdx])
        else:
            evtUxiDict['uxi'][key]=uxiDict[key][uxiIdx]

    ####
    #work on the uxi data stuff here. Ideally w/ function calls.
    ####
    for fN, frame in zip(['A','B'],[evtUxiDict['uxi']['frameA'], evtUxiDict['uxi']['frameB']]):
        #common mode (unbonded):get mean of two rows of unbonded pixels as possible estimate of the common mode
        evtUxiDict['uxi']['cmUnb_median_'+fN] = getUnbonded(frame)
        evtUxiDict['uxi']['cmUnb_mean_'+fN] = getUnbonded(frame, useMed=False)

        #common mode (get position of zero photon peak using a gaussian fit)
        fitResult = getZeroPeakFit(frame)
        for key in fitResult.keys():
            fkey = 'cmFit_'+key+'_'+fN
            evtUxiDict['uxi'][fkey] = fitResult[key]
            #if key.find('error')>=0:
            #    evtUxiDict['uxi'][fkey] = fitResult[key]
            #if not key.find('estimate')>=0:
            #    evtUxiDict['uxi'][fkey] = fitResult[key]
    
    smldata.event(evtUxiDict)

    #print evttime, evtfid, secIdx,' == ',uxiDict['lcls_ts_fids'][secIdx]
    #DEBUG: quit after first event.
    #break

    ##acqiris functions into own utility file (working on waveform), call here
    #if 'acq01' in userData.keys():
    #    acqData = userData['acq01']


print 'rank %d on %s is finished'%(ds.rank, hostname)
####
#save the config info
####
smldata.save(configDict)
#smldata.save()
