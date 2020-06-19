import re
import sys
import os
import time
import numpy as np
from collections import namedtuple
import tables
import requests
from datetime import datetime
from krtc import KerberosTicket
try:
    #python3
    from urllib.parse import urlparse
except:
    #python 2
    from urlparse import urlparse

abspath=(os.path.abspath(os.path.dirname(__file__))).replace('/examples','')
sys.path.append(abspath)
sys.path.append(abspath+'/smalldata_tools')

from utilities import dictToHdf5
from GaussFit import GaussFit

uxi_shape = (1024, 512)
uxiDarkRuns = [32,36,57,71,72,80,86,109]
uxiHistogramEdges = np.arange(-150,200)

def getUxiDict(run):
    ##
    # get timestamp range for xtc file (???)
    ##
    xtcBeginTime=None
    xtcEndTime=None
    args_url = "https://pswww.slac.stanford.edu/ws-kerb/lgbk/"
    ws_url = args_url + "/lgbk/{0}/ws/runs".format('xcsx31116')
    krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    r = requests.get(ws_url, headers=krbheaders)
    runDict = (r.json()["value"])
    #print(json.dumps(r.json()["value"], sort_keys=True, indent=4))
    xtcEndTime=-1
    for thisRun in runDict:
        if thisRun['num']==int(run):
            #print(thisRun['num'], thisRun['begin_time'], thisRun['end_time'])
            xtcBeginTimeStr = (thisRun['begin_time'].split('+')[0]).replace('T',' ')
            xtcEndTimeStr = (thisRun['end_time'].split('+')[0]).replace('T',' ')
            xtcBeginTime = int((datetime.strptime(xtcBeginTimeStr, '%Y-%m-%d %H:%M:%S')-datetime.utcfromtimestamp(0)).total_seconds())
            xtcEndTime = int((datetime.strptime(xtcEndTimeStr, '%Y-%m-%d %H:%M:%S')-datetime.utcfromtimestamp(0)).total_seconds())

    if xtcEndTime<0:
        print('something went wrong with getting the xtc time from the run database')
        return None

    ##
    # put some UXI stuff here.
    ##
    uxiPath='/reg/d/psdm/xcs/xcsx31116/usrdaq/uxi/'
    #start with actual run first.
    fname='uxi_x311_%04d.dat'%int(run)
    fnamePre='uxi_x311_%04d.dat'%(int(run)-1)
    if int(run) == 19: #yes, this is special....
        fnamePre='uxi_x311_%04d.dat'%(int(run)-2)
    fnamePost='uxi_x311_%04d.dat'%(int(run)+1)

    print('look for uxi: ',uxiPath+fname)
    minUxi=None; maxUxi=None
    print('xtc begin: ',xtcBeginTime, xtcEndTime)

    uxiDict={}
    try:
        configDict, minUxi, maxUxi=fillUxiDict(uxiPath+fname, xtcBeginTime, xtcEndTime, uxiDict)
    except:
        pass

    #print 'DEBUG times: ',xtcBeginTime,minUxi
    if minUxi is None:
        configDict, tmp, tmp2 = fillUxiDict(uxiPath+fnamePre, xtcBeginTime, xtcEndTime, uxiDict)
    if xtcBeginTime<minUxi and os.path.isfile(uxiPath+fnamePre):
        print('checking the run before....')
        tmp1, tmp2, tmp3 = fillUxiDict(uxiPath+fnamePre, xtcBeginTime, xtcEndTime, uxiDict)

    if xtcEndTime>maxUxi and os.path.isfile(uxiPath+fnamePost):
        print('checking the run after....')
        tmp1, tmp2, tmp3 = fillUxiDict(uxiPath+fnamePost, xtcBeginTime, xtcEndTime, uxiDict)
    
    if uxiDict=={}:
        return uxiDict, None
    fidmask = int('0x1ffff',16)
    uxiSecs=[]; uxiNsecs=[]; uxiFids=[]
    for tsec, tnsec, tfid in zip(uxiDict['lcls_ts_secs'], uxiDict['lcls_ts_necs'], uxiDict['lcls_ts_high']):
        uxiSecs.append(np.int64(tsec))
        uxiNsecs.append(np.int64(tnsec))
        uxiFids.append(int(tfid)&int(fidmask))

    uxiDict['lcls_ts_secs'] = np.array(uxiSecs)
    uxiDict['lcls_ts_necs'] = np.array(uxiNsecs)
    uxiDict['lcls_ts_fids'] = np.array(uxiFids)
        
    return uxiDict, configDict


def fillUxiDict(fname, xtcBeginTime, xtcEndTime, uxiDict, debug=False):
    if not os.path.isfile(fname): return
    dataDictTS = {}
    #read the uxi file, get dict of frame fids&evttimes & images & config.
    t_preUxiRead=time.time()
    dataDictTS['uxi'] = get_uxi_timestamps(fname)    
    
    print('UXI file ',fname, dataDictTS['uxi']['lcls_ts_secs'])
    minUxi = int(min(dataDictTS['uxi']['lcls_ts_secs']))
    maxUxi = int(max(dataDictTS['uxi']['lcls_ts_secs']))
    if debug: print('min/max time: ',minUxi, maxUxi)
    
    configDict=None
    #if uxi event time stamps lie in uxi range, get all uxi data
    if minUxi<=xtcEndTime and maxUxi>=xtcBeginTime:
        inXtc = np.array([(int(ts) >= xtcBeginTime and int(ts) <= xtcEndTime) for ts in dataDictTS['uxi']['lcls_ts_secs']])
        uxiDictAll,configDict = read_uxi(fname, returnConfig=True)    
        for key in uxiDictAll.keys():
            uxiDictAllValue = uxiDictAll[key]
            if key.find('frame')<0 and key.find('lcls_ts')<0:
                if uxiDictAllValue[0].find('.')<0:
                    uxiDictAllValue = np.array(uxiDictAllValue).astype(int)
                else:
                    uxiDictAllValue = np.array(uxiDictAllValue).astype(float)
            if key not in uxiDict.keys():
                uxiDict[key] = uxiDictAllValue#np.array((uxiDictAll[key]))[inXtc]
            else:
                uxiDict[key] = np.append(uxiDict[key], np.array(uxiDictAllValue)[inXtc], axis=0)
                #uxiDict[key] = np.append(uxiDict[key], np.array((uxiDictAll[key]))[inXtc], axis=0)
        if debug: print('check lengths: ',len(uxiDictAll['lcls_ts_secs']),' inXtc: ',inXtc.sum(),' out dict ',len(uxiDict['lcls_ts_secs']))
    return configDict, minUxi, maxUxi

def read_uxi(fname, returnConfig=False):
    header = re.compile("#\s+Frame\s+header\s+info:\s+(?P<meta>.*)")
    config_pat = re.compile("#\s+Config:\s+(?P<conf>.*)")
    pat = re.compile("\s+")
    config = None
    Metadata = None
    MetadataList = []
    configDict={}
    with open(fname, 'r') as infile:
        dataDict={}
        #nline=0 #DEBUG
        while True:
            line = infile.readline()
            if not line:
                break
            #nline+=1 #DEBUG
            #if nline>10: #DEBUG
            #             break
            if line.startswith("#"):
                config_match = config_pat.match(line)
                header_match = header.match(line)
                if header_match:
                    Metadata = namedtuple('Metadata', header_match.group("meta"))
                elif config_match:
                    configDict = eval(config_match.group("conf"))
                    print(configDict)
                continue
            frameA = np.reshape(np.load(infile), uxi_shape)
            frameB = np.reshape(np.load(infile), uxi_shape)
            metadata = Metadata._make(pat.split(line.rstrip()))

            if 'frameA' not in dataDict.keys():
                dataDict['frameA'] = [frameA]
            else:
                dataDict['frameA'].append(frameA)
            if 'frameB' not in dataDict.keys():
                dataDict['frameB'] = [frameB]
            else:
                dataDict['frameB'].append(frameB)

            MetadataDict = metadata._asdict()
            for key in MetadataDict.keys():
                if key not in dataDict.keys():
                    dataDict[key] = [MetadataDict[key]]
                else:
                    dataDict[key].append(MetadataDict[key])

            #print(metadata)
            #print(frameA)
            #print(frameB)

        if returnConfig:
            return dataDict, configDict
        else:
            return dataDict



def uxi_pedestal(fname):
    header = re.compile("#\s+Frame\s+header\s+info:\s+(?P<meta>.*)")
    config_pat = re.compile("#\s+Config:\s+(?P<conf>.*)")
    pat = re.compile("\s+")
    config = None
    Metadata = None
    MetadataList = []
    with open(fname, 'r') as infile:
        dataDict={}
        #nline=0 #DEBUG
        while True:
            line = infile.readline()
            if not line:
                break
            if line.startswith("#"):
                config_match = config_pat.match(line)
                header_match = header.match(line)
                if header_match:
                    Metadata = namedtuple('Metadata', header_match.group("meta"))
                elif config_match:
                    configDict = eval(config_match.group("conf"))
                    print(configDict)
                continue
            frameA = np.load(infile)
            frameB = np.load(infile)
            metadata = Metadata._make(pat.split(line.rstrip()))

            frameA = np.reshape(np.array(frameA), uxi_shape )
            frameB = np.reshape(np.array(frameB), uxi_shape )

            if 'frameA' not in dataDict.keys():
                dataDict['frameA'] = [frameA]
            else:
                dataDict['frameA'].append(frameA)
            if 'frameB' not in dataDict.keys():
                dataDict['frameB'] = [frameB]
            else:
                dataDict['frameB'].append(frameB)

            #print(metadata)
            #print(frameA)
            #print(frameB)

    runStr = fname.split('_')[-1].replace('.dat','')
    print('fname ',fname, runStr)
    run = int(runStr)
    fnameOut='/reg/d/psdm/xcs/xcsx31116/calib/uxi/Run%04d_pedestals.h5'%run

    outDict={}
    #outDict['meanA'] = np.reshape(np.array(dataDict['frameA'])[1:].mean(axis=0), uxi_shape)
    #outDict['medA']  = np.reshape(np.median(np.array(dataDict['frameA'])[1:],axis=0), uxi_shape)
    #outDict['meanB'] = np.reshape(np.array(dataDict['frameB'])[1:].mean(axis=0), uxi_shape)
    #outDict['medB']  = np.reshape(np.median(np.array(dataDict['frameB'])[1:],axis=0), uxi_shape)
    outDict['meanA'] = np.array(dataDict['frameA'])[1:].mean(axis=0)
    outDict['medA']  = np.median(np.array(dataDict['frameA'])[1:],axis=0)
    outDict['meanB'] = np.array(dataDict['frameB'])[1:].mean(axis=0)
    outDict['medB']  = np.median(np.array(dataDict['frameB'])[1:],axis=0)
    outDict['rawA']  = np.array(dataDict['frameA'])
    outDict['rawB']  = np.array(dataDict['frameB'])

    dictToHdf5(fnameOut,outDict)
    

def get_uxi_timestamps(fname):
    header = re.compile("#\s+Frame\s+header\s+info:\s+(?P<meta>.*)")
    config_pat = re.compile("#\s+Config:\s+(?P<conf>.*)")
    pat = re.compile("\s+")
    config = None
    Metadata = None
    MetadataList = []
    with open(fname, 'r') as infile:
        dataDict={}
        #nline=0 #DEBUG
        while True:
            line = infile.readline()
            if not line:
                break
            if line.startswith("#"):
                config_match = config_pat.match(line)
                header_match = header.match(line)
                if header_match:
                    Metadata = namedtuple('Metadata', header_match.group("meta"))
                elif config_match:
                    continue
                    #config = eval(config_match.group("conf"))
                    #print("Config: %s"%config)
                continue
            frameA = np.load(infile)
            frameB = np.load(infile)
            metadata = Metadata._make(pat.split(line.rstrip()))

            MetadataDict = metadata._asdict()
            for key in MetadataDict.keys():
                if key not in dataDict.keys():
                    dataDict[key] = [MetadataDict[key]]
                else:
                    dataDict[key].append(MetadataDict[key])

        return dataDict

def getDarks(run, useMed=True):
    iDark=-1
    for d in uxiDarkRuns:
        if d > run: continue
        iDark=d
    if iDark<0: iDark=uxiDarkRuns[0]

    dat = tables.open_file('/reg/d/psdm/xcs/xcsx31116/calib/uxi/Run%04d_pedestals.h5'%iDark)
    if useMed: return iDark, dat.root.medA,dat.root.medB
    else: return iDark,dat.root.meanA,dat.root.meanB

def getUnbonded(frame, useMed=True):
    if frame.shape != uxi_shape:
        print('frame does not have right shape, expected ',uxi_shape,', and got ',frame.shape)
        return
    unbData = frame[:,255:257]
    if useMed:
        return np.nanmedian(unbData)
    else:
        return np.nanmean(unbData)

#FIX ME
#this is actually a fit of the low intensity part of the spectrum, NOT the zero photon (or first) peak.
def getZeroPeakFit(frame, histogramEdges=None):
    if histogramEdges is None:
        histogramEdges = uxiHistogramEdges
    his = np.histogram(frame.flatten(), histogramEdges)
    #print his[0]
    fitResult = GaussFit(his[0])
    #print fitResult
    #print 'DEBUG: ',fitResult['mean_estimate'],(fitResult['mean_estimate']!=np.nan),np.isnan(fitResult['mean_estimate'])
    meanVal = fitResult['mean']
    if np.isnan(fitResult['mean_estimate']):
        return fitResult
    if fitResult['mean_estimate'] < his[1].shape:
        fitResult['mean_estimate']=float(his[1][fitResult['mean_estimate']])
    else:
        fitResult['mean_estimate']=np.nan
    if int(meanVal)+(his[1][1]-his[1][0])*(meanVal-int(meanVal)) < his[1].shape:
        fitResult['mean']=float(his[1][int(meanVal)+(his[1][1]-his[1][0])*(meanVal-int(meanVal))])
    else:
        fitResult['mean']=np.nan
    #print 'mean in x: ',fitResult
    return fitResult
