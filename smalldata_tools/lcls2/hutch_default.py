import numpy as np
#import psana
from smalldata_tools.SmallDataDefaultDetector import *
from smalldata_tools.epicsarchive import EpicsArchive
from datetime import datetime


def defaultDetectors(hutch, run=None, env=None):
    if hutch.lower()=='det':
        dets = detDetectors()
    elif hutch.lower()=='dia':
        dets = diaDetectors()
    elif hutch.lower()=='tmo':
        dets = tmoDetectors(run)
    elif hutch.lower()=='rix':
        dets = rixDetectors(run)
    elif hutch.lower()=='ued':
        dets = uedDetectors(run)
    else:
        dets = []
    detsInRun= [ det for det in dets if det.inRun() ]
    return detsInRun


def tmoDetectors(run, beamCodes=[[162],[91]]):
    dets=[]
    dets.append(scanDetector('scan', run))
    dets.append(genlcls2Detector('timing',run))
    dets.append(lcls2_lightStatus(beamCodes,run))
    dets.append(genlcls2Detector('gmd',run))
    dets.append(genlcls2Detector('xgmd',run))
    dets.append(genlcls2Detector('ebeam',run))
    dets.append(genlcls2Detector('pcav',run))
    dets.append(fimfexDetector('tmo_fim0',run))
    dets.append(fimfexDetector('tmo_fim1',run))
    dets.append(ttlcls2Detector('tmoopal2',run, saveTraces=True))
    dets.append(lcls2_epicsDetector(PVlist=['MR1K4_pitch', 'MR2K4_pitch'],run=run))
    return dets


def rixDetectors(run, beamCodes=[[-136],[77]]):
    dets=[]
    dets.append(scanDetector('scan', run))
    dets.append(genlcls2Detector('timing',run))
    dets.append(lcls2_lightStatus(beamCodes,run))
    dets.append(genlcls2Detector('gmd',run))
    dets.append(genlcls2Detector('xgmd',run))
    dets.append(genlcls2Detector('ebeam',run))
    dets.append(genlcls2Detector('pcav',run))
    dets.append(fimfexDetector('rix_fim0',run))
    dets.append(fimfexDetector('rix_fim1',run))
    dets.append(fimfexDetector('rix_fim2',run))
    dets.append(fimfexDetector('crix_w8',run))
    dets.append(genlcls2Detector('mono_encoder',run))

    dets.append(ttlcls2Detector('atmopal',run, saveTraces=True))
    #check a RIX scan to figure out controlDetector.
    return dets

def lcls2_detOnceData(det, data, ts, noarch):
    if not noarch:
        # If we've missed something, look for it in the archiver.
        # ts is our timestamp.
        data = addArchiverData(det, data, ts)
    return data