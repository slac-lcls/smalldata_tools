import numpy as np
#import psana
from smalldata_tools.lcls2.default_detectors import *
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
    detsInRun= [ det for det in dets if det.in_run() ]
    return detsInRun


def tmoDetectors(run, beamCodes=[[162], [91]]):
    dets=[]
    dets.append(scanDetector('scan', run))
    dets.append(genericDetector('timing',run))
    dets.append(lightStatus(beamCodes,run))
    dets.append(genericDetector('gmd',run))
    dets.append(genericDetector('xgmd',run))
    dets.append(genericDetector('ebeam',run))
    dets.append(genericDetector('pcav',run))
    dets.append(fimfexDetector('tmo_fim0',run))
    dets.append(fimfexDetector('tmo_fim1',run))
    dets.append(ttDetector('tmoopal2',run, saveTraces=True))
    dets.append(epicsDetector(PVlist=['MR1K4_pitch', 'MR2K4_pitch'], run=run))
    return dets


def rixDetectors(run, beamCodes=[[-136], [77]]):
    dets=[]
    # dets.append(scanDetector('scan', run))
    dets.append(genericDetector('timing', run))
    dets.append(lightStatus(beamCodes, run))
    dets.append(genericDetector('gmd', run))
    dets.append(genericDetector('xgmd', run))
    dets.append(genericDetector('ebeam', run))
    dets.append(genericDetector('pcav', run))
    dets.append(fimfexDetector('rix_fim0', run))
    dets.append(fimfexDetector('rix_fim1', run))
    dets.append(fimfexDetector('rix_fim2', run))
    dets.append(fimfexDetector('crix_w8', run))
    dets.append(genericDetector('mono_encoder', run))
    dets.append(ttDetector('atmopal', run, saveTraces=True))
    # dets.append(epicsDetector(PVlist=[], run=run))
    #check a RIX scan to figure out controlDetector.
    return dets