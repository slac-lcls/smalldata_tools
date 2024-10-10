import numpy as np
from enum import Enum
from smalldata_tools.lcls2.default_detectors import *
from datetime import datetime


class BeamDestination(Enum):
    """Current reference:
    https://confluence.slac.stanford.edu/display/~carolina/Destinations+SC+agreement+MPS+and+Timing
    """

    LASER_Standby = 0
    SC_DIAG0 = 1
    SC_BSYD = 2
    SC_HXR = 3
    SC_SXR = 4
    SC_LESA = 5
    MPS_RESERVED_Laser_Heater_Shutter = 12
    MPS_RESERVED_Mechanical_Shutter = 13
    MPS_RESERVED = 14
    BSYD_STANDBY_TRIGGERS = 15


def defaultDetectors(hutch, run=None, env=None):
    if hutch.lower() == "det":
        dets = detDetectors()
    elif hutch.lower() == "dia":
        dets = diaDetectors()
    elif hutch.lower() == "tmo":
        dets = tmoDetectors(run)
    elif hutch.lower() == "rix":
        dets = rixDetectors(run)
    elif hutch.lower() == "ued":
        dets = uedDetectors(run)
    else:
        dets = []
    detsInRun = [det for det in dets if det.in_run()]
    return detsInRun


def tmoDetectors(run, beamCodes=[[162], [91]]):
    dets = []
    dets.append(scanDetector("scan", run))
    dets.append(genericDetector("timing", run))
    dets.append(lightStatus(beamCodes, run))
    dets.append(genericDetector("gmd", run))
    dets.append(genericDetector("xgmd", run))
    dets.append(genericDetector("ebeam", run))
    dets.append(genericDetector("pcav", run))
    dets.append(fimfexDetector("tmo_fim0", run))
    dets.append(fimfexDetector("tmo_fim1", run))
    dets.append(ttDetector("tmoopal2", run, saveTraces=True))
    dets.append(epicsDetector(PVlist=["MR1K4_pitch", "MR2K4_pitch"], run=run))
    return dets


# def rixDetectors(run, beam_destination=BeamDestination.SC_SXR, laser_codes=[-272]):
def rixDetectors(run, beam_destination=BeamDestination.SC_SXR, laser_codes=[-284]):
    dets = []
    dets.append(scanDetector('scan', run))
    dets.append(genericDetector("timing", run))
    dets.append(lightStatus(beam_destination, laser_codes, run))
    dets.append(genericDetector("gmd", run))
    dets.append(genericDetector("xgmd", run))
    dets.append(genericDetector("ebeam", run))
    dets.append(genericDetector("pcav", run))
    dets.append(fimfexDetector("rix_fim0", run))
    dets.append(fimfexDetector("rix_fim1", run))
    dets.append(fimfexDetector("rix_fim2", run))
    dets.append(fimfexDetector("crix_w8", run))
    dets.append(interpolatedEncoder("mono_encoder", run))
    dets.append(ttDetector("atmopal", run, saveTraces=True))
    # dets.append(epicsDetector(PVlist=[], run=run))
    return dets
