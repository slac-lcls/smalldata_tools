import sys
import psana

from pathlib import Path

smdata_tools_path = Path(__file__).parent.parent
sys.path.insert(0, str(smdata_tools_path))

from smalldata_tools.lcls2 import default_detectors, hutch_default
from smalldata_tools.lcls2.DetObject import DetObject

from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, imageFunc

from smalldata_tools.lcls2.cube.event_screener import *


def detectors(run: psana.psexp.run.Run):
    dets = []

    full_roi = {"thresADU": None, "writeArea": True, "calcPars": False, "ROI": None}

    # epix_quad
    epix_quad = DetObject("epixquad", run)
    roi = {"thresADU": 70, "writeArea": True, "calcPars": False, "ROI": None}
    roi_func = ROIFunc(**roi)

    img_func = imageFunc()
    roi_func.addFunc(img_func)

    epix_quad.addFunc(roi_func)

    dets = [epix_quad]
    return dets


def screener(run):
    return None
