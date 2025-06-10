import sys
import psana

from pathlib import Path

smdata_tools_path = Path(__file__).parent.parent
sys.path.insert(0, str(smdata_tools_path))

from smalldata_tools.lcls2 import default_detectors, hutch_default
from smalldata_tools.lcls2.DetObject import DetObject

from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, projectionFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.detector_corrections import PolynomialCurveCorrection
from smalldata_tools.ana_funcs.sparsifyFunc import unsparsifyFunc
from smalldata_tools.ana_funcs.waveformFunc import WfIntegration
from smalldata_tools.ana_funcs.waveformFunc import hsdsplitFunc, hsdROIFunc

from smalldata_tools.lcls2.cube.event_screener import *


def detectors(run: psana.psexp.run.Run):
    dets = []

    full_roi = {"thresADU": None, "writeArea": True, "calcPars": False, "ROI": None}

    # andor_dir
    andor_dir = DetObject("andor_dir", run)
    func = ROIFunc(**full_roi)
    andor_dir.addFunc(func)

    # andor_vls
    andor_vls = DetObject("andor_vls", run)
    func = ROIFunc(**full_roi)
    andor_vls.addFunc(func)

    # axis_svls: Dropplet + curve correction + projection
    axis_svls = DetObject("axis_svls", run)
    droplet = dropletFunc(
        threshold=3,
        thresADU=0,  # discard droplet whose total ADU is below this value
        useRms=False,
    )
    unsparsify = unsparsifyFunc()
    poly_corr = PolynomialCurveCorrection(
        polynomial_coefficients=[1.29626379e-6, -1.33220754e-2, 4.70900708e1],
        axis=1,
        method="roll",
    )
    projection = projectionFunc(axis=1)
    poly_corr.addFunc(projection)
    unsparsify.addFunc(poly_corr)
    droplet.addFunc(unsparsify)
    axis_svls.addFunc(droplet)

    # FIMs
    fim0 = DetObject("rix_fim0", run)
    func = ROIFunc(**full_roi)
    fim0.addFunc(func)

    fim1 = DetObject("rix_fim1", run)
    func = ROIFunc(**full_roi)
    fim1.addFunc(func)

    crix_w8 = DetObject("crix_w8",run)
    func = ROIFunc(**full_roi)
    crix_w8.addFunc(func)

    # HSD
    hsd_rois = {
        "hsd_1": [0, -1],
        "hsd_2": [0, -1],
    }
    hsd = DetObject("hsd", run)
    hsdsplit = hsdsplitFunc(writeHsd=False)
    for hsd_name, hsd_roi in hsd_rois.items():
        hsd_roi = {"name": hsd_name + "__ROI", "writeArea": True, "ROI": hsd_roi}
        func = hsdROIFunc(**hsd_roi)
        hsdsplit.addFunc(func)
    hsd.addFunc(hsdsplit)

    dets = [fim0, fim1, crix_w8, hsd, andor_dir, andor_vls, axis_svls]
    return dets


def screener(run):
    ddets = hutch_default.rixDetectors(run)
    lightstatus_det = [d for d in ddets if d.name == "lightStatus"][0]

    # Individual filters
    laser_on = BoolFilter(
        lightstatus_det, "laser", expected_state=True, label="laser_on"
    )

    laser_off = BoolFilter(
        lightstatus_det, "laser", expected_state=False, label="laser_off"
    )

    xray_on = BoolFilter(
        lightstatus_det,
        "xray",
        expected_state=True,
    )

    xray_off = BoolFilter(
        lightstatus_det, "xray", expected_state=False, label="dropped_shots"
    )

    # Combine individual filters into composite filters for different states
    screener_on = CompositeFilter([laser_on, xray_on])
    screener_off = CompositeFilter([laser_off, xray_on])

    # Combine all screeners into an OR screener to be run by the cube
    event_screener = CompositeFilter([laser_on, laser_off], require_all=False)
    return event_screener
