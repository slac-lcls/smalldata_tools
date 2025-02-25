#!/usr/bin/env python
import psana
import numpy as np
import holoviews as hv

hv.extension("bokeh")
import panel as pn

pn.extension()
import os
import argparse
import sys
import logging
import requests
from pathlib import Path
from requests.auth import HTTPBasicAuth
import socket
from typing import Optional, Union
import mimetypes

try:
    basestring
except NameError:
    basestring = str
fpath = os.path.dirname(os.path.abspath(__file__))
fpathup = "/".join(fpath.split("/")[:-1])
try:
    fpath = os.environ.get("MYDIR", fpathup).replace("/arp_scripts", "")
except:
    fpath = fpathup
sys.path.append(fpath)
from smalldata_tools.utilities import rebin
from smalldata_tools.utilities import postRunTable
from smalldata_tools.utilities import getElogBasicAuth
from smalldata_tools.utilities import postElogMsg
from smalldata_tools.utilities import getRunsWithTag
from summaries.PedestalPlot import gainSwitching
from summaries.PedestalPlot import statusDict
from summaries.PedestalPlot import ped_rms_histograms
from summaries.PedestalPlot import plotPedImgs

# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument(
    "--run", help="run", type=str, default=os.environ.get("RUN_NUM", "")
)
parser.add_argument(
    "--experiment",
    help="experiment name",
    type=str,
    default=os.environ.get("EXPERIMENT", ""),
)
parser.add_argument(
    "--nopostElog",
    help="do not post plot &status to elog",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--pedImgs", help="make images of pedestals", action="store_true", default=False
)
parser.add_argument(
    "--pedDiffImgs",
    help="make images of first 10 subtracted images ",
    action="store_true",
    default=False,
)
parser.add_argument("--url", default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
args = parser.parse_args()
logger.debug("Args to be used for pedestal plots: {0}".format(args))

nosave_elog = args.nopostElog
make_ped_imgs = args.pedImgs
make_ped_data_imgs = args.pedDiffImgs
expname = args.experiment
run = int(args.run)

SIT_PSDM_DATA = Path(os.environ.get("SIT_PSDM_DATA"))

####
# helper functions
####

def statusStats(det_name, printme=False, request_run=None):
    print(det_name)
    print(run)
    status_type = "camera"  # see if different detector types here are different?
    if isinstance(run, int):
        det = psana.Detector(det_name)
        if request_run:
            statusmask = det.mask(request_run, status=True)
            status = det.status(request_run)
        else:
            statusmask = det.mask(run, status=True)
            status = det.status(run)
        if det.is_epix10ka_any():
            status_type = "epix10k"
        elif det.is_jungfrau():
            status_type = "jungfrau"
        elif det.is_epix100a():
            status_type = "epix100"
        elif det.is_cspad() or det.is_cspad2x2():
            status_type = "cspad"
    else:
        det = run.Detector(det_name)
        detrawid = det.raw._uniqueid
        runnum = run.runnum
        from psana.pscalib.calib.MDBWebUtils import calib_constants
        if request_run:
            status = calib_constants(detrawid, exp=expname, ctype="status", run=request_run)[0]
            statusmask = calib_constants(detrawid, exp=expname, ctype="mask", run=request_run)[0]
        else:
            status = calib_constants(detrawid, exp=expname, ctype="status", run=runnum)[0]
            statusmask = calib_constants(detrawid, exp=expname, ctype="mask", run=runnum)[0]
        if det.calibconst['pedestals'][1]['dettype']=='epix10ka':
            status_type = "epix10k"
        elif det.calibconst['pedestals'][1]['dettype']=='jungfrau':
            status_type = "jungfrau"
        elif det.calibconst['pedestals'][1]['dettype']=='epix100':
            status_type = "epix100"
        #this is not how to get status or mask...
        print(statusmask)
        print(status)    
            
    status_unp = np.array(
        [
            (
                np.unpackbits(tstatus.flatten().astype(np.uint8), bitorder="little")
            ).reshape([int(tstatus.flatten().shape[0]), 8])
            for tstatus in status
        ]
    )

    statusStatDict = {}
    for istatus, statusk in enumerate(statusDict[status_type]):
        if status_type in gainSwitching:
            statusStatDict[statusk] = int(
                (status_unp.sum(axis=0)[:, istatus] > 0).sum()
            )
        else:
            statusStatDict[statusk] = int((status_unp[:, istatus] > 0).sum())
    statusStatDict["total_masked"] = int(
        statusmask.flatten().shape[0] - statusmask.sum()
    )
    if status_type in gainSwitching:
        print(status.shape)
        for icycle in range(status.shape[0]):
            statusStatDict["cycle%i" % icycle] = int((status[icycle] > 0).sum())

    if printme:
        for k, v in statusStatDict.items():
            print(k, v)
    return statusStatDict


def getKerberosAuthHeaders() -> dict: ...



def postBadPixMsg(
    detectors: list,
    exp: str,
    run: int,
    *,
    tag: str = "SUMMARY_BAD_PIX",
    title: str = "Detector Bad Pixel Info",
) -> None:
    """Post bad pixel data for a given detector and run to the eLog.

    Parameters
    ----------
    detectors (str) Names of detectors to pull bad pixel data for.
    exp (str) Experiment name.
    run (int) Run number. Pulls data for this run and all previous DARK runs.
    tag (str) Optional. Tag for the bad pixel summary posts.
    title (str) Optional. Title for bad pixel summary posts.
    """
    http_auth: HTTPBasicAuth = getElogBasicAuth(exp=exp)

    dark_runs: list = getRunsWithTag(exp=exp, tag="DARK", http_auth=http_auth)

    if dark_runs:
        dark_runs = [dr for dr in dark_runs if dr <= run]

        table_header: str = (
            '<thead><tr><th colspan="3">'
            f"<center>{title}</center>"
            "</th></tr></thead>"
        )
        table_body: str = (
            "<tbody><tr>"
            "<td><b><center>Detector</center></b></td>"
            "<td><b><center>Number of bad pixels</center></b></td>"
            "<td><b><center>Difference vs previous DARK</center></b></td></tr>"
        )

        for det_name in detectors:
            bad_pix: list = []
            for dr in dark_runs:
                try:
                    stat_dict: dict = statusStats(det_name, request_run=dr)
                    bad_pix.append(stat_dict["total_masked"])
                except TypeError as err:
                    # `statusStats` may throw TypeError if detector not in run
                    logger.debug(
                        f"Is {det_name} not present in DARK run {dr}?\n" f"ERROR: {err}"
                    )
                except KeyError as err:
                    # `statusStats` may throw KeyError if detector not in run
                    logger.debug(
                        f"Is {det_name} not present in DARK run {dr}?\n" f"ERROR: {err}"
                    )

            # Report current DARK run bad pix and the delta vs previous DARK run
            curr_bad_pix = bad_pix[-1]
            if len(bad_pix) > 1:
                diff_bad_pix = bad_pix[-1] - bad_pix[-2]
            else:
                diff_bad_pix = "-"
            det_entry: str = (
                f"<tr><td><center>{det_name}</center></td>"
                f"<td><center>{curr_bad_pix}</center></td>"
                f"<td><center>{diff_bad_pix}</center></td></tr>"
            )
            table_body += det_entry
        table_body += "</tbody>"
        msg: str = f'<table border="1">{table_header}{table_body}</table>'
    else:
        msg: str = "No DARK runs or cannot communicate with eLog."
        logger.debug(msg)

    #DEBUG - turn this back on before deploying!!!
    #postElogMsg(exp=exp, msg=msg, run=run, tag=tag, title=title)


def plotDataImgs(expname, run, det_name, nCycles, plotInfo=None):
    import smalldata_tools.SmallDataAna_psana as sda

    anaps = sda.SmallDataAna_psana(expname, run, plotWith=None)

    gspecI = pn.GridSpec(
        sizing_mode="stretch_both", max_width=900, name="Pedestal data - subtracted"
    )
    iwidth = 3
    iheight = 3
    gspecI[0, 0 : (iwidth * 3)] = pn.Row(
        "# Pedestal data - subtracted - Run %04d" % run
    )

    print(expname, run, det_name)
    for i in range(min(5, nCycles)):
        common_mode = None
        if det_name.find("Epix") >= 0:
            common_mode = 80
        anaps.AvImage(det_name, common_mode=common_mode, nSkip=1200 * i, numEvts=10)

        retp = anaps.plotAvImage(returnIt=True, plotWith=None)
        imageDim = hv.Dimension(
            ("ped_subtr", "in keV"),
            range=(np.nanpercentile(retp, 0.1), np.nanpercentile(retp, 99.9)),
        )
        if plotInfo is not None:
            xrange, yrange, xDim, yDim = plotInfo
        if xrange is not None:
            xrange = [0, retp.shape[0]]
        if yrange is not None:
            yrange = [0, retp.shape[1]]
        timg = hv.Image(
            retp,
            bounds=(xrange[0], yrange[0], xrange[1], yrange[1]),
            kdims=[xDim, yDim],
            vdims=[imageDim],
            label="Image, Cycle %d" % i,
        ).options(colorbar=True, aspect="equal", cmap="rainbow")
        row = 1 + iheight * int(i * 0.5)
        col1 = iwidth * (i % 2)
        col2 = iwidth * (i % 2 + 1)
        gspecI[row, col1:col2] = pn.Column(timg)

    return gspecI


def allPlots(
    det_name,
    run,
    make_ped_imgs=False,
    make_ped_data_imgs=False,
    tabs=None,
    detImgMaxSize=400,
):
    print("Working on plots for ", det_name)
    det = psana.Detector(det_name)
    peds = det.pedestals(run)
    try:
        peds_pre = det.pedestals(run - 1)
        diffPeds = peds - peds_pre
    except:
        diffPeds = None
        diffPedsImgs = None
    noise = det.rms(run)
    runnum = run
    
    xDim = hv.Dimension(("x", "x in micron"))
    yDim = hv.Dimension(("y", "y in micron"))
    try:
        xcoords, ycoords = det.coords_xy(run)
        xrange = (np.nanmin(xcoords), np.nanmax(xcoords))
        yrange = (np.nanmin(ycoords), np.nanmax(ycoords))
    except:
        if len(noise.shape) == 2:
            xmax = noise.shape[0]
            ymax = noise.shape[1]
        else:
            xmax = noise[0].shape[0]
            ymax = noise[0].shape[1]
        xrange = (0, xmax)
        yrange = (0, ymax)
        xDim = hv.Dimension(("x", "x in pixel"))
        yDim = hv.Dimension(("y", "y in pixel"))

    plotInfo = [xrange, yrange, xDim, yDim]

    nCycles = 1
    if len(peds.shape) >= 3 and det_name.find("CsPad") < 0:
        nCycles = peds.shape[0]
    if nCycles > 5:
        nCycles = 5

    #make images for pedestals&noise
    if nCycles == 1:
        pedImgs = [ det.image(runnum, peds) ]
        noiseImgs = [ det.image(runnum, noise) ]
        if diffPeds is not None:
            diffPedsImgs = det.image(runnum, diffPeds)
    else:
        pedImgs = [ det.image(runnum, thisPed) for thisPed in peds ]
        noiseImgs = [ det.image(runnum, thisNoise) for thisNoise in noise ]
        if diffPeds is not None:
            diffPedsImgs = [ det.image(runnum, thisDiff) for thisDiff in diffPeds ]

    pedHists, noiseHists, diffHists = ped_rms_histograms(
        nCycles, peds, noise, diffPeds, det_name
    )
    gspecH = pn.GridSpec(
        sizing_mode="stretch_width", max_width=500, name="Histogram - %s" % det_name
    )
    gspecH[0, 0:8] = pn.Row("# Pedestals&RMS Histograms - Run %04d" % (runnum))
    gspecH[1:4, 0:8] = pn.Column(hv.Overlay(pedHists))
    gspecH[4:7, 0:8] = pn.Column(hv.Overlay(noiseHists))
    if diffHists is not None:
        gspecH[7:10, 0:8] = pn.Column(hv.Overlay(diffHists))
    if tabs is None:
        tabs = pn.Tabs(gspecH)
        # this is for debugging.....
        # return tabs
    else:
        tabs.append(gspecH)

    if make_ped_imgs:
        pedImgs, rmsImgs, diffImgs = plotPedImgs(
            nCycles,
            runnum,
            peds,
            noise,
            pedImgs,
            noiseImgs,
            diffPeds,
            diffPedsImgs,
            detImgMaxSize=detImgMaxSize,
            plotInfo=plotInfo,
        )

        gspec = pn.GridSpec(
            sizing_mode="stretch_both", max_width=1000, name="Det Imgs - %s" % det_name
        )
        iwidth = 3
        iheight = 3
        gspec[0, 0 : (iwidth * 3)] = pn.Row("# Pedestals&RMS - Run %04d" % (runnum))
        if nCycles == 1:
            gspec[1 : (1 * iheight) + 1, 0:iwidth] = pn.Column(pedImgs[0])
            gspec[(1 * iheight) + 1 : (2 * iheight) + 1, 0:iwidth] = pn.Column(
                rmsImgs[0]
            )
            if len(diffImgs) == len(pedImgs):
                gspec[(2 * iheight) + 1 : (3 * iheight) + 1, 0:iwidth] = pn.Column(
                    diffImgs[0]
                )
        else:
            for i in range(nCycles):
                gspec[
                    (i * iheight + 1) : ((i + 1) * iheight + 1),
                    (iwidth * 0) : (iwidth * 1),
                ] = pn.Column(pedImgs[i])
                gspec[
                    (i * iheight + 1) : ((i + 1) * iheight + 1),
                    (iwidth * 1) : (iwidth * 2),
                ] = pn.Column(rmsImgs[i])
                if len(diffImgs) == len(pedImgs):
                    gspec[
                        (i * iheight + 1) : ((i + 1) * iheight + 1),
                        (iwidth * 2) : (iwidth * 3),
                    ] = pn.Column(diffImgs[i])
        tabs.append(gspec)

    if make_ped_data_imgs:
        gspecI = plotDataImgs(
            det.env.experiment(), runnum, det.name.__str__(), nCycles, plotInfo=plotInfo
        )
        tabs.append(gspecI)

    return tabs


def plotPedestals(
    expname="mfxc00118",
    run=364,
    nosave_elog=False,
    make_ped_imgs=False,
    make_ped_data_imgs=False,
    detImgMaxSize=400,
):
    ds_name = "exp={}:run={}:smd".format(expname, run)
    ds = psana.DataSource(ds_name)

    det_names = [
        dn[0]
        for dn in psana.DetNames()
        if dn[0].find("Jungfrau") >= 0
        or dn[0].find("Epix") >= 0
        or dn[0].find("Cspad") >= 0
        or dn[0].find("Uxi") >= 0
    ]
    aliases = [
        dn[1]
        for dn in psana.DetNames()
        if dn[0].find("Jungfrau") >= 0
        or dn[0].find("Epix") >= 0
        or dn[0].find("Cspad") >= 0
        or dn[0].find("Uxi") >= 0
    ]
    runnum = run

    tabs = None
    runTableData = {}
    for det_name, alias in zip(det_names, aliases):
        # print(det_name, alias)
        this_det_name = alias
        if alias == "":
            this_det_name = det_name
        tabs = allPlots(
            this_det_name,
            run,
            make_ped_imgs=make_ped_imgs,
            make_ped_data_imgs=make_ped_data_imgs,
            tabs=tabs,
            detImgMaxSize=detImgMaxSize,
        )
        print(this_det_name)
        statusDict = statusStats(this_det_name, run)
        for k, v in statusDict.items():
            runTableData[f"Pixel Status/{this_det_name}_n_{k}"] = v
        print("runTableData:")
        print(runTableData)

    postBadPixMsg(detectors=sorted(det_names, reverse=True), exp=expname, run=run)
    if not nosave_elog:
        elogDir = (
            Path(SIT_PSDM_DATA)
            / expname[:3]
            / expname
            / f"stats/summary/Pedestals/Pedestals_Run{runnum:03d}"
        )

        import os

        if not os.path.isdir(elogDir):
            os.makedirs(elogDir)
        print("Made Directory to save data:", elogDir)
        tabs.save(("%s/report.html" % elogDir))

        postRunTable(runTableData, args.experiment, args.run)

    return tabs


tabs = plotPedestals(
    expname, run, make_ped_imgs=True, make_ped_data_imgs=False, nosave_elog=nosave_elog
)
