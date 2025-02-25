#!/usr/bin/env python
########################################################################
## load tools and packages required for data handling and plotting
########################################################################
import panel as pn
import h5py
import os
import argparse
import logging
import requests
import numpy as np
from requests.auth import HTTPBasicAuth
import holoviews as hv
from holoviews import dim

hv.extension("bokeh")
pn.extension()
import sys

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
from smalldata_tools.lcls1.SmallDataAna_psana import SmallDataAna_psana as sdaps
from smalldata_tools.utilities import image_from_dxy
from smalldata_tools.utilities import rebin
from smalldata_tools.utilities import evtt2Rt
from smalldata_tools.utilities import postRunTable

def makeRunTableData(ana, ipmUpDim, ipmDownDim, Filter, scanName):
    n162 = ana.getVar("evr/code_162").sum()
    ana.addCut("evr/code_162", -0.5, 0.5, "xon")
    ana.addCut("evr/code_137", 0.5, 1.5, "xon")
    nOff = ana.getFilter("xon").shape[0] - ana.getFilter("xon").sum()
    # data to be posted to the run table if so requested.
    runtable_data = {"N dropped Shots": int(nOff), "N BYKIK 162": int(n162)}
    if scanName != "":
        runtable_data["scanName"] = scanName

    ipmUpVar = ana.getVar(ipmUpDim.name, useFilter=Filter)
    ipmDownVar = ana.getVar(ipmDownDim.name, useFilter=Filter)
    ipmUpP = np.nanpercentile(ipmUpVar, [25, 50, 75])
    ipmDownP = np.nanpercentile(ipmDownVar, [25, 50, 75])
    runtable_data["%s_1qt" % (ipmUpDim.name.replace("/", "__"))] = ipmUpP[0]
    runtable_data["%s_med" % (ipmUpDim.name.replace("/", "__"))] = ipmUpP[1]
    runtable_data["%s_3qt" % (ipmUpDim.name.replace("/", "__"))] = ipmUpP[2]
    runtable_data["%s_1qt" % (ipmDownDim.name.replace("/", "__"))] = ipmDownP[0]
    runtable_data["%s_med" % (ipmDownDim.name.replace("/", "__"))] = ipmDownP[1]
    runtable_data["%s_3qt" % (ipmDownDim.name.replace("/", "__"))] = ipmDownP[2]

    return runtable_data


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
parser.add_argument("--stn", help="hutch station", type=int, default=0)
parser.add_argument("--nevents", help="number of events", type=int, default=1e9)
parser.add_argument(
    "--directory",
    help="directory to read files from (def <exp>/hdf5/smalldata)",
    default=None,
)
parser.add_argument(
    "--postElog", help="post plot to elog", action="store_true", default=True
)
parser.add_argument(
    "--postStats",
    help="post summary numbers to run tables",
    action="store_true",
    default=False,
)
# parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-kerb/lgbk/")
parser.add_argument("--url", default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
args = parser.parse_args()
logger.debug("Args to be used for data quality plots: {0}".format(args))


##############################################
## Setup Global parameters and run numbers ###
##############################################
save_elog = args.postElog
detImgMaxSize = 500  # max dimension of image.
expname = args.experiment
run = int(args.run)

if int(os.environ.get("RUN_NUM", "-1")) > 0:
    requests.post(
        os.environ["JID_UPDATE_COUNTERS"],
        json=[{"key": "<b>BeamlineSummary Plots </b>", "value": "Started"}],
    )

######################################
### load data for the chosen run  ####
######################################
# get the ana & anaps objects (from smalldata_tools
if args.directory is not None:
    anaps = sdaps(expname, run, dirname=args.directory)
else:
    anaps = sdaps(expname, run)

ana = anaps.sda  #

## Defining initial selection (laser-on events)
iniFilter = "initial"
ana.addCut("lightStatus/xray", 0.5, 1.5, iniFilter)
ana.addCut("lightStatus/laser", 0.5, 1.5, iniFilter)

### Get data & define axis title&ranges.

ipmUpDim = hv.Dimension(("ipm4/sum", "ipm4 Sum"))
ipmDownDim = hv.Dimension(("ipm5/sum", "ipm5 Sum"))
scatterDim = hv.Dimension(("epix10k2M/ROI_0_sum", "epix10k2M intensity"))
eventTimeDim = hv.Dimension(("eventTimeR", "relative event time"))
l3eDim = hv.Dimension(("l3e", "L3 Energy"))

scanVar = ana.getScanName()
try:
    scanDim = hv.Dimension(("scan/%s" % scanVar, "%s" % scanVar))
except:
    scanDim = None
nevtsDim = hv.Dimension(("nevents", "N events / scan point"))
nevtsLxtDim = hv.Dimension(("neventslxt", "N events / lxt"))

# timing vars.
lxtDim = hv.Dimension(("epics/lxt", "lxt"))

ipmUpVar = ana.getVar(ipmUpDim.name, useFilter=iniFilter)
ipmDownVar = ana.getVar(ipmDownDim.name, useFilter=iniFilter)
stepVar = ana.getVar("scan/varStep", useFilter=iniFilter)
l3eVar = ana.getVar("ebeam/L3_energy", useFilter=iniFilter)
eventTimeRaw = ana.getVar("event_time", useFilter=iniFilter)
eventTime = (eventTimeRaw >> 32).astype(float) + ((eventTimeRaw << 32) >> 32).astype(
    float
) * 1e-9
eventTimeR = eventTime - eventTime[0]

eventTimeRMed = [
    np.nanmedian(eventTimeR[i * 120 : i * 120 + 120])
    for i in range(int(eventTimeR.shape[0] / 120))
]
ipmUpMed = [
    np.nanmedian(ipmUpVar[i * 120 : i * 120 + 120])
    for i in range(int(eventTimeR.shape[0] / 120))
]
ipmDownMed = [
    np.nanmedian(ipmDownVar[i * 120 : i * 120 + 120])
    for i in range(int(eventTimeR.shape[0] / 120))
]

try:
    azav = ana.getVar("epix10k2M/azav_azav", useFilter=iniFilter)
    azav_sum = np.nanmean(azav, axis=0)
    azav_peak = np.argmax(azav_sum)
    if len(azav.shape) > 2:
        azav = np.nanmean(azav, axis=1)
    scatterVar = np.nanmean(
        azav[:, max(0, azav_peak - 50) : min(azav.shape[1], azav_peak + 50)], axis=1
    )
    if len(scatterVar.shape) > 1:
        scatterVar = np.nanmean(scatterVar, axis=1)
except:
    scatterVar = None

### Scan Variable

try:
    isStepScan = np.nanmax(stepVar) > 0
    scanVarBins = np.bincount(stepVar, weights=scatterVar)
    scanNsteps = np.bincount(stepVar)
except:
    isStepScan = False

### Fast delay stage

lxt_fast_his = None
try:
    lxt_fast = ana.getVar("enc/lasDelay", useFilter=iniFilter)
    print(np.nanstd(lxt_fast))
    if lxt_fast is not None and np.nanstd(lxt_fast) < 1e-4:
        lxt_fast_his = np.histogram(
            lxt_fast,
            np.linspace(
                np.nanpercentile(lxt_fast, 1), np.nanpercentile(lxt_fast, 99), 100
            ),
        )
except:
    pass

# droppled sthots.
ana.addCut("lightStatus/xray", -0.5, 0.5, "off")
ana.addCut("evr/code_137", -0.5, 0.5, "hxroff")
if ana.getFilter("hxroff").sum() > ana.getFilter("off").sum():
    offFilter = "hxroff"
else:
    offFilter = "off"
nOff = ana.getFilter(offFilter).sum()

# plots.
ipmUpTime = hv.HexTiles(
    (
        eventTimeR[ipmUpVar < np.nanpercentile(ipmUpVar, 99)],
        ipmUpVar[ipmUpVar < np.nanpercentile(ipmUpVar, 99)],
    ),
    kdims=[eventTimeDim, ipmUpDim],
).opts(cmap="Blues")
ipmUpTimeMed = hv.Points(
    (eventTimeRMed, ipmUpMed), kdims=[eventTimeDim, ipmUpDim], label=ipmUpDim.label
).options(color="r")
ipmDownTimeMed = hv.Points(
    (eventTimeRMed, ipmDownMed), kdims=[eventTimeDim, ipmUpDim], label=ipmDownDim.label
).options(color="m")

ipmTimeLayout = ipmUpTime * ipmUpTimeMed * ipmDownTimeMed

treeSel = l3eVar > np.nanpercentile(l3eVar, 1)
treePlot = hv.HexTiles((l3eVar[treeSel], ipmUpVar[treeSel]), kdims=[l3eDim, ipmUpDim])

ipmPlot = hv.HexTiles((ipmUpVar, ipmDownVar), kdims=[ipmUpDim, ipmDownDim])
ipmLayout = ipmPlot.hist(dimension=[ipmUpDim.name, ipmDownDim.name])
if scatterVar is not None:
    ipmscatterPlot = hv.HexTiles(
        (scatterVar, ipmDownVar), kdims=[scatterDim, ipmDownDim]
    )

stepPlot = None
if isStepScan:
    try:
        stepPlot = hv.Points(
            (scanVarBins / scanNsteps, scanNsteps), kdims=[scanDim, nevtsDim]
        )
    except:
        print(
            "Failed to make stepPlot",
            np.nanmax(stepVar),
            scanNsteps.shape,
            scanVarBins.shape,
        )
        print("DataPoints: ", scanVarBins / scanNsteps)
        print("Dimensions:", scanDim, nevtsDim)
if lxt_fast_his is not None:
    lxtPlot = hv.Points(
        (0.5 * (lxt_fast_his[1][:-1] + lxt_fast_his[1][1:]), lxt_fast_his[0]),
        kdims=[lxtDim, nevtsLxtDim],
    )
else:
    lxtPlot = None

gspec = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="Data Quality - Run %d" % run
)
# gspec[0,0:8] = pn.Row('# Data Quality Plot - Run %04d'%run)
gspec[0:2, 0:8] = pn.Column(ipmTimeLayout)
gspec[2:5, 0:4] = pn.Column(ipmLayout)
gspec[2:5, 4:8] = pn.Column(treePlot)

gspecS = pn.GridSpec(sizing_mode="stretch_both", max_width=700, name="Scan&Scatter")
# gspec[0,0:8] = pn.Row('# Data Quality Plot - Run %04d'%run)
if scatterVar is not None:
    gspecS[0:4, 0:8] = pn.Column(ipmscatterPlot)
maxRow = 4
if stepPlot is not None:
    gspecS[maxRow, 0:8] = pn.Row("## Scan Variable")
    gspecS[maxRow + 1 : maxRow + 3, 0:8] = pn.Column(stepPlot)
    maxRow = 7
if lxtPlot is not None:
    gspecS[maxRow, 0:8] = pn.Row("## Laser - xray Timing")
    gspecS[maxRow + 1 : maxRow + 3, 0:8] = pn.Column(lxtPlot)

# Detector stuff.

detImgs = []
detGrids = []
for detImgName in ana.Keys("Sums"):
    image = ana.fh5.get_node("/%s" % detImgName).read()
    if len(image.shape) > 2:
        if detImgName.find("135") < 0:
            detName = (
                detImgName.replace("Sums/", "")
                .replace("_calib", "")
                .replace("_dropped", "")
                .replace("_square", "")
                .replace("_max", "")
                .replace("_thresADU1", "")
                .replace("_thresADU5", "")
                .replace("_skipFirst", "")
            )
            ix = ana.fh5.get_node("/UserDataCfg/%s/ix" % detName).read()
            iy = ana.fh5.get_node("/UserDataCfg/%s/iy" % detName).read()
            image = image_from_dxy(image, ix, iy)
        else:
            # somehow the epix10k135 has the wrong shape....
            image = image[0]
            # image = image.squeeze()
    if max(image.shape[0], image.shape[1]) > detImgMaxSize:
        rebinFactor = float(detImgMaxSize) / max(image.shape[0], image.shape[1])
        imageR = rebin(
            image,
            [int(image.shape[0] * rebinFactor), int(image.shape[1] * rebinFactor)],
        ) / (ana.getVar("fiducials").shape[0])
    else:
        imageR = image / (ana.getVar("fiducials").shape[0])
    # imgArrays.append(imageR/ana.getVar('fiducials').shape[0])
    imgDim = hv.Dimension(
        ("image", detImgName.replace("Sums/", "").replace("_calib_img", " Mean Image")),
        range=(np.nanpercentile(imageR, 1), np.nanpercentile(imageR, 99.0)),
    )
    detImgs.append(
        hv.Image(imageR, vdims=[imgDim], name=imgDim.label).options(
            colorbar=True, cmap="rainbow"
        )
    )

    detGrid = pn.GridSpec(
        sizing_mode="stretch_both", max_width=700, name=detImgName.replace("Sums/", "")
    )
    detGrid[0, 0] = pn.Row(detImgs[-1])
    detGrids.append(detGrid)

if nOff > 100:
    for detImgName in ana.Keys("Sums"):
        detName = (
            detImgName.replace("_calib", "").replace("_img", "").replace("Sums/", "")
        )
        try:
            common_mode = 0
            if detName.find("epix10k"):
                common_mode = 80
            anaps.AvImage(
                detName,
                useFilter=offFilter,
                numEvts=min(1000, nOff),
                common_mode=common_mode,
            )
        except:
            print("failed to get off shot data for detector %s" % detName)
            continue
        avData = anaps.getAvImage(detName)[1]
        try:
            image = anaps.__dict__[detName].det.image(run, avData)
        except:
            print("failed to make image for detector %s" % detName)
            continue
        if max(image.shape[0], image.shape[1]) > detImgMaxSize:
            rebinFactor = float(detImgMaxSize) / max(image.shape[0], image.shape[1])
            imageR = rebin(
                image,
                [int(image.shape[0] * rebinFactor), int(image.shape[1] * rebinFactor)],
            )
        else:
            imageR = image
        imgOffDim = hv.Dimension(
            (
                "image_off",
                detImgName.replace("Sums/", "").replace(
                    "_calib_img", " Mean Image Off"
                ),
            ),
            range=(np.nanpercentile(imageR, 1), np.nanpercentile(imageR, 99.0)),
        )
        detImgs.append(
            hv.Image(imageR, vdims=[imgOffDim], name=imgOffDim.label).options(
                colorbar=True, cmap="rainbow"
            )
        )

        detGrid = pn.GridSpec(
            sizing_mode="stretch_both",
            max_width=700,
            name="%s, dropped shots" % detName,
        )
        detGrid[0, 0] = pn.Row(detImgs[-1])
        detGrids.append(detGrid)

tabs = pn.Tabs(gspec)
tabs.append(gspecS)
for detGrid in detGrids:
    tabs.append(detGrid)

if int(os.environ.get("RUN_NUM", "-1")) > 0:
    requests.post(
        os.environ["JID_UPDATE_COUNTERS"],
        json=[{"key": "<b>BeamlineSummary Plots </b>", "value": "Done"}],
    )

SIT_PSDM_DATA = os.getenv("SIT_PSDM_DATA")
elogDir = f"{SIT_PSDM_DATA}/{expname[:3]}/{expname}/stats/summary/BeamlineSummary/BeamlineSummary_Run{run:04d}"

if save_elog:
    from summaries.summary_utils import prepareHtmlReport

    pageTitleFormat = "BeamlineSummary/BeamlineSummary_Run{run:04d}"
    prepareHtmlReport(tabs, expname, run, pageTitleFormat)

    if int(os.environ.get("RUN_NUM", "-1")) > 0:
        requests.post(
            os.environ["JID_UPDATE_COUNTERS"],
            json=[{"key": "<b>BeamlineSummary Plots </b>", "value": "Posted"}],
        )

if args.postStats:
    if scanVar == "":
        encDelay = ana.getVar("enc/lasDelay")
        delta_encDelay = np.nanmax(encDelay) - np.nanmin(encDelay)
        if delta_encDelay > 0.5:
            scanVar = "delay"
    elif scanVar.find("lxt"):
        scanVar = "delay"
    runtable_data = makeRunTableData(ana, ipmUpDim, ipmDownDim, iniFilter, scanVar)
    postRunTable(runtable_data)
