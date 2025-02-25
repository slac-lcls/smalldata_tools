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
from pathlib import Path
from requests.auth import HTTPBasicAuth
import holoviews as hv
from holoviews import dim

hv.extension("bokeh")
pn.extension()
import sys

fpath = os.path.dirname(os.path.abspath(__file__))
fpathup = "/".join(fpath.split("/")[:-1])
try:
    fpath = os.environ.get("MYDIR", fpathup).replace("/arp_scripts", "")
except:
    fpath = fpathup
sys.path.append(fpath)
from smalldata_tools.lcls1.SmallDataAna_psana import SmallDataAna_psana as sdaps
from smalldata_tools.utilities import rebin
from smalldata_tools.utilities import evtt2Rt
from smalldata_tools.utilities import postRunTable

def makeRunTableData(ana, ipmUpDim, ipmDownDim, Filter):
    n162 = ana.getVar("evr/code_162").sum()
    ana.addCut("evr/code_162", -0.5, 0.5, "xon")
    ana.addCut("evr/code_137", 0.5, 1.5, "xon")
    nOff = ana.getFilter("xon").shape[0] - ana.getFilter("xon").sum()
    # data to be posted to the run table if so requested.
    runtable_data = {"N dropped Shots": int(nOff), "N BYKIK 162": int(n162)}

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


logging.basicConfig(level=logging.DEBUG)
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

SIT_PSDM_DATA = Path(os.environ.get("SIT_PSDM_DATA"))

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

ana = anaps.sda  # ?? what is the ana object?
## can add some info about the ana and anaps objects here, where to look it up

## Defining basic filter for beam ON events for data loading
iniFilter = "initial"
ana.addCut("lightStatus/xray", -1, 1.5, iniFilter)
ana.addCut("lightStatus/laser", -1, 1.5, iniFilter)

## define shot time for each event starting from the start of the run
eventTimeR = evtt2Rt(ana.getVar("event_time", useFilter=iniFilter))
eventTimeRMed = [
    np.nanmedian(eventTimeR[i * 120 : i * 120 + 120])
    for i in range(int(eventTimeR.shape[0] / 120))
]  # median value of the event time per second ... why are we doing this??

### define data, data label, and load
ipm2Dim = hv.Dimension(("ipm2/sum", "ipm2 Sum"))
ipm3Dim = hv.Dimension(("ipm3/sum", "ipm3 Sum"))
eventTimeDim = hv.Dimension(("eventTimeR", "relative event time (s)"))
l3eDim = hv.Dimension(("ebeam/L3_energy", "L3 Energy"))

## Define point detector values to be plotted
ipm2Var = ana.getVar(ipm2Dim.name, useFilter=iniFilter)
ipm3Var = ana.getVar(ipm3Dim.name, useFilter=iniFilter)
## obtain median values for IPM2 and IPM3 for 'average trace'
ipm2Med = [
    np.nanmedian(ipm2Var[i * 120 : i * 120 + 120])
    for i in range(int(eventTimeR.shape[0] / 120))
]
ipm3Med = [
    np.nanmedian(ipm3Var[i * 120 : i * 120 + 120])
    for i in range(int(eventTimeR.shape[0] / 120))
]

## Define the plot contents
ipm2vTime = hv.Scatter(
    (
        eventTimeR[ipm2Var < np.nanpercentile(ipm2Var, 99)],
        ipm2Var[ipm2Var < np.nanpercentile(ipm2Var, 99)],
    ),
    kdims=[eventTimeDim, ipm2Dim],
).opts(color="g", alpha=0.2)
ipm2vTimeMed = hv.Curve(
    (eventTimeRMed, ipm2Med), kdims=[eventTimeDim, ipm2Dim], label=ipm2Dim.label
).options(color="k")
ipm3vTimeMed = hv.Curve(
    (eventTimeRMed, ipm3Med), kdims=[eventTimeDim, ipm3Dim], label=ipm3Dim.label
).options(color="r")

# create a multi trace plot that contains the 3 datasets
ipmTimeLayout = ipm2vTime * ipm2vTimeMed * ipm3vTimeMed

## IPM2 vs L3E plot (or as we call the TREE plot)
try:
    l3eVar = ana.getVar(l3eDim.name, useFilter=iniFilter)
    treeSel = l3eVar > np.nanpercentile(l3eVar, 1)
    treePlot = hv.HexTiles(
        (l3eVar[treeSel], ipm2Var[treeSel]), kdims=[l3eDim, ipm2Dim]
    ).opts(cmap="Blues", min_count=0)
except:
    treePlot = None

#########################################
## defining the first Tab
#########################################
tabIPM = pn.GridSpec(
    sizing_mode="stretch_both",
    max_width=800,
    name="Data Quality - Run %d, Beamline" % run,
)
tabIPM[0:2, 0:8] = pn.Column(ipmTimeLayout)
tabIPM[2:4, 0:8] = pn.Column(treePlot)

tabs = pn.Tabs(tabIPM)  # IPM tab

############################################
## TT data tab
############################################

ttPosDim = hv.Dimension(("tt/FLTPOS", "tt position (pixel)"))
ttAmpDim = hv.Dimension(("tt/AMPL", "tt amplitude (%)"))
if ttPosDim.name in ana.Keys():
    ttPos = ana.getVar(ttPosDim.name, useFilter=iniFilter)
    ttAmp = ana.getVar(ttAmpDim.name, useFilter=iniFilter)
    ttPosvTime = hv.HexTiles(
        (eventTimeR, ttPos), kdims=[eventTimeDim, ttPosDim]
    ).options(cmap="summer")
    ttAmpvIpm2 = hv.Scatter((ipm2Var, ttAmp), kdims=[ipm2Dim, ttAmpDim]).options(
        color="g", alpha=0.3
    )

    tabTT = pn.GridSpec(sizing_mode="stretch_both", max_width=800, name="Time Tool")
    tabTT[0:2, 0:6] = pn.Column(ttPosvTime)
    tabTT[2:5, 0:6] = pn.Column(ttAmpvIpm2)
    tabs.append(tabTT)

##############################################################
## User Diode tabs with time plots and correlation plots
##############################################################

diodeUDim = hv.Dimension(("diodeU/channels", "ipmU"))
if diodeUDim.name in ana.Keys():
    diodeU = ana.getVar(diodeUDim.name, useFilter=iniFilter)
    diodeUMed = np.array(
        [
            np.nanmedian(diodeU[i * 120 : i * 120 + 120, :], axis=0)
            for i in range(int(eventTimeR.shape[0] / 120))
        ]
    )
    diodeU0vTime = hv.Scatter(
        (eventTimeR, diodeU[:, 0]), kdims=[eventTimeDim, diodeUDim]
    ).options(color="g", alpha=0.3)
    diodeU0vTimeMed = hv.Curve(
        (eventTimeRMed, diodeUMed[:, 0]), kdims=[eventTimeDim, diodeUDim]
    ).options(color="k")
    diodeU1vTime = hv.Scatter(
        (eventTimeR, diodeU[:, 1]), kdims=[eventTimeDim, diodeUDim]
    ).options(color="g", alpha=0.3)
    diodeU1vTimeMed = hv.Curve(
        (eventTimeRMed, diodeUMed[:, 1]), kdims=[eventTimeDim, diodeUDim]
    ).options(color="k")
    diodeU2vTime = hv.Scatter(
        (eventTimeR, diodeU[:, 2]), kdims=[eventTimeDim, diodeUDim]
    ).options(color="g", alpha=0.3)
    diodeU2vTimeMed = hv.Curve(
        (eventTimeRMed, diodeUMed[:, 2]), kdims=[eventTimeDim, diodeUDim]
    ).options(color="k")
    diodeU3vTime = hv.Scatter(
        (eventTimeR, diodeU[:, 3]), kdims=[eventTimeDim, diodeUDim]
    ).options(color="g", alpha=0.3)
    diodeU3vTimeMed = hv.Curve(
        (eventTimeRMed, diodeUMed[:, 3]), kdims=[eventTimeDim, diodeUDim]
    ).options(color="k")

    diodeU0vIpm2 = hv.Scatter(
        (ipm2Var, diodeU[:, 0]), kdims=[ipm2Dim, diodeUDim]
    ).options(color="g", alpha=0.3)
    diodeU1vIpm2 = hv.Scatter(
        (ipm2Var, diodeU[:, 1]), kdims=[ipm2Dim, diodeUDim]
    ).options(color="g", alpha=0.3)
    diodeU2vIpm2 = hv.Scatter(
        (ipm2Var, diodeU[:, 2]), kdims=[ipm2Dim, diodeUDim]
    ).options(color="g", alpha=0.3)
    diodeU3vIpm2 = hv.Scatter(
        (ipm2Var, diodeU[:, 3]), kdims=[ipm2Dim, diodeUDim]
    ).options(color="g", alpha=0.3)

    ## defining subplots
    tabDiodeU = pn.GridSpec(sizing_mode="stretch_both", max_width=1200, name="IPM User")

    # time plots for the 4 user IPM channels
    tabDiodeU[0:1, 0:4] = pn.Column(diodeU0vTime * diodeU0vTimeMed)
    tabDiodeU[1:2, 0:4] = pn.Column(diodeU1vTime * diodeU1vTimeMed)
    tabDiodeU[2:3, 0:4] = pn.Column(diodeU2vTime * diodeU2vTimeMed)
    tabDiodeU[3:4, 0:4] = pn.Column(diodeU3vTime * diodeU3vTimeMed)

    # correlation plots for the 4 user IPM channels
    tabDiodeU[4:6, 0:2] = pn.Column(diodeU0vIpm2)
    tabDiodeU[4:6, 2:4] = pn.Column(diodeU1vIpm2)
    tabDiodeU[6:8, 0:2] = pn.Column(diodeU2vIpm2)
    tabDiodeU[6:8, 2:4] = pn.Column(diodeU3vIpm2)
    tabs.append(tabDiodeU)


diode2Dim = hv.Dimension(("diode2/channels", "PIM3"))
if diode2Dim.name in ana.Keys():
    diode2 = ana.getVar(diode2Dim.name, useFilter=iniFilter)
    diode2Med = np.array(
        [
            np.nanmedian(diode2[i * 120 : i * 120 + 120, :], axis=0)
            for i in range(int(eventTimeR.shape[0] / 120))
        ]
    )

    diode20vTime = hv.Scatter(
        (eventTimeR, diode2[:, 0]), kdims=[eventTimeDim, diode2Dim]
    ).options(color="g", alpha=0.3)
    diode20vTimeMed = hv.Curve(
        (eventTimeRMed, diode2Med[:, 0]), kdims=[eventTimeDim, diode2Dim]
    ).options(color="k")
    diode21vTime = hv.Scatter(
        (eventTimeR, diode2[:, 1]), kdims=[eventTimeDim, diode2Dim]
    ).options(color="g", alpha=0.3)
    diode21vTimeMed = hv.Curve(
        (eventTimeRMed, diode2Med[:, 1]), kdims=[eventTimeDim, diode2Dim]
    ).options(color="k")
    diode22vTime = hv.Scatter(
        (eventTimeR, diode2[:, 2]), kdims=[eventTimeDim, diode2Dim]
    ).options(color="g", alpha=0.3)
    diode22vTimeMed = hv.Curve(
        (eventTimeRMed, diode2Med[:, 2]), kdims=[eventTimeDim, diode2Dim]
    ).options(color="k")
    diode23vTime = hv.Scatter(
        (eventTimeR, diode2[:, 3]), kdims=[eventTimeDim, diode2Dim]
    ).options(color="g", alpha=0.3)
    diode23vTimeMed = hv.Curve(
        (eventTimeRMed, diode2Med[:, 3]), kdims=[eventTimeDim, diode2Dim]
    ).options(color="k")

    diode20vIpm2 = hv.Scatter(
        (ipm2Var, diode2[:, 0]), kdims=[ipm2Dim, diode2Dim]
    ).options(color="g", alpha=0.3)
    diode21vIpm2 = hv.Scatter(
        (ipm2Var, diode2[:, 1]), kdims=[ipm2Dim, diode2Dim]
    ).options(color="g", alpha=0.3)
    diode22vIpm2 = hv.Scatter(
        (ipm2Var, diode2[:, 2]), kdims=[ipm2Dim, diode2Dim]
    ).options(color="g", alpha=0.3)
    diode23vIpm2 = hv.Scatter(
        (ipm2Var, diode2[:, 3]), kdims=[ipm2Dim, diode2Dim]
    ).options(color="g", alpha=0.3)

    ## defining subplots
    tabDiode2 = pn.GridSpec(sizing_mode="stretch_both", max_width=1200, name="PIM3")

    # time plots for the 4 PIM3 channels
    tabDiode2[0:1, 0:4] = pn.Column(diode20vTime * diode20vTimeMed)
    tabDiode2[1:2, 0:4] = pn.Column(diode21vTime * diode21vTimeMed)
    tabDiode2[2:3, 0:4] = pn.Column(diode22vTime * diode22vTimeMed)
    tabDiode2[3:4, 0:4] = pn.Column(diode23vTime * diode23vTimeMed)

    # correlation plots for the 4 PIM3 channels
    tabDiode2[4:6, 0:2] = pn.Column(diode20vIpm2)
    tabDiode2[4:6, 2:4] = pn.Column(diode21vIpm2)
    tabDiode2[6:8, 0:2] = pn.Column(diode22vIpm2)
    tabDiode2[6:8, 2:4] = pn.Column(diode23vIpm2)
    tabs.append(tabDiode2)

if ana.getFilter("hxroff").sum() > ana.getFilter("off").sum():
    offFilter = "hxroff"
else:
    offFilter = "off"
nOff = ana.getFilter(
    offFilter
).sum()  # total number of OFF shots, if this number is larger than 100, will display a tab for the average Dark.


####################################################
## Area Detector Example Images as the next tabs
####################################################
#
# under development, should be shared code.
#

##################################
## save the html file
##################################

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
    print("posting to the run tables - ipm values.")
    runtable_data = makeRunTableData(ana, ipm2Dim, ipm3Dim, iniFilter)
    postRunTable(runtable_data)
