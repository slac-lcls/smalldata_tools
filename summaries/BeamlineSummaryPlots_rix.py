#!/usr/bin/env python
########################################################################
## load tools and packages required for data handling and plotting
########################################################################

import tables
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

fpath = os.path.dirname(os.path.abspath(__file__))
fpathup = "/".join(fpath.split("/")[:-1])
try:
    fpath = os.environ.get("MYDIR", fpathup).replace("/arp_scripts", "")
except:
    fpath = fpathup
sys.path.append(fpath)
from smalldata_tools.lcls2.SmallDataAna_psana import SmallDataAna_psana as sdaps
from smalldata_tools.utilities import rebin


def click_policy(plot, element):
    plot.state.legend.click_policy = "hide"


## function that chops the 64 bit time integer into soemthing a bit more realistic
def evtt2Rt(event_time):
    evtt0 = event_time >> 32
    evtt1 = (event_time << 32) >> 32
    evtt_sec = evtt0.astype(float)
    evtt_ns = evtt1.astype(float) * 1e-9
    Rt = evtt_sec + evtt_ns
    Rt = Rt - Rt[0]
    return Rt


def postRunTable(runtable_data):
    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(args.experiment)
    print("URL:", ws_url)
    user = args.experiment[:3] + "opr"
    elogPostFile = "/cds/home/opr/%s/forElogPost.txt" % user
    hostname = socket.gethostname()
    if hostname.find("sdf") >= 0:
        elogPostFile = "/sdf/group/lcls/ds/tools/forElogPost.txt"
    with open(elogPostFile, "r") as reader:
        answer = reader.readline()
    r = requests.post(
        ws_url,
        params={"run_num": args.run},
        json=runtable_data,
        auth=HTTPBasicAuth(args.experiment[:3] + "opr", answer[:-1]),
    )
    # we might need to use this for non=current expetiments. Currently does not work in ARP
    # krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    # r = requests.post(ws_url, headers=krbheaders, params={"run_num": args.run}, json=runtable_data)
    print(r)


def makeRunTableData(dat):
    xray = dat.lightStatus.xray.read()
    laser = dat.lightStatus.laser.read()
    on = xray.astype(bool) & laser.astype(bool)
    off = xray.astype(bool) & ~laser.astype(bool)
    dropped = ~xray.astype(bool)
    total = xray.shape[0]

    runtable_data = {
        "N dropped Shots": int(dropped.sum()),
        "N total": int(total),
        "N laser on": int(on.sum()),
        "N laser off": int(off.sum()),
    }

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
dat = tables.open_file(
    "/reg/d/psdm/%s/%s/hdf5/smalldata/%s_Run%04d.h5"
    % (expname[:3], expname, expname, run)
).root

xray = dat.lightStatus.xray.read()
laser = dat.lightStatus.laser.read()
xoff = ~xray.astype(bool)
on = xray.astype(bool) & laser.astype(bool)
off = xray.astype(bool) & ~laser.astype(bool)
# iniFilter='initial'
# ana.addCut('lightStatus/xray',0.5,1.5,iniFilter)
# ana.addCut('lightStatus/laser',0.5,1.5,iniFilter)

andor_roi = slice(1050, 1150)
andor_bg_roi = slice(400, 700)
fim_roi = slice(102, 120)
fim_bg_roi = slice(0, 50)
fim2_roi = slice(55, 70)
fim2_bg_roi = slice(0, 20)
# rawfims = [d for d in dir(dat.rix_fim0_raw) if d[0]!='_']

fimDim = hv.Dimension(("fim0", "fim"))
fim1Dim = hv.Dimension(("fim1", "fim1"))
fim2Dim = hv.Dimension(("fim2", "fim2"))
scatterDim = hv.Dimension(("andor", "andor intensity"))
eventTimeDim = hv.Dimension(("eventTimeR", "relative event time"))
l3eDim = hv.Dimension(("l3e", "L3 Energy"))
ttposDim = hv.Dimension(("ttpos", "tt fitted position"), range=(0, 1000))
ttamplDim = hv.Dimension(("ttampl", "tt amplitude"))

# fim0All = -1.*dat.det_rix_fim0.full_fimSum_fimSum.read().sum(axis=1)
fim0 = -1.0 * dat.det_rix_fim0.full_sigbkg_sum.read()[:, 4:].sum(axis=1)
fim1 = -1.0 * dat.det_rix_fim1.full_sigbkg_sum.read()[:, 5]

eventTimeRaw = dat.timestamp.read()
eventTime = (eventTimeRaw >> 32).astype(float) + ((eventTimeRaw << 32) >> 32).astype(
    float
) * 1e-9
eventTimeR = eventTime - eventTime[0]

l3eVar = dat.ebeam.ebeamL3Energy.read()

i0Var = fim0
i0Dim = fimDim

eventTimeRMed = [
    np.nanmedian(eventTimeR[i * 120 : i * 120 + 120])
    for i in range(int(eventTimeR.shape[0] / 120))
]
i0Med = [
    np.nanmedian(i0Var[i * 120 : i * 120 + 120])
    for i in range(int(eventTimeR.shape[0] / 120))
]

# this should be the andor when present
# scatterVar = dat.gmd.avgIntensity.read()
try:
    scatterVar = (
        dat.andor_dir.ROI_sig_sum.read() - dat.andor_dir.ROI_bkg_sum.read() / 550 * 400
    )
except:
    scatterVar = None

i0Time = hv.HexTiles(
    (
        eventTimeR[i0Var < np.nanpercentile(i0Var, 99)],
        i0Var[i0Var < np.nanpercentile(i0Var, 99)],
    ),
    kdims=[eventTimeDim, i0Dim],
).opts(cmap="Blues")
i0TimeMed = hv.Points(
    (eventTimeRMed, i0Med), kdims=[eventTimeDim, i0Dim], label=i0Dim.label
).options(color="r")

ipmTimeLayout = i0Time * i0TimeMed

treeSel = l3eVar > np.nanpercentile(l3eVar, 1)
treePlot = hv.HexTiles((l3eVar[treeSel], i0Var[treeSel]), kdims=[l3eDim, i0Dim])

try:
    andor_trace = np.nanmean(dat.andor_dir.full_area, axis=0)
    andorDim = hv.Dimension(("andor_dir", "AndOR dir"))
    pixelDim = hv.Dimension(("pixel", "pixel"))
    andorPlot = hv.Curve(
        (np.arange(andor_trace.shape[0]), andor_trace),
        kdims=[pixelDim],
        vdims=[andorDim],
    )
except:
    andorPlot = None

ipmPlot = hv.HexTiles((fim0, fim1), kdims=[fimDim, fim1Dim])
ipmLayout = ipmPlot.hist(dimension=[fimDim.name, fim1Dim.name])
if scatterVar is not None:
    ipmscatterPlot = hv.HexTiles((scatterVar, i0Var), kdims=[scatterDim, i0Dim])
fim1Plot = None
if (fim1 is not None) and (scatterVar is not None):
    fim1Plot = hv.HexTiles((scatterVar, fim1), kdims=[scatterDim, fim1Dim])

try:
    ttpos = dat.tt.fltpos.read()
    ttampl = dat.tt.ampl.read()
    if eventTimeR[on].shape[0] < 5000:
        ttposPlot = hv.Scatter(
            (eventTimeR[on], ttpos[on]), kdims=[eventTimeDim, ttposDim]
        )
        ttamplPlot = hv.Scatter((ttampl[on], i0Var[on]), kdims=[ttamplDim, i0Dim])
        ttamplposPlot = hv.Scatter((ttpos[on], ttampl[on]), kdims=[ttposDim, ttamplDim])
    else:
        ttposPlot = hv.HexTiles(
            (eventTimeR[on], ttpos[on]), kdims=[eventTimeDim, ttposDim]
        )
        ttamplPlot = hv.HexTiles((ttampl[on], i0Var[on]), kdims=[ttamplDim, i0Dim])
        ttamplposPlot = hv.HexTiles(
            (ttpos[on], ttampl[on]), kdims=[ttposDim, ttamplDim]
        )
except:
    ttpos = None

fimDims = {
    "ch%d" % i: hv.Dimension(("fim_ch%d" % i, "channel %d" % i))
    for i in np.arange(dat.det_rix_fim0.full_area.shape[1])
}
fimDims["legend"] = hv.Dimension(("fim_legend", "legend"))
fimTimeDim = hv.Dimension(("fim_time", "fim points"))
fimsums = {
    "fim%d" % ifim: getattr(dat, "det_rix_fim%d" % ifim).full_area.read().sum(axis=0)
    for ifim in np.arange(3)
}


fimCurves = {
    "channel_%d"
    % i: hv.NdOverlay(
        {
            fn: hv.Curve(
                (np.arange(256), fimsums[fn][i]),
                kdims=[fimTimeDim],
                vdims=[fimDims["ch%d" % i]],
            ).opts(show_legend=False)
            for fn in fimsums
        }
    )
    for i in np.arange(fimsums["fim0"].shape[0])
}
fimCurves["legend"] = hv.NdOverlay(
    {
        fn: hv.Curve(
            (np.arange(256), np.zeros(256)),
            kdims=[fimTimeDim],
            vdims=[fimDims["legend"]],
        ).opts(hooks=[click_policy])
        for fn in fimsums
    }
)

gspecFim = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="Summed FIM traces - Run %d" % run
)
for k, plot in fimCurves.items():
    try:
        irow = int(int(k[-1]) / 3)
        icol = int(k[-1]) % 3
    except:
        irow = 2
        icol = 2
    gspecFim[irow * 3 : (irow + 1) * 3, icol * 3 : (icol + 1) * 3] = pn.Column(
        plot.opts(legend_position="bottom_right")
    )

gspec = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="Data Quality - Run %d" % run
)
# gspec[0,0:8] = pn.Row('# Data Quality Plot - Run %04d'%run)
gspec[0:2, 0:8] = pn.Column(ipmTimeLayout)
gspec[2:6, 0:4] = pn.Column(ipmLayout)
gspec[2:6, 4:8] = pn.Column(treePlot)

gspecS = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="Andor, Scan&Scatter"
)
# gspec[0,0:8] = pn.Row('# Data Quality Plot - Run %04d'%run)
maxRow = 0
if andorPlot is not None:
    gspecS[0:2, 0:8] = pn.Column(andorPlot)
    maxRow += 2
if scatterVar is not None:
    gspecS[maxRow : maxRow + 4, 0:8] = pn.Column(ipmscatterPlot)
    maxRow += 4
if fim1Plot is not None:
    # gspecS[maxRow,0:8] = pn.Row('## Scan Variable')
    gspecS[maxRow : maxRow + 4, 0:8] = pn.Column(fim1Plot)
    maxRow += 4
# if stepPlot is not None:
#    gspecS[maxRow,0:8] = pn.Row('## Scan Variable')
#    gspecS[maxRow+1:maxRow+3,0:8] = pn.Column(stepPlot)
#    maxRow=7
# if lxtPlot is not None:
#    gspecS[maxRow,0:8] = pn.Row('## Laser - xray Timing')
#    gspecS[maxRow+1:maxRow+3,0:8] = pn.Column(lxtPlot)

tabs = pn.Tabs(gspec)
tabs.append(gspecFim)
# if maxRow>0:
#    tabs.append(gspecS)

if ttpos is not None:
    gspecT = pn.GridSpec(sizing_mode="stretch_both", max_width=700, name="Timetool")
    gspecT[0:4, 0:8] = pn.Column(ttposPlot)
    gspecT[4:8, 0:4] = pn.Column(ttamplposPlot)
    gspecT[4:8, 4:8] = pn.Column(ttamplPlot)
    tabs.append(gspecT)
# for detGrid in detGrids:

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
    runtable_data = makeRunTableData(dat)
    postRunTable(runtable_data)
