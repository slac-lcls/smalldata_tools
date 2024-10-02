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
from smalldata_tools.lcls1.SmallDataAna_psana import SmallDataAna_psana as sdaps
from smalldata_tools.utilities import image_from_dxy
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
    # xray = dat.lightStatus.xray.read()
    # data was 'broken' and we did not reprocess yet.
    xray = dat.evr.core_137.read()
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

if args.directory is None:
    args.directory = "/sdf/data/lcls/ds/%s/%s/hdf5/smalldata" % (
        args.experiment[:3],
        args.experiment,
    )
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
######################################
### load data for the chosen run  ####
######################################
# get the ana & anaps objects (from smalldata_tools
if args.directory is not None:
    anaps = sdaps(expname, run, dirname=args.directory)
else:
    anaps = sdaps(expname, run)

ana = anaps.sda  #

dat = tables.open_file("%s/%s_Run%04d.h5" % (args.directory, expname, run)).root

xray = dat.lightStatus.xray.read()
laser = dat.lightStatus.laser.read()
xoff = ~xray.astype(bool)
on = xray.astype(bool) & laser.astype(bool)
off = xray.astype(bool) & ~laser.astype(bool)
# iniFilter='initial'
# ana.addCut('lightStatus/xray',0.5,1.5,iniFilter)
# ana.addCut('lightStatus/laser',0.5,1.5,iniFilter)

gdetDim = hv.Dimension(("gdet", "gdet"))
dg2Dim = hv.Dimension(("dg2", "dg2 bmmon"))
scatterDim = hv.Dimension(("jungfrau", "jungfrau intensity"))
eventTimeDim = hv.Dimension(("eventTimeR", "relative event time"))
ttposDim = hv.Dimension(("ttpos", "tt fitted position"), range=(0, 1000))
ttamplDim = hv.Dimension(("ttampl", "tt amplitude"))
ttfwhmDim = hv.Dimension(("ttfwhm", "tt step width"))

# fim0All = -1.*dat.det_rix_fim0.full_fimSum_fimSum.read().sum(axis=1)
gdet = dat.gas_detector.f_11_ENRC.read()
try:
    dg2ipm = dat.ipm_dg2.sum.read()
except:
    dg2ipm = dat.ipm_hfx_dg2.sum.read()

eventTimeRaw = dat.event_time.read()
eventTime = (eventTimeRaw >> 32).astype(float) + ((eventTimeRaw << 32) >> 32).astype(
    float
) * 1e-9
eventTimeR = eventTime - eventTime[0]

i0Var = gdet
i0Dim = dg2Dim

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
    scatterVar = dat.jungfrau4M.Full_thres_sum.read()
except:
    scatterVar = None

##droppled sthots.
ana.addCut("lightStatus/xray", -0.5, 0.5, "off")
ana.addCut("evr/code_137", -0.5, 0.5, "hxroff")
offFilter = "hxroff"
# if ana.getFilter('hxroff').sum() >  ana.getFilter('off').sum():
#    offFilter = 'hxroff'
# else:
#    offFilter = 'off'
nOff = ana.getFilter(offFilter).sum()

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

try:
    j4m_trace = np.nanmean(dat.jungfrau4M.azav_azav.read(), axis=0).squeeze()
    j4m_qbins = dat.UserDataCfg.jungfrau4M.azav__azav_qbins.read()[1:]
    j4mDim = hv.Dimension(("j4m_dir", "jungfrau img"))
    pixelDim = hv.Dimension(("q", "q [angstrom^-1]"))
    j4mPlot = hv.Curve(
        (np.arange(j4m_qbins), j4m_trace), kdims=[pixelDim], vdims=[j4mDim]
    )
except:
    j4mPlot = None

ipmPlot = hv.HexTiles((gdet, dg2ipm), kdims=[gdetDim, dg2Dim])
ipmLayout = ipmPlot.hist(dimension=[gdetDim.name, dg2Dim.name])
if scatterVar is not None:
    ipmscatterPlot = hv.HexTiles((scatterVar, i0Var), kdims=[scatterDim, i0Dim])

try:
    ttpos = dat.tt.FLTPOS.read()
    ttampl = dat.tt.AMPL.read()
    ttfwhm = dat.tt.FLTPOSFWHM.read()
    if eventTimeR[on].shape[0] < 5000:
        ttposPlot = hv.Scatter(
            (eventTimeR[on], ttpos[on]), kdims=[eventTimeDim, ttposDim]
        )
        ttamplPlot = hv.Scatter((ttampl[on], i0Var[on]), kdims=[ttamplDim, i0Dim])
        ttamplposPlot = hv.Scatter((ttpos[on], ttampl[on]), kdims=[ttposDim, ttamplDim])
        ttfwhmposPlot = hv.Scatter((ttpos[on], ttfwhm[on]), kdims=[ttposDim, ttfwhmDim])
    else:
        ttposPlot = hv.HexTiles(
            (eventTimeR[on], ttpos[on]), kdims=[eventTimeDim, ttposDim]
        )
        ttamplPlot = hv.HexTiles((ttampl[on], i0Var[on]), kdims=[ttamplDim, i0Dim])
        ttamplposPlot = hv.HexTiles(
            (ttpos[on], ttampl[on]), kdims=[ttposDim, ttamplDim]
        )
        ttfwhmposPlot = hv.HexTiles(
            (ttpos[on], ttfwhm[on]), kdims=[ttposDim, ttfwhmDim]
        )
except:
    ttpos = None


# Detector stuff.

detImgs = []
detGrids = []
for detImgName in ana.Keys("Sums"):
    if detImgName.find("Acqiris") >= 0:
        continue
    image = ana.fh5.get_node("/%s" % detImgName).read()
    if len(image.shape) > 2:
        if detImgName.find("135") < 0:
            # detName = detImgName.replace('Sums/','').replace('_calib','')
            detName = detImgName.replace("Sums/", "").split("_")[0]
            ix = ana.fh5.get_node("/UserDataCfg/%s/ix" % detName).read()
            iy = ana.fh5.get_node("/UserDataCfg/%s/iy" % detName).read()
            image = image_from_dxy(image, ix, iy)
        else:
            # somehow the epix10k135 has the wrong shape....
            image = image[0]
            # image = image.squeeze()
    if max(image.shape[0], image.shape[1]) > detImgMaxSize:
        rebinFactor = float(detImgMaxSize) / max(image.shape[0], image.shape[1])
        # print('rebinning: ',rebinFactor, image.shape, detImgName)
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

    detGrid = pn.GridSpec(max_width=700, name=detImgName.replace("Sums/", ""))
    detGrid[0, 0] = pn.Row(detImgs[-1])
    detGrids.append(detGrid)

if nOff > 100:
    detNames = []
    for detImgName in ana.Keys("Sums"):
        detName = (
            detImgName.replace("_calib", "")
            .replace("_img", "")
            .replace("_dropped", "")
            .replace("_square", "")
            .replace("Sums/", "")
        )
        if detName in detNames:
            continue
        detNames.append(detName)
        try:
            common_mode = getattr(dat.UserDataCfg, detName).common_mode.read()
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
            sizing_mode="fixed", max_width=700, name="%s, dropped shots" % detName
        )
        detGrid[0, 0] = pn.Row(detImgs[-1])
        detGrids.append(detGrid)


gspec = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="Data Quality - Run %d" % run
)
# gspec[0,0:8] = pn.Row('# Data Quality Plot - Run %04d'%run)
gspec[0:2, 0:8] = pn.Column(ipmTimeLayout)
gspec[2:6, 0:4] = pn.Column(ipmLayout)

gspecS = pn.GridSpec(
    sizing_mode="stretch_both", max_width=700, name="jungfrau4M, Scan&Scatter"
)
# gspec[0,0:8] = pn.Row('# Data Quality Plot - Run %04d'%run)
maxRow = 0
if j4mPlot is not None:
    gspecS[0:2, 0:8] = pn.Column(j4mPlot)
    maxRow += 2
if scatterVar is not None:
    gspecS[maxRow : maxRow + 4, 0:8] = pn.Column(ipmscatterPlot)
    maxRow += 4
# if stepPlot is not None:
#    gspecS[maxRow,0:8] = pn.Row('## Scan Variable')
#    gspecS[maxRow+1:maxRow+3,0:8] = pn.Column(stepPlot)
#    maxRow=7
# if lxtPlot is not None:
#    gspecS[maxRow,0:8] = pn.Row('## Laser - xray Timing')
#    gspecS[maxRow+1:maxRow+3,0:8] = pn.Column(lxtPlot)


tabs = pn.Tabs(gspec)
# if maxRow>0:
#    tabs.append(gspecS)

if ttpos is not None:
    gspecT = pn.GridSpec(sizing_mode="stretch_both", max_width=700, name="Timetool")
    gspecT[0:4, 0:8] = pn.Column(ttposPlot)
    gspecT[4:8, 0:4] = pn.Column(ttamplposPlot)
    gspecT[4:8, 4:8] = pn.Column(ttamplPlot)
    gspecT[8:12, 0:4] = pn.Column(ttfwhmposPlot)
    tabs.append(gspecT)

for detGrid in detGrids:
    tabs.append(detGrid)

##################################
## save the html file
##################################

from pathlib import Path

SIT_PSDM_DATA = Path(os.environ.get("SIT_PSDM_DATA"))
runnum = int(args.run)
elogDir = (
    Path(SIT_PSDM_DATA)
    / expname[:3]
    / expname
    / f"stats/summary/BeamlineSummary/BeamlineSummary_Run{runnum:03d}"
)
if save_elog:
    import os

    if not os.path.isdir(elogDir):
        os.makedirs(elogDir)
    print("Made Directory to save data:", elogDir)
    # gspec.save(('%s/report.html'%elogDir))
    tabs.save(("%s/report.html" % elogDir))

if int(os.environ.get("RUN_NUM", "-1")) > 0:
    requests.post(
        os.environ["JID_UPDATE_COUNTERS"],
        json=[{"key": "<b>BeamlineSummary Plots </b>", "value": "Posted"}],
    )

if args.postStats:
    print("posting to the run tables - ipm values.")
    runtable_data = makeRunTableData(dat)
    postRunTable(runtable_data)
