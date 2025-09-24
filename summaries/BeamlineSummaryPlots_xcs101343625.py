#!/usr/bin/env python
########################################################################
## load tools and packages required for data handling and plotting
########################################################################
import panel as pn
import h5py
import pandas as pd
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

###########***********Yanwen Plots********################

#read in packages
import numpy as np
import h5py as h5
from matplotlib import pyplot as plt
import time
import psana
from tqdm import tqdm
import sys
import tables

from pathlib import Path

try:
    #read in detector
    def get_smd(exp,run):
        h5dir = Path('/sdf/data/lcls/ds/xcs/{}/hdf5/smalldata'.format(expname))
        fn = '{}_Run{:04d}.h5'.format(expname,run)
        fn = h5dir / fn
        f = h5.File(fn, 'r')
        print (fn)
        data = tables.open_file(fn).root
        return data


    data = get_smd(expname, run)

    #check the fit
    res = 0.34 #eV spectrometer resolution

    nspec = len(data.Alvium_spec.spectrum_proj_h)
    x = np.arange(data.Alvium_spec.spectrum_proj_h.shape[1])
    E = x*res

    nplot = 10 # pick 10 random plots to see fitting
    plot_nums = np.random.randint(low = 0, high = nspec-1, size = nplot)


    plot_list = []
    for i in plot_nums:
        # Get the data for this specific spectrum
        intensity_data = data.Alvium_spec.spectrum_proj_h[i, :]
        fit_data = data.Alvium_spec.spectrum_spectrum_fit[i, :]
    
        # Create the title string for this subplot
        fwhm = data.Alvium_spec.spectrum_sigma[i] * 2.35 * res
        r2_spec = data.Alvium_spec.spectrum_R2_h[i]
        r2_beam = data.Alvium_spec.spectrum_R2_v[i]
        title = f'FWHM: {fwhm:.2f} eV, R2 spec: {r2_spec:.2f}, R2 beam: {r2_beam:.2f}'
        
        # Create the Scatter element for the raw data
        # We pass data as a tuple: (x_values, y_values)
        scatter_plot = hv.Scatter((E, intensity_data), kdims='Energy', vdims='Intensity')
        
        # Create the Curve element for the fit data
        curve_plot = hv.Curve((E, fit_data), kdims='Energy', vdims='Fit')
        
        # Overlay the two elements using the '*' operator
        # and apply the title to this specific combined plot
        overlay = (scatter_plot * curve_plot).opts(title=title)
    
        # Add the finished subplot to our list
        plot_list.append(overlay)


        # --- 3. Arrange the Plots into a Layout and Apply Global Styles ---

        # Create a Layout object from our list of plots
        # and arrange it into 2 columns
    final_layout = hv.Layout(plot_list).cols(2)

    # Apply global styling options to the entire layout
    # HoloViews is smart enough to apply Curve styles only to Curves, etc.
    final_checksum_plot = final_layout.opts(
        hv.opts.Curve(color='red', width=450, height=300, 
                      xlabel='Δ E (eV)', ylabel='intensity (a.u.)'),
        hv.opts.Scatter(color='black', size=3),
        hv.opts.Layout(shared_axes=False) # Options for the final layout
    )


    #histogram statistics
    # ===================================================================
    # Part 1: Convert the Grid of Subplots
    # ===================================================================
    plot_list = []
    for i in tqdm(plot_nums):
        # Get the data for this specific spectrum
        intensity_data = data.Alvium_spec.spectrum_proj_h[i, :]
        fit_data = data.Alvium_spec.spectrum_spectrum_fit[i, :]
        
        # Calculate values for the title
        fwhm = data.Alvium_spec.spectrum_sigma[i] * 2.35 * res
        r2_spec = data.Alvium_spec.spectrum_R2_h[i]
        r2_beam = data.Alvium_spec.spectrum_R2_v[i]
        title = f'FWHM: {fwhm:.2f} eV, R2 spec: {r2_spec:.2f}, R2 beam: {r2_beam:.2f}'

        print("fwhm",fwhm)

        # ===================================================================
        # Part 2: Convert the Histogram
        # ===================================================================
        
        # --- Key Difference: Pre-bin the data with NumPy ---
        # plt.hist does this automatically; hv.Histogram expects binned data.
        counts1, edges1 = np.histogram(
            data.Alvium_spec.spectrum_R2_h, 
            bins=100, 
            range=(0, 1)
        )
            
        counts2, edges2 = np.histogram(
            data.Alvium_spec.spectrum_R2_v, 
            bins=100, 
            range=(0, 1)
        )
            
        # Convert the EArray to a NumPy array by slicing it
        bandwidth_numpy_array = data.Alvium_spec.spectrum_sigma[:]*res*2.35
        
        counts3, edges3 = np.histogram(
            bandwidth_numpy_array, 
            bins=50, 
            range=(1, 50)
            )

        r2_h_numpy_array = data.Alvium_spec.spectrum_R2_h[:]
        ipm4_numpy_array = data.ipm4.sum[:]
        
        # Create the HoloViews Histogram object
        # We pass it a tuple of (bin_edges, bin_counts)
        histogram_plot1 = hv.Histogram((edges1, counts1), 
                                       kdims=['R^2 spectrum'], 
                                       vdims=['shots'])
    
        histogram_plot2 = hv.Histogram((edges2, counts2), 
                                       kdims=['R^2 beam'], 
                                       vdims=['shots'])
        
        histogram_plot3 = hv.Histogram((edges3, counts3), 
                                       kdims=['Bandwidth (eV)'], 
                                       vdims=['shots'])

        scatter_plot_ipm = hv.Scatter((r2_h_numpy_array, ipm4_numpy_array), kdims='R^2 spectrum', vdims='i0')

        
        # Apply styling to the histogram
        histogram_plot1.opts(
            width=250, 
            height=300, 
            title="R^2 Spectrum",
            tools=['hover'] # Add hover tool for interactivity
        )

        histogram_plot2.opts(
            width=250, 
            height=300, 
            title="R^2 Beam",
            tools=['hover'] # Add hover tool for interactivity
        )

        histogram_plot3.opts(
            width=250, 
            height=300, 
            title="Distribution of Bandwidth (eV)",
            tools=['hover'] # Add hover tool for interactivity
        )
            
        final_histo_layout = histogram_plot1 + histogram_plot2 + histogram_plot3 + scatter_plot_ipm
        
        
        # Display the final plot
        #Create a tab in the summaries to save the plots in
        gspecCS = pn.GridSpec(sizing_mode="stretch_both", max_width=700, name="CheckSpectrum")
        gspecCS[0:2, 0:8] = pn.Column(final_checksum_plot)
        gspecCS                                                                                                                                    
        print("created check spectrum tab on run summaries")
        

        gspecHisto = pn.GridSpec(sizing_mode="stretch_both", max_width=700, name="CheckSpectrumHisto")
        gspecHisto[0:2, 0:8] = pn.Column(final_histo_layout)
        gspecHisto
        print("created check spectrum histogram tab on run summaries")

except:
    print("No alvium spec in this run")

##################################end Yanwen plots##################
    
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
try:
    tabs.append(gspecCS)
    tabs.append(gspecHisto)
except:
    print("No tab for checksum alvium spec for this run")
    
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
