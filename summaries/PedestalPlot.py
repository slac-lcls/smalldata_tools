#!/usr/bin/env python
import psana
import numpy as np
import holoviews as hv

hv.extension("bokeh")
import panel as pn

pn.extension()
#import os
#import sys
import logging
from pathlib import Path
import mimetypes

from smalldata_tools.utilities import postElogMsg
from smalldata_tools.utilities import rebin

####
# detector info to be used for status#s to elog
####
gainSwitching = ["jungfrau", "epix10k"]

statusDict = {}
# https://confluence.slac.stanford.edu/display/PSDMInternal/Pixel+status+in+data
statusDict["pnccd"] = {
    "rms_high": 0x1,
    "adu_high_frac": 0x2,
    "adu_low_frac": 0x4,
    "rms_low": 0x8,
    "adu_av_high": 0x10,
    "adu_av_low": 0x20,
    "adu_max_high": 0x40,
    "adu_min_low": 0x80,
}

# https://confluence.slac.stanford.edu/display/PSDM/Jungfrau+bad+gain+mode+switch+status#Jungfraubadgainmodeswitchstatus-pixelstatusstatistics
statusDict["jungfrau"] = {
    "rms_high": 0x1,
    "rms_low": 0x2,
    "adu_high_frac": 0x4,
    "adu_low_frac": 0x8,
    "adu_av_high": 0x10,
    "adu_av_low": 0x20,
    "bad_switch": 0x40,
}

# from logfile for epix10k2M using ana-4.0.48
statusDict["epix10k"] = {
    "rms_high": 0x1,
    "rms_low": 0x2,
    "adu_high_frac": 0x4,
    "adu_low_frac": 0x8,
    "adu_av_high": 0x10,
    "adu_av_low": 0x20,
}

# epix100 for epix100 using ana-4.0.48 for xpptut15. run 260
statusDict["epix100"] = {
    "rms_high": 0x1,
    "rms_low": 0x2,
    "adu_high_frac": 0x4,
    "adu_low_frac": 0x8,
    "adu_av_high": 0x10,
    "adu_av_low": 0x20,
}

# cspad from cs140 using ana-4.0.48 for xpptut15. run 201
statusDict["cspad"] = {
    "rms_high": 0x1,
    "rms_low": 0x2,
    "adu_high_frac": 0x4,
    "adu_low_frac": 0x8,
    "adu_av_high": 0x10,
    "adu_av_low": 0x20,
}

# (other detectors)
statusDict["camera"] = {
    "rms_high": 0x1,
    "rms_low": 0x2,
    "adu_high_frac": 0x4,
    "adu_low_frac": 0x8,
    "adu_av_high": 0x10,
    "adu_av_low": 0x20,
}


####
# helper functions
####
def ped_rms_histograms(nCycles, peds, noise, diff, alias=""):
    min5Ped = 1e6
    max95Ped = -1e6
    max95Noise = -1e6

    for i in range(nCycles):
        if nCycles > 1:
            thisPed = peds[i]
            thisNoise = noise[i]
        else:
            thisPed = peds
            thisNoise = noise

        if np.nanpercentile(thisNoise, 95) > max95Noise:
            max95Noise = np.nanpercentile(thisNoise, 95)
        if np.nanpercentile(thisPed, 95) > max95Ped:
            max95Ped = np.nanpercentile(thisPed, 95)
        if np.nanpercentile(thisPed, 5) < min5Ped:
            min5Ped = np.nanpercentile(thisPed, 5)

    # for the UXI, most pixels read back as 0. Widen the range to look for decent min/max values.
    if max95Ped == min5Ped:
        for i in range(nCycles):
            if nCycles > 1:
                thisPed = peds[i]
                thisNoise = noise[i]
            else:
                thisPed = peds
                thisNoise = noise
            if np.nanpercentile(thisNoise, 99.9) > max95Noise:
                max95Noise = np.nanpercentile(thisNoise, 99.9)
            if np.nanpercentile(thisPed, 99.9) > max95Ped:
                max95Ped = np.nanpercentile(thisPed, 99.9)
            if np.nanpercentile(thisPed, 0.1) < min5Ped:
                min5Ped = np.nanpercentile(thisPed, 0.1)

    max95Noise *= 1.25
    max95Ped *= 1.25
    min5Ped *= 0.9

    pedBins = np.linspace(min5Ped, max95Ped, 200)
    pedHistograms = []
    pedBinDim = hv.Dimension(
        ("ped_bin%s" % alias, "pedestal in ADU"), range=(min5Ped, max95Ped)
    )
    noiseBins = np.linspace(0, max95Noise, 200)
    noiseHistograms = []
    noiseBinDim = hv.Dimension(
        ("rms_bin%s" % alias, "noise in ADU"), range=(0, max95Noise)
    )

    noiseMax = 0
    pedMax = 0
    for i in range(max(1, nCycles)):
        if nCycles > 1:
            thisNoise = noise[i]
            thisPed = peds[i]
        else:
            thisPed = peds
            thisNoise = noise

        pedHistograms.append(np.histogram(thisPed.flatten(), pedBins))
        noiseHistograms.append(np.histogram(thisNoise.flatten(), noiseBins))
        if pedMax < np.nanmax(pedHistograms[-1][0]):
            pedMax = np.nanmax(pedHistograms[-1][0])
        if noiseMax < np.nanmax(noiseHistograms[-1][0]):
            noiseMax = np.nanmax(noiseHistograms[-1][0])
    pedMax *= 1.1
    noiseMax *= 1.1
    evtsPDim = hv.Dimension(
        ("evtsP%s" % alias, "N pixels / pedestal"), range=(0, pedMax)
    )
    evtsDim = hv.Dimension(("evtsN%s" % alias, "N pixels / noise"), range=(0, noiseMax))

    pedHists = []
    noiseHists = []
    i = 0
    for pedH, noiseH in zip(pedHistograms, noiseHistograms):
        pedHists.append(
            hv.Points(
                (0.5 * (pedBins[1:] + pedBins[:-1]), pedH[0]),
                label="Cycle %d" % i,
                kdims=[pedBinDim, evtsPDim],
            )
            * hv.Curve((0.5 * (pedBins[1:] + pedBins[:-1]), pedH[0]))
        )
        noiseHists.append(
            hv.Points(
                (0.5 * (noiseBins[1:] + noiseBins[:-1]), noiseH[0]),
                label="Cycle %d" % i,
                kdims=[noiseBinDim, evtsDim],
            )
            * hv.Curve((0.5 * (noiseBins[1:] + noiseBins[:-1]), noiseH[0]))
        )
        i += 1

    if diff is None:
        return pedHists, noiseHists, None

    min5Diff = 1e6
    max95Diff = -1e6

    for i in range(nCycles):
        if nCycles > 1:
            thisDiff = diff[i]
        else:
            thisDiff = diff

        if np.nanpercentile(thisDiff, 98) > max95Diff:
            max95Diff = np.nanpercentile(thisDiff, 98)
        if np.nanpercentile(thisDiff, 2) < min5Diff:
            min5Diff = np.nanpercentile(thisDiff, 2)

    if max95Diff >= 0:
        max95Diff *= 1.25
    else:
        max95Diff *= 0.8
    if min5Diff >= 0:
        min5Diff *= 0.8
    else:
        min5Diff *= 1.25

    diffBins = np.linspace(min5Diff, max95Diff, 200)
    diffHistograms = []
    diffBinDim = hv.Dimension(
        ("diff_bin%s" % alias, "diff in ADU"), range=(min5Diff, max95Diff)
    )

    diffMax = 0
    for i in range(max(1, nCycles)):
        if nCycles > 1:
            thisDiff = diff[i]
        else:
            thisDiff = diff

        diffHistograms.append(np.histogram(thisDiff.flatten(), diffBins))
        if diffMax < np.nanmax(diffHistograms[-1][0]):
            diffMax = np.nanmax(diffHistograms[-1][0])
    diffMax *= 1.1
    evtsDDim = hv.Dimension(("evtsD%s" % alias, "N pixels / diff"), range=(0, diffMax))

    diffHists = []
    i = 0
    for diffH in diffHistograms:
        diffHists.append(
            hv.Points(
                (0.5 * (diffBins[1:] + diffBins[:-1]), diffH[0]),
                label="Cycle %d" % i,
                kdims=[diffBinDim, evtsDDim],
            )
            * hv.Curve((0.5 * (diffBins[1:] + diffBins[:-1]), diffH[0]))
        )
        i += 1

    return pedHists, noiseHists, diffHists


def plotPedImgs(
    nCycles,
    run,
    peds,
    noise,
    pedImgs,
    noiseImgs,
    diffPeds=None,
    diffPedsImgs=None,
    detImgMaxSize=500,
    plotInfo=None,
):
    hpedImgs = []
    hrmsImgs = []
    hdiffImgs = []

    for i in range(nCycles):
        pedImg = pedImgs[i]
        rmsImg = noiseImgs[i]
        if nCycles > 1:
            thisPed = peds[i]
            thisNoise = noise[i]
            tpedDim = hv.Dimension(
                ("ped_%d" % i, "pedestal in ADU, cycle %d" % i),
                range=(np.nanpercentile(thisPed, 0.1), np.nanpercentile(thisPed, 99.9)),
            )
            trmsDim = hv.Dimension(
                ("rms_%d" % i, "noise in ADU, cycle %d" % i),
                range=(
                    np.nanpercentile(thisNoise, 0.1),
                    np.nanpercentile(thisNoise, 99.9),
                ),
            )
            if diffPeds is not None:
                try:
                    thisDiff = diffPeds[i]
                    diffImg = diffPedsImgs[i]
                    tdiffDim = hv.Dimension(
                        ("ped_%d" % i, "delta pedestal in ADU, cycle %d" % i),
                        range=(
                            np.nanpercentile(thisDiff, 0.1),
                            np.nanpercentile(thisDiff, 99.9),
                        ),
                    )
                except:
                    pass
        else:
            thisPed = peds
            thisNoise = noise
            tpedDim = hv.Dimension(
                ("ped", "pedestal in ADU"),
                range=(np.nanpercentile(thisPed, 0.1), np.nanpercentile(thisPed, 99.9)),
            )
            trmsDim = hv.Dimension(
                ("rms", "noise in ADU"),
                range=(
                    np.nanpercentile(thisNoise, 0.1),
                    np.nanpercentile(thisNoise, 99.9),
                ),
            )
            if diffPeds is not None:
                try:
                    thisDiff = diffPeds
                    diffImg = diffPedsImgs
                    tdiffDim = hv.Dimension(
                        ("ped_%d" % i, "delta pedestal in ADU"),
                        range=(
                            np.nanpercentile(thisDiff, 0.1),
                            np.nanpercentile(thisDiff, 99.9),
                        ),
                    )
                except:
                    pass

        if pedImg is None:
            pedImg = thisPed
            rmsImg = thisNoise
            if peds_pre is not None:
                diffImg = thisDiff
        if max(pedImg.shape[0], pedImg.shape[1]) > detImgMaxSize:
            rebinFactor = float(detImgMaxSize) / max(pedImg.shape[0], pedImg.shape[1])
            pedImg = rebin(
                pedImg,
                [
                    int(pedImg.shape[0] * rebinFactor),
                    int(pedImg.shape[1] * rebinFactor),
                ],
            )
            rmsImg = rebin(
                rmsImg,
                [
                    int(rmsImg.shape[0] * rebinFactor),
                    int(rmsImg.shape[1] * rebinFactor),
                ],
            )
            if diffPeds is not None:
                diffImg = rebin(
                    diffImg,
                    [
                        int(diffImg.shape[0] * rebinFactor),
                        int(diffImg.shape[1] * rebinFactor),
                    ],
                )

        if plotInfo is not None:
            xrange, yrange, xDim, yDim = plotInfo
        if xrange is not None:
            xrange = [0, pedImg.shape[0]]
        if yrange is not None:
            yrange = [0, pedImg.shape[1]]
        hpedImgs.append(
            hv.Image(
                pedImg,
                bounds=(xrange[0], yrange[0], xrange[1], yrange[1]),
                kdims=[xDim, yDim],
                vdims=[tpedDim],
                label="Pedestal, Cycle %d" % i,
            ).options(colorbar=True, aspect="equal", cmap="rainbow")
        )
        hrmsImgs.append(
            hv.Image(
                rmsImg,
                bounds=(xrange[0], yrange[0], xrange[1], yrange[1]),
                kdims=[xDim, yDim],
                vdims=[trmsDim],
                label="Rms, Cycle %d" % i,
            ).options(colorbar=True, aspect="equal", cmap="rainbow")
        )
        if diffPeds is not None:
            hdiffImgs.append(
                hv.Image(
                    diffImg,
                    bounds=(xrange[0], yrange[0], xrange[1], yrange[1]),
                    kdims=[xDim, yDim],
                    vdims=[tdiffDim],
                    label="Diff current-prevous, Cycle %d" % i,
                ).options(colorbar=True, aspect="equal", cmap="rainbow")
            )

    return hpedImgs, hrmsImgs, hdiffImgs


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

