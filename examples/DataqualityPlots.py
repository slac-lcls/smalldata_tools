import tables
import panel as pn
import h5py
import numpy as np
import holoviews as hv
from holoviews import dim
hv.extension('bokeh')
pn.extension()
import sys
sys.path.append('/reg/d/psdm/xcs/xcslu7818/results/smalldata_tools')
from smalldata_tools.SmallDataAna_psana import SmallDataAna_psana as sdaps
from smalldata_tools.utilities import rebin

## Loading the data (single run)

save_elog = True
detImgMaxSize = 500 #max dimension of image.
expname = 'xcslu7818'
if len(sys.argv)<2:
    print('We need a run#!')
    sys.exit()
run = int(sys.argv[1])
anaps = sdaps(expname,run)
ana = anaps.sda
#ana = sda(expname,run)

## Defining initial selection (laser-on events)
iniFilter='initial'
ana.addCut('lightStatus/xray',0.5,1.5,iniFilter)
ana.addCut('lightStatus/laser',0.5,1.5,iniFilter)

### Get data & define axis title&ranges.

ipmUpDim = hv.Dimension(('ipm4/sum','ipm4 Sum'))
ipmDownDim = hv.Dimension(('ipm5/sum','ipm5 Sum'))
scatterDim = hv.Dimension(('epix10k2M/ROI_0_sum','epix10k2M intensity'))
eventTimeDim = hv.Dimension(('eventTimeR','relative event time'))
l3eDim = hv.Dimension(('l3e','L3 Energy'))

scanVar = ana.getScanName()
try:
    scanDim = hv.Dimension(('scan/%s'%scanVar,'%s'%scanVar))
except:
    scanDim = None
nevtsDim = hv.Dimension(('nevents','N events / scan point'))
nevtsLxtDim = hv.Dimension(('neventslxt','N events / lxt'))

#timing vars.
lxtDim = hv.Dimension(('epics/lxt','lxt'))

ipmUpVar = ana.getVar(ipmUpDim.name,useFilter=iniFilter)
ipmDownVar = ana.getVar(ipmDownDim.name,useFilter=iniFilter)
stepVar = ana.getVar('scan/varStep',useFilter=iniFilter)
l3eVar = ana.getVar('ebeam/L3_energy',useFilter=iniFilter)
eventTimeRaw = ana.getVar('event_time',useFilter=iniFilter)
eventTime = (eventTimeRaw>>32).astype(float)+((eventTimeRaw<<32)>>32).astype(float)*1e-9
eventTimeR = eventTime-eventTime[0]

eventTimeRMed = [np.nanmedian(eventTimeR[i*120:i*120+120]) for i in range(int(eventTimeR.shape[0]/120))]
ipmUpMed =  [np.nanmedian(ipmUpVar[i*120:i*120+120]) for i in range(int(eventTimeR.shape[0]/120))]
ipmDownMed =  [np.nanmedian(ipmDownVar[i*120:i*120+120]) for i in range(int(eventTimeR.shape[0]/120))]

try:
    azav = ana.getVar('epix10k2M/azav_azav',useFilter=iniFilter)
    azav_sum = np.nanmean(azav, axis=0)
    azav_peak = np.argmax(azav_sum)
    scatterVar = np.nanmean(azav[:,max(0,azav_peak-50):min(azav.shape[1],azav_peak+50)], axis=1)
    if len(scatterVar.shape)>1:
        scatterVar = np.nanmean(scatterVar,axis=1)
except:
    scatterVar = None

### Scan Variable

isStepScan = np.nanmax(stepVar)>0
scanVarBins = np.bincount(stepVar,weights=scatterVar)
scanNsteps = np.bincount(stepVar)

### Fast delay stage 

lxt_fast_his = None
try:
    lxt_fast = ana.getVar('enc/lasDelay',useFilter=iniFilter)
    print(np.nanstd(lxt_fast))
    if lxt_fast is not None and np.nanstd(lxt_fast)<1e-4:
        lxt_fast_his = np.histogram(lxt_fast, np.linspace(np.nanpercentile(lxt_fast,1), np.nanpercentile(lxt_fast,99),100))
except:
    pass

#droppled sthots.
ana.addCut('lightStatus/xray',-0.5,0.5,'off')
ana.addCut('evr/code_137',-0.5,0.5,'hxroff')
if ana.getFilter('hxroff').sum() >  ana.getFilter('off').sum():
    offFilter = 'hxroff'
else:
    offFilter = 'off'    
nOff = ana.getFilter(offFilter).sum()

#plots.
ipmUpTime = hv.HexTiles((eventTimeR[ipmUpVar<np.nanpercentile(ipmUpVar,99)],
                         ipmUpVar[ipmUpVar<np.nanpercentile(ipmUpVar,99)]),
                        kdims=[eventTimeDim, ipmUpDim]).\
                        opts(cmap='Blues')
ipmUpTimeMed = hv.Points((eventTimeRMed, ipmUpMed), kdims=[eventTimeDim,ipmUpDim],label=ipmUpDim.label).\
    options(color='r')
ipmDownTimeMed = hv.Points((eventTimeRMed, ipmDownMed), kdims=[eventTimeDim,ipmUpDim],label=ipmDownDim.label).\
    options(color='m')
    
ipmTimeLayout = ipmUpTime*ipmUpTimeMed*ipmDownTimeMed

treeSel = (l3eVar>np.nanpercentile(l3eVar,1))
treePlot = hv.HexTiles((l3eVar[treeSel],ipmUpVar[treeSel]),kdims=[l3eDim,ipmUpDim])

ipmPlot = hv.HexTiles((ipmUpVar, ipmDownVar), kdims=[ipmUpDim, ipmDownDim])
ipmLayout = ipmPlot.hist(dimension=[ipmUpDim.name,ipmDownDim.name])
if scatterVar is not None:
    ipmscatterPlot = hv.HexTiles((scatterVar, ipmDownVar), kdims=[scatterDim, ipmDownDim])

stepPlot = None
if isStepScan:
    try:
        stepPlot = hv.Points((scanVarBins/scanNsteps, scanNsteps), kdims=[scanDim,nevtsDim])
    except:
        print('Failed to make stepPlot', np.nanmax(stepVar), scanNsteps.shape, scanVarBins.shape)
        print('DataPoints: ',scanVarBins/scanNsteps)
        print('Dimensions:', scanDim, nevtsDim)
if lxt_fast_his is not None:
    lxtPlot = hv.Points( (0.5*(lxt_fast_his[1][:-1]+lxt_fast_his[1][1:]), lxt_fast_his[0]), \
                             kdims=[lxtDim,nevtsLxtDim])
else:
    lxtPlot = None

gspec = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name='Data Quality - Run %d'%run)
#gspec[0,0:8] = pn.Row('# Data Quality Plot - Run %04d'%run)
gspec[0:2,0:8] = pn.Column(ipmTimeLayout)
gspec[2:5,0:4] = pn.Column(ipmLayout)
gspec[2:5,4:8] = pn.Column(treePlot)

gspecS = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name='Scan&Scatter')
#gspec[0,0:8] = pn.Row('# Data Quality Plot - Run %04d'%run)
if scatterVar is not None:
    gspecS[0:4,0:8] = pn.Column(ipmscatterPlot)
maxRow=4
if stepPlot is not None:
    gspecS[maxRow,0:8] = pn.Row('## Scan Variable')
    gspecS[maxRow+1:maxRow+3,0:8] = pn.Column(stepPlot)
    maxRow=7
if lxtPlot is not None:
    gspecS[maxRow,0:8] = pn.Row('## Laser - xray Timing')
    gspecS[maxRow+1:maxRow+3,0:8] = pn.Column(lxtPlot)

# Detector stuff. 

detImgs=[]
detGrids=[]
for detImgName in ana.Keys('Sums'):
    image = ana.fh5.get_node('/%s'%detImgName)
    if max(image.shape[0], image.shape[1])>detImgMaxSize:
        rebinFactor = float(detImgMaxSize)/max(image.shape[0],image.shape[1])
        imageR = rebin(image, [int(image.shape[0]*rebinFactor), int(image.shape[1]*rebinFactor)])/ana.getVar('fiducials').shape[0]
    else:
        imageR = image/ana.getVar('fiducials').shape[0]
    #imgArrays.append(imageR/ana.getVar('fiducials').shape[0])
    imgDim = hv.Dimension(('image',detImgName.replace('Sums/','').replace('_calib_img',' Mean Image')),
                                    range=(np.nanpercentile(imageR,1), np.nanpercentile(imageR,99.)))
    detImgs.append(hv.Image(imageR, vdims=[imgDim], name=imgDim.label).options(colorbar=True, cmap='rainbow'))
        
    detGrid = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name=detImgName.replace('Sums/',''))
    detGrid[0,0] = pn.Row(detImgs[-1])
    detGrids.append(detGrid)

if nOff>100:
    for detImgName in ana.Keys('Sums'):
        detName = detImgName.replace('_calib_img','').replace('Sums/','')
        try:
            anaps.AvImage(detName, useFilter=offFilter, numEvts=min(1000, nOff))
        except:
            print('failed to get off shot data for detector %s'%detName)
            continue
        avData = anaps.getAvImage(detName)[1]
        try:
            image = anaps.__dict__[detName].det.image(run, avData)
        except:
            print('failed to make image for detector %s'%detName)
            continue
        if max(image.shape[0], image.shape[1])>detImgMaxSize:
            rebinFactor = float(detImgMaxSize)/max(image.shape[0],image.shape[1])
            imageR = rebin(image, [int(image.shape[0]*rebinFactor), int(image.shape[1]*rebinFactor)])
        else:
            imageR = image
        imgOffDim = hv.Dimension(('image_off',detImgName.replace('Sums/','').replace('_calib_img',' Mean Image Off')),
                                    range=(np.nanpercentile(imageR,1), np.nanpercentile(imageR,99.)))
        detImgs.append(hv.Image(imageR, vdims=[imgOffDim], name=imgOffDim.label).options(colorbar=True, cmap='rainbow'))
        
        detGrid = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name='%s, dropped shots'%detName)
        detGrid[0,0] = pn.Row(detImgs[-1])
        detGrids.append(detGrid)       

tabs = pn.Tabs(gspec)
tabs.append(gspecS)
for detGrid in detGrids:
    tabs.append(detGrid)

elogDir = '/reg/d/psdm/%s/%s/stats/summary/DataQuality_Run%03d'%\
                (expname[0:3],expname,run)
if save_elog:
    import os
    if not os.path.isdir(elogDir):
        os.makedirs(elogDir)
    print('Made Directory to save data:', elogDir)
    #gspec.save(('%s/report.html'%elogDir)) 
    tabs.save(('%s/report.html'%elogDir))
