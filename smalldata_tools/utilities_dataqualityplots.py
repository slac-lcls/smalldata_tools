#!/usr/bin/env python
########################################################################
## load tools and packages required for data handling and plotting
########################################################################
import panel as pn
import h5py
import os
import requests
import numpy as np
from requests.auth import HTTPBasicAuth
import holoviews as hv
import sys
try:
    basestring
except NameError:
    basestring = str
fpath=os.path.dirname(os.path.abspath(__file__))
fpathup = '/'.join(fpath.split('/')[:-1])
try:
    fpath = os.environ.get('MYDIR', fpathup).replace('/arp_scripts','')
except:
    fpath = fpathup
sys.path.append(fpath)
from smalldata_tools.SmallDataAna_psana import SmallDataAna_psana as sdaps
from smalldata_tools.utilities import image_from_dxy
from smalldata_tools.utilities import rebin

#sys.path.append('/reg/data/drpsrcf/mec/meclu9418/scratch/smalldata_tools/')
from smalldata_tools.utilities_plotting import hv_image
import smalldata_tools.SmallDataAna as sda
from smalldata_tools.DropletAna import droplets

## function that chops the 64 bit time integer into something a bit more useful
def evtt2Rt(event_time):
    evtt0 = event_time>>32
    evtt1 = (event_time<<32)>>32
    evtt_sec = evtt0.astype(float)
    evtt_ns = evtt1.astype(float)*1e-9
    Rt = evtt_sec + evtt_ns
    Rt = Rt-Rt[0]
    return Rt

def postRunTable(runtable_data, experiment, run):
    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(experiment)
    print('URL:',ws_url)
    user=experiment[:3]+'opr'
    with open('/cds/home/opr/%s/forElogPost.txt'%user,'r') as reader:
        answer = reader.readline()
    r = requests.post(ws_url, params={"run_num": run}, json=runtable_data, \
                      auth=HTTPBasicAuth(experiment[:3]+'opr', answer[:-1]))
    #we might need to use this for non=current expetiments. Currently does not work in ARP
    #krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    #r = requests.post(ws_url, headers=krbheaders, params={"run_num": args.run}, json=runtable_data)
    print(r)

def makeRunTableData(ana, scanName, Filter = None, varNames=[]):
    n162 = ana.getVar('evr/code_162').sum()
    ana.addCut('evr/code_162',-0.5,0.5,'xon')
    ana.addCut('evr/code_137',0.5,1.5,'xon')
    nOff = ana.getFilter('xon').shape[0]-ana.getFilter('xon').sum()
    #data to be posted to the run table if so requested.
    runtable_data = {"N dropped Shots":int(nOff),
                     "N BYKIK 162":int(n162)}
    if scanName != '':
        runtable_data['scanName'] = scanName

    for var in varNames:
        if Filter is not None:
            thisVar = ana.getVar(var,useFilter=Filter)
        else:
            thisVar = ana.getVar(var)
            thisVarP = np.nanpercentile(thisVar,[25,50,75,99.9])
            runtable_data["%s_1qt"%(var.replace('/','__'))]=thisVarP[0]
            runtable_data["%s_med"%(var.replace('/','__'))]=thisVarP[1]
            runtable_data["%s_3qt"%(var.replace('/','__'))]=thisVarP[2]
            runtable_data["%s_99p9"%(var.replace('/','__'))]=thisVarP[3]
            runtable_data["%s_max"%(var.replace('/','__'))]=np.nanmax(thisVar)
    
    return runtable_data


def ndropletPlot(sda):
    if len(sda.Keys('nDroplets'))==0:
        return
    scatterPlots=[]
    for nDropname in sda.Keys('nDroplets'):
        nDrop = sda.getVar(nDropname)
        nDrop995 = int(np.nanpercentile(nDrop, 99.5))
        nDropHis = np.histogram(nDrop, np.linspace(0, max(10,nDrop995),
                                           min(100,int(max(10,nDrop995)))))
        nDropDim = hv.Dimension(('nDrop%s'%nDropname.split('/')[1],'N droplet'))
        #nDropEvt = hv.Dimension(('nADUEvt%s'%nDropname.split('/')[1],'n events / droplet'))
        #scatterPlot = hv.Scatter((nDropHis[1][1:], nDropHis[0]),kdims=[nDropDim, nDropEvt])
        scatterPlot = hv.Histogram(nDropHis, kdims=[nDropDim])
        scatterPlots.append(scatterPlot)
    return scatterPlots

def specdropletPlot(sda):
    specPlots = []
    for spectrumname in sda.Keys('histogram'):
        spectrum = sdaR.getVar(spectrumname).mean(axis=0)
        specDim = hv.Dimension(('spec%s'%spectrumname.split('/')[1],'droplet ADU'))
        nADUEvt = hv.Dimension(('nADUEvt%s'%spectrumname.split('/')[1],'droplet / ADU'))
        nADUEvtLog = hv.Dimension(('nADUEvtLog%s'%spectrumname.split('/')[1],'log(droplet) / ADU'))

        binName = 'UserDataCfg/'+spectrumname.split('/')[0]+'/'+spectrumname.split('/')[1].replace('histogram','bins')
        binName = binName.replace('droplet','droplet_')
        bins = sdaR.getVar(binName)
        specPlot = hv.Scatter((bins[:-1], spectrum[:len(bins)-1]),kdims=[specDim, nADUEvt])
        specLogPlot = hv.Scatter((bins[:-1], np.log(spectrum[:len(bins)-1])),kdims=[specDim, nADUEvtLog])
        specPlots.append(specPlot)  
        specPlots.append(specLogPlot)  
    return specPlots


def dropletPics(sda, detname, dropname):
#sda = sda.SmallDataAna('meclq8515',324)#, dirname='/reg/d/psdm/xcs/xcsl2516/res/smallData/output_scratch')
##sda = sda.SmallDataAna('meclq8515',266)#, dirname='/reg/d/psdm/xcs/xcsl2516/res/smallData/output_scratch')
    #drop = droplets(sda.fh5,detName='epix100a',dropName='droplet', plotWith='None')#holoviews')
    try:
        drop = droplets(sda.fh5,detName=detname,dropName=dropname, plotWith='None')#holoviews')
        drop.fillDropArrays(only_XYADU=True)
        drop.flattenDropArray()
    except:
        return None
    if len(sda.Keys('nDroplets'))==0:
        return None
    ndropplot =  ndropletPlot(sda)
    drop.set_plotWith('holoviews')
    hnpix = drop.plotnpixel()
    specplot = drop.plotSpectrum(aduRange=[60,500], plotLog=False)

    gspec = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name='Droplets')
    gspec[0:3,0:8] = pn.Column((ndropplot[0]+hnpix).options(shared_axes=False))
    gspec[3:5,0:7] = pn.Column(specplot)

    #check if we save enough droplets for this plot to make sense.
    nDropSaved = sda._fields['%s/%s_%s'%(detname,dropname,drop.xName)][0][1]
    nDrops = sda.getVar('%s/%s_nDroplets'%(detname,dropname))
    nDropPercentile = np.nanpercentile(nDrops,[50,95,99,99.9])
    if nDropPercentile[0]>nDropSaved:
        gspec[5:6,0:7] = pn.Column('# Less than 50% of events have all droplets saved - no XY plot ')
        return gspec

    aduXYDim = hv.Dimension(('adu','ADU'), range=(90,1000))
    #hXY = drop.plotXY(aduXYDim.range, npix=1)
    hXY = drop.plotXY(aduXYDim.range)
    xDim = hv.Dimension(('xDrop','x'), range=(0, drop.detSize[0]))
    yDim = hv.Dimension(('yDrop','y'), range=(0, drop.detSize[1]))
    boundsXY=(0,0,drop.detSize[0],drop.detSize[1])
    imkwargsXY = {'kdims':[xDim, yDim]}
    imkwargsXY['bounds']=boundsXY
    imkwargsXY['cmap']='gray_r'
    imkwargsXY['imgHigh']=np.nanmax(hXY)
    hvimg = hv_image(hXY, **imkwargsXY)
 
    if nDropPercentile[1]>nDropSaved:
        gspec[5:6,0:7] = pn.Column('# 5 or more % of events have not all droplets saved ')
        gspec[6:11,0:7] = pn.Column(hvimg)
    else:
        gspec[5:10,0:7] = pn.Column(hvimg)

    return gspec





def detplots_fromSums(ana, anaps=None):
    # Detector stuff. 
    detImgs=[]
    detGrids=[]
    detNames=[]
    for detImgName in ana.Keys('Sums'):
        image = ana.fh5.get_node('/%s'%detImgName).read()
        if len(image.shape)>2:
            if detImgName.find('135')<0:
                detName = detImgName.replace('Sums/','').split('_')[0]
                if detName in detNames:
                    continue
                detNames.append(detName)
                ix = ana.fh5.get_node('/UserDataCfg/%s/ix'%detName).read()
                iy = ana.fh5.get_node('/UserDataCfg/%s/iy'%detName).read()
                image = image_from_dxy(image, ix, iy)
            else:
                #somehow the epix10k135 has the wrong shape....
                image = image[0]
                #image = image.squeeze()
        if max(image.shape[0], image.shape[1])>detImgMaxSize:
            rebinFactor = float(detImgMaxSize)/max(image.shape[0],image.shape[1])
            imageR = rebin(image, [int(image.shape[0]*rebinFactor), int(image.shape[1]*rebinFactor)])/(ana.getVar('fiducials').shape[0])
        else:
            imageR = image/(ana.getVar('fiducials').shape[0])
        #imgArrays.append(imageR/ana.getVar('fiducials').shape[0])
        imgDim = hv.Dimension(('image',detImgName.replace('Sums/','').replace('_calib_img',' Mean Image')),
                                        range=(np.nanpercentile(imageR,1), np.nanpercentile(imageR,99.)))
        detImgs.append(hv.Image(imageR, vdims=[imgDim], name=imgDim.label).options(colorbar=True, cmap='rainbow'))
        
        detGrid = pn.GridSpec(sizing_mode='stretch_both', max_width=700, name=detImgName.replace('Sums/',''))
        detGrid[0,0] = pn.Row(detImgs[-1])
        detGrids.append(detGrid)

    if nOff>100 and anaps is not None:
        for detImgName in ana.Keys('Sums'):
            detName = detImgName.replace('_calib','').replace('_img','').replace('Sums/','')
            try:
                common_mode=0
                if detName.find('epix10k'): common_mode=80
                anaps.AvImage(detName, useFilter=offFilter, numEvts=min(1000, nOff), common_mode=common_mode)
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
    return detGrid
