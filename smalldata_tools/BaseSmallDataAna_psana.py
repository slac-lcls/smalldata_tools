import numpy as np
try:
    from pylab import ginput
except:
    print('could not import ginput')
    pass
try:
    from matplotlib import pyplot as plt
except:
    print('could not import pyplot')
    pass
from matplotlib import gridspec
from matplotlib import path
from itertools import count as icount
import os
import socket
import holoviews as hv
import bokeh.plotting as bp
from matplotlib import gridspec
import logging
import time
import requests
try:
    basestring
except NameError:
    basestring = str
try:
    raw_input
except NameError:
    raw_input = input


import smalldata_tools.SmallDataAna as sda

from smalldata_tools.DetObject import DetObject, DetObjectClass
from smalldata_tools.SmallDataUtils import getUserData
from smalldata_tools.utilities import printR
from smalldata_tools.utilities import addToHdf5
from smalldata_tools.utilities import rename_reduceRandomVar
from smalldata_tools.utilities_plotting import plotImageBokeh
from smalldata_tools.utilities_plotting import hv_image
from smalldata_tools.utilities_plotting import hv_image_ctl
from smalldata_tools.utilities_plotting import hv_3dimage
from smalldata_tools.utilities_FitCenter import FindFitCenter
from smalldata_tools.utilities_FitCenter import fitCircle
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, spectrumFunc, projectionFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning
from smalldata_tools.ana_funcs.sparsifyFunc import sparsifyFunc

import smalldata_tools.cube.cube_mpi_fun as mpi_fun

from mpi4py import MPI
import h5py
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class BaseSmallDataAna_psana(object):
    """ Base class for SmallDataAna_psana. Encapsulate all the analysis functions common to lcls-I
    and lcls-II.
    psana specific things (data loading, etc) must be implemented in the respective lcls-I and 
    lcls-II sub-classes. See SmallDataAna_psana.py and SmallDataAna_psana_lcls2.py
    """
    
    def __init__(self, expname='', run=-1, dirname='', filename='', plotWith='matplotlib'):
        self.run = int(run)
        self.expname = expname
        self.hutch = expname[0:3]
        self.plotWith = plotWith
        self.sda = None

        ws_url = "https://pswww.slac.stanford.edu/ws/lgbk"
        if self.hutch == 'cxi':
            print('Will assume the first CXI station, if this is wrong, make sure to add  -e <expname> on commandline')
        resp = requests.get(ws_url + "/lgbk/ws/activeexperiment_for_instrument_station", {"instrument_name": self.hutch.upper(), 
                                                                                          "station": 0})
        try:
            currExpname = resp.json().get("value", {}).get("name")
            print('Current experiment for %s is %s'%(self.hutch, currExpname))
            rundoc = requests.get(ws_url + "/lgbk/" + expname  + "/ws/current_run").json()["value"]
            if not rundoc:
                logger.error("Invalid response from server")
            lastRun = int(rundoc['num'])
            #the tutorial data has missing run numbers, breaking this check
            if self.run > lastRun and self.expname.find('tut')<0:
                printR(rank, 'experiment %s does only have %d runs, requested %d'%(expname, lastRun, self.run))
                return None
        except:
            lastRun = -1
        
        self.isLive = False
        if self.run == lastRun:
            end_time = rundoc.get('end_time', None)
            if end_time is None:
                self.isLive = True
        
        printR(rank, '\nTry to make SmallDataAna using dirname %s, for exp %s and run %s.'%(dirname,expname,run))
        try:
            printR(rank, 'Setting up SmallData ana from anaps.')
            self.sda = sda.SmallDataAna(
                expname,
                self.run,
                dirname=dirname,
                filename=filename,
                plotWith=plotWith)
        except:
            printR(rank, 'Failed, set anaps.sda to None.')
            self.sda = None
        self.calibhisto={}
        self.jobsIds = []
        self.commonModeStrings=['raw','pedSub','unb','hist','histAmiLike','median','medianNoNorm',
                                'medianSN','median45','cm47','cm71','cm72','cm10','cm145','cm146',
                                'cm147','cm110','calib','cm80','cm81']
        return

    def commonModeStr(self, common_mode=0):
        if common_mode<0:
            return 'raw_'
        elif common_mode==5:
            return 'unb_'
        elif common_mode==4 or common_mode==1:
            return 'hist_'
        elif common_mode==6:
            return 'median_'
        elif common_mode==30:
            return 'calib_'
        elif common_mode==34:
            return 'histAmiLike_'
        elif common_mode==36:
            return 'medianNoNorm_'
        elif common_mode==45:
            return 'median45_'
        elif common_mode==46:
            return 'medianSN_'
        elif common_mode==47:
            return 'cm47_'
        elif common_mode==71:
            return 'cm71_'
        elif common_mode==72:
            return 'cm72_'
        elif common_mode==80:
            return 'cm80_'
        elif common_mode==81:
            return 'cm81_'
        elif common_mode==10:
            return 'cm10_'
        elif common_mode==105:
            return 'cm105_'
        elif common_mode==145:
            return 'cm145_'
        elif common_mode==146:
            return 'cm146_'
        elif common_mode==147:
            return 'cm147_'
        elif common_mode==110:
            return 'cm110_'
        elif common_mode==0:
            return 'pedSub_'
        else:
            return ''


    def plotVar(self, plotvar, numBins=[100], useFilter=False, limits=[1,99],fig=None,asHist=False):
        self.sda.plotVar(plotvar=plotvar, numBins=numBins, useFilter=useFilter, limits=limits,fig=fig,asHist=asHist)
        return
    

    def plotScan(self, ttCorr=False, sig='diodeU/channels', sigROI=[], i0='ipm3/sum', numBins=100):
        self.sda.plotScan(ttCorr=ttCorr, sig=sig, sigROI=sigROI, i0=i0, numBins=numBins)
        return


    def _getDetName(self, detname=None):
        detNameList = ['cs','Cs','epix','Epix','opal','Opal','zyla','Zyla',
                       'jungfrau','Jungfrau','gige','Camera','icarus','rayonix']
        #look for detector
        aliases=[]
        for key in self.Keys():
            if key.alias()!='':
                kname = key.alias()
            else:
                kname = key.src().__repr__()
            for detName in detNameList:
                if kname.find(detName)>=0 and kname.find('Epics')<0:
                    aliases.append(kname)
        if len(aliases)==1:
            return aliases[0]
        elif detname is not None:
            alias_match = [alias for alias in aliases if alias.find(detname)>=0]
            return alias_match
        else:
            return aliases
    

    def getAvImage(self,detname=None, imgName=None):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('AvImg')!=0:
                continue
            if key.find('_mask_')>=0:
                continue
            if key.find('_azint_')>=0:
                continue
            if imgName is not None and key.find(imgName)<0:
                continue
            if detname is not None and key.find(detname)<0:
                continue
            avImages.append(key)
        if len(avImages)==0:
            print('please create the AvImage first!')
            return
        elif len(avImages)>1:
            print('we have the following options: ',avImages)
            avImage=raw_input('type the name of the AvImage to use:')
        else:
            avImage=avImages[0]
        detname = self._getDetName_from_AvImage(avImage)
        img = self.__dict__[avImage]
        return detname, img, avImage

    
    def getAzInt(self,detname=None):
        azInts=[]
        for key in self.__dict__.keys():
            if key.find('_mask_')>=0:
                continue
            print('key ',key)
            if key.find('azint')>=0:
                if detname is not None and key.find(detname)>=0:
                    azInts.append(key)
                elif detname is None:
                    azInts.append(key)
        if len(azInts)==0:
            print('please create an azimuthal integral first!')
            return
        elif len(azInts)>1:
            print('we have the following options: ',azInts)
            avInt=raw_input('type the name of the AvImage to use:')
        else:
            avInt=azInts[0]
        detname = self._getDetName_from_AvImage(avInt.replace('_azint_',''))
        values = self.__dict__[avInt]
        #azavName = 'azav'
        #bins = self.__dict__[detname].__dict__[azavName+'_q']
        return detname, values, avInt

    
    def _getDetName_from_AvImage(self,avimage):
        detname=''
        avimage=avimage.replace('std_','')
        avimage=avimage.replace('AvImg_','')
        for thisCmString in self.commonModeStrings:
            avimage=avimage.replace(thisCmString+'_','')
        dns = avimage.split('_')
        for ddns in dns:
            if ddns.find('thres')<0 and ddns.find('Filter')<0:
                detname+=ddns;detname+='_'
        if detname[-1]=='_':
            detname = detname[:-1]
        return detname
    
    
    def plotAvImage(self, detname=None, 
                    imgName=None, 
                    use_mask=False, 
                    ROI=[], limits=[5,99.5], 
                    returnIt=False, 
                    plotWith=None,
                    debugPlot=-1, 
                    inImage={}):
        if inImage!={}:
            img=inImage['image']
            if 'mask' in inImage:
                mask=inImage['mask']
            else:
                mask=np.ones_like(img)
            detname=inImage['detname']
            avImage='input Image'
            if 'name' in inImage:
                avImage=inImage['name']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)
        else:
            detname, img, avImage = self.getAvImage(detname=detname, imgName=imgName)

        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
            img = (img*mask)

        plotMax = np.nanpercentile(img, limits[1])
        plotMin = np.nanpercentile(img, limits[0])
        print(f'plot {avImage} using the {limits[0]} {limits[1]} percentiles as plot min/max: ({plotMin}, {plotMax})')

        if len(img.shape)>2:
            image = self.__dict__[detname].det.image(self.run, img)
            if image is None:
                if self.__dict__[detname].det.dettype != DetObjectClass.Icarus:
                    image = img.squeeze()
                else:
                    imgNew = img[0]
                    for iFrame, frame in enumerate(img):
                        if iFrame>0:
                            imgNew = np.append(imgNew, frame, axis=1)
                    image = imgNew
        else:
            image = img

        #have issues with Jupyter notebook io size.
        if plotWith is None:
            plotWith=self.plotWith
        if plotWith=='bokeh_notebook':
            plot_title="%s in Run %d"%(avImage, self.run)
            if debugPlot>0:
                img = image[:debugPlot, :debugPlot]
            else:
                img = image
            layout, p, im = plotImageBokeh(img, 
                                           plot_title=plot_title, 
                                           plotMaxP=np.nanpercentile(img,99),
                                           plotMinP=np.nanpercentile(img,1), 
                                           plotWidth=700, 
                                           plotHeight=700)
            bp.output_notebook()
            bp.show(layout)
        else:
            fig=plt.figure(figsize=(10,6))
            if ROI!=[]:
                gs=gridspec.GridSpec(1,2,width_ratios=[2,1])        
                im1 = plt.subplot(gs[1]).imshow(img[ROI[0][0],ROI[1][0]:ROI[1][1],ROI[2][0]:ROI[2][1]],
                                                clim=[plotMin,plotMax],
                                                nterpolation='None')
                cbar1 = plt.colorbar(im1)
            else:
                gs=gridspec.GridSpec(1,2,width_ratios=[99,1])        
            im0 = plt.subplot(gs[0]).imshow(image, clim=[plotMin,plotMax], interpolation='None')
            cbar1 = plt.colorbar(im0)
            plt.title(avImage)
            plt.show()

        if returnIt:
            return image
        return
    
    
    def saveAvImage(self,detname=None,dirname=''):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('_mask_')>=0:
                continue
            if key.find('_azint_')>=0:
                continue
            if detname is not None and key.find(detname)<0:
                continue
            if key.find('AvImg')==0:
                if detname is not None and key.find(detname)>=0:
                    avImages.append(key)
                elif detname is None:
                    avImages.append(key)

        detname, img, avImage = self.getAvImage(detname=None)

        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True
        data_to_save={}
        data_to_save['ped'] = self.__dict__[detname].ped
        data_to_save['rms'] = self.__dict__[detname].rms
        data_to_save['mask'] = self.__dict__[detname].mask
        data_to_save['cmask'] = self.__dict__[detname].cmask
        if self.__dict__[detname].x is not None:
            data_to_save['x'] = self.__dict__[detname].x
            data_to_save['y'] = self.__dict__[detname].y
        for rawImg in avImages:
            data_to_save[rawImg] = self.__dict__[rawImg]
            if needsGeo:
                data_to_save[rawImg+'_image'] = self.__dict__[detname].det.image(self.run,self.__dict__[rawImg])
        if dirname=='':
            dirname=self.sda.dirname
        if len(avImages)<=0:
            return
        fname = '%s%s_Run%03d_%s.h5'%(dirname,self.expname,self.run,avImages[0])
        print('now save information to file: ',fname)
        imgFile = h5py.File(fname, "w")
        for key in data_to_save:
            addToHdf5(imgFile, key, data_to_save[key])
        imgFile.close()

        
    #bokeh, use BoxSelectTool to get selected region
    #https://stackoverflow.com/questions/34164587/get-selected-data-contained-within-box-select-tool-in-bokeh
    def SelectRegion(self,detname=None, limits=[5,99.5]):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('AvImg')==0:
                if detname is not None and key.find(detname)>=0:
                    avImages.append(key)
                elif detname is None:
                    avImages.append(key)
        if len(avImages)==0:
            print('please create the AvImage first!')
            return
        elif len(avImages)>1:
            print('we have the following options: ',avImages)
            avImage=raw_input('type the name of the AvImage to use:')
        else:
            avImage=avImages[0]
        detname = self._getDetName_from_AvImage(avImage)
        img = self.__dict__[avImage]

        plotMax = np.nanpercentile(img, limits[1])
        plotMin = np.nanpercentile(img, limits[0])
        print('plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        fig=plt.figure(figsize=(10,6))
        gs=gridspec.GridSpec(1,2,width_ratios=[2,1])

        needsGeo=False
        if self.__dict__[detname].det.dettype==DetObjectClass.Jungfrau:
            if self.__dict__[detname].ped[0].shape != self.__dict__[detname].imgShape:
                needsGeo=True
        elif self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if needsGeo:
            image = self.__dict__[detname].det.image(self.run, img)
        else:
            image = img.squeeze()

        plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')
        plt.pause(0.0001)

        ix = self.__dict__[detname+'_ix']
        iy = self.__dict__[detname+'_iy']

        happy = False
        while not happy:
            p =np.array(ginput(2))
            p_axis1=[int(p[:,1].min()),int(p[:,1].max())+1]
            p_axis0=[int(p[:,0].min()),int(p[:,0].max())+1]
            print('points:',p)
            plt.subplot(gs[1]).imshow(image[p_axis1[0]:p_axis1[1],p_axis0[0]:p_axis0[1]],clim=[plotMin,plotMax],interpolation='None')
            plt.pause(0.0001)
            if needsGeo:
                mask_roi=np.zeros_like(image)
                mask_roi[p_axis1[0]:p_axis1[1],p_axis0[0]:p_axis0[1]]=1
                mask_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(ix,iy)] )
            if raw_input("Happy with this selection:\n") in ["y","Y"]:
                happy = True
                    
        if not needsGeo:
            print('ROI: [[%i,%i], [%i,%i]]'%(p[:,1].min(),p[:,1].max(),p[:,0].min(),p[:,0].max()))
            return

        if len(mask_nda.shape)>2:
            for itile,tile in enumerate(mask_nda):
                if tile.sum()>0:
                    ax0 = np.arange(0,tile.sum(axis=0).shape[0])[tile.sum(axis=0)>0]
                    ax1 = np.arange(0,tile.sum(axis=1).shape[0])[tile.sum(axis=1)>0]
                    print('ROI: [[%i,%i], [%i,%i], [%i,%i]]'%(itile,itile+1,ax1.min(),ax1.max(),ax0.min(),ax0.max()))
                    fig=plt.figure(figsize=(6,6))
                    plt.imshow(img[itile,ax1.min():ax1.max(),ax0.min():ax0.max()],interpolation='none')
                    plt.pause(0.0001)
                    plt.show()
        else:
            tile=mask_nda
            if tile.sum()>0:
                ax0 = np.arange(0,tile.sum(axis=0).shape[0])[tile.sum(axis=0)>0]
                ax1 = np.arange(0,tile.sum(axis=1).shape[0])[tile.sum(axis=1)>0]
                print('ROI: [[%i,%i], [%i,%i]]'%(ax1.min(),ax1.max(),ax0.min(),ax0.max()))
        return
    
    
    def FitCircleAuto(self, detname=None, plotRes=True, forceMask=False, inImage={},inParams={}):
        if inImage!={}:
            img=inImage['image']
            if 'mask' in inImage:
                mask=inImage['mask']
            else:
                mask=np.ones_like(img)
            detname=inImage['detname']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)

        else:
            detname, img, avImage = self.getAvImage(detname=None)
            try:
               mask = self.__dict__[detname].det.cmask.astype(bool)
            except:
               mask = self.__dict__[detname].det.mask(self.run, calib=True, status=True).astype(bool)
            nPixRaw=1.
            for idim in mask.shape:
                nPixRaw=nPixRaw*idim
            print('check calib mask: ',mask.sum(),nPixRaw,(mask.sum()/nPixRaw)>0.5,(mask.sum().astype(float)/nPixRaw))
            if (mask.sum()/nPixRaw)<0.5 and not forceMask:
                mask=~mask
            try:
                maskgeo = self.__dict__[detname].det.mask_geo(self.run).astype(bool)
                mask = mask*maskgeo
            except:
                pass
        #apply mask to image.
        img = (img*mask)

        needsGeo=False
        if self.__dict__[detname].ped is not None and self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if needsGeo:
            image = self.__dict__[detname].det.image(self.run, img)
            mask =  self.__dict__[detname].det.image(self.run, mask)
            #limits =  [[self.__dict__[detname+'_x'].min(),self.__dict__[detname+'_x'].max()],
            #           [self.__dict__[detname+'_y'].min(),self.__dict__[detname+'_y'].max()]]
        else:
            image = img
        limits=[[0,image.shape[0]],[0,image.shape[1]]]

        #plot image in grayscale
        combRes, ringInfo, arSparse = FindFitCenter(image, mask, inParams=inParams)

        if not plotRes:
            return combRes
        #plot image in grayscale
        plt.figure(figsize=[12,12])
        plt.imshow(image, interpolation='none', cmap='gray',clim=[0,np.nanpercentile(img.flatten(),99.5)])
        #plot sparse images in blue.
        plt.plot(arSparse.col, arSparse.row,markersize=5,color='#ff9900',marker='.',linestyle='None')
        if combRes==-1:
            return -1
        greens=['#666600','#669900','#66cc00','#66cc99','#6699cc','#6633ff']
        print('center ',combRes['xCen'],combRes['yCen'])
        for ir,thisRingInfo,rFit in zip(icount(),ringInfo,combRes['R']):
           #plot ransac selected data in green <-- consider different greens for first circles.
            plt.plot(arSparse.col[thisRingInfo['pointsInCircle']],
                     arSparse.row[thisRingInfo['pointsInCircle']],
                     marker='.',
                     color=greens[ir],
                     linestyle='None', 
                     markersize=4)
            circle1 = plt.Circle((combRes['xCen'], combRes['yCen']),
                                 rFit,color=greens[ir],
                                 fill=False,
                                 linestyle='dashed',linewidth=1)
            plt.gca().add_artist(circle1)
        #circle1 = plt.Circle((combRes['xCen'],combRes['yCen']),ringInfo[0]['rInput'],color='r',fill=False,linestyle='dashed')
        #plt.gca().add_artist(circle1)
        plt.plot([combRes['xCen'],combRes['xCen']],[combRes['yCen']-20,combRes['yCen']+20],color='m')
        plt.plot([combRes['xCen']-20,combRes['xCen']+20],[combRes['yCen'],combRes['yCen']],color='m')
        plt.xlim(limits[0])
        plt.ylim(limits[1])
        plt.show()
        if needsGeo:
            geo = self.__dict__[detname].det.geometry(self.run)
            d=160000.
            ix0,iy0 = geo.point_coord_indexes(p_um=(0,0))
            ixd,iyd = geo.point_coord_indexes(p_um=(d,d))
            combRes['yCen'] = d*(combRes['yCen']-ix0)/(ixd-ix0)
            combRes['xCen'] = d*(combRes['xCen']-iy0)/(iyd-iy0)
            print('center Final mid',combRes['xCen'],combRes['yCen'])
            helpVar =  combRes['xCen']
            combRes['xCen'] = combRes['yCen']
            combRes['yCen'] = helpVar
            for r in combRes['R']:
                print('aradii: ',r,d*r/(ixd-ix0))
        else:
            for r in combRes['R']:
                print('aradii: ',r)
        print('center Final ',combRes['xCen'],combRes['yCen'])
        return combRes
    
    
    def FitCircleMouse(self, detname=None, use_mask=False, limits=[5,99.5], use_mask_local=False):
        self.FitCircle(detname=detname, 
                       use_mouse=True, 
                       use_mask=use_mask, 
                       use_mask_local=use_mask_local, 
                       limits=[5,99.5])
        return
    

    def FitCircleThreshold(self, 
                           detname=None, 
                           use_mask=False, 
                           use_mask_local=False, 
                           limits=[5,99.5],
                           plotIt=False, 
                           thresIn=None):
        self.FitCircle(detname=detname, 
                       use_mouse=False, 
                       use_mask=use_mask, 
                       use_mask_local=use_mask_local, 
                       limits=[5,99.5], 
                       thresIn=thresIn, 
                       plotIt=plotIt)
        return
    
    
    def FitCircle(self, **kwargs):
        """
        function to fit (single) circles in images
        options are to fit point selected by mouse or by threshold
        
        Parameters
        ----------
        detname: name of detector we have created an average image for
        use_mouse: select points to be fitted by mouse (if false, use threshold)
        use_mask: use the mask of detector as saved in calibration data, default is False
        use_mask_local: use a locally defined mask, default is False
        limits: range of z-axisof plot in percent of image data, default is [5, 99.5]
        plotIt: show plots of results & for selection of parameters. Default is True. Turn off for running this in batch
        thresIn: supply a minimum threshold for points selected for fit in percentile of image values
        singleTile: select a single tile of image (e.g. useful for Icarus detector), default -1 (use whole detector)
        """
        detname = kwargs.pop("detname",None)
        use_mouse = kwargs.pop("use_mouse",None)
        use_mask = kwargs.pop("use_mask",False)
        use_mask_local = kwargs.pop("use_mask_local",False)
        plotIt = kwargs.pop("plotIt",True)
        thresIn = kwargs.pop("thresIn",None)
        limits = kwargs.pop("limits",[5, 99.5])
        singleTile = kwargs.pop("singleTile",-1)

        try: 
            detname, img, avImage = self.getAvImage(detname=detname)
        except:
            detname, img, avImage = self.getAvImage(detname=None)

        if use_mouse is None or (not isinstance(use_mouse, bool) and not isinstance(use_mouse, basestring)):
            #set pplotting to yes
            plotIt=True
            if raw_input("Select Circle Points by Mouse?:\n") in ["y","Y"]:
                use_mouse = True
            else:
                use_mouse = False

        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
            img = (img*mask)

        if use_mask_local and self.__dict__.has_key('_mask_'+avImage):
                mask = self.__dict__['_mask_'+avImage]
                img = (img*mask)
        elif use_mask_local and not use_mouse:
            maskName = raw_input('no local mask defined for %s, please enter file name '%avImage)
                
            if maskName!='' and os.path.isfile(maskName):
                locmask=np.loadtxt(maskName)
                if self.__dict__[detname].x is not None and locmask.shape != self.__dict__[detname].x.shape:
                    if locmask.shape[1] == 2:
                        locmask = locmask.transpose(1,0)
                    locmask = locmask.reshape(self.__dict__[detname].x.shape)
                self.__dict__['_mask_'+avImage] = locmask
                img = (img*locmask)

        plotMax = np.nanpercentile(img[img!=0], 99.5)
        plotMin = np.nanpercentile(img[img!=0], 5)
        printR(rank, 'plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        needsGeo=False
        if self.__dict__[detname].ped is not None and self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        x = self.__dict__[detname+'_x']
        y = self.__dict__[detname+'_y']
        if x is None:
            x = np.arange(0, image.shape[1])
            y = np.arange(0, image.shape[0])
            x,y = np.meshgrid(x,y)
        extent=[x.min(), x.max(), y.min(), y.max()]

        if needsGeo:
            if self.__dict__[detname].det.dettype==DetObjectClass.Icarus:
                if singleTile<0:
                    print('we need to select a tile for the icarus')
                    return
                elif singleTile>=img.shape[0]:
                    print('requested tile %d, detector %s only has %d tiles'%(singleTile, detname, image.shape[0]))
                    return
                    
                image = img[singleTile]
                x = x[singleTile]
                y = y[singleTile]
                extent=[x.min(), x.max(), y.min(), y.max()]
                print('icarus image: ',img.shape, image.shape, extent)
            else:
                image = self.__dict__[detname].det.image(self.run, img)
        else:
            image = img

        happy = False
        if use_mouse:
            while not happy:
                fig=plt.figure(figsize=(10,10))
                if needsGeo:
                    plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
                else:
                    plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None')
                points=ginput(n=0)
                parr=np.array(points)
                if parr.shape[0]==1:
                    print('[x,y]: [%f, %f]'%(parr[0][0],parr[0][1]))
                    happy=True
                    break
                res = fitCircle(parr[:,0],parr[:,1])
                #draw the circle. now need to redraw the whole image thanks to the conda matplotlib
                fig=plt.figure(figsize=(10,10))
                if needsGeo:
                    plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
                else:
                    plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None')
                circle = plt.Circle((res['xCen'],res['yCen']),res['R'],color='b',fill=False)
                plt.gca().add_artist(circle)
                plt.plot([res['xCen'],res['xCen']],[y.min(),y.max()],'r')
                plt.plot([x.min(),x.max()],[res['yCen'],res['yCen']],'r')
                print('[x,y] (micron): [%f, %f] R (in mm): %f '%(res['xCen'],res['yCen'],res['R']/1000.))
                if raw_input("Happy with this selection:\n") in ["y","Y"]:
                    happy = True

        else:
            while not happy:
                if thresIn is not None:
                    thres = thresIn
                    thresP = np.nanpercentile(img[img!=0], thres)
                    happy = True
                else:                    
                    if not plotIt:
                        print('this is not going to work, you either need to specify a threshold or require plots')
                        return
                    fig=plt.figure(figsize=(10,10))
                    if needsGeo:
                        plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
                    else:
                        plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None')
                    thres = float(raw_input("min percentile % of selected points:\n"))

                    thresP = np.nanpercentile(img[img!=0], thres)
                    print('thresP',thresP)
                    imageThres=image.copy()
                    imageThres[image>thresP]=1
                    imageThres[image<thresP]=0

                    fig=plt.figure(figsize=(5,5))
                    plt.imshow(imageThres,clim=[-0.1,1.1])
                    if raw_input("Happy with this threshold (y/n):\n") in ["y","Y"]:
                        happy=True

            if singleTile<0:
                res = fitCircle(x.flatten()[img.flatten()>thresP],y.flatten()[img.flatten()>thresP])
            else:
                if len(x.shape)==len(img.shape) and len(x.shape)>2:
                    res = fitCircle(x[singleTile].flatten()[img[singleTile].flatten()>thresP],
                                    y[singleTile].flatten()[img[singleTile].flatten()>thresP])
                else:
                    res = fitCircle(x.flatten()[img[singleTile].flatten()>thresP],
                                    y.flatten()[img[singleTile].flatten()>thresP])
            circleM = plt.Circle((res['xCen'],res['yCen']),res['R'],color='b',fill=False)
            print('[x,y] (micron): [%f, %f] R (in mm): %f '%(res['xCen'],res['yCen'],res['R']/1000.))
            if plotIt:
                fig=plt.figure(figsize=(10,10))
                if needsGeo:
                    plt.imshow(np.rot90(image),extent=extent,clim=[plotMin,plotMax],interpolation='None')
                else:
                    plt.imshow(image,extent=extent,clim=[plotMin,plotMax],interpolation='None')
                plt.gca().add_artist(circleM)
                plt.plot([res['xCen'],res['xCen']],[y.min(),y.max()],'r')
                plt.plot([x.min(),x.max()],[res['yCen'],res['yCen']],'r')
                plt.show()
        return
    
            
    def MakeMask(self, detname=None, limits=[5,99.5], singleTile=-1, extMask=None):
        detname, img, avImage = self.getAvImage(detname=None)

        plotMax = np.nanpercentile(img, limits[1])
        plotMin = np.nanpercentile(img, limits[0])
        print('plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if self.__dict__[detname].det.dettype==DetObjectClass.Icarus:
            needsGeo=False
            if singleTile<0 or singleTile>= img.shape[0]:
                image = img.sum(axis=0)
            else:
                image = img[singleTile]
        elif needsGeo:
            image = self.__dict__[detname].det.image(self.run, img)
            if image is None and self.__dict__[detname].ped.shape[1]==1:
                image = img.squeeze()
        else:
            image = img

        det = self.__dict__[detname].det
        x = self.__dict__[detname+'_x']
        y = self.__dict__[detname+'_y']
        if x is None:
            xVec = np.arange(0, image.shape[1])
            yVec = np.arange(0, image.shape[0])
            x, y = np.meshgrid(xVec, yVec)
            self.__dict__[detname+'_x'] = x
            self.__dict__[detname+'_y'] = y
        ix = self.__dict__[detname+'_ix']
        iy = self.__dict__[detname+'_iy']
        extent=[x.min(), x.max(), y.min(), y.max()]
        #print('DEBUG: extent(x,y min,max)',extent)
        if self.__dict__[detname].det.dettype==DetObjectClass.Icarus:
            if singleTile<0 or singleTile>= img.shape[0]:
                x=x[0]
                y=y[0]
            else:
                x=x[singleTile]
                y=y[singleTile]
        
        mask=[]
        mask_r_nda=None
        select=extMask is None  # Short-circuit the loop if given an external mask.
        while select:
            fig=plt.figure(figsize=(12,10))
            gs=gridspec.GridSpec(1,2,width_ratios=[2,1])
            #print("rectangle(r-click, R-enter), circle(c), polygon(p), dark(d), noise(n) or edgepixels(e)?:")
            #needs to be pixel coordinates for rectable selection to work.
            plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')

            #edgepixels(e)?:\n")
            #polygon(p), 
            #rectangle(r-click,
            #rectangle(R-enter), 

            #circle(c), 
            #dark(d), 
            #noise(n) or 

            shape = raw_input("rectangle(r-click, R-enter), circle(c), polygon(p), dark(d), noise(n) or edgepixels(e)?:\n")
            #shape = raw_input()
            #this definitely works for the rayonix...
            if shape=='r':
                print('select two corners: ')
                p =np.array(ginput(2))
                mask_roi=np.zeros_like(image)
                p=p.astype(int)
                mask_roi[p[:,1].min():p[:,1].max(),p[:,0].min():p[:,0].max()]=1
                if needsGeo:
                    mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(ix,iy)] )
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    mask_r_nda = mask_roi
                    plt.subplot(gs[1]).imshow(mask_r_nda)
                    if self.__dict__[detname].det.dettype==DetObjectClass.Icarus:
                        maskTuple=[mask_r_nda for itile in range(self.__dict__[detname].ped.shape[0])]
                        mask_r_nda = np.array(maskTuple)
                print('mask from rectangle (shape):',mask_r_nda.shape)
            elif shape=='R':
                print('coordinates to select: ')
                htot = raw_input("horizontal,vertical h1 h2 v1 v2 ?\n")
                h = htot.split(' ');h1=float(h[0]);h2=float(h[1]);v1=float(h[2]);v2=float(h[3]);
                mask_roi=np.zeros_like(image)
                mask_roi[int(h1):int(h2),int(v1):int(v2)]=1
                if needsGeo:
                    mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(ix,iy)] )
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    mask_r_nda = mask_roi
                    plt.subplot(gs[1]).imshow(mask_r_nda)
                print('mask from rectangle (shape):',mask_r_nda.shape)
            elif shape=='c':
                fig=plt.figure(figsize=(12,10))
                gs=gridspec.GridSpec(1,2,width_ratios=[2,1])
                if det.dettype==DetObjectClass.Icarus:
                    plt.subplot(gs[0]).imshow(image,
                                              clim=[plotMin,plotMax],
                                              interpolation='None',
                                              extent=(y.min(),y.max(),x.min(),x.max()))
                else:
                    plt.subplot(gs[0]).imshow(np.rot90(image),
                                              clim=[plotMin,plotMax],
                                              interpolation='None',
                                              extent=(x.min(),x.max(),y.min(),y.max()))
                if raw_input("Select center by mouse?\n") in ["y","Y"]:
                    c=ginput(1)
                    cx=c[0][0];cy=c[0][1]
                    print('Coordinates of selected center: ',cx,' ',cy)
                else:
                    ctot = raw_input("center (x y)?\n")
                    c = ctot.split(' ');cx=float(c[0]);cy=float(c[1]);
                if raw_input("Select outer radius by mouse?\n") in ["y","Y"]: 
                    r=ginput(1)
                    rox=r[0][0];roy=r[0][1]
                    print('outer radius point: ',r[0])
                    ro=np.sqrt((rox-cx)**2+(roy-cy)**2)
                    if raw_input("Select inner radius by mouse (for donut-shaped mask)?\n") in ["y","Y"]:
                        r=ginput(1)
                        rix=r[0][0];riy=r[0][1]
                        print('inner radius point: ',r[0])
                        ri=np.sqrt((rix-cx)**2+(riy-cy)**2)
                    else:
                        ri=0
                    print('radii: ',ro,' ',ri)
                else:
                    rtot = raw_input("radii (r_outer r_inner)?\n")
                    r = rtot.split(' ');ro=float(r[0]);ri=max(0.,float(r[1]));        
                mask_router_nda = np.array( [(ix-cx)**2+(iy-cy)**2<ro**2 for ix, iy in zip(x,y)] )
                mask_rinner_nda = np.array( [(ix-cx)**2+(iy-cy)**2<ri**2 for ix, iy in zip(x,y)] )
                mask_r_nda = mask_router_nda&~mask_rinner_nda
                print('mask from circle (shape):',mask_r_nda.shape)
                if raw_input("Invert [y/n]? (n/no inversion: masked pixels will get rejected)?\n") in ["y","Y"]:
                    mask_r_nda = ~mask_r_nda

                if needsGeo:
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    plt.subplot(gs[1]).imshow(mask_r_nda)
            elif shape=='p':
                fig=plt.figure(figsize=(12,10))
                gs=gridspec.GridSpec(1,2,width_ratios=[2,1])
                if needsGeo:
                    plt.subplot(gs[0]).imshow(np.rot90(image),
                                              clim=[plotMin,plotMax],
                                              interpolation='None',
                                              extent=(x.min(),x.max(),y.min(),y.max()))
                elif det.dettype==DetObjectClass.Epix:
                    plt.subplot(gs[0]).imshow(np.rot90(image),
                                              clim=[plotMin,plotMax],
                                              interpolation='None',
                                              extent=(x.min(),x.max(),y.min(),y.max()))
                elif det.dettype==DetObjectClass.Rayonix:
                    plt.subplot(gs[0]).imshow(image,
                                              clim=[plotMin,plotMax],
                                              interpolation='None',
                                              extent=(x.min(),x.max(),y.min(),y.max()))
                elif det.dettype==DetObjectClass.Icarus:
                    plt.subplot(gs[0]).imshow(np.rot90(image),
                                              clim=[plotMin,plotMax],
                                              interpolation='None',
                                              extent=(x.min(),x.max(),y.min(),y.max()))
                else:
                    plt.subplot(gs[0]).imshow(image,
                                              clim=[plotMin,plotMax],
                                              interpolation='None',
                                              extent=(x.min(),x.max(),y.min(),y.max()))
                nPoints = int(raw_input("Number of Points (-1 until middle mouse click)?\n"))
                p=np.array(ginput(nPoints))
                print(p)
                mpath=path.Path(p)
                all_p = np.array([ (ix,iy) for ix,iy in zip(x.flatten(),y.flatten()) ] )
                mask_r_nda = np.array([mpath.contains_points(all_p)]).reshape(x.shape)
                if needsGeo:
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    plt.subplot(gs[1]).imshow(mask_r_nda)
                if self.__dict__[detname].det.dettype==DetObjectClass.Icarus:
                    maskTuple=[mask_r_nda for itile in range(self.__dict__[detname].ped.shape[0])]
                    mask_r_nda = np.array(maskTuple)
                print('mask from polygon (shape):',mask_r_nda.shape)
            elif shape=='d' or shape=='n':
                figDark=plt.figure(figsize=(12,10))
                gsPed=gridspec.GridSpec(1,2,width_ratios=[1,1])
                if shape=='d':
                    pedResult = det.pedestals(self.run)
                    if det.dettype in [DetObjectClass.Jungfrau, DetObjectClass.Epix10k,
                                       DetObjectClass.Epix10k2M, DetObjectClass.Epix10k2M_quad]:
                        pedResult = pedResult[0]
                    hstPed = np.histogram(pedResult.flatten(), np.arange(0, pedResult.max()*1.05))
                else:
                    pedResult = det.rms(self.run)
                    if det.dettype in [DetObjectClass.Jungfrau, DetObjectClass.Epix10k,
                                       DetObjectClass.Epix10k2M, DetObjectClass.Epix10k2M_quad]:
                        pedResult = pedResult[0]
                    hstPed = np.histogram(pedResult.flatten(), np.arange(0, pedResult.max()*1.05,0.05))
                if needsGeo:
                    pedResultImg = det.image(self.run,pedResult)
                else:
                    pedResultImg = pedResult.copy()
                plt.subplot(gsPed[0]).imshow(pedResultImg,clim=[np.nanpercentile(pedResult,1),np.nanpercentile(pedResult,99)])
                plt.subplot(gsPed[1]).plot(hstPed[1][:-1],np.log(hstPed[0]),'o')
                ctot=raw_input("Enter allowed pedestal range (min max)")
                c = ctot.split(' ');pedMin=float(c[0]);pedMax=float(c[1]);
                mask_r_nda=np.zeros_like(pedResult)
                mask_r_nda[pedResult<pedMin]=1
                mask_r_nda[pedResult>pedMax]=1
                mask_r_nda = (mask_r_nda.astype(bool)).astype(int)
                if needsGeo:
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda))
                else:
                    plt.subplot(gs[1]).imshow(mask_r_nda)
            if shape=='e':
                ctot=raw_input("Enter number of edge rows that should be masked:")
                try:
                    nEdge = int(ctot)
                except:
                    print('please enter an integer')
                    continue
                mask_r_nda = np.zeros_like(det.coords_x(self.run))
                if needsGeo:
                    for tile in mask_r_nda:
                        tile[0:nEdge,:]=1
                        tile[tile.shape[0]-nEdge:,:]=1
                        tile[:,0:nEdge]=1
                        tile[:,tile.shape[1]-nEdge:]=1
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda.astype(bool)))
                else:
                    tile=mask_r_nda
                    tile[0:nEdge,:]=1
                    tile[nEdge:-1,:]=1
                    tile[:,0:nEdge]=1
                    tile[:,nEdge:-1]=1
                    plt.subplot(gs[1]).imshow(mask_r_nda.astype(bool))

            if shape=='cen' or shape=='center':
                ctot=raw_input("Enter number of center rows that should be masked:")
                try:
                    nEdge = int(ctot)
                except:
                    print('please enter an integer')
                    continue
                mask_r_nda = np.zeros_like(det.coords_x(self.run))
                if needsGeo:
                    for tile in mask_r_nda:
                        tile[int(tile.shape[0]/2-nEdge):int(tile.shape[0]/2+nEdge),:]=1
                        tile[:,int(tile.shape[1]/2-nEdge):int(tile.shape[1]/2+nEdge)]=1
                    plt.subplot(gs[1]).imshow(det.image(self.run,mask_r_nda.astype(bool)))
                else:
                    tile=mask_r_nda
                    tile[int(tile.shape[0]/2-nEdge):int(tile.shape[0]/2+nEdge),:]=1
                    tile[:,int(tile.shape[1]/2-nEdge):int(tile.shape[1]/2+nEdge)]=1
                    plt.subplot(gs[1]).imshow(mask_r_nda.astype(bool))

            if mask_r_nda is not None:
                print('created a mask....',len(mask))
                mask.append(mask_r_nda.astype(bool).copy())
                countMask=1
                totmask_nm1 = mask[0]
                for thismask in mask[1:-1]:
                    totmask_nm1 = np.logical_or(totmask_nm1,thismask)
                    countMask+=1
                if len(mask)>1:
                    totmask = np.logical_or(totmask_nm1,mask[-1])
                    countMask+=1
                else:
                    totmask = totmask_nm1
            #print('DEBUG: ',mask_r_nda.shape, totmask_nm1.shape, x.shape)
            print('masked in this step: ',np.ones_like(self.__dict__[detname].x)[mask_r_nda.astype(bool)].sum())
            print('masked up to this step: ',np.ones_like(self.__dict__[detname].x)[totmask_nm1].sum())
            print('masked tot: ',np.ones_like(self.__dict__[detname].x)[totmask].sum())

            if len(mask)>1:
                fig=plt.figure(figsize=(15,9))
                gs2=gridspec.GridSpec(1,2,width_ratios=[1,1])
                plt.show()
                image_mask = img.copy(); image_mask[totmask]=0;
                #image_mask_nm1 = img.copy(); image_mask_nm1[totmask_nm1]=0;
                #print('DEBUG: ',(img-image_mask_nm1).sum(), (img-image_mask).sum())
                plt.subplot(gs2[0]).imshow(image,clim=[plotMin,plotMax])
                if needsGeo:
                    plt.subplot(gs2[1]).imshow(det.image(self.run,image_mask),clim=[plotMin,plotMax])
                else:
                    plt.subplot(gs2[1]).imshow(image_mask,clim=[plotMin,plotMax])
            else:
                fig=plt.figure(figsize=(12,10))
                plt.show()
                image_mask = img.copy(); image_mask[totmask]=0;
                if needsGeo:
                    plt.imshow(det.image(self.run,image_mask),clim=[plotMin,plotMax])
                else:
                    if (self.__dict__[detname].det.dettype==DetObjectClass.Icarus):
                        if singleTile<0 or singleTile>= img.shape[0]:
                            image = image_mask.sum(axis=0)
                        else:
                            image = [singleTile]
                        plt.imshow(image,clim=[plotMin,plotMax])
                    else:
                        plt.imshow(image_mask,clim=[plotMin,plotMax])
                

            if raw_input("Add this mask?\n") in ["n","N"]:
                mask = mask[:-1]

            if len(mask)>0:
                #remake the image with mask up to here.
                totmask = mask[0]
                for thismask in mask[1:]:
                    totmask = np.logical_or(totmask,thismask)
                image_mask = img.copy(); image_mask[totmask]=0;
                if needsGeo:
                    image = self.__dict__[detname].det.image(self.run, image_mask)
                else:
                    image = image_mask

            if raw_input("Done?\n") in ["y","Y"]:
                select = False

        #end of mask creating loop

        if extMask is not None:
            if extMask.shape != self.__dict__[detname].rawShape:
                print("extMask is the wrong shape: %s actual, %s desired." % (
                      extMask.shape, self.__dict__[detname].rawShape))
                return None
            if extMask.dtype != np.dtype('int64'):
                print("extMask is the wrong type: %s actual, %s desired." % (
                      extMask.shape, np.dtype('int64')))
                return None
            mask = [extMask]

        if len(mask)==0:
            return
        totmask = mask[0]
        for thismask in mask[1:]:
            totmask = np.logical_or(totmask,thismask)

        if extMask is None:
            prompt = "Invert [y/n]? (n/no inversion: masked pixels will get rejected)?\n"
        else:
            prompt = "Does the ext. mask have 0 for masked pixels [y/n] (n: non-zero pixels are rejected)?\n"
        if raw_input(prompt) in ["y","Y"]:
            totmask = (totmask.astype(bool)).astype(int)
        else:
            totmask = (~(totmask.astype(bool))).astype(int)

        self.__dict__['_mask_'+avImage]=totmask

        if  det.dettype == DetObjectClass.CsPad:
            mask=totmask.reshape(2,185*388).transpose(1,0)
        elif det.dettype == DetObjectClass.CsPad2M:
            mask=totmask.reshape(32*185,388)
        elif det.dettype == DetObjectClass.Epix:
            mask=totmask.reshape(704,768)
        else:
            mask=totmask
        #2x2 save as 71780 lines, 2 entries
        #cspad save as 5920 lines, 388 entries
        if raw_input("Save to calibdir?\n") in ["y","Y"]:
            srcStr=det.source.__str__().replace('Source("DetInfo(','').replace(')")','')
            calibDir=det.env.calibDir()
            if det.dettype==DetObjectClass.CsPad:
                dirname='/%s/CsPad2x2::CalibV1/%s/'%(calibDir,srcStr)
            elif det.dettype==DetObjectClass.CsPad2M:
                dirname='%s/CsPad::CalibV1/%s/'%(calibDir,srcStr)        
            elif det.dettype==DetObjectClass.Epix:
                dirname='%s/Epix100a::CalibV1/%s/'%(calibDir,srcStr)
            elif det.dettype==DetObjectClass.Rayonix:
                dirname='%s/Camera::CalibV1/%s/'%(calibDir,srcStr)        
            elif det.dettype==DetObjectClass.Epix10k:
                dirname='%s/Epix10ka::CalibV1/%s/'%(calibDir,srcStr)
            elif det.dettype==DetObjectClass.Epix10k2M:
                dirname='%s/Epix10ka2M::CalibV1/%s/'%(calibDir,srcStr)
  
            dirname+='pixel_mask/'
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            fname='%s-end.data'%self.run
            print('save mask in %s as %s '%(dirname,fname))
            #np.savetxt(dirname+fname,mask)
            det.save_txtnda(dirname+fname,mask.astype(float), fmt='%d',addmetad=True)
        elif raw_input("Save to local?\n") in ["y","Y"]:
            if len(mask.shape)<=2:
                np.savetxt('Mask_%s_%s_Run%03d.data'%(avImage,self.expname,int(self.run)), mask)
            elif len(mask.shape)==3:
                np.savetxt('Mask_%s_%s_Run%03d.data'%(avImage,self.expname,int(self.run)),
                           mask.reshape(mask.shape[0]*mask.shape[1], mask.shape[2]))
                
        return mask
    

    def addAzInt(self, 
                 detname=None, 
                 phiBins=1, 
                 qBin=0.01, 
                 eBeam=9.5, 
                 center=None, 
                 dis_to_sam=None, 
                 name='azav', 
                 Pplane=1,
                 userMask=None,
                 tx=0, 
                 ty=0, 
                 geomCorr=True, 
                 polCorr=True, 
                 inImage={}):
        if inImage!={}:
            img=inImage['image']
            if 'mask' in inImage:
                mask=inImage['mask']
            else:
                mask=np.ones_like(img)
            detname=inImage['detname']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)
        else:
            detname, img, avImage = self.getAvImage(detname=detname)

        if dis_to_sam==None:
            dis_to_sam=float(raw_input('please enter the detector distance'))
        if center==None or len(center)!=2:
            centerString=raw_input('please enter the coordinates of the beam center as c1,c2 or [c1,c2]:')
            center=[float(centerString.replace('[','').replace(']','').split(',')[0]),
                    float(centerString.replace('[','').replace(']','').split(',')[1])]

        azav = azimuthalBinning(center=center, 
                                dis_to_sam=dis_to_sam,  
                                phiBins=phiBins, 
                                eBeam=eBeam, 
                                Pplane=Pplane, 
                                userMask=userMask, 
                                qbin=qBin, 
                                tx=tx, 
                                ty=ty, 
                                geomCorr=geomCorr, 
                                polCorr=polCorr, 
                                name=name)
        getattr(self,detname).addFunc(azav)

    
    def getAzAvs(self,detname=None):
        if detname is None:
            detname, img, avImage = self.getAvImage(detname=None)   
            if detname is None:
                return
        
        azintArray = [ 
            getattr(getattr(self,detname),key) for key in getattr(self, detname).__dict__.keys() \
            if isinstance(getattr(getattr(self, detname),key), azimuthalBinning) 
        ]
        
        azintNames = [ 
            key for key in getattr(self,detname).__dict__.keys() if isinstance(getattr(getattr(self, detname),key), azimuthalBinning) 
        ]
        return azintNames, azintArray

    
    def AzInt(self, 
              detname=None, 
              use_mask=False, 
              use_mask_local=False, 
              plotIt=False, 
              azintName=None, 
              inImage={}, 
              imgName=None):
        avImage=None
        if inImage!={}:
            img = inImage['image']
            if 'mask' in inImage:
                mask = inImage['mask']
            else:
                mask = np.ones_like(img)
            avImage='input Image'
            if 'name' in inImage:
                avImage=inImage['name']
            detname = inImage['detname']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)
        else:
            detname, img, avImage = self.getAvImage(detname=detname, imgName=None)
        mask=np.ones_like(img)
        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
        if use_mask_local:
            if avImage is None:
                detname_dump, img_dump, avImage = self.getAvImage(detname=detname)
            if self.__dict__.has_key('_mask_'+avImage):
                mask = self.__dict__['_mask_'+avImage]
            else:
                print('no local mask defined for ',avImage)

        img = (img*mask)
        azIntNames,azIntegrations = self.getAzAvs(detname)
        if len(azIntegrations)>1:
            if azintName is None:
                print('we have the following options: ',azIntNames)
                azintName=raw_input('type the name of the Azimuthal integral to use:')
            azIntegs = [ iiaz for iName, iiaz in zip(azIntNames,azIntegrations) if iName==azintName]
            azInteg = azIntegs[0]
        else:
            azInteg = azIntegrations[0]
            azintName = azIntNames[0]
        azintValues = azInteg.doCake(img).squeeze()
        self.__dict__['_azint_'+azintName] = azintValues
        if plotIt:
            fig=plt.figure(figsize=(8,5))
            if len(azintValues.shape)==1:
                qVals = getattr(getattr(getattr(self,detname), azintName),'q')
                plt.plot(qVals,azintValues,'o')
            elif len(azintValues.shape)==2:
                plt.imshow(azintValues,aspect='auto',interpolation='none')
                plt.colorbar()
        else:
            return azintValues

    def plotAzInt(self, detname=None, azintName=None):
        detname, img, avImage = self.getAvImage(detname=detname)
        azIntNames,azIntegrations = self.getAzAvs(detname)
        azint=''
        if len(azIntegrations)>1:
            if azintName is None:
                print('we have the following options: ',azIntNames)
                azintName=raw_input('type the name of the Azimuthal integral to use:')
            azIntegs = [ iiaz for iName, iiaz in zip(azIntNames,azIntegrations) if iName==azintName]
            azInteg = azIntegs[0]
            azint = azintName
        else:
            azInteg = azIntegrations[0]
            azint = azIntNames[0]
        if azint=='':
            print('did not find azimuthal integral asked for')
            return
        else:
            if ('_azint_'+azint) not in self.__dict__.keys():
                print('did not find azint ',azint,', all keys are: ',self.__dict__.keys())
        try:
            azintValues = self.__dict__['_azint_'+azint]
            fig=plt.figure(figsize=(8,5))
            if len(azintValues.shape)==1:
                qVals = getattr(getattr(getattr(self,detname), azint),'q')
                print(azintValues.shape, qVals.shape)
                plt.plot(qVals,azintValues,'o')
            elif len(azintValues.shape)==2:
                plt.imshow(azintValues,aspect='auto',interpolation='none')
                plt.colorbar()
        except:
            pass

    def AzInt_centerVar(self, 
                        detname=None, 
                        use_mask=False, 
                        center=None, 
                        data=None, 
                        varCenter=110., 
                        zoom=None, 
                        qBin=0.001, 
                        phiBins=13, 
                        dis_to_sam=1000., 
                        inImage={}, 
                        imgName=None, 
                        plotGrid=0):
        if inImage!={}:
            img=inImage['image']
            if 'mask' in inImage:
                mask=inImage['mask']
            else:
                mask=np.ones_like(img)
            detname=inImage['detname']
            if not detname in self.__dict__.keys():
                self.addDetInfo(detname=detname)
        else:
            detname, img, avImage = self.getAvImage(detname=detname, imgName=imgName)
        if use_mask:
            mask = self.__dict__[detname].det.mask_calib(self.run)
            img = (img*mask)

        if center==None or len(center)!=2:
            centerString=raw_input('please enter the coordinates of the beam center as c1,c2 or [c1,c2]:')
            center=[float(centerString.replace('[','').replace(']','').split(',')[0]),
                    float(centerString.replace('[','').replace(']','').split(',')[1])]
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, 
                      center=[center[0],center[1]], dis_to_sam=dis_to_sam, 
                      name='c00', inImage=inImage)
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, 
                      center=[center[0]-varCenter,center[1]], dis_to_sam=dis_to_sam, 
                      name='cm10', inImage=inImage)
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, 
                      center=[center[0]+varCenter,center[1]], dis_to_sam=dis_to_sam, 
                      name='cp10', inImage=inImage)
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, 
                      center=[center[0],center[1]-varCenter], dis_to_sam=dis_to_sam, 
                      name='cm01', inImage=inImage)
        self.addAzInt(detname=detname, phiBins=phiBins, qBin=qBin, eBeam=9.5, 
                      center=[center[0],center[1]+varCenter], dis_to_sam=dis_to_sam, 
                      name='cp01', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='c00', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='cm10', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='cp10', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='cm01', inImage=inImage)
        self.AzInt(detname=detname, use_mask=use_mask, azintName='cp01', inImage=inImage)

        azintValues_c00 = self.__dict__['_azint_c00']
        azintValues_cm10 = self.__dict__['_azint_cm10']
        azintValues_cp10 = self.__dict__['_azint_cp10']
        azintValues_cm01 = self.__dict__['_azint_cm01']
        azintValues_cp01 = self.__dict__['_azint_cp01']
        
        try:
            fig=plt.figure(figsize=(10,6))
            from matplotlib import gridspec
            p33_11 = plt.subplot2grid((3,3),(1,1))
            p33_10 = plt.subplot2grid((3,3),(1,0))
            p33_12 = plt.subplot2grid((3,3),(1,2))
            p33_01 = plt.subplot2grid((3,3),(0,1))
            p33_21 = plt.subplot2grid((3,3),(2,1))
        except:
            print('failed with getting plot ready')
            return

        while 1:
            ymin=0;ymax=azintValues_c00.shape[1]-1
            maxVal = np.nanpercentile(azintValues_c00, 99.99)
            if zoom is not None:
                if len(zoom)==1:
                    maxVal = np.nanpercentile(azintValues_c00, zoom[0])
                    zoom = None
                elif len(zoom)==3:
                    maxVal = np.nanpercentile(azintValues_c00, zoom[2])
                else:
                    if isinstance(zoom, list) and isinstance(zoom[0], int) and zoom[0]>=0\
                    and zoom[1]>=0 and zoom[0]!=zoom[1]:
                        ymin=zoom[0]
                        ymax=zoom[1]
                    else:
                        yString=raw_input('please enter x-boundaries of the zoomed figure as c1,c2 or [c1,c2]:')
                        ymin=int(yString.replace('[','').replace(']','').split(',')[0])
                        ymax=int(yString.replace('[','').replace(']','').split(',')[1])
            print('we will plot from bin %d to %d '%(ymin, ymax))

            try:
                p33_11.imshow(azintValues_c00[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                p33_10.imshow(azintValues_cm10[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                p33_12.imshow(azintValues_cp10[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                p33_01.imshow(azintValues_cm01[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                p33_21.imshow(azintValues_cp01[:,ymin:ymax],aspect='auto',interpolation='none',vmax=maxVal)
                maxVal_c00 =  np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_c00[:,ymin:ymax])])
                maxVal_cm10 = np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_cm10[:,ymin:ymax])])
                maxVal_cm01 = np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_cm01[:,ymin:ymax])])
                maxVal_cp10 = np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_cp10[:,ymin:ymax])])
                maxVal_cp01 = np.array([[ip,np.nanargmax(phiBin)] for ip,phiBin in enumerate(azintValues_cp01[:,ymin:ymax])])
                try:
                    p33_11.plot(maxVal_c00[:,1], maxVal_c00[:,0],'r+')
                    p33_10.plot(maxVal_cm10[:,1], maxVal_cm10[:,0],'r+')
                    p33_12.plot(maxVal_cp10[:,1], maxVal_cp10[:,0],'r+')
                    p33_01.plot(maxVal_cm01[:,1], maxVal_cm01[:,0],'r+')
                    p33_21.plot(maxVal_cp01[:,1], maxVal_cp01[:,0],'r+')
                    if plotGrid==0:
                        cMaxVal_c00 = np.nanargmax(azintValues_c00[:,ymin:ymax].sum(axis=0))
                        p33_11.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                        p33_10.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                        p33_12.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                        p33_01.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                        p33_21.plot([cMaxVal_c00,cMaxVal_c00],[0, azintValues_c00.shape[0]],'m')
                    else:
                        for il,lineVal in enumerate(np.linspace(ymin,ymax,int(plotGrid))):
                            lw=1
                            if plotGrid>15 and il%5>0:
                                lw=0.5
                            p33_11.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                            p33_10.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                            p33_12.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                            p33_01.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                            p33_21.plot([lineVal,lineVal],[0, azintValues_c00.shape[0]],'b',linewidth=lw,linestyle='dotted')
                except:
                    print('failed at plotting')

                if ymin==0 and ymax==azintValues_c00.shape[1]-1:
                    break
                if raw_input("done? (y/n):\n") in ["y","Y"]:
                    break
            except:
                print('failed in try-except.')
                break

                
    def SelectRegionDroplet(self, detname=None, limits=[5,99.5]):
        avImages=[]
        for key in self.__dict__.keys():
            if key.find('AvImg')==0:
                if detname is not None and key.find(detname)>=0:
                    avImages.append(key)
                elif detname is None:
                    avImages.append(key)
        if len(avImages)==0:
            print('please create the AvImage first!')
            return
        elif len(avImages)>1:
            print('we have the following options: ',avImages)
            avImage=raw_input('type the name of the AvImage to use:')
        else:
            avImage=avImages[0]
        detname = self._getDetName_from_AvImage(avImage)
        img = self.__dict__[avImage]
        print(img.shape)

        plotMax = np.nanpercentile(img, limits[1])
        plotMin = np.nanpercentile(img, limits[0])
        print('plot %s using the %g/%g percentiles as plot min/max: (%g, %g)'%(avImage,limits[0],limits[1],plotMin,plotMax))

        needsGeo=False
        if self.__dict__[detname].ped.shape != self.__dict__[detname].imgShape:
            needsGeo=True

        if needsGeo:
            image = self.__dict__[detname].det.image(self.run, img)
        else:
            image = img

        det = self.__dict__[detname].det
        x = self.__dict__[detname+'_x']
        y = self.__dict__[detname+'_y']
        ix = self.__dict__[detname+'_ix']
        iy = self.__dict__[detname+'_iy']
        extent=[x.min(), x.max(), y.min(), y.max()]

        fig=plt.figure(figsize=(10,6))
        from matplotlib import gridspec
        gs=gridspec.GridSpec(1,2,width_ratios=[1,1])
        
        mask=None
        mask_r_nda=None

        plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')

        print('select two corners: ')
        p =np.array(ginput(2))
        mask_roi=np.zeros_like(image)
        mask_roi[p[:,1].min():p[:,1].max(),p[:,0].min():p[:,0].max()]=1
        if needsGeo:
            mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(ix,iy)] )
        else:
            mask_r_nda = mask_roi

        if mask_r_nda is not None:
            #print('created a mask....')
            if mask is None:
                mask = mask_r_nda.astype(bool).copy()
            else:
                mask = np.logical_or(mask,mask_r_nda)
        print('masked: ',np.ones_like(x)[mask.astype(bool)].sum())

        image_mask = img.copy(); image_mask[~mask]=0;
        if needsGeo:
            plt.subplot(gs[1]).imshow(det.image(self.run,image_mask),clim=[plotMin,plotMax])
        else:
            plt.subplot(gs[1]).imshow(image_mask,clim=[plotMin,plotMax])

        if  det.dettype == DetObjectClass.CsPad:
            mask=mask.reshape(2,185*388).transpose(1,0)
        elif det.dettype == DetObjectClass.CsPad2M:
            mask=mask.reshape(32*185,388)
        elif det.dettype == DetObjectClass.Epix:
            mask=mask.reshape(704,768)

        #2x2 save as 71780 lines, 2 entries
        #cspad save as 5920 lines, 388 entries
        if raw_input("Save to calibdir?\n") in ["y","Y"]:
            mask = (~(mask.astype(bool))).astype(int)
            srcStr=det.source.__str__().replace('Source("DetInfo(','').replace(')")','')
            if det.dettype==DetObjectClass.CsPad:
                dirname='/reg/d/psdm/%s/%s/calib/CsPad2x2::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)
            elif det.dettype==DetObjectClass.CsPad2M:
                dirname='/reg/d/psdm/%s/%s/calib/CsPad::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            elif det.dettype==DetObjectClass.Epix:
                dirname='/reg/d/psdm/%s/%s/calib/Epix100a::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            elif det.dettype==DetObjectClass.Epix10k:
                dirname='/reg/d/psdm/%s/%s/calib/Epix10ka::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            elif det.dettype==DetObjectClass.Epix10k2M:
                dirname='/reg/d/psdm/%s/%s/calib/Epix10ka2M::CalibV1/%s/pixel_mask/'%(self.sda.expname[:3],self.sda.expname,srcStr)        
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            fname='%s-end.data'%self.run
            det.save_txtnda(dirname+fname,mask, fmt='%d',addmetad=True)
        elif raw_input("Save to local?\n") in ["y","Y"]:
            mask = (~(mask.astype(bool))).astype(int)
            np.savetxt('%s_mask_run%s.data'%(self.sda.expname,self.run),mask)
        return mask


