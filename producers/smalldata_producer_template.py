#!/usr/bin/env python

import numpy as np
import psana
import time
from datetime import datetime
begin_job_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
import argparse
import socket
import os
import logging 
import requests
import sys
from glob import glob
from PIL import Image
from requests.auth import HTTPBasicAuth

########################################################## 
##
## User Input start --> 
##
########################################################## 
##########################################################
# functions for run dependant parameters
##########################################################

# 1) REGIONS OF INTEREST
def getROIs(run):
    """ Set parameter for ROI analysis. Set writeArea to True to write the full ROI in the h5 file.
    See roi_rebin.py for more info
    """
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run>0:
        roi_dict = {}
        roi_dict['ROIs'] = [ [[1,2], [157,487], [294,598]] ] # can define more than one ROI
        roi_dict['writeArea'] = True
        roi_dict['thresADU'] = None
        ret_dict['jungfrau1M'] = roi_dict
    return ret_dict


# 2) AZIMUTHAL INTEGRATION
def getAzIntParams(run):
    """ Parameters for azimuthal integration
    See azimuthalBinning.py for more info
    """
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run>0:
        az_dict = {'eBeam': 18.0} # keV
        az_dict['center'] = [87526.79161840, 92773.3296889500] # um
        az_dict['dis_to_sam'] = 80. # mm
        az_dict['tx'] = 0 # deg
        az_dict['ty'] = 0 # deg
        ret_dict['jungfrau1M'] = az_dict
    return ret_dict

def getAzIntPyFAIParams(run):
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run>0:
        az_dict = {}
        pix_size = 176e-6
        ai_kwargs = {'dist':1, 'poni1':960*pix_size, 'poni2':960*pix_size}
        az_dict['ai_kwargs'] = ai_kwargs
        az_dict['npts'] = 512
        az_dict['int_units'] = '2th_deg'
        az_dict['return2d'] = False
        ret_dict['Rayonix'] = az_dict
    return ret_dict


# 3) PHOTON COUNTING AND DROPLET
# Photon
def getPhotonParams(run):
    """ Parameters for droplet algorithm
    See photons.py for more info
    """
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run>0:
        photon_dict = {}
        photon_dict['ADU_per_photon'] = 9.5
        photon_dict['thresADU'] = 0.8 # fraction of ADU_per_photon
        ret_dict['jungfrau1M'] = photon_dict
    return ret_dict


# Droplet algorithm
def getDropletParams(run):
    """ Parameters for droplet algorithm
    See droplet.py for more info
    """
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run>0:
        droplet_dict = {
            'theshold': 5,
            'thresholdLow':5,
            'thresADU':60,
            'useRms':True,
            'nData': 1e5 # int or None: sparsify arg
        }
        ret_dict['epix_1'] = droplet_dict
    return ret_dict


# Droplet to photon algorithm (greedy guess)
def getDroplet2Photons(run):
    """ Set parameter for droplet2photon analysis. The analysis uses two functions, each with their
    own dictionary of argument:
        1. droplet: see droplet.py for args documentation
        2. photonize droplets: see droplet2Photons.py for args documentation
    """
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run>0:
        d2p_dict = {}
        # droplet args
        d2p_dict['droplet'] = {
            'theshold': 10,
            'thresholdLow':3,
            'thresADU':10,
            'useRms':True
        }
        # droplet2Photons args
        d2p_dict['d2p'] = {
            'aduspphot': 162,
            'mask': np.load('/reg/d/psdm/xpp/xppx49520/results/haoyuan/' +
                             'check_data/haoyuan_epix1_roi.npy'),
            'cputime': True
        }
        d2p_dict['nData'] = 3e4

        ret_dict['epix_alc1'] = d2p_dict
    return ret_dict


# 4) WAVEFORM ANALYSIS (SVD, peak finding)
def getSvdParams(run):
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run>0:
        svd_dict = {}
        svd_dict['basis_file'] = None
        svd_dict['n_pulse'] = 1
        svd_dict['delay'] = None
        svd_dict['return_reconstructed'] = True
        ret_dict['acq_0'] = svd_dict
    return ret_dict


# 5) autocorrelation
def getAutocorrParams(run):
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run>0:
        autocorr_dict = {}
#         autocorr_dict['ROIs'] = [ [[100,200], [100,200]] ] # can define more than one ROI
        autocorr_dict['mask'] = '/cds/home/e/espov/dataAna/mask_epix.npy' # path to mask saved as a npy file. Can define multiple mask if 3D.
        autocorr_dict['thresADU'] = [72.,1e6]
        autocorr_dict['save_range'] = [70,50] # range to save around the autcorrelation center
        autocorr_dict['save_lineout'] = True # save horiz / vert lineout through center instead of autocorr array
        ret_dict['epix_2'] = autocorr_dict
    return ret_dict



##########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
#epicsPV = ['slit_s1_hw'] 
epicsPV = []
#fix timetool calibration if necessary
#ttCalib=[0.,2.,0.]
ttCalib=[]
#ttCalib=[1.860828, -0.002950]
#decide which analog input to save & give them nice names
#aioParams=[[1],['laser']]
aioParams=[]
########################################################## 
##
## <-- User Input end
##
##########################################################


# DEFINE DETECTOR AND ADD ANALYSIS FUNCTIONS
def define_dets(run):
    detnames = ['jungfrau1M', 'epix_alc1'] # add detector here
    dets = []
    
    # Load DetObjectFunc parameters (if defined)
    try:
        ROIs = getROIs(run)
    except Exception as e:
        print(f'Can\'t instantiate ROI args: {e}')
        ROIs = []
    try:
        az = getAzIntParams(run)
    except Exception as e:
        print(f'Can\'t instantiate azimuthalBinning args: {e}')
        az = []
    try:
        az_pyfai = getAzIntPyFAIParams(run)
    except Exception as e:
        print(f'Can\'t instantiate AzIntPyFAI args: {e}')
        az_pyfai = []
    try:
        phot = getPhotonParams(run)
    except Exception as e:
        print(f'Can\'t instantiate Photon args: {e}')
        phot = []
    try:
        drop = getDropletParams(run)
    except Exception as e:
        print(f'Can\'t instantiate Droplet args: {e}')
        drop = []
    try:
        drop2phot = getDroplet2Photons(run)
    except Exception as e:
        print(f'Can\'t instantiate Droplet2Photons args: {e}')
        drop2phot = []
    try:
        auto = getAutocorrParams(run)
    except Exception as e:
        print(f'Can\'t instantiate Autocorrelation args: {e}')
        auto = []
    try:
        svd = getSvdParams(run)
    except Exception as e:
        print(f'Can\'t instantiate SVD args: {e}')
        svd = []
        
    # Define detectors and their associated DetObjectFuncs
    for detname in detnames:
        havedet = checkDet(ds.env(), detname)
        # Common mode
        if havedet:
            if detname=='': 
                # change here to specify common mode for detname if desired. Else default is used
                common_mode=0
            else:
                common_mode=None
            det = DetObject(detname ,ds.env(), int(run), common_mode=common_mode)
            
            # Analysis functions
            # ROIs:
            if detname in ROIs:
                for iROI,ROI in enumerate(ROIs[detname]['ROIs']):
                    det.addFunc(ROIFunc(name='ROI_%d'%iROI,
                                        ROI=ROI,
                                        writeArea=ROIs[detname]['writeArea'],
                                        thresADU=ROIs[detname]['thresADU']))
            # Azimuthal binning
            if detname in az:
                det.addFunc(azimuthalBinning(**az[detname]))
            if detname in az_pyfai:
                det.addFunc(azav_pyfai(**az_pyfai[detname]))
            # Photon count
            if detname in phot:
                det.addFunc(photonFunc(**phot[detname]))
            # Droplet algo
            if detname in drop:
                if nData in drop[detame]:
                    nData = drop[detname].pop('nData')
                else:
                    nData = None
                func = dropletFunc(**drop[detname])
                func.addFunc(roi.sparsifyFunc(nData=nData))
                det.addFunc(func)
            # Droplet to photons
            if detname in drop2phot:
                if 'nData' in drop2phot[detname]:
                    nData = drop2phot[detname].pop('nData')
                else:
                    nData = None
                # getp droplet dict
                droplet_dict = drop2phot[detname]['droplet']
                #get droplet2Photon dict
                d2p_dict = drop2phot[detname]['d2p']
                dropfunc = dropletFunc(**droplet_dict)
                drop2phot = droplet2Photons(**d2p_dict)
                # add sparsify to put photon coord to file
                sparsify = sparsifyFunc(nData=nData)
                # assemble function stack: droplet -> photon -> sparsify
                drop2phot.addFunc(sparsify)
                dropfunc.addFunc(drop2phot)
                det.addFunc(dropfunc)
            # Autocorrelation
            if detname in auto:
                det.addFunc(Autocorrelation(**auto[detname]))
            # SVD waveform analysis
            if detname in svd:
                det.addFunc(svdFit(**svd[detname]))

            det.storeSum(sumAlgo='calib')
            #det.storeSum(sumAlgo='calib_img')
            dets.append(det)
    return dets