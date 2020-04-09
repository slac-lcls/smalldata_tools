from smalldata_tools.SmallDataUtils import setParameter, defaultDetectors, detData
from smalldata_tools.SmallDataDefaultDetector import ttRawDetector, wave8Detector, epicsDetector

########################################################## 
##
## Run dependent functions to set data reduction params
##
########################################################## 
def getNmaxDrop(run):
    if isinstance(run,basestring):
        run=int(run)

    if run >= 10:
        return 2000
    else:
        return 400

def getAzIntParams(run):
    if isinstance(run,basestring):
        run=int(run)
        
    ret_dict = {'eBeam': 9.5}
    ret_dict['cspad_center'] = [87526.79161840, 92773.3296889500]
    ret_dict['cspad_dis_to_sam'] = 80.
    return ret_dict

def getROI(run):
    if isinstance(run,basestring):
        run=int(run)
    if run <=6:
        return [ [[0,1], [1,74], [312,381]],
                 [[8,9], [8,89], [218,303]] ]
    else:
        return [ [[0,1], [1,74], [312,381]],
                 [[8,9], [8,89], [218,303]] ]

########################################################## 
##
## Definition of area detectors & functions
##
########################################################## 
def addROIDet(detname, run, env):
    det = None
    ROI = getROI(int(run))
    if checkDet(env, detname):
        det = DetObject(detname ,env, int(run))#, name='Rowland')
        for iROI,ROI in enumerate(ROI):
            det.addFunc(ROIFunc(ROI=ROI, name='ROI_%d'%iROI))
    return det

def addAzIntDet(detname, run, env):
    det = None
    azIntParams = getAzIntParams(run)
    haveDet = checkDet(ds.env(), detname)
    if haveDet:
        cspad = DetObject(detname , env, int(run))#, name=detname)

    cspad.azav_eBeam=azIntParams['eBeam']
    if azIntParams.has_key('cspad_center'):
        try:
            azav = azimuthalBinning(center=azIntParams['cspad_center'], azIntParams['cspad_dis_to_sam'], phiBins=11, Pplane=0)
            cspad.addFunc(azav)
        except:
            pass

    cspad.storeSum(sumAlgo='calib')
    cspad.storeSum(sumAlgo='square')
    return det

def addDropletDet(detname, run, env):
    det = None
    nMaxDrop = getNmaxDrop(run)
    #create detector object. needs run for calibration data
    #common mode: 46 is the "original" to psana method 7(?)
    #row & column correction on data w/ photons&neighbors removed.
    det = DetObject(epixname , env, int(run), common_mode=46)

    #two threshold droplet finding.
    #for data w/ > 1 photon energy this is the only thing that will work.
    #Tends to add photons together into single droplet if occupancy
    #is not low, might need photonizing step to get single photon positions
    droplet = droplet(threshold=10., thresholdLow=3., thresADU=0.,name='droplet')
    #droplet.add_aduHist() 
    droplet.addDropletSave(maxDroplets=nDrop)
    det.addFund(droplet)

    #now add photon algorithms. Only works for single energy photon data
    # ADU_per_photon: expected ADU for photon of expected energy
    # thresADU: fraction of photon energy in photon candidate
                #(2 neighboring pixel)
    #retImg: 0 (def): only histogram of 0,1,2,3,...,24 photons/pixel is returned
    #        1 return Nphot, x, y arrays
    #        2 store image using photons /event
    if (int(run)==444):
        det.addFunc(photon(ADU_per_photon=165, thresADU=0.9, retImg=2, nphotMax=200))

    return det

#
# add all detectors to a list.
#
def addDets(env, run):
    dets=[]
    #use a 'with' construct here?
    #det = None

    dets.append(addROIDet('roiDet'), run, env)
    dets.append(addDropletDet('dropletDet'), run, env)
    return dets

########################################################## 
##
## Run Independent Functions for small detectors
##
########################################################## 

#
# here you would add PVs that you expect to change during the run
# you can use the alias as it also appears in ami
# there will be a reading on each event, but this is NOT properly
# timestamped: readings will only go as fast as EPICS and contain
# extra uncertainty as to when the DAQ read the PV
#
# an example are temperatures read by controllers, or scan motors 
# without access to the interpolator jitter/lag of motion could affect
# your experiment
def addEpicsPV_shotbyshot():
    return []


# 
# here you can add timetool calibration parameters of up to pol2 form
# to override the calibration during data taking
# passing an empty list will default back to the original calibration
def ttCalib():
    #ttCalib=[0.,2.,0.]
    ttCalib=[]
    return ttCalib


#
# decide which analog input to save & give them nice names
def aioParams():
    #aioParams=[[1],['laser']]
    aioParams=[]
    return aioParams

def addTimetoolTraces():
    return False;
