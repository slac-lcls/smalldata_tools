# importing generic python modules
import numpy as np
#import psana
from SmallDataDefaultDetector import *

def defaultDetectors(hutch):
    if hutch.lower()=='amo':
        dets = amoDetectors()
    elif hutch.lower()=='sxr':
        dets = sxrDetectors()
    elif hutch.lower()=='xpp':
        dets = xppDetectors()
    elif hutch.lower()=='xcs':
        dets = xcsDetectors()
    elif hutch.lower()=='mfx':
        dets = mfxDetectors()
    elif hutch.lower()=='cxi':
        dets = cxiDetectors()
    elif hutch.lower()=='mec':
        dets = mecDetectors()
    else:
        dets = []
    detsInRun= [ det for det in dets if det.inRun() ]
    #print 'found %d detectors '%len(detsInRun)
    return detsInRun

def amoDetectors():
    return []

def sxrDetectors():
    return []

def xppDetectors(beamCodes=[[162],[91]]):
    dets=[]
    dets.append(lightStatus(codes=beamCodes))
    dets.append(epicsDetector(PVlist=['att_T', 'att_T3rd', 'slit_s1_hw', 'slit_s1_vw', 'slit_s2_hw', 'slit_s2_vw', 'slit_s3_hw', 'slit_s3_vw', 'slit_s4_hw', 'slit_s4_vw', 'lxt_vitara', 'lxt', 'lxt_ttc', 'lxe', 'ccm_E', 'lom_E', 'lom_EC', 'gon_v', 'gon_h', 'gon_r', 'gon_x', 'gon_y', 'gon_z', 'gon_roll', 'gon_pitch', 'gon_kappa_eta', 'gon_kappa_kappa', 'gon_kappa_phi', 'gon_kappa_samx','gon_kappa_samy', 'gon_kappa_samz', 'gon_sam_phi', 'gon_sam_z', 'robot_x', 'robot_y', 'robot_z', 'robot_rx', 'robot_ry', 'robot_rz', 'robot_azi', 'robot_ele', 'robot_rad', 'las_comp_wp', 'las_opa_wp']))
    dets.append(controlDetector())
    dets.append(ipmDetector('NH2-SB1-IPM-01','ipm1'))
    dets.append(ipmDetector('NH2-SB1-IPM-02','ipm1c'))
    dets.append(ipmDetector('XppMon_Pim0','lombpm'))
    dets.append(ipmDetector('XppMon_Pim1','lomdiode'))
    dets.append(ipmDetector('XppSb2_Ipm','ipm2'))
    dets.append(ipmDetector('XppSb3_Ipm','ipm3'))
    dets.append(ipmDetector('XppSb3_Pim','diode2'))
    dets.append(ipmDetector('XppSb4_Pim','diode2'))
    dets.append(ipmDetector('XppEnds_Ipm0','diodeU'))
    dets.append(aiDetector('XPP-AIN-01','ai'))
    dets.append(adcDetector('adc','adc'))
    dets.append(ttDetector(baseName='XPP:TIMETOOL:'))
    #dets.append(ttDetector(baseName='TTSPEC:'))
    dets.append(bmmonDetector('HX2-BEAMMON-01','ipm_hx2'))
    try:
        dets.append(encoderDetector('usbencoder','enc'))
    except:
        print 'did not add encoder detector'
        pass
    dets.append(damageDetector())
    setParameter(dets, dets, detName='damage')
    return dets

def xcsDetectors(beamCodes=[[162],[81]]):
    dets=[]
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(ipmDetector('XCS-IPM-02','ipm2'))
    dets.append(ipmDetector('XCS-IPM-03','ipm3'))
    dets.append(ipmDetector('XCS-IPM-04','ipm4'))
    dets.append(ipmDetector('XCS-IPM-05','ipm5'))
    dets.append(ipmDetector('XCS-DIO-03','dio3'))
    dets.append(ipmDetector('XCS-DIO-05','dio5'))
    dets.append(ipmDetector('XCS-IPM-mono','diodeMono'))
    dets.append(ipmDetector('XCS-IPM-gon','diodeGon'))
    #dets.append(ipmDetector('XCS-IPM-lam','diodeLadm'))
    #dets.append(ipmDetector('XcsUsrIpm01','diodeU'))
    dets.append(bmmonDetector('HX2-BEAMMON-01','ipm_hx2'))
    dets.append(aiDetector('XCS-AIN-01','ai'))
    dets.append(epicsDetector(PVlist=['att_transmission', 'att_transmission_3rd_h', 'ccm_E', 'lom_E', 'DIFF_phis', 'DIFF_th', 'DIFF_tth', 'DIFF_xs', 'DIFF_ys', 'DIFF_zs', 'DIFF_x', 'DIFF_y', 'DIFF_chis','DIFF_dety','ladm_theta','LAM_Z','LAM_X1','LAM_X2','LAM_Y1','LAM_Y2','LAM_DET_Y','LAM_DET_X']))
    dets.append(ttDetector(baseName='XCS:TIMETOOL:'))
    try:
        dets.append(encoderDetector('XCS-USB-ENCODER-01','enc'))
    except:
        print 'did not add encoder detector'
        pass
    dets.append(damageDetector())
    setParameter(dets, dets, detName='damage')
    return dets

def mfxDetectors():
    return []

def cxiDetectors():
    return []

def mecDetectors():
    return []


def detData(detList, evt):
    #mpiDataSource Issue? cpo, tjlane
    #if one of the keys here contains an epty dict (say: no user PV that are in the data)
    #then it will not be saved. Also any keys _after_ that dir will be lost!
    #test: misspell userPV in producer file and see that the scan directory disappears....
    data={}
    for det in detList:
        try:
            data[det.name] = det.data(evt)
        except:
            #print 'could not get data in this event for detector ',det.name
            pass
    return data

def setParameter(detList, Params, detName='tt'):
    for det in detList:
        if det.name==detName:
            det.setPars(Params)

#userData functions
def getUserData(det):
    det_dict= {}
    try:
        userData_keys = [ key for key in det.evt.__dict__.keys() if key.find('write_')>=0 ]
        #print 'DEBUG SmallDataUtils: getting keys: ',userData_keys
        for key in userData_keys:
            #print 'DEBUG SmallDataUtils: getting data for key: ',key
            if isinstance(det.evt[key], np.ndarray):
                if isinstance(det.evt[key], np.ma.masked_array):
                    data = det.evt[key].data
                    if np.issubdtype(data.dtype, np.integer):
                        data[det.evt[key].mask]=0
                    else:
                        data[det.evt[key].mask]=np.nan
                    #print eventNr,' is ndarray ',key
                    det_dict[key.replace('write_','')] = data
                else:
                    det_dict[key.replace('write_','')] = det.evt[key]
            else:
                det_dict[key.replace('write_','')] = np.array(det.evt[key])
    except:
        pass
    return det_dict

def getUserEnvData(det):
    env_dict= {}
    try:
        envData_keys = [ key for key in det.evt.__dict__.keys() if key.find('env_')>=0 ]
        for key in envData_keys:
            if isinstance(det.evt[key], np.ndarray):
                if isinstance(det.evt[key], np.ma.masked_array):
                    data = det.evt[key].data
                    data[det.evt[key].mask]=np.nan
                    env_dict[key.replace('env_','')] = data
                else:
                    env_dict[key.replace('env_','')] = det.evt[key]
            else:
                env_dict[key.replace('env_','')] = np.array(det.evt[key])
    except:
        pass
    return env_dict

#get userData configuration data
def getCfgOutput(det):
    cfgDict={}
    #baseName=det._name+'_'
    baseName=''
    for key in det.__dict__.keys():
        if key == 'mask_ROI' or key == 'mask_ROI_shape' or key == 'bankMasks':
            continue
        if isinstance(det[key], list) or isinstance(det[key], np.ndarray):
            if isinstance(det[key], list):
                cfgDict[baseName+key] = np.array(det[key])
            elif isinstance(det[key], np.ndarray):
                cfgDict[baseName+key] = np.array(det[key])
        elif isinstance(det[key], float) or isinstance(det[key], int):
            cfgDict[baseName+key] = np.array([det[key]])
    #now add ROI boundaries (so we can use mask, rms from main det object later....)
    for ROI in det.getROIs():
        cfgDict[baseName+ROI.name+'_bound'] = np.array(ROI.bound)

    for drops in det.getDroplets():
        for aduHist in drops.aduHists:
            cfgDict[baseName+aduHist.name+'_bins'] = np.array(aduHist.bins)
        
    return cfgDict

