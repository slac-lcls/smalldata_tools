# importing generic python modules
import numpy as np
#import psana
from smalldata_tools.SmallDataDefaultDetector import *
from smalldata_tools.epicsarchive import EpicsArchive
from datetime import datetime

def defaultDetectors(hutch, run=None, env=None):
    if hutch.lower()=='amo':
        dets = amoDetectors()
    elif hutch.lower()=='sxr':
        dets = sxrDetectors()
    elif hutch.lower()=='xpp':
        dets = xppDetectors(env=env)
    elif hutch.lower()=='xcs':
        dets = xcsDetectors(env=env)
    elif hutch.lower()=='mfx':
        dets = mfxDetectors(env=env)
    elif hutch.lower()=='cxi':
        dets = cxiDetectors(env=env)
    elif hutch.lower()=='mec':
        dets = mecDetectors(env=env)
    elif hutch.lower()=='det':
        dets = detDetectors()
    elif hutch.lower()=='dia':
        dets = diaDetectors()
    elif hutch.lower()=='tmo':
        dets = tmoDetectors(run)
    elif hutch.lower()=='rix':
        dets = rixDetectors(run)
    elif hutch.lower()=='ued':
        dets = uedDetectors(run)
    else:
        dets = []
    detsInRun= [ det for det in dets if det.inRun() ]
    #print('found %d detectors '%len(detsInRun))
    return detsInRun

def amoDetectors():
    return []

def tmoDetectors(run, beamCodes=[[162],[91]]):
    dets=[]
    dets.append(scanDetector('scan', run))
    dets.append(genlcls2Detector('timing',run))
    dets.append(lcls2_lightStatus(beamCodes,run))
    dets.append(genlcls2Detector('gmd',run))
    dets.append(genlcls2Detector('xgmd',run))
    dets.append(genlcls2Detector('ebeam',run))
    dets.append(genlcls2Detector('pcav',run))
    dets.append(fimfexDetector('tmo_fim0',run))
    dets.append(fimfexDetector('tmo_fim1',run))
    dets.append(ttlcls2Detector('tmoopal2',run, saveTraces=True))
    dets.append(lcls2_epicsDetector(PVlist=['MR1K4_pitch', 'MR2K4_pitch'],run=run))

    return dets

def rixDetectors(run, beamCodes=[[-136],[77]]):
    dets=[]
    dets.append(scanDetector('scan', run))
    dets.append(genlcls2Detector('timing',run))
    dets.append(lcls2_lightStatus(beamCodes,run))
    dets.append(genlcls2Detector('gmd',run))
    dets.append(genlcls2Detector('xgmd',run))
    dets.append(genlcls2Detector('ebeam',run))
    dets.append(genlcls2Detector('pcav',run))
    dets.append(fimfexDetector('rix_fim0',run))
    dets.append(fimfexDetector('rix_fim1',run))
    dets.append(fimfexDetector('rix_fim2',run))
    dets.append(genlcls2Detector('mono_encoder',run))

    dets.append(ttlcls2Detector('atmopal',run, saveTraces=True))
    #check a RIX scan to figure out controlDetector.
    return dets

def uedDetectors(run, beamCodes=[[162],[91]]):
    dets=[]
    dets.append(scanDetector('scan', run))
    dets.append(genlcls2Detector('timing',run))
    dets.append(lcls2_lightStatus(beamCodes,run))
    dets.append(lcls2_epicsDetector(PVlist=['MOTR_AS01_MC05_CH1'],run=run))

    return dets

def sxrDetectors(beamCodes=[[162],[91]]):
    dets=[]
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(epicsDetector(PVlist=['GATT:FEE1:310:P_DES','GATT:FEE1:310:R_ACT','GMD_ACQ_RAW']))
    dets.append(ttDetector())
    dets.append(GMDDetector())
    return dets

def xppDetectors(beamCodes=[[-137],[91]], env=None):
    dets=[]
    dets.append(lightStatus(codes=beamCodes))
    dets.append(epicsDetector(PVlist=['att_T', 'att_T3rd', 'slit_s1_hw', 'slit_s1_vw', 'slit_s2_hw', 'slit_s2_vw', 'slit_s3_hw', 'slit_s3_vw', 'slit_s4_hw', 'slit_s4_vw', 'lxt_vitara', 'lxt', 'lxt_ttc', 'lxe', 'ccm_E', 'lom_E', 'lom_EC', 'gon_v', 'gon_h', 'gon_r', 'gon_x', 'gon_y', 'gon_z', 'gon_roll', 'gon_pitch', 'gon_kappa_eta', 'gon_kappa_kappa', 'gon_kappa_phi', 'gon_kappa_samx','gon_kappa_samy', 'gon_kappa_samz', 'gon_sam_phi', 'gon_sam_z', 'robot_x', 'robot_y', 'robot_z', 'robot_rx', 'robot_ry', 'robot_rz', 'robot_azi', 'robot_ele', 'robot_rad', 'las_comp_wp', 'las_opa_wp', 'las_drift_correction']))
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
    dets.append(xtcavDetector('xtcav','xtcav'))
    try:
        dets.append(adcDetector('adc','adc'))
    except:
        print('did not find slow Adc detector')
    dets.append(ttDetector(env=env))
    dets.append(bmmonDetector('HX2-SB1-BMMON','ipm_hx2'))
    dets.append(bmmonDetector('XPP-SB2-BMMON','ipm2'))
    dets.append(bmmonDetector('XPP-SB3-BMMON','ipm3'))
    dets.append(feeBldDetector('FEE-SPEC0','feeBld'))
    dets.append(encoderDetector('XPP-USB-ENCODER-02','enc'))
    try:
        dets.append(encoderDetector('usbencoder','enc'))
    except:
        try:
            print('did not find encoder detector with alias, look for DAQ name')
            dets.append(encoderDetector('XppEndstation.0:USDUSB.0','enc'))
        except:
            print('did not add encoder detector')
            pass
    dets.append(encoderDetector('XPP-USB-ENCODER-01','lom_enc'))
    dets.append(l3tDetector())
    dets.append(damageDetector())
    setParameter(dets, dets, detName='damage')
    return dets

def xcsDetectors(beamCodes=[[-137],[89]], env=None):
    dets=[]
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(ipmDetector('XCS-IPM-02','ipm2'))
    dets.append(ipmDetector('XCS-IPM-03','ipm3'))
    dets.append(ipmDetector('XCS-IPM-04','ipm4'))
    dets.append(ipmDetector('XCS-IPM-05','ipm5'))
    dets.append(ipmDetector('XCS-DIO-03','dio3'))
    dets.append(ipmDetector('XCS-DIO-04','dio4'))
    dets.append(ipmDetector('XCS-DIO-05','dio5'))
    dets.append(ipmDetector('XCS-IPM-mono','diodeMono'))
    dets.append(ipmDetector('XCS-IPM-gon','diodeGon'))
    dets.append(bmmonDetector('HX2-SB1-BMMON','ipm_hx2'))
    dets.append(bmmonDetector('XCS-SND-DIO','snd_dio',savePos=False))
    dets.append(bmmonDetector('XCS-SB1-BMMON','ipm4'))
    dets.append(bmmonDetector('XCS-SB2-BMMON','ipm5'))
    dets.append(bmmonDetector('HFX-DG2-BMMON','ipm2'))
    dets.append(aiDetector('XCS-AIN-01','ai'))
    dets.append(epicsDetector(PVlist=['att_T', 'att_T3rd', 'ccm_E', 'lom_E', 'diff_phis', 'diff_th', 'diff_tth', 'diff_xs', 'diff_ys', 'diff_zs', 'diff_x', 'diff_y', 'diff_chis','diff_dety','ladm_theta','lam_z','lam_x1','lam_x2','lam_y1','lam_y2','lam_det_y','lam_det_x','las_comp_wp', 'las_opa_wp', 'las_drift_correction', 'lxt_vitara', 'lxt', 'lxt_ttc', 'lxe' ]))
    dets.append(ttDetector(env=env))
    dets.append(feeBldDetector('FEE-SPEC0','feeBld'))
    dets.append(xtcavDetector('xtcav','xtcav'))
    try:
        dets.append(encoderDetector('XRT-USB-ENCODER-01','xrt_enc'))
        dets.append(encoderDetector('XCS-USB-ENCODER-01','enc'))
    except:
        print('did not add encoder detector')
        pass
    dets.append(l3tDetector())
    dets.append(damageDetector())
    setParameter(dets, dets, detName='damage')
    return dets

def mfxDetectors(beamCodes=[[-137],[204]], env=None):
    dets=[]
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(feeBldDetector('FEE-SPEC0','feeBld'))
    dets.append(aiDetector('MFX-AIN-01','ai')) 
    dets.append(ttDetector(env=env))
    dets.append(bmmonDetector('MFX-DG1-BMMON','ipm_dg1'))
    dets.append(bmmonDetector('MFX-DG2-BMMON','ipm_dg2'))
    dets.append(damageDetector())
    dets.append(xtcavDetector('xtcav','xtcav'))
    try:
        dets.append(encoderDetector('XRT-USB-ENCODER-01','xrt_enc'))
        dets.append(encoderDetector('MFX-USB-ENCODER-01','enc'))
    except:
        pass
    dets.append(epicsDetector(PVlist=['atten_trans1','atten_trans3',
                                      'fee_Attenuator_transmission',
                                      'lens_energy',
                                      'BeamMonitor_target','Dg1Ipm_target',
                                      'lxt', 'txt']))
    return dets

def cxiDetectors(beamCodes=[[-137],[184]], env=None):
    dets=[]
    #dets.append(lightStatus(codes=beamCodes, evrName='evr0'))
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(feeBldDetector('FEE-SPEC0','feeBld'))
    dets.append(epicsDetector(PVlist=['dia_trans1', 'dia_trans3', 'atc_trans1', 'atc_trans3', 'laser_time_target', 'lxt', 'las_uv_waveplate', 'las_evo_waveplate', 'las_drift_correction']))
    dets.append(bmmonDetector('HX2-SB1-BMMON','ipm_hx2'))
    dets.append(bmmonDetector('HFX-DG2-BMMON','ipm_hfx_dg2'))
    dets.append(bmmonDetector('CXI-DG2-BMMON','ipm_dg2'))
    dets.append(bmmonDetector('CXI-DG3-BMMON','ipm_dg3'))
    dets.append(bmmonDetector('CXI-USR-DIO','usr_dio',savePos=False))
    dets.append(ttDetector(env=env))
    dets.append(xtcavDetector('xtcav','xtcav'))
    try:
        dets.append(impDetector('Sc1Imp'))
    except:
        pass
    dets.append(damageDetector())
    try:
        dets.append(encoderDetector('XRT-USB-ENCODER-01','xrt_enc'))
    except:
        pass
    try:
        dets.append(encoderDetector('CXI-USB-ENCODER-01','enc'))
    except:
        pass
    try:
        dets.append(encoderDetector('CxiEndstation.0:USDUSB.0','KbEncoder'))
    except:
        pass
    return dets

def mecDetectors(beamCodes=[[-137],[-182]], env=None):
    dets=[]
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(ipmDetector('MEC-XT2-PIM-02','pim2')) 
    dets.append(ipmDetector('MEC-XT2-PIM-03','pim3')) 
    dets.append(ipmDetector('MEC-TCTR-DI-01','pips-diode-air'))
    dets.append(ipmDetector('MEC-LAS-EM-01','dio_las'))
    dets.append(ipmDetector('MEC-XT2-IPM-02','ipm2'))
    dets.append(ipmDetector('MEC-XT2-IPM-03','ipm3'))
    dets.append(bmmonDetector('MEC-XT2-BMMON-02','xt2_ipm2'))
    dets.append(bmmonDetector('MEC-XT2-BMMON-03','xt2_ipm3'))
    dets.append(ttDetector(env=env))
    dets.append(aiDetector('MEC-AIN-01','ai')) 
    dets.append(feeBldDetector('FEE-SPEC0','feeBld'))
    dets.append(xtcavDetector('xtcav','xtcav'))
    dets.append(damageDetector())
    setParameter(dets, dets, detName='damage') 
    dets.append(epicsDetector(PVlist=['belens_z',
                                      'MEC:NOTE:LAS:NSDELAY',
                                      'MEC:NS1:MMS:01.RBV',
                                      'MEC:NS1:MMS:02.RBV',
                                      'MEC:LAS:MMN:30.RBV',
                                      'MEC:LAS:MMN:29.RBV',
                                      'MEC:PPL:MMS:09.RBV',
                                      'MEC:USR:MMS:17.RBV',
                                      'MEC:LAS:DDG:05:aDelaySI',
                                      'MEC:LAS:DDG:05:eDelaySI',
                                      'SIOC:SYS0:LM00:AO627',
                                      'MEC:TC1:GCC:01:PMON',
                                      'MEC:TC1:GPI:01:PMON',
                                      'MEC:TC1:VGC:01:OPN_OK',
                                      'MEC:TC1:VGC:02:OPN_OK']))
    dets.append(damageDetector())
    try:
        dets.append(encoderDetector('XRT-USB-ENCODER-01','xrt_enc'))
    except:
        pass
    return dets

def diaDetectors():
    dets=[]
    dets.append(feeBldDetector('FEE-SPEC0','feeBld'))
    return dets

def detDetectors():
    return []

def detData(detList, evt):
    #mpiDataSource Issue? cpo, tjlane
    #if one of the keys here contains an empty dict (say: no user PV that are in the data)
    #then it will not be saved. Also any keys _after_ that dir will be lost!
    #test: misspell userPV in producer file and see that the scan directory disappears....
    data={}
    for det in detList:
        try:
            data[det.name] = det.data(evt)
        except:
            #print('could not get data in this event for detector ',det.name)
            pass
    return data

def addArchiverData(det, data, ts):
    try:
        arch = EpicsArchive()
    except Exception:
        print('Failed to create the EPICS archiver')
        return data
    for (i, pv) in enumerate(det.missingPV):
        try:
            (t, v) = arch.get_points(pv, start=ts, end=ts+30, two_lists=True, raw=True)
        except Exception:
            continue # Not in the archiver!
        #
        # OK, t[0] < ts, but t[1] > ts, *if* it exists.
        # So, let's take t[1] if it exists, and t[0] if it's "recent-ish".
        #
        if len(t) >= 2:
            t = t[1]
            v = v[1]
        elif len(t) == 1 and abs(t[0]-ts) < 30:
            t = t[0]
            v = v[0]
        else:
            continue
        al = det.missing[i]
        det.addPV(al, pv)
        data['epicsOnce'][al] = v
    return data

def detOnceData(det, evt, noarch):
    data = detData([det], evt)
    if not noarch:
        # If we've missed something, look for it in the archiver.
        # Look at the timestamp of evt and convert to linux epoch.
        ts = evt.get(psana.EventId).time()[0]
        data = addArchiverData(det, data, ts)
    return data

def lcls2_detOnceData(det, data, ts, noarch):
    if not noarch:
        # If we've missed something, look for it in the archiver.
        # ts is our timestamp.
        data = addArchiverData(det, data, ts)
    return data

def setParameter(detList, Params, detName='tt'):
    for det in detList:
        if det.name==detName:
            det.setPars(Params)
###
#userData functions
###
def getUserData(det):
    """return dictionary with event-based user data from input detector, apply mask here."""
    det_dict= {}
    try:
        userData_keys = [ key for key in det.evt.__dict__.keys() if key.find('_write_')>=0 ]
        #print('DEBUG SmallDataUtils: getting keys: ',userData_keys)
        for key in userData_keys:
            #print('DEBUG SmallDataUtils: getting data for key: ',key)
            if isinstance(getattr(det.evt, key),dict):
                for kkey in getattr(det.evt,key).keys():
                    #print('DEBUG SmallDataUtils: getting data for key, kkey: ',key, kkey)
                    fieldName = ('%s_%s'%(key,kkey)).replace('_write_','')
                    #print 'fieldname ',fieldName
                    if isinstance(getattr(det.evt,key)[kkey], tuple):
                        det_dict[fieldName] = np.array(getattr(det.evt,key)[kkey])
                    else:
                        det_dict[fieldName] = getattr(det.evt,key)[kkey]
                    #print 'local det_dict ',det_dict
            elif isinstance(det.evt.__dict__[key], np.ndarray):
                if isinstance(det.evt.__dict__[key], np.ma.masked_array):
                    data = det.evt.__dict__[key].data
                    if np.issubdtype(data.dtype, np.integer):
                        data[det.evt.__dict__[key].mask]=0
                    else:
                        data[det.evt.__dict__[key].mask]=np.nan
                    #print(eventNr,' is ndarray ',key)
                    det_dict[key.replace('_write_','')] = data
                else:
                    det_dict[key.replace('_write_','')] = det.evt.__dict__[key]
            else:
                det_dict[key.replace('_write_','')] = np.array(det.evt.__dict__[key])
    except:
        pass
    newDict = {}
    for key in det_dict:
        if key.find('ragged')>=0:
            newDict['ragged_%s'%(key.replace('ragged_',''))] = det_dict[key]
        elif key.find('var')>=0:
            fname = key.split('_')[-1]
            dname = key.replace(f'_{fname}','').replace('_var','')
            if 'var_%s'%dname in newDict:
                newDict[f'var_{dname}'][fname] = det_dict[key]
            else:
                newDict[f'var_{dname}'] = {fname : det_dict[key]}


        else:
            newDict['%s'%(key)] = det_dict[key]
    det_dict = newDict
    #print 'det_dict ',det_dict.keys()
    return det_dict

def getUserEnvData(det):
    """return environment data for a detectors as a dictionary (temperatures,....)"""
    env_dict= {}
    try:
        envData_keys = [ key for key in det.evt.__dict__.keys() if key.find('env_')>=0 ]
        for key in envData_keys:
            if isinstance(det.evt.__dict__[key], np.ndarray):
                if isinstance(det.evt.__dict__[key], np.ma.masked_array):
                    data = det.evt.__dict__[key].data
                    data[det.evt.__dict__[key].mask]=np.nan
                    env_dict[key.replace('env_','')] = data
                else:
                    env_dict[key.replace('env_','')] = det.evt.__dict__[key]
            else:
                env_dict[key.replace('env_','')] = np.array(det.evt.__dict__[key])
    except:
        pass
    return env_dict

def getCfgOutput(det):
    """return the configuration data for the user data (parameters for feature extraction) as a dict"""
    cfgDict={}
    #baseName=det._name+'_'
    baseName=''
    for key in det.__dict__.keys():
        if key == 'mask_ROI' or key == 'mask_ROI_shape' or key == 'bankMasks' or key.find('__')==0:
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
        if 'ped' in ROI.__dict__.keys() and ROI.ped is not None:
            cfgDict[baseName+ROI.name+'_ped'] = np.array(ROI.ped)
        if 'mask' in ROI.__dict__.keys() and ROI.mask is not None:
            cfgDict[baseName+ROI.name+'_mask'] = np.array(ROI.mask)

    for drops in det.getDroplets():
        for aduHist in drops.aduHists:
            cfgDict[baseName+aduHist.name+'_bins'] = np.array(aduHist.bins)
            if aduHist.ROI is not None:
                cfgDict[baseName+aduHist.name+'_ROI'] = np.array(aduHist.ROI)
        for dropSave in drops.dropletSaves:
            cfgDict[baseName+dropSave.name+'_thresADU'] = np.array(dropSave.thresADU)
            if dropSave.ROI is not None:
                cfgDict[baseName+dropSave.name+'_ROI'] = np.array(dropSave.ROI)
    return cfgDict
