import numpy as np
from smalldata_tools.lcls1.default_detectors import *
from datetime import datetime


def defaultDetectors(hutch, run=None, env=None):
    if hutch.lower() == "amo":
        dets = amoDetectors()
    elif hutch.lower() == "sxr":
        dets = sxrDetectors()
    elif hutch.lower() == "xpp":
        dets = xppDetectors(env=env)
    elif hutch.lower() == "xcs":
        dets = xcsDetectors(env=env)
    elif hutch.lower() == "mfx":
        dets = mfxDetectors(env=env)
    elif hutch.lower() == "cxi":
        dets = cxiDetectors(env=env)
    elif hutch.lower() == "mec":
        dets = mecDetectors(env=env)
    elif hutch.lower() == "det":
        dets = detDetectors()
    elif hutch.lower() == "dia":
        dets = diaDetectors()
    elif hutch.lower() == "tmo":
        dets = tmoDetectors(run)
    elif hutch.lower() == "rix":
        dets = rixDetectors(run)
    elif hutch.lower() == "ued":
        dets = uedDetectors(run)
    else:
        dets = []
    detsInRun = [det for det in dets if det.in_run()]
    # print('found %d detectors '%len(detsInRun))
    return detsInRun


def amoDetectors():
    return []


def uedDetectors(run, beamCodes=[[162], [91]]):
    dets = []
    dets.append(scanDetector("scan", run))
    dets.append(genlcls2Detector("timing", run))
    dets.append(lcls2_lightStatus(beamCodes, run))
    dets.append(lcls2_epicsDetector(PVlist=["MOTR_AS01_MC05_CH1"], run=run))

    return dets


def sxrDetectors(beamCodes=[[162], [91]]):
    dets = []
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(
        epicsDetector(
            PVlist=["GATT:FEE1:310:P_DES", "GATT:FEE1:310:R_ACT", "GMD_ACQ_RAW"]
        )
    )
    dets.append(ttDetector())
    dets.append(GMDDetector())
    return dets


def xppDetectors(beamCodes=[[-137], [91]], env=None):
    dets = []
    dets.append(lightStatus(codes=beamCodes))
    dets.append(
        epicsDetector(
            PVlist=[
                "att_T",
                "att_T3rd",
                "slit_s1_hw",
                "slit_s1_vw",
                "slit_s2_hw",
                "slit_s2_vw",
                "slit_s3_hw",
                "slit_s3_vw",
                "slit_s4_hw",
                "slit_s4_vw",
                "lxt_vitara",
                "lxt",
                "lxt_ttc",
                "lxe",
                "ccm_E",
                "lom_E",
                "lom_EC",
                "gon_v",
                "gon_h",
                "gon_r",
                "gon_x",
                "gon_y",
                "gon_z",
                "gon_roll",
                "gon_pitch",
                "gon_kappa_eta",
                "gon_kappa_kappa",
                "gon_kappa_phi",
                "gon_kappa_samx",
                "gon_kappa_samy",
                "gon_kappa_samz",
                "gon_sam_phi",
                "gon_sam_z",
                "robot_x",
                "robot_y",
                "robot_z",
                "robot_rx",
                "robot_ry",
                "robot_rz",
                "robot_azi",
                "robot_ele",
                "robot_rad",
                "las_comp_wp",
                "las_opa_wp",
                "las_drift_correction",
            ]
        )
    )
    dets.append(controlDetector())
    dets.append(ipmDetector("NH2-SB1-IPM-01", "ipm1"))
    dets.append(ipmDetector("NH2-SB1-IPM-02", "ipm1c"))
    dets.append(ipmDetector("XppMon_Pim0", "lombpm"))
    dets.append(ipmDetector("XppMon_Pim1", "lomdiode"))
    dets.append(ipmDetector("XppSb2_Ipm", "ipm2"))
    dets.append(ipmDetector("XppSb3_Ipm", "ipm3"))
    dets.append(ipmDetector("XppSb3_Pim", "diode2"))
    dets.append(ipmDetector("XppSb4_Pim", "diode2"))
    dets.append(ipmDetector("XppEnds_Ipm0", "diodeU"))
    dets.append(aiDetector("XPP-AIN-01", "ai"))
    dets.append(xtcavDetector("xtcav", "xtcav"))
    try:
        dets.append(adcDetector("adc", "adc"))
    except:
        print("did not find slow Adc detector")
    dets.append(ttDetector(env=env))
    dets.append(bmmonDetector("HX2-SB1-BMMON", "ipm_hx2"))
    dets.append(bmmonDetector("XPP-SB2-BMMON", "ipm2"))
    dets.append(bmmonDetector("XPP-SB3-BMMON", "ipm3"))
    dets.append(feeBldDetector("FEE-SPEC0", "feeBld"))
    dets.append(encoderDetector("XPP-USB-ENCODER-02", "enc"))
    try:
        dets.append(encoderDetector("usbencoder", "enc"))
    except:
        try:
            print("did not find encoder detector with alias, look for DAQ name")
            dets.append(encoderDetector("XppEndstation.0:USDUSB.0", "enc"))
        except:
            print("did not add encoder detector")
            pass
    dets.append(encoderDetector("XPP-USB-ENCODER-01", "lom_enc"))
    dets.append(l3tDetector())
    dets.append(damageDetector())
    setParameter(dets, dets, detName="damage")
    return dets


def xcsDetectors(beamCodes=[[-137], [89]], env=None):
    dets = []
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(ipmDetector("XCS-IPM-02", "ipm2"))
    dets.append(ipmDetector("XCS-IPM-03", "ipm3"))
    dets.append(ipmDetector("XCS-IPM-04", "ipm4"))
    dets.append(ipmDetector("XCS-IPM-05", "ipm5"))
    dets.append(ipmDetector("XCS-DIO-03", "dio3"))
    dets.append(ipmDetector("XCS-DIO-04", "dio4"))
    dets.append(ipmDetector("XCS-DIO-05", "dio5"))
    dets.append(ipmDetector("XCS-IPM-mono", "diodeMono"))
    dets.append(ipmDetector("XCS-IPM-gon", "diodeGon"))
    dets.append(bmmonDetector("HX2-SB1-BMMON", "ipm_hx2"))
    dets.append(bmmonDetector("XCS-SND-DIO", "snd_dio", savePos=False))
    dets.append(bmmonDetector("XCS-SB1-BMMON", "ipm4"))
    dets.append(bmmonDetector("XCS-SB2-BMMON", "ipm5"))
    dets.append(bmmonDetector("HFX-DG2-BMMON", "ipm2"))
    dets.append(aiDetector("XCS-AIN-01", "ai"))
    dets.append(
        epicsDetector(
            PVlist=[
                "att_T",
                "att_T3rd",
                "ccm_E",
                "lom_E",
                "diff_phis",
                "diff_th",
                "diff_tth",
                "diff_xs",
                "diff_ys",
                "diff_zs",
                "diff_x",
                "diff_y",
                "diff_chis",
                "diff_dety",
                "ladm_theta",
                "lam_z",
                "lam_x1",
                "lam_x2",
                "lam_y1",
                "lam_y2",
                "lam_det_y",
                "lam_det_x",
                "las_comp_wp",
                "las_opa_wp",
                "las_drift_correction",
                "lxt_vitara",
                "lxt",
                "lxt_ttc",
                "lxe",
            ]
        )
    )
    dets.append(ttDetector(env=env))
    dets.append(feeBldDetector("FEE-SPEC0", "feeBld"))
    dets.append(xtcavDetector("xtcav", "xtcav"))
    try:
        dets.append(encoderDetector("XRT-USB-ENCODER-01", "xrt_enc"))
        dets.append(encoderDetector("XCS-USB-ENCODER-01", "enc"))
    except:
        print("did not add encoder detector")
        pass
    dets.append(l3tDetector())
    dets.append(damageDetector())
    setParameter(dets, dets, detName="damage")
    return dets


def mfxDetectors(beamCodes=[[-137], [204]], env=None):
    dets = []
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(feeBldDetector("FEE-SPEC0", "feeBld"))
    dets.append(aiDetector("MFX-AIN-01", "ai"))
    dets.append(ttDetector(env=env))
    dets.append(bmmonDetector("MFX-DG1-BMMON", "ipm_dg1"))
    dets.append(bmmonDetector("MFX-DG2-BMMON", "ipm_dg2"))
    dets.append(damageDetector())
    dets.append(xtcavDetector("xtcav", "xtcav"))
    try:
        dets.append(encoderDetector("XRT-USB-ENCODER-01", "xrt_enc"))
        dets.append(encoderDetector("MFX-USB-ENCODER-01", "enc"))
    except:
        pass
    dets.append(
        epicsDetector(
            PVlist=[
                "atten_trans1",
                "atten_trans3",
                "fee_Attenuator_transmission",
                "lens_energy",
                "lens_z",
                "BeamMonitor_target",
                "Dg1Ipm_target",
                "lxt",
                "lxt_fs_tgt",
                "lxt_fast",
                "lxt_fast_mot",
                "txt",
                "txt_mot",
            ]
        )
    )
    return dets


def cxiDetectors(beamCodes=[[-137], [184]], env=None):
    dets = []
    # dets.append(lightStatus(codes=beamCodes, evrName='evr0'))
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(feeBldDetector("FEE-SPEC0", "feeBld"))
    dets.append(
        epicsDetector(
            PVlist=[
                "dia_trans1",
                "dia_trans3",
                "atc_trans1",
                "atc_trans3",
                "laser_time_target",
                "lxt",
                "las_uv_waveplate",
                "las_evo_waveplate",
                "las_drift_correction",
            ]
        )
    )
    dets.append(bmmonDetector("HX2-SB1-BMMON", "ipm_hx2"))
    dets.append(bmmonDetector("HFX-DG2-BMMON", "ipm_hfx_dg2"))
    dets.append(bmmonDetector("CXI-DG2-BMMON", "ipm_dg2"))
    dets.append(bmmonDetector("CXI-DG3-BMMON", "ipm_dg3"))
    dets.append(bmmonDetector("CXI-USR-DIO", "usr_dio", savePos=False))
    dets.append(ttDetector(env=env))
    dets.append(xtcavDetector("xtcav", "xtcav"))
    try:
        dets.append(impDetector("Sc1Imp"))
    except:
        pass
    dets.append(damageDetector())
    try:
        dets.append(encoderDetector("XRT-USB-ENCODER-01", "xrt_enc"))
    except:
        pass
    try:
        dets.append(encoderDetector("CXI-USB-ENCODER-01", "enc"))
    except:
        pass
    try:
        dets.append(encoderDetector("CxiEndstation.0:USDUSB.0", "KbEncoder"))
    except:
        pass
    return dets


def mecDetectors(beamCodes=[[-137], [-182]], env=None):
    dets = []
    dets.append(lightStatus(codes=beamCodes))
    dets.append(controlDetector())
    dets.append(ipmDetector("MEC-XT2-PIM-02", "pim2"))
    dets.append(ipmDetector("MEC-XT2-PIM-03", "pim3"))
    dets.append(ipmDetector("MEC-TCTR-DI-01", "pips-diode-air"))
    dets.append(ipmDetector("MEC-LAS-EM-01", "dio_las"))
    dets.append(ipmDetector("MEC-XT2-IPM-02", "ipm2"))
    dets.append(ipmDetector("MEC-XT2-IPM-03", "ipm3"))
    dets.append(bmmonDetector("MEC-XT2-BMMON-02", "xt2_ipm2"))
    dets.append(bmmonDetector("MEC-XT2-BMMON-03", "xt2_ipm3"))
    dets.append(ttDetector(env=env))
    dets.append(aiDetector("MEC-AIN-01", "ai"))
    dets.append(feeBldDetector("FEE-SPEC0", "feeBld"))
    dets.append(xtcavDetector("xtcav", "xtcav"))
    dets.append(damageDetector())
    setParameter(dets, dets, detName="damage")
    dets.append(
        epicsDetector(
            PVlist=[
                "belens_z",
                "MEC:NOTE:LAS:NSDELAY",
                "MEC:NS1:MMS:01.RBV",
                "MEC:NS1:MMS:02.RBV",
                "MEC:LAS:MMN:30.RBV",
                "MEC:LAS:MMN:29.RBV",
                "MEC:PPL:MMS:09.RBV",
                "MEC:USR:MMS:17.RBV",
                "MEC:LAS:DDG:05:aDelaySI",
                "MEC:LAS:DDG:05:eDelaySI",
                "SIOC:SYS0:ML00:AO627",
                "MEC:TC1:GCC:01:PMON",
                "MEC:TC1:GPI:01:PMON",
                "MEC:TC1:VGC:01:OPN_OK",
                "MEC:TC1:VGC:02:OPN_OK",
            ]
        )
    )
    dets.append(damageDetector())
    try:
        dets.append(encoderDetector("XRT-USB-ENCODER-01", "xrt_enc"))
    except:
        pass
    return dets


def diaDetectors():
    dets = []
    dets.append(feeBldDetector("FEE-SPEC0", "feeBld"))
    return dets


def detDetectors():
    return []


###
# userData functions
###


def getCfgOutput(det):
    """Return the configuration data for the user data (parameters for feature extraction) as a dict"""
    cfgDict = {}
    # baseName=det._name+'_'
    baseName = ""
    for key in det.__dict__.keys():
        if (
            key == "mask_ROI"
            or key == "mask_ROI_shape"
            or key == "bankMasks"
            or key.find("__") == 0
        ):
            continue
        if isinstance(det[key], list) or isinstance(det[key], np.ndarray):
            if isinstance(det[key], list):
                cfgDict[baseName + key] = np.array(det[key])
            elif isinstance(det[key], np.ndarray):
                cfgDict[baseName + key] = np.array(det[key])
        elif isinstance(det[key], float) or isinstance(det[key], int):
            cfgDict[baseName + key] = np.array([det[key]])
    # now add ROI boundaries (so we can use mask, rms from main det object later....)
    for ROI in det.getROIs():
        cfgDict[baseName + ROI.name + "_bound"] = np.array(ROI.bound)
        if "ped" in ROI.__dict__.keys() and ROI.ped is not None:
            cfgDict[baseName + ROI.name + "_ped"] = np.array(ROI.ped)
        if "mask" in ROI.__dict__.keys() and ROI.mask is not None:
            cfgDict[baseName + ROI.name + "_mask"] = np.array(ROI.mask)

    for drops in det.getDroplets():
        for aduHist in drops.aduHists:
            cfgDict[baseName + aduHist.name + "_bins"] = np.array(aduHist.bins)
            if aduHist.ROI is not None:
                cfgDict[baseName + aduHist.name + "_ROI"] = np.array(aduHist.ROI)
        for dropSave in drops.dropletSaves:
            cfgDict[baseName + dropSave.name + "_thresADU"] = np.array(
                dropSave.thresADU
            )
            if dropSave.ROI is not None:
                cfgDict[baseName + dropSave.name + "_ROI"] = np.array(dropSave.ROI)
    return cfgDict
