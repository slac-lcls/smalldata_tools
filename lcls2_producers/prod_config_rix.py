import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding 
# list empty.
detectors = ['hsd', 'rix_fim0', 'rix_fim1', 'crix_w8', 'c_piranha']
# integrating_detectors = []
integrating_detectors = ['andor_dir', 'andor_vls', 'andor_norm', 'archon']
# Comment: the first integrating detector will set the sub-sampling of all
# integrating detectors.
slow_detectors = []  # NOT IMPLEMENTED


def getROIs(run):
    ret_dict = {}

    full_roi = {
        'thresADU' : None,
        'writeArea' : True,
        'calcPars' : False,
        'ROI' : None
    }

    if run>0:
        # Save the full ANDOR.
        ret_dict['andor_dir'] = [full_roi]
        ret_dict['andor_vls'] = [full_roi]
        ret_dict['andor_norm'] = [full_roi]
        ret_dict['c_piranha'] = [full_roi]
        
        # and the FIM waveforms
        ret_dict['rix_fim0'] = [full_roi]
        ret_dict['rix_fim1'] = [full_roi]
        # ret_dict['crix_w8'] = [full_roi]
        # ret_dict['qrix_w8'] = [full_roi]

        # hsd
        # hsd_dict = {}
        # hsd_dict['hsd_0'] = [0,-1]
        # hsd_dict['hsd_1'] = [0,-1]
        # hsd_dict['hsd_2'] = [0,-1]
        # hsd_dict['hsd_3'] = [0,-1]
        # ret_dict['hsd'] = hsd_dict

    return ret_dict


def getDetImages(run):
    ret_dict = {}

    if run>0:
        ret_dict['archon'] = {}
    return ret_dict


def get_wf_hitfinder(run):
    ret_dict = {}

    if run>0:
        andor_vls = {
            'threshold' : 4,
            'treshold_max' : 3500,
            'use_rms' : True,
            'bkg_roi' : (500, 800)
        }
        ret_dict['andor_vls'] = andor_vls
    return ret_dict


def get_wf_integrate(run):
    ret_dict = {}

    if run>0:
        ret_dict['rix_fim0'] = {
            "sig_roi" : slice(102,108),
            "bkg_roi" : slice(0,80),
            "negative_signal" : True
        }

        ret_dict['rix_fim1'] = {
            "sig_roi" : slice(116,122),
            "bkg_roi" : slice(0,80),
            "negative_signal" : True
        }

        ret_dict['crix_w8'] = {
            "sig_roi" : slice(69,76),
            "bkg_roi" : slice(0,50),
            "negative_signal" : True
        }
    return ret_dict


##########################################################
# run independent parameters
##########################################################
# These lists are either PV names, aliases, or tuples with both.
#epicsPV = ['las_fs14_controller_time']
#epicsOncePV = ['m0c0_vset', ('TMO:PRO2:MPOD:01:M2:C3:VoltageMeasure', 'MyAlias'),
#               'IM4K4:PPM:SPM:VOLT_RBV', "FOO:BAR:BAZ", ("X:Y:Z", "MCBTest"), "A:B:C"]
#epicsOncePV = [('GDET:FEE1:241:ENRC', "MyTest"), 'GDET:FEE1:242:ENRC', "FOO:BAR:BAZ"]
epicsPV = []
epicsPV = [
    "EM2K0:XGMD:ETM:01:Reading",
    "EM2K0:XGMD:ETM:02:Reading",
    "EM2K0:XGMD:HPS:milliJoulesPerPulse",
    "EM2K0:XGMD:HPS:AvgPulseIntensityy",
]
epicsOncePV = []


##########################################################
# psplot config
##########################################################

import psplot

def get_psplot_configs(run):
    configs = {}

    if run>0:
        andor_vls = {
            'callback' : psplot.SpectrumScan,
            'data_fields' : {
                'data' : 'andor_vls/unaligned_hitfinder',
                'norm' : 'andor_vls/unaligned_norm_det_rix_fim1_sum_wfintegrate',
                'count' : 'andor_vls/unaligned_norm_count',
                'scan' : 'andor_vls/unaligned_norm_scan_sum_mono_ev',
            },
            'bins' : np.linspace(520, 540, 30),
            'spectrum_range' : (900, 1050),
            'lineouts_idx' : (50, 70, 100)
        }
        configs['andor_vls'] = andor_vls
    return configs