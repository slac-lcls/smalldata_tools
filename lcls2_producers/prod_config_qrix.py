import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding
# list empty.
detectors = ["hsd", "rix_fim0", "rix_fim1", "qrix_w8", "q_piranha", "q_atmopal"]


def get_intg(run):
    """
    Returns
    -------
    intg_main (str):  This detector ill be passed to the psana datasource. It should be the SLOWEST of
                all integrating detectors in the data
    intg_addl (list of str): The detectors in this list will be analyzed as integrating detectors. It is
                             important that the readout frequency of these detectors is commensurate and
                             in-phase with intg_main.
    """
    run = int(run)
    intg_main = None
    intg_addl = []

    # if run > 0:
    #     intg_main = "archon"
    #     intg_addl = []

    return intg_main, intg_addl


slow_detectors = []  # NOT IMPLEMENTED


def getROIs(run):
    ret_dict = {}

    full_roi = {"thresADU": None, "writeArea": True, "calcPars": False, "ROI": None}

    if run > 0:
        # Save the full DET.
        ret_dict["q_piranha"] = [full_roi]
        ret_dict["q_atmopal"] = [full_roi]
        ret_dict["archon"] = [full_roi]

        # and the FIM waveforms
        ret_dict["rix_fim0"] = [full_roi]
        ret_dict["rix_fim1"] = [full_roi]
        ret_dict["qrix_w8"] = [full_roi]

        # hsd channels used by qRIXS
        hsd_dict = {}
        hsd_dict["hsd_0"] = [0, -1]
        hsd_dict["hsd_3"] = [0, -1]
        ret_dict["hsd"] = hsd_dict

    return ret_dict


def get_droplet2photon(run):
    ret_dict = {}

    if run < 0:
        # ##### Archon setup #####
        d2p_dict = {}
        # droplet args
        d2p_dict["droplet"] = {
            "threshold": 65,
            # 'thresholdLow': 0,
            "thresADU": 0,  # discard droplet whose total ADU is below this value
            "useRms": False,
        }

        # droplet2Photons args
        d2p_dict["d2p"] = {
            "aduspphot": 400,
            # 'roi_mask': np.load('path_to_mask.npy'),
            "cputime": True,
        }
        d2p_dict["nData"] = None
        d2p_dict["get_photon_img"] = False

        ret_dict["archon"] = d2p_dict

    return ret_dict


def get_wf_integrate(run):
    ret_dict = {}

    if run > 0:
        ret_dict["rix_fim0"] = {
            "sig_roi": slice(104, 112),
            "bkg_roi": slice(0, 80),
            "negative_signal": True,
        }

        ret_dict["rix_fim1"] = {
            "sig_roi": slice(120, 127),
            "bkg_roi": slice(0, 80),
            "negative_signal": True,
        }

        ret_dict["crix_w8"] = {
            "sig_roi": slice(76, 95),
            "bkg_roi": slice(0, 50),
            "negative_signal": True,
        }
    return ret_dict


def get_wf_svd(run):
    ret_dict = {}

    if run > 0:
        det_kwargs = {}
        w8s = ["rix_fim0", "rix_fim1", "crix_w8"]
        for w8 in w8s:
            # one basis file per channel
            run_basis = 146
            path = "/sdf/data/lcls/ds/rix/rixx1011723/hdf5/smalldata/svd_basis"
            basis_files = [
                f"{path}/wave_basis_{w8}_ch{ch}_r{str(run_basis).zfill(4)}.h5"
                for ch in range(8)
            ]
            det_kwargs[w8] = []
            for basis_file in basis_files:
                kwargs = {
                    "n_components": 3,
                    "mode": "both",
                    "return_reconstructed": True,
                    "basis_file": basis_file,
                }
                det_kwargs[w8].append(kwargs)
    return det_kwargs


##########################################################
# run independent parameters
##########################################################
# These lists are either PV names, aliases, or tuples with both.
# epicsPV = ['las_fs14_controller_time']
# epicsOncePV = ['m0c0_vset', ('TMO:PRO2:MPOD:01:M2:C3:VoltageMeasure', 'MyAlias'),
#               'IM4K4:PPM:SPM:VOLT_RBV', "FOO:BAR:BAZ", ("X:Y:Z", "MCBTest"), "A:B:C"]
# epicsOncePV = [('GDET:FEE1:241:ENRC', "MyTest"), 'GDET:FEE1:242:ENRC', "FOO:BAR:BAZ"]
epicsPV = []
epicsPV = [
    "EM2K0:XGMD:ETM:01:Reading",
    "EM2K0:XGMD:ETM:02:Reading",
    "EM2K0:XGMD:HPS:milliJoulesPerPulse",
    "EM2K0:XGMD:HPS:AvgPulseIntensity",
]
epicsOncePV = []


##########################################################
# psplot config
##########################################################

import psplot


def get_psplot_configs(run):
    configs = {}

    if run > 0:
        andor_vls = {
            "callback": psplot.SpectrumScan,
            "data_fields": {
                "data": "andor_vls/unaligned_hitfinder",
                "norm": "andor_vls/unaligned_norm_det_rix_fim1_sum_wfintegrate",
                "count": "andor_vls/unaligned_norm_count",
                "scan": "andor_vls/unaligned_norm_scan_sum_mono_ev",
            },
            "bins": np.linspace(520, 540, 30),
            "spectrum_range": (900, 1050),
            "lineouts_idx": (50, 70, 100),
        }
        configs["andor_vls"] = andor_vls
    return configs
