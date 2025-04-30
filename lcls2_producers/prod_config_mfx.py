import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding
# list empty.
# detectors = ['jungfrau','epix100']
detectors = []
integrating_detectors = []


# NOTE TO MFX: Do not modify get_intg. This does not apply to MFX but is required for
# smalldata
# Comment: the first integrating detector will set the sub-sampling of all
# integrating detectors.
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
    if run > 0:
        intg_main = ""
        intg_addl = []
    return intg_main, intg_addl


def getROIs(run):
    ret_dict = {}

    jungfrau_roi = {"thresADU": None, "writeArea": True, "calcPars": False, "ROI": None}
    epix100_roi = {"thresADU": None, "writeArea": True, "calcPars": False, "ROI": None}

    if run > 0:
        # ret_dict['jungfrau'] = [jungfrau_roi]
        # ret_dict['epix100'] = [epix100_roi]
        ...
    return ret_dict


def get_droplet2photon(run):
    ret_dict = {}

    if run > 0:
        # ##### epixquad (1010667 setup) #####
        d2p_dict = {}
        # droplet args
        d2p_dict["droplet"] = {
            #'threshold': 12.5,
            "threshold": 50,
            "thresholdLow": 12.5,
            #'threshold': 70,
            #'thresholdLow': 25,
            "thresADU": 70,  # discard droplet whose total ADU is below this value
            "useRms": False,
        }

        # droplet2Photons args
        d2p_dict["d2p"] = {
            "aduspphot": 20,
            # 'roi_mask': np.load('path_to_mask.npy'),
            "cputime": True,
        }
        d2p_dict["nData"] = None
        d2p_dict["get_photon_img"] = False

        ret_dict["epix100"] = d2p_dict

    return ret_dict


##########################################################
# run independent parameters
##########################################################
# These lists are either PV names, aliases, or tuples with both.
# epicsPV = ['las_fs14_controller_time']
# epicsOncePV = ['m0c0_vset', ('TMO:PRO2:MPOD:01:M2:C3:VoltageMeasure', 'MyAlias'),
#               'IM4K4:PPM:SPM:VOLT_RBV', "FOO:BAR:BAZ", ("X:Y:Z", "MCBTest"), "A:B:C"]
epicsPV = []
epicsOncePV = []


##########################################################
# psplot config
##########################################################

import psplot
