import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding
# list empty.
detectors = ["epix100_0", "epix100_1", "alvium"]  # , 'qadc1']
# detectors = []
integrating_detectors = []


def getROIs(run):
    ret_dict = {}

    jungfrau_roi = {"thresADU": None, "writeArea": True, "calcPars": False, "ROI": None}
    alv_roi = {"writeArea": True, "calcPars": False, "ROI": None, "thresADU": None}

    epix100_0_roi = {
        "thresADU": None,
        "writeArea": True,
        "calcPars": False,
        "ROI": [[0, 300], [0, 710]],
    }

    epix100_1_roi = {
        "thresADU": None,
        "writeArea": True,
        "calcPars": False,
        "ROI": [[0, 300], [0, 710]],
    }

    if run > 0:
        # ret_dict['jungfrau'] = [jungfrau_roi]
        ret_dict["epix100_0"] = [epix100_0_roi]
        ret_dict["epix100_1"] = [epix100_1_roi]
        ret_dict["alvium"] = [alv_roi]
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

        # ret_dict["epix100"] = d2p_dict

    return ret_dict


def get_droplet(run):
    ret_dict = {}
    if run > 0:
        dfunc_dict = {}
        dfunc_dict["threshold"] = 110.0 / 32
        # dfunc_dict['thresholdLow'] = 10./8
        dfunc_dict["thresADU"] = 20.0  # 4.
        dfunc_dict["useRms"] = False
        ret_dict["epix100"] = dfunc_dict
    return ret_dict


def get_azav(run):
    ret_dict = {}
    if run > 0:
        az_dict = {"eBeam": 9.8}  # Photon energy in keV
        az_dict["center"] = [-1107, 918]  # Position in um
        az_dict["dis_to_sam"] = 92.51  # try to overwrite dis_to_sam
        az_dict["tx"] = 0  # Tilt in x in degrees
        az_dict["ty"] = 0  # Tilt in y in degrees
        az_dict["phiBins"] = 11  # Number of phi bins for azint
        az_dict["qbin"] = 0.025  # Bin width in q for az int
        az_dict["userMask"] = None
        # ret_dict["jungfrau"] = az_dict

    return ret_dict


def get_sum_algos(run):
    """Specify detector sum algorithms for specific detectors or `all` detectors.

    Provide a list of algorithms for each detector `detname`, or provide a list
    under the `all` key which will be applied to each detector.

    Possible algorithms:
    - calib: Sum calib.
    - calib_square: Square of image.
    - calib_thresADU1 (or 5, 10 etc): Apply threshold using the number after `ADU`
    - calib_max: Maximum projection instead of sum.
    """
    ret_dict = {}
    if run > 0:
        ret_dict["all"] = [
            "calib",
            "calib_square",
            "calib_thresADU1",
        ]
        ret_dict["jungfrau"] = ["calib_thresADU5", "calib_max"]
        # ret_dict["epix100_0"] = []
        # ret_dict["epix100_1"] = []
    return ret_dict


##########################################################
# run independent parameters
##########################################################
# Will be taken from the archiver
# This list is PV names or tuples (PV name, alias)
epicsPV = []

# Will be taken from the XTC (epicsArch file)
# NOTE: Data appears to be shot-to-shot in the HDF5 but is NOT! This for convenience only
# These lists are either PV names, aliases, or tuples (PV name, alias)
# epicsArchFilePV = ['las_fs14_controller_time']
# epicsOncePV = ['m0c0_vset', ('TMO:PRO2:MPOD:01:M2:C3:VoltageMeasure', 'MyAlias'),
#               'IM4K4:PPM:SPM:VOLT_RBV', "FOO:BAR:BAZ", ("X:Y:Z", "MCBTest"), "A:B:C"]
epicsArchFilePV = []
epicsOncePV = []


##########################################################
# psplot config
##########################################################

import psplot
