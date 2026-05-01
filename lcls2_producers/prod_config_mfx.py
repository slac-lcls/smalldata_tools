detectors = ["jungfrau", "epix100_0", "qadc"]


def getROIs(run):
    ret_dict = {}

    if run > 0:

        # Create list of dicts for epix100_0
        roi_dicts = []
        # Create a dict for this set
        roi_dict = {}
        roi_dict["ROIs"] = [[[530, 680], [270, 350]]]

        roi_dict["writeArea"] = True

        roi_dict["thresADU"] = 3.5

        roi_dict["calcPars"] = True

        roi_dicts.append(roi_dict)

        # Add list of dicts for epix100_0 to total dictionary
        ret_dict["epix100_0"] = roi_dicts

    return ret_dict


def get_droplet(run):
    if isinstance(run, str):
        run = int(run)
    ret_dict = {}
    if run > 0:
        droplet_dict = {}

        droplet_dict["threshold"] = 3.4375

        droplet_dict["thresADU"] = 3.5

        droplet_dict["useRms"] = False

        droplet_dict["name"] = "droplet"

        ret_dict["epix100_0"] = droplet_dict

    return ret_dict


def get_azav(run):
    if isinstance(run, str):
        run = int(run)
    ret_dict = {}
    if run > 0:
        az_dict = {}

        az_dict["eBeam"] = 7.1

        az_dict["center"] = [0, 0]

        az_dict["dis_to_sam"] = 150

        az_dict["phiBins"] = 11

        az_dict["qbin"] = 0.02

        az_dict["tx"] = 0

        az_dict["ty"] = 0

        ret_dict["jungfrau"] = az_dict

    return ret_dict


def get_sum_algos(run):
    ret_dict = {}
    if run > 0:

        ret_dict["epix100_0"] = ["calib", "calib_max"]

        ret_dict["jungfrau"] = ["calib", "calib_max"]

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
