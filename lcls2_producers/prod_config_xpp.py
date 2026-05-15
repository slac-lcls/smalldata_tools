import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding
# list empty.
detectors = [ 'alvium_1','jungfrau1M']

def getROIs(run):

    # full_roi = {"thresADU": None, "writeArea": True, "calcPars": False, "ROI": None}
    # sum_thresh = {"thresADU": [8.0,9.0], "writeArea": True, "ROI": [ [[0,1], [0, 300], [200,700]],
                             # [[1,2], [0, 200], [200,700]] ] }
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if (run>= 0):
        roi_dict = {}
        roi_dict['ROI'] = None # can define more than one ROI
        roi_dict['writeArea'] = True
        ret_dict['alvium_1'] = [roi_dict]
    return ret_dict


##########################################################
# run independent parameters
##########################################################
# These lists are either PV names, aliases, or tuples with both.
# epicsPV = ['las_fs14_controller_time']
# epicsOncePV = ['m0c0_vset', ('TMO:PRO2:MPOD:01:M2:C3:VoltageMeasure', 'MyAlias'),
#               'IM4K4:PPM:SPM:VOLT_RBV', "FOO:BAR:BAZ", ("X:Y:Z", "MCBTest"), "A:B:C"]
# epicsOncePV = [('GDET:FEE1:241:ENRC', "MyTest"), 'GDET:FEE1:242:ENRC', "FOO:BAR:BAZ"]
epicsPV = []
epicsOncePV = []


##########################################################
# psplot config
##########################################################
