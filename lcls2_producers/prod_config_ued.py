import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding 
# list empty.
detectors = ['andor']
#detectors = ['epixquad', 'andor']
# integrating_detectors = []
integrating_detectors = []
# Comment: the first integrating detector will set the sub-sampling of all
# integrating detectors.


def getROIs(run):
    ret_dict = {}

    full_roi = {
        'thresADU' : None,
        'writeArea' : True,
        'calcPars' : False,
        'ROI' : None
    }

    if run>0:
        # Save the full ANDOR and QUAD.
        ret_dict['epixquad'] = [full_roi]
        ret_dict['andor'] = [full_roi]
    return ret_dict


##########################################################
# run independent parameters
##########################################################
# These lists are either PV names, aliases, or tuples with both.
#epicsPV = ['las_fs14_controller_time']
#epicsOncePV = ['m0c0_vset', ('TMO:PRO2:MPOD:01:M2:C3:VoltageMeasure', 'MyAlias'),
#               'IM4K4:PPM:SPM:VOLT_RBV', "FOO:BAR:BAZ", ("X:Y:Z", "MCBTest"), "A:B:C"]
epicsPV = [
]
epicsOncePV = []


##########################################################
# psplot config
##########################################################

import psplot
