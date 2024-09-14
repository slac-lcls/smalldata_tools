import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding 
# list empty.
detectors = ['hsd', 'rix_fim0', 'rix_fim1', 'crix_w8', 'c_piranha']
integrating_detectors = []
# integrating_detectors = ['andor_dir', 'andor_vls', 'andor_norm']
# Comment: the first integrating detector will set the sub-sampling of all
# integrating detectors.
slow_detectors = [] # To do


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
        # ret_dict['andor_dir'] = [full_roi]
        # ret_dict['andor_vls'] = [full_roi]
        # ret_dict['andor_norm'] = [full_roi]
        # ret_dict['c_piranha'] = [full_roi]
        
        # and the FIM waveforms
        ret_dict['rix_fim0'] = [full_roi]
        ret_dict['rix_fim1'] = [full_roi]
        ret_dict['crix_w8'] = [full_roi]
    return ret_dict


def getFIMs(run):
    ret_dict = {}

    ret_dict['rix_fim0'] = {
        "sigROI" : slice(102,108),
        "bkgROI" : slice(0,80)
    }

    ret_dict['rix_fim1'] = {
        "sigROI" : slice(116,122),
        "bkgROI" : slice(0,80)
    }

    ret_dict['crix_w8'] = {
        "sigROI" : slice(69,76),
        "bkgROI" : slice(0,50)
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
epicsOncePV = []