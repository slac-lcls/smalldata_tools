import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding
# list empty.
detnames = ['alvium_spot',
            'alvium_nf_ff',
            'Zyla_xray_east',
            'alvium_tt',
            'epix100a_2_bxrts',
            'epix100a_1_xtcs',
            'Epix10kaQuad2',
            'Epix10kaQuad3']

# 1) REGIONS OF INTEREST
def getROIs(run):
    """ Set parameter for ROI analysis. Set writeArea to True to write the full ROI in the h5 file.
    See roi_rebin.py for more info
    """
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run > 0:
        # Zyla_xray_east
        roi_dict = {
            'name':'ROI_0',
            'writeArea': True,
            'thresADU': None,
            'ROI':[[157,487], [294,598]]
        }
        ret_dict['Zyla_xray_east'] = [ roi_dict ]

        # epix100a_1_xtcs
        roi0_dict = {'name':'ROI_0',
                    'writeArea': True,
                    'thresADU': None,
                    'ROI':[[357,707], [0,768]],
                    'pj': [{'name':'pj1', 'axis':0, 'thresADU': 0.0 }]}
        roi1_dict = {'name':'ROI_1',
                    'writeArea': True,
                    'thresADU': 7.0,
                    'ROI':[[357,707], [0,768]],
                    'pj': [{'name':'pjt', 'axis':0, 'thresADU': 7.0 }]}
        ret_dict['epix100a_1_xtcs'] = [ roi0_dict, roi1_dict ]

        # epix100a_2_bxrts
        roi0_dict = {'name':'ROI_0',
                     'writeArea': True,
                     'thresADU': None,
                    'ROI':[[0,707], [360,370]],
                    'pj': [{'name':'pj21', 'axis':1, 'thresADU': 0.0 }]}
        roi1_dict = {'name':'ROI_1',
                     'writeArea': True,
                     'thresADU': None,
                     'ROI':[[0,707], [310,320]],
                     'pj': [{'name':'pj22', 'axis':1, 'thresADU': 0.0 }]}
        roi2_dict = {'name':'ROI_2',
                     'writeArea': True,
                     'thresADU': 1,
                     'ROI':[[0,707], [360,370]],
                     'pj': [{'name':'pj2t', 'axis':1, 'thresADU': 1.0 }]}
        ret_dict['epix100a_2_bxrts'] = [ roi0_dict, roi1_dict, roi2_dict ]
        
    return ret_dict


def get_mectt(run):
    if isinstance(run,str):
        run=int(run)
    tt_dict = {'bkgdet':'alvium_tt_calib', #needs to be the name of the Sums dataset
               'bkgfile':'/sdf/data/lcls/ds/mec/mec100860324/hdf5/smalldata/mec100860324_Run0165.h5',
               'sm':3,
               'tilt':-4.7, 
               'roi_s': np.array([[200,1450],[800,1900]]),
               'roi_n': np.array([[100,1400],[200,400]]),
               'sigmaTT': 18,
               'ttRaw': True,
               'debug': False
    }
    ret_dict = {'alvium_tt': [ tt_dict ]}
    return ret_dict


# pyFAI function. 
def getAzIntPyFAIParams(run):
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}

    if run > 0:
        az_dict = {}
        #pix_size = 100e-6
        #ai_kwargs = {
        #    "dist": 0.12022,  
        #    "poni1": 0.1,  
        #    "poni2": 0.05,  
        #}
        az_dict["return2d"] = False
        # az_dict["ai_kwargs"] = ai_kwargs
        az_dict["poni_file"] = (
            "/sdf/data/lcls/ds/mec/mec100860324/results/njh/Diffraction/run175_rot.poni"
        )
        az_dict["npts"] = 700
        #az_dict["npts_az"] = 11
        az_dict["int_units"] = "2th_deg"
        az_dict["threshold"] = 2.5
        az_dict["polarization_factor"] = -1
        ret_dict["Epix10kaQuad3"] = [ az_dict ]

    return ret_dict


def getDetSums(run):
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    
    ret_dict['Epix10kaQuad2'] = ['calib', 'calib_thresADU1', 'calib_dropped']
    ret_dict['Epix10kaQuad3'] = ['calib', 'calib_thresADU1', 'calib_dropped']
    ret_dict['alvium_tt'] = ['calib', 'calib_thresADU10']
    # add the maximum-image for calibration runs.
    if run > 0 :
        ret_dict['Epix10kaQuad2'].append('calib_max')
        ret_dict['Epix10kaQuad3'].append('calib_max')
    
    return ret_dict


##########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
# These lists are either PV names, aliases, or tuples with both.
epicsOncePV = []
epicsPV = []
#fix timetool calibration if necessary
ttCalib=[]

