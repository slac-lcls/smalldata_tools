import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding
# list empty.
detnames = ['epix_2',
            'epix10k2M']

# 1) REGIONS OF INTEREST
roi_detname: str = "epix_2"
azint_detname: str = "epix10k2M"

# 1) REGIONS OF INTEREST                                                                              
def getROIs(run):
    """ Set parameter for ROI analysis. Set writeArea to True to write the full ROI in the h5 file.   
    See roi_rebin.py for more info                                                                    
    """
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run > 0:
        roi_dict = {}
        roi_dict['ROI'] = [[200,500],[300,580]] 
        roi_dict['writeArea'] = False
        roi_dict['thresADU'] = 2.5 
        ret_dict[roi_detname] = [ roi_dict ]
    return ret_dict


def getAzIntParams(run):
    """Define parameters for Azimuthal integration without PyFAI."""
    if isinstance(run, str):
        run = int(run)

    ret_dict = {}
    if run > 0:
                                                                   
        az_dict = {'eBeam': 6} # Photon energy in keV                 
        az_dict['center'] = [-628.983255, -461.137890]# Position in um                
        az_dict['dis_to_sam'] = 60.0 # Position in mm                          
        az_dict['tx'] = 0 #0 # Tilt in x in degrees                                                 
        az_dict['ty'] = 0 #0 # Tilt in y in degrees                                                 
        az_dict['phiBins'] = 11 # Number of phi bins for azint                                      
        az_dict['qbin'] = 0.025 # Bin width in q for az int
        ret_dict[azint_detname] = [ az_dict ]

    return ret_dict

def getDetSums(run):
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    
    ret_dict['epix_2']=['calib','calib_thresADU1',
                        'calib_dropped','calib_dropped_square']
    ret_dict['epix10k2M']=['calib','calib_thresADU1',
                           'calib_dropped','calib_dropped_square']
    
    return ret_dict

#########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
# These lists are either PV names, aliases, or tuples with both.
epicsOncePV = []
epicsPV = ['las_Counter_Readback_ns', ("SP1L0:DCCM:MMS:TH1.RBV", 'dccm_th1'), ("SP1L0:DCCM:MMS:TH2.RBV", 'dccm_th2'), ("SP1L0:DCCM:MMS:TH1", 'dccm_th1_setpoint'), ("SP1L0:DCCM:MMS:TH2", 'dccm_th2_setpoint'),("SP1L0:DCCM:energy:OphydReadback",'dccm_E'),("SP1L0:DCCM:energy:OphydSetpoint",'dccm_E_setpoint')]
#fix timetool calibration if necessary
#ttCalib=[0.,2.,0.]
ttCalib=[]
#ttCalib=[1.860828, -0.002950]
