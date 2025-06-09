import numpy as np

#This is the XCS configuration file that will automatically generate h5 files from the smd_producer.py code. 

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding
# list empty.
detnames = ['epix_2']

# 1) REGIONS OF INTEREST
roi_detname: str = "epix_2"


# 1) REGIONS OF INTEREST                                                                              
def getROIs(run):
    """ Set parameter for ROI analysis. Set "writeArea=True" to write the defined ROI to the h5 file.    
    See roi_rebin.py for more info                                                                    
    """
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run > 0:
        roi_dict = {}
        roi_dict['ROI'] = [[200,500],[300,580]] #[[x1,x2], [y1,y2]] 
        roi_dict['writeArea'] = True
        roi_dict['thresADU'] = 0.0
        ret_dict[roi_detname] = [ roi_dict ]
    return ret_dict


def getDetSums(run):
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    
    ret_dict['epix_2']=['calib','calib_thresADU1',
                        'calib_dropped','calib_dropped_square']
    
    return ret_dict



#########################################################
# run independent parameters 
##########################################################
#aliases for experiment specific PVs go here
# These lists are either PV names, aliases, or tuples with both.
epicsOncePV = []
epicsPV = []
#fix timetool calibration if necessary
ttCalib=[]

