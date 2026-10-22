import numpy as np

# These lists are needed, do not delete them
# If no detector in a given category, leave the corresponding
# list empty.
detnames = []

# define functions here
# 1) AZIMUTHAL INTEGRATION
def getAzIntLCLSGeomParams(run):
    if isinstance(run,str):
        run=int(run)
    ret_dict = {}
    if run>0:
        az_dict = {}
        az_dict["poni_file"] = 'path/to/poni/calib.poni'
        az_dict["npts"] = 256
        az_dict["npts_az"] = 36
        az_dict["int_units"] = 'q_A^-1'
        az_dict["return2d"] = False
        az_dict["azav_kwargs"] = {}
        #ret_dict["detname"] = [az_dict]
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

