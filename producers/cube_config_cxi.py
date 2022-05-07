import numpy as np

# Custom bins
# currently used for timing scans w/ timetool correction where you might want binning !- step size.
def binBoundaries(run):
    if isinstance(run,str):
        run=int(run)
    if run>0:
        # return np.arange(-5.,50.,0.2)
        return np.arange(-5.,10.,1)
    return None

# ##### FILTERS ##########################################################################
# filters to apply to the data
# format: list of [det (field), low, high, name]
# 'filter1' is the standard name and will not be added to the h5 file name.
filters = [
    ['lightStatus/xray',0.5,1.5,'filter1'],
    # ['ipm_dg1/sum',1e3,4.5e4,'filter1'],
    ['ipm_dg2/sum',50,125000.,'filter1'],
    # ['tt/FLTPOSFWHM',60,230,'filter1'],
    # ['tt/AMPL',0.005,0.19,'filter1']
    # ['evr/code_41',0.5,1.5,'custom']
]

# Laser on/off
laser = True
use_tt = False

# ##### DETECTORS ########################################################################
# List detectors to be cubed. Area detector have additional options such as threshold
# Full images will always be saved. 
# More area detector can be defined following the same syntax, and adding them to varList
detDict = {'source':'jungfrau4M', 
           'full':1, 
           'image':0, 
           'thresADU':1.2, 
           'common_mode':30,
           'maxHis':20}


# make list of all variables to be added to the cube
varList = ['ipm_dg2/sum','ipm_dg3/sum', 'ipm_dg2/peaks','ipm_dg3/peaks', 'jungfrau4M/azav_azav', detDict]
#varList = [detDict]
#varList = ['ipm_dg2/sum', detDict]

# ##### LOGBOOK SUMMARIES #################################################################
# histogram configuration. Usually does not need to be changed
# field: destination in smd, list: [low,high,n] or None (then default to some percentile)
hist_list = {
    # 'ipm_dg1/sum': [0, 5e4, 70],
    'ipm_dg2/sum': [-1000, 1.75e5, 175],
    'ipm_dg3/sum': [-5e4, 0, 125],
    'tt/FLTPOS_PS': [-0.5, 0.5, 70],
    'tt/AMPL': [0, 0.2, 70],
    'tt/FLTPOSFWHM': [0, 300, 70]
}

# save as tiff file (ingore unless MEC)
save_tiff = False
