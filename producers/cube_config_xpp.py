import numpy as np

# Custom bins
def binBoundaries(run):
    if isinstance(run,str):
        run=int(run)
    if run>0:
        # return np.arange(-5.,50.,0.2)
        return np.arange(-5.,10.,1)
    return None

# Multi-dimensional binning
# Additional variables to bin along can be passed as a dictionary. Use similar run logic
# as for main bin definition
def get_addBinVars(run):
    if isinstance(run,str):
        run=int(run)
    addBinVars = None
    if run==128:
        addBinVars = {'ipm2/sum': np.linspace(0,4e4,4)}
    return addBinVars


# filters to apply to the data
# format: list of [det (field), low, high, name]
# 'filter1' is the standard name and will not be added to the h5 file name.
filters = [
    ['lightStatus/xray',0.5,1.5,'filter1'],
    #['ipm2/sum',1e3,4.5e4,'filter1'],
    # ['ipm3/sum',200,3e3,'filter1'],
    #['tt/FLTPOSFWHM',60,230,'filter1'],
    #['tt/AMPL',0.005,0.19,'filter1']
    # ['evr/code_41',0.5,1.5,'custom']
]

# Laser on/off
laser = True
use_tt = True

# List detectors to be cubed. Area detector have additional options such as threshold
# For now only full image work. TODO: Add photon maps. And then any detObjectFunc
# Detectors should be added to varList
detDict = {'source':'jungfrau1M', 
           'full':1, 
           'image':1, 
           'thresADU':5, 
           'common_mode':None}

varList = ['ipm2/sum','ipm3/sum','diodeU/channels', detDict]


# histogram configuration. Usually does not need to be changed
# field: destination in smd, list: [low,high,n] or None (then default to some percentile)
hist_list = {
    'ipm2/sum': [0, 5e4, 70],
    'ipm3/sum': [0, 5e3, 70],
    'tt/FLTPOS_PS': [-0.5, 0.5, 70],
    'tt/AMPL': [0, 0.2, 70],
    'tt/FLTPOSFWHM': [0, 300, 70]
}


# save as tiff file (ingore unless MEC)
save_tiff = False
