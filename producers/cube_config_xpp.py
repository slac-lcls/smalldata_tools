import numpy as np

# custom bins
def binBoundaries(run):
    if isinstance(run,str):
        run=int(run)
    if run>0:
        return np.arange(-5.,15.,0.05)
    return None


# filters to apply ont the data
# format: list of [det (field), low, high, name]
# 'filter1' is the standard name and will not be added to the h5 file name.
filters = [
    ['lightStatus/xray',0.5,1.5,'filter1'],
    # ['ipm2/sum',3e2,6e4,'filter1'],
    ['ipm3/sum',200,3e3,'filter1'],
    ['tt/FLTPOSFWHM',60,230,'filter1'],
    ['tt/AMPL',0.005,0.19,'filter1']
#     ['evr/code_41',0.5,1.5,'custom']
]

# Laser on/off
laser = True

# List detectors to be cubed. Area detector have additional options such as threshold
# For now only full image work. TODO: Add photon maps. And then any detObjectFunc
# Detectors should be added to varList
detDict = {'source':'jungfrau1M', 
           'full':1, 
           'image':1, 
           'thresADU':6.5, 
           'common_mode':0}

varList = ['ipm2/sum','ipm3/sum','diodeU/channels', detDict]


# histogram configuration. Usually does not need to be changed (not yet implemented)
# field: destination in smd, list: [low,high,n] or [n]
hist_list = {
    'ipm2/sum': [],
    'ipm3/sum': [],
    'tt/FLTPOS_PS': [],
    'tt/AMPL': [],
    'tt/FLTPOSFWHM': []
}


# save as tiff file (ingore unless MEC)
save_tiff = False