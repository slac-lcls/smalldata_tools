import numpy as np

# custom bins
def binBoundaries(run):
    if isinstance(run,str):
        run=int(run)
    if run>0:
        return np.arange(-5.,15.,0.05)
    return None


# filters to apply ont the data
# format: list of [det, low, high, name]
# 'filter1' is the standard name that will not be added to the h5 file name.
def get_filters(run):
    if isinstance(run,str):
        run = int(run)
    if run > 0:
        filters = [
            ['lightStatus/xray',0.5,1.5,'filter1'],
            ['ipm4/sum',2e3,3e4,'filter1'],
            ['tt/FLTPOSFWHM',30,175,'filter1'],
            ['tt/FLTPOS_PS',.05,.75,'filter1'],
            ['tt/AMPL',.02,.23,'filter1']
        ]
    return filters

# Laser on/off
"""
Set to laser=True if the laser sequence is being used, set to False if the laser sequence is not being used
"""
laser = False
use_tt = True

# List detectors to be cubed. Area detector have additional options such as threshold
# For now only full image work. TODO: Add photon maps. And then any detObjectFunc
# Detectors should be added to varList
detDict = {'source':'epix_1',
           'full':1,
           'image':1,
           'thresADU':6.5,
           'common_mode':None}

varList = ['ipm4/sum','ipm5/sum', detDict]


# histogram configuration. Usually does not need to be changed (not yet implemented)
# field: destination in smd, list: [low,high,n] or [n]
hist_list = {
    'ipm4/sum': [0, 5e4, 70],
    'ipm5/sum': [0, 5e3, 70],
    'tt/FLTPOS_PS': [0,1,50],
    'tt/AMPL': [0,.25,50],
    'tt/FLTPOSFWHM': [0,300,60]
}

# save as tiff file (ingore unless MEC)
save_tiff = False
