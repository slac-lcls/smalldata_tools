import numpy as np

# Custom bins
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
    ['ipm_dg1/sum',1e3,4.5e4,'filter1'],
    # ['ipm3/sum',200,3e3,'filter1'],
    # ['tt/FLTPOSFWHM',60,230,'filter1'],
    # ['tt/AMPL',0.005,0.19,'filter1']
    # ['evr/code_41',0.5,1.5,'custom']
]

# Laser on/off
laser = False
use_tt = False

# ##### DETECTORS ########################################################################
# List detectors to be cubed. Area detector have additional options such as threshold
# Full images will always be saved. 
# More area detector can be defined following the same syntax, and adding them to varList
detDict = {'source':'Rayonix', 
           'full':1, 
           'image':1, 
           'thresADU':-1e5, 
           'common_mode':0}

pix_size = 176e-6
func_kwargs = {
    'name': 'azav_pyfai', # must be the name of the smalldata_tools ana_func
    'ai_kwargs': {'dist':1, 'poni1':960*pix_size, 'poni2':960*pix_size},
    # 'poni_file': '<poni_file_path>' # alternative to ai_kwargs if a pyfai poni file exists
    'npts': 512,
    'int_units' : '2th_deg',
    'return2d' : False
}
det_proc = [func_kwargs]
detDict['det_proc'] = det_proc

# make list of all variables to be added to the cube
varList = ['ipm_dg1/sum','ipm_dg2/sum', detDict]




# ##### LOGBOOK SUMMARIES #################################################################
# histogram configuration. Usually does not need to be changed
# field: destination in smd, list: [low,high,n] or None (then default to some percentile)
hist_list = {
    'ipm_dg1/sum': [0, 5e4, 70],
    'ipm_dg2/sum': [0, 5e3, 70],
    # 'tt/FLTPOS_PS': [-0.5, 0.5, 70],
    # 'tt/AMPL': [0, 0.2, 70],
    # 'tt/FLTPOSFWHM': [0, 300, 70]
}


# save as tiff file (ingore unless MEC)
save_tiff = False
