import h5py
import numpy as np
import matplotlib.pyplot as plt
import logging
from scipy.optimize import curve_fit
from argparse import ArgumentParser
import json
import os
import sys
from bokeh.plotting import figure, show
from bokeh.models import Span, Legend, LegendItem, ColorBar, LinearColorMapper
from bokeh.io import output_notebook
import panel as pn

# Constants

# Hutches that support jet tracking
HUTCHES = [
    'cxi',
    'mfx'
]

DETS = [
    'DsaCsPad',
    'DsbCsPad',
    'DscCsPad'
]

# The possible intensity monitors to use, prefer f_2.. since
# this is after the solid attenuators (maybe still true?)
# TODO: find wave8 key and add to I_MONS
I_MONS = [
    'f_11_ENRC',
    'f_12_ENRC',
    'f_21_ENRC',
    'f_22_ENRC',
    'f_63_ENRC',
    'f_64_ENRC'
]

I_MON_DEFAULT = 'f_21_ENRC' 

# base path to experiment data
DATA_PATH = '/reg/d/psdm/'

# extension to small data hdf5 files
HDF5_EXT = '/hdf5/smalldata/'

# Calib directory to write results to
CAL_DIR = '/calib/'
CAL_FILE = 'jt_cal'

# dataset keys
GDET_KEY = 'gas_detector'
AZAV_KEY = 'azav_azav'

# i0 treshold
I0_THRESH = 0.1

# Guassian Fitting
LINE_FIT_POINTS = 50
Q_RANGE = 5

# Lower bound on azav values
AZAV_LOW_THRESH = 0.1

# Logger
f = '%(asctime)s - %(levelname)s - %(filename)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.DEBUG, format=f)
logger = logging.getLogger(__name__)


def gaussian(x, a, mean, std, m, b):
    """Equation for gaussian with linear component/offset"""
    return (a * np.exp(-((x - mean) / 2 / std) ** 2)) + (m * x + b)

def fit_line(array, fit_points=LINE_FIT_POINTS):
    """Fit the line from edges of array"""
    azav_len = len(ave_azav)
    x0 = fit_points / 2
    x1 = azav_len - (fit_points / 2)
    y0 = np.mean(ave_azav[:fit_points])
    y1 = np.mean(ave_azav[azav_len - fit_points:])
    m, b = np.polyfit((x0, x1), (y0, y1), 1)

    return m, b

def peak_lr(array_data, threshold=0.05, bins=100):
    """Find max of normal distribution from histogram, search right and left until
    population falls below threshold
    """
    hist, edges = np.histogram(array_data, bins=bins)
    
    # Find peak information
    peak_val = hist.max()
    peak_idx = np.where(hist == peak_val)[0][0]
    peak_bin = edges[peak_idx]
    
    #search right
    right = np.argmax(hist[peak_idx:] < threshold * peak_val)
    right += peak_idx
    
    # search left
    left_array = hist[:peak_idx]
    left = peak_idx - np.argmax(left_array[::-1] < threshold * peak_val)
    
    return hist, edges, peak_idx, left, right

def det_azav(f, idxs_use):
    """Get the average azimuthally q binned array from used events"""
    # Find the detector term
    shared_det = set(DETS).intersection(f.keys())
    if not bool(shared_det):
        raise ValueError('smalldata file does not contain known detector')

    # Assuming they didn't have two different detector keys...
    det = shared_det.pop()
    logger.debug(f'Getting azimuthal values for detector: {det}')
    det_group = f[det]
    azav = np.array(det_group[AZAV_KEY])
    azav_use = azav[idxs_use]
    ave_azav = (azav_use.sum(axis=0) / len(azav_use))[0]
    logger.debug('Found azimuthal average arrays for all used events')

    return azav_use, ave_azav

def get_peak_azav_bin(ave_azav):
    """Fit the gaussian with linear offset and get the peak 
    and integrated intensity around peak
    """
    logger.debug('Fitting Guassian of average azimuthal binned array')
    azav_len = len(ave_azav)
    # Fit the line
    m, b = fit_line(ave_azav)

    # Estimate mean and std
    x = np.arange(azav_len)
    mean = sum(x * ave_azav) / sum(ave_azav)
    std = np.sqrt(sum((x - mean) ** 2 / azav_len))

    # Fit Gaussian and get center and integrated intensity
    popt, _ = curve_fit(gaussian, x, ave_azav, p0=[max(ave_azav), mean, std, m, b])
    peak = int(round(popt[1]))
    #bin_start = peak - delta_bin
    #bin_end = peak + delta_bin
    #integrated_intensity = round(ave_azav[rad_start:rad_end].sum(axis=0) / median_intensity, 2)
    #logger.debug(f'Peak found at {peak}, with integrated intensity of {integrated_intensity}') 

    return peak

def peak_in_bins(azav_use, peak, i0_data, left_bin, right_bin):
    """Get the peak in the bins we're using"""
    peak_vals = np.array([azav[0][low:hi].sum(axis=0) for azav in azav_use])
    peak_idxs = np.where(peak_vals > low_limit)
    peak_vals = peak_vals[peak_idxs]
    ratio = peak_vals.mean() / i_dat.mean()

    return peak_vals, peak_idxs, ratio

def fit_limits(i0_data, azav_intensity_vals, i0_edges):
    m, b = np.polyfit(i0_data, azav_intensity_vals, 1)
    x = np.arange(i0_edges[0], i0_edges[-1], 0.1)
    sigma = np.std(azav_intensity_vals)
    y = x * m + b
    
    return x, y, m, b, sigma

def write_html(cal_dir, results, run):
    """Write the data to the calib directory"""
    if not os.path.exists(cal_dir):
        logger.debug('calibration directory does not exist, creating')
        os.makedirs(cal_dir)

    logger.debug(f'Writing results to {cal_dir}')
    with open(''.join([cal_dir, CAL_FILE, '_', run]), 'w') as output:
        json.dump(results, output, indent=4, sort_keys=True)

####### HTML file figures ##########

def i0_fig(intensity_monitor, hist, edges, peak_idx, left, right):
    """General histogram plotter with peak location and left/right limits plotted"""
    fig = figure(
        title=f'Used Intensity Distribution for {intensity_monitor}. Low/Hi: {round(edges[left], 2)}/{round(edges[right], 2)}',
        x_axis_label='Intensity Values',
        y_axis_label='Counts'
    )
    fig.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:])
    left_line = Span(location=edges[left], dimension='height', line_color='black')
    right_line = Span(location=edges[right], dimension='height', line_color='black')
    peak_line = Span(location=edges[peak_idx], dimension='height', line_color='red')
    fig.renderers.extend([left_line, right_line, peak_line])
    
    return fig

def ave_azav_fig(ave_azav, peak_bin, bin_start, bin_end):
    """Generate the azav fig for html file"""
    x_vals = np.arange(len(ave_azav))
    fig = figure(
        title=f'Average Azimuthal Binned Array: Center - {peak_bin}, min/max - {bin_start}/{bin_end}',
        x_axis_label='Bins',
        y_axis_label='Intensity',
    )

    peak_line = Span(location=peak_bin, dimension='height', line_color='green', line_width=2)
    lower_line = Span(location=bin_start, dimension='height', line_color='black')
    upper_line = Span(location=bin_end, dimension='height', line_color='black')
    ave_azav_curve = fig.scatter(x_vals, ave_azav)
    fig.renderers.extend([peak_line, lower_line, upper_line])

    azav_legend = Legend(items=[
        LegendItem(label='Azimuthal Average', renderers=[ave_azav_curve])
    ])
    fig.add_layout(azav_legend)

    return fig
   
def azav_intensity_fig(azav_intensity_vals, low, bins=100):
    hist, edges = np.histogram(azav_intensity_vals, bins=bins)
    fig = figure(
        title=f'Azav Integrated Intensities. Low Cut Val: {round(low, 2)}',
        x_axis_label='Integrated Intensities',
        y_axis_label='Counts'
    )
    fig.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:])
    low_line = Span(location=low, dimension='height', line_color='red')
    fig.renderers.extend([low_line])
    
    return fig
    
def intensity_vs_peak_fig(i0_data, azav_intensity_vals, x, y, slope, intercept, sigma, threshold=1):
    """Simple plot of intensity vs peak value"""
    fig = figure(
        title=f'Peak value vs Intensity. Slope = {round(slope, 2)}, Intercept = {round(intercept, 2)}',
        x_axis_label='Intensity Monitor Value',
        y_axis_label='Peak Values'
    )
    fig.x_range.range_padding = fig.y_range.range_padding = 0
    h, y_edge, x_edge = np.histogram2d(azav_intensity_vals, i0_data, bins=100)
    fig.image(image=[h], x=x_edge[0], y=y_edge[0], dh=y_edge[-1]-y_edge[0], dw=x_edge[-1]-x_edge[0], palette="Spectral11")
    color_mapper = LinearColorMapper(palette="Spectral11", low=h.min(), high=h.max())
    color_bar = ColorBar(color_mapper=color_mapper, location=(0,0))
    fig.add_layout(color_bar, 'right')
    fig.line(x, y, color='red')
    fig.line(x, y - threshold * sigma, color='orange')
    fig.line(x, y + threshold * sigma, color='orange')
    
    return fig

if __name__ == '__main__':
 	# If run through arp make os vars default
    exp = os.environ.get('EXPERIMENT', 'cxilr6716')
    run = os.environ.get('RUN_NUM', '139')
    parser = ArgumentParser()
    parser.add_argument('--exp', type=str, default=exp)
    parser.add_argument('--run', type=str, default=run)
    parser.add_argument('--imon', type=str, default=I_MON_DEFAULT)
    parser.add_argument('--i0_thresh', type=float, default=I0_THRESH)
    parser.add_argument('--det', type=str)
    parser.add_argument('--q_range', type=int, default=Q_RANGE)
    parser.add_argument('--azav_thresh', type=float, default=AZAV_LOW_THRESH)
    parser.add_argument('--results_dir', type=str)
    parser.add_argument('--cal_dir', type=str)
    args = parser.parse_args()

    if not args.exp or not args.run:
        logger.debug('You have not specified an experiment or run, exiting')
        sys.exit()

    logger.debug(f'Starting Jet Tracking Calibration for Experiment {args.exp}, Run {args.run}')
    hutch = exp[:3]
    exp_path = f'{DATA_PATH}{hutch}/{exp}'
    results_path = f'{exp_path}/stats/summary/Run_{run}/jt_cal'
    sd_file = f'{exp_path}{HDF5_EXT}{exp}_Run{run}.h5'
    cal_path = f'{exp_path}/calib/jt_cal/'

    f = h5py.File(sd_file, 'r')
    logger.debug(f'Loaded smalldata file {sd_file} with keys {list(f.keys())}')

    # This will hold all the plots for the html file
    gspec = pn.GridSpec(sizing_mode='stretch_both')

    # Filter on i0 data
    i0_data = np.array(f[GDET_KEY][args.imon])
    i0_data = np.nan_to_num(i0_data, copy=False)
    i0_hist, i0_edges, i0_peak_idx, i0_left_idx, i0_right_idx = peak_lr(i0_data)
    i0_low = i0_edges[i0_left_idx]
    i0_high = i0_edges[i0_right_idx]
    i0_idxs = np.where((i0_data > i0_low) & (i0_data < i0_high))
    i0_data_use = i0_data[i0_idxs]

    # Generate Plot for i0 cut
    gspec[0:3, 0:3] = i0_fig(args.imon, i0_hist, i0_edges, i0_peak_idx, i0_left_idx, i0_right_idx)
    
    # Check if we can find the detector if not passed
    det = args.det
    if not det:
        det = set(DETS).intersection(f.keys())
        if not bool(det):
            logger.debug('Found no common detector and detector was not passed in, exiting')
            sys.exit()
        det = det.pop()
        logger.debug(f'You have not specified a detector, found {det} in dataset.  Will gather azimuthal data')

    # Get azimuthally binned data
    det_group = f[det]
    azav_data = np.array(det_group[AZAV_KEY])

    logger.debug('Filtering azav data on i0 cut and getting average azav array')
    # Only use data where we have i0 in threshold
    azav_data_filter = azav_data[i0_idxs]

    # Get the average azav array from the data filtered on i0
    ave_azav = (azav_data_filter.sum(axis=0) / len(azav_data_filter))[0]

    # Get the peak and range of the average azav curve
    peak_azav_bin = get_peak_azav_bin(ave_azav)
    left_azav_bin = peak_azav_bin - args.q_range
    right_azav_bin = peak_azav_bin + args.q_range

    logger.debug(f'Found azav peak at {peak_azav_bin}, with left bin at {left_azav_bin} and right bin at {right_azav_bin}')

    # Generate plot for average azimuthal array after i0 cut
    gspec[4:6, 0:3] = ave_azav_fig(ave_azav, peak_azav_bin, left_azav_bin, right_azav_bin)
    
    # Getting the intensity integration values for the azavs being used
    azav_intensity_vals = np.array([azav[0][left_azav_bin:right_azav_bin].sum(axis=0) for azav in azav_data_filter])
    low_cut = args.azav_thresh * np.mean(azav_intensity_vals)

    # Remove low shots
    azav_idxs = np.where(azav_intensity_vals > low_cut)

    # Generate plot for the azimuthal intensity cut
    gspec[7:9, 0:3] = azav_intensity_fig(azav_intensity_vals, low_cut) 

    # Get i0 and azav intensity data after filtering on both signals
    i0_data_use = i0_data_use[azav_idxs]
    azav_intensity_vals = azav_intensity_vals[azav_idxs]

    # Fit the line through the used points
    logger.debug('Fitting the distribution of i0 vs azav intensity values')
    x, y, m, b, sigma = fit_limits(i0_data_use, azav_intensity_vals, i0_edges)
    logger.debug(f'Found slope of fit: {round(m, 2)} with intercept of {round(b, 2)}')

    # 2d histogram with limits plotted
    gspec[10:12, 0:3] = intensity_vs_peak_fig(i0_data_use, azav_intensity_vals, x, y, m, b, sigma)

    results_dir = args.results_dir
    if not results_dir:
        results_dir = f'{results_path}/results.html'

    if not os.path.exists(results_dir):
        logger.debug('results directory does not exist, creating')
        os.makedirs(results_dir)

    logger.debug(f'Saving html file to {results_dir}')
    gspec.save(f'{results_path}/results.html')

    results = {
        'i0_low': i0_low,
        'i0_high': i0_high,
        'peak_azav_bin': peak_azav_bin,
        'left_azav_bin': left_azav_bin,
        'right_azav_bin': right_azav_bin,
        'slope_fit': m,
        'intercept': b,
        'sigma': sigma
    }
    cal_dir = args.cal_dir
    if not cal_dir:
        cal_dir = f'{cal_path}/Run_{run}'

    write_html(cal_dir, results, run)
