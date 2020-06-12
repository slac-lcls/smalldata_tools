#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import logging
from scipy.optimize import curve_fit
from argparse import ArgumentParser
import json
import os
import sys
from bokeh.plotting import figure
from bokeh.models import Span, Legend, LegendItem
import panel as pn

# CONSTANTS
# It'd be nice to do this in a full stack-ish way with a DB

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

# Guassian Fitting
LINE_FIT_POINTS = 50
RADIAL_RANGE = 5

# Logger
f = '%(asctime)s - %(levelname)s - %(filename)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.DEBUG, format=f)
logger = logging.getLogger(__name__)

# Math helper functions
def gaussian(x, a, mean, std, m, b):
	"""Equation for gaussian with linear component/offset"""
	return (a * np.exp(-((x - mean) / 2 / std) ** 2)) + (m * x + b)

def fit_line(ave_azav):
	"""Fit the line from edges of array"""
	azav_len = len(ave_azav)
	x0 = LINE_FIT_POINTS / 2
	x1 = azav_len - (LINE_FIT_POINTS / 2)
	y0 = np.mean(ave_azav[0:LINE_FIT_POINTS])
	y1 = np.mean(ave_azav[azav_len - LINE_FIT_POINTS:])
	m, b = np.polyfit((x0, x1), (y0, y1), 1)

	return m, b

def calc_intensity_bounds(f, i_mon, low, hi):
	"""Calculate lower and upper bounds to be used in jet tracking
	and generate indices of events to be used in calibration

	Parameters
	----------
	f: h5py File Object
		The smalldata file to be analyzed
	i_mon: str
		The name of the intensity monitor we'll use
	low: float
		The number of sigmas below mean to include
	hi: float
		The number of sigmas above mean to include

	Returns
	-------
	i_low: float
		lower bound for intensity data
	i_hi: float
		upper bound for intensity data
	i_dat: np.ndarray
		The intensity data
	"""
	# Set intensity monitor get data and set nan to 0.0
	if GDET_KEY not in f.keys():
		raise KeyError('Could not find intensity monitor key')

	logger.debug('Getting intensity data for {i_mon}')
	i_dat = np.array(f[GDET_KEY][i_mon])
	i_dat = np.nan_to_num(i_dat, copy=False)

	# Get upper and lower bounds
	mean = i_dat.mean()
	std = np.std(i_dat)

	i_low = round(mean - low * std, 2)
	i_hi = round(mean + hi * std, 2)
	logger.debug(f'Lower:Upper bounds for intensity monitor: {i_low}:{i_hi}')

	return i_dat, i_low, i_hi


class JTCal:
	def __init__(self, exp, run):
		self.hutch = exp[:3]  # meh, string parsing...
		self.exp = exp
		self.run = run
		self._i_dat = None
		self._idxs_use = None
		self._i_mon = None
		self._i_low = None
		self._i_hi = None
		self._azav_use = None
		self._ave_azav = None
		self._peak = None
		self._intensity = None
		self._rad_start = None
		self._rad_end = None
		self._idat_ave = None
		if not os.path.isfile(self.sd_file):
			raise OSError('hdf5 file does not exist')
		self.f = h5py.File(self.sd_file, 'r')

	@property
	def hutch(self):
		"""Current hutch set to get calib data for"""
		return self._hutch

	@hutch.setter
	def hutch(self, hutch):
		"""Set the hutch you want to use"""
		if hutch.lower() not in HUTCHES:
			raise ValueError('Passed Hutch not in list of jet tracking hutches')
		
		self._hutch = hutch

	@property
	def run(self):
		"""Run number to analyze as string"""
		return self._run

	@run.setter
	def run(self, run):
		"""Get the current run number being analyzed"""
		if not isinstance(run, str):
			raise TypeError('Run must be a string')

		self._run = run

	@property
	def exp(self):
		"""Get the experiment name we're looking at"""
		return self._exp

	@exp.setter
	def exp(self, exp):
		"""Set the experiment to analyze"""
		if not isinstance(exp, str):
			raise TypeError('Experiment must be a string')

		self._exp = exp

	@property
	def exp_path(self):
		"""Path to experiment directory"""
		return f'{DATA_PATH}{self.hutch}/{self.exp}'

	@property
	def sd_file(self):
		"""Full path and file name for smalldata file"""
		return f'{self.exp_path}{HDF5_EXT}{self.exp}_Run{self.run}.h5'

	@property
	def f(self):
		"""This is the h5py File object"""
		return self._f

	@f.setter
	def f(self, f):
		"""Set the h5py File Object"""
		if not isinstance(f, h5py._hl.files.File):
			raise TypeError('You must pass an h5py File object')

		self._f = f

	@property
	def i_dat(self):
		"""Intensity data being used for given monitor (pmt, wave8)"""
		return self._i_dat

	@i_dat.setter
	def i_dat(self, i_dat):
		"""Set the intensity data, must be np.array"""
		if not isinstance(i_dat, np.ndarray):
			raise TypeError('Intensity data must be an np.array')

		self._i_dat = i_dat

	@property
	def ave_azav(self):
		"""The average binned azimuthal average from all events"""
		return self._ave_azav

	@ave_azav.setter
	def ave_azav(self, ave_azav):
		"""Set the average azav array for all used events"""
		if not isinstance(ave_azav, np.ndarray):
			raise TypeError('Average azimuthal binned array must by np.array')

		self._ave_azav = ave_azav

	@property
	def idxs_use(self):
		"""The indices of each event we are going to use in analysis"""
		return self._idxs_use

	@idxs_use.setter
	def idxs_use(self, idxs_use):
		"""Set the indices of events to use"""
		if not isinstance(idxs_use, tuple):
			raise TypeError('Indices must be a tuple')

		self._idxs_use = idxs_use

	@property
	def i_mon(self):
		"""The intensity monitor being used"""
		return self._i_mon

	@i_mon.setter
	def i_mon(self, i_mon):
		"""Set the intensity monitor"""
		if i_mon not in I_MONS:
			raise ValueError('Passed Intensity Monitor is not supported')

		self._i_mon = i_mon

	@property
	def peak(self):
		"""The peak location in terms of q-bins of gaussian fit of ave_azav"""
		return self._peak

	@peak.setter
	def peak(self, peak):
		"""Set the peak location of guassian fit of ave_azav"""
		if not isinstance(peak, int):
			raise TypeError('Peak q bin of ave azav must be int')

		self._peak = peak

	@property
	def intensity(self):
		"""Intensity of average azav q binned array with some range around peak"""
		return self._intensity

	@intensity.setter
	def intensity(self, intensity):
		"""Set the intensity of the ave azav array around peak"""
		if not isinstance(intensity, float):
			raise TypeError('Inensity must be a float')

		self._intensity = intensity

	@property
	def i_low(self):
		"""Lower bound of intensity monitor used to find events"""
		return self._i_low

	@i_low.setter
	def i_low(self, i_low):
		"""Set lower bound of intensity monitor"""
		if not isinstance(i_low, float):
			raise ValueError('Lower bound of intensity monitor must be a float')

		self._i_low = i_low

	@property
	def i_hi(self):
		"""Upper bound of intensity monitor used to find events"""
		return self._i_hi

	@i_hi.setter
	def i_hi(self, i_hi):
		"""Set the upper bound of the intensity monitor"""
		if not isinstance(i_hi, float):
			raise ValueError('Upper bound of intensity monitor must be a float')

		self._i_hi = i_hi

	@property
	def results(self):
		"""Get results from calibration calculations"""
		res_dict = {
			'bin_peak': self.peak,
            'intensity': self.intensity,
			'bin_peak_min': self.rad_start,
			'bin_peak_max': self.rad_end
		}
		
		return res_dict

	@property
	def rad_start(self):
		"""Starting value for intensity integration on azav array"""
		return self._rad_start

	@rad_start.setter
	def rad_start(self, rad_start):
		self._rad_start = rad_start

	@property
	def rad_end(self):
		"""Ending value for intensity integration on azav array"""
		return self._rad_end

	@rad_end.setter
	def rad_end(self, rad_end):
		self._rad_end = rad_end

	@property
	def i_dat_mean(self):
		return self.i_dat.mean()

	@property
	def azav_use(self):
		"""Get the azimuthal avarage arrays of used indices"""
		return self._azav_use

	@azav_use.setter
	def azav_use(self, azav_use):
		self._azav_use = azav_use

	def calc_i_dat_events(self, i_dat, i_low, i_hi):
		"""Get the indices of events to use and the intensity data from those events"""
		self.idxs_use = np.where((i_dat >= i_low) & (i_dat <= i_hi))
		self.i_dat = i_dat[self.idxs_use]

	def det_azav(self):
		"""Get the average azimuthally q binned array from used events"""
		# Find the detector term
		shared_det = set(DETS).intersection(self.f.keys())
		if not bool(shared_det):
			raise ValueError('smalldata file does not contain known detector')

		# Assuming they didn't have two different detector keys...
		det = shared_det.pop()
		logger.debug(f'Getting azimuthal values for detector: {det}')
		det_group = self.f[det]
		azav = np.array(det_group[AZAV_KEY])
		self.azav_use = azav[self.idxs_use]
		self.ave_azav = (self.azav_use.sum(axis=0) / len(self.azav_use))[0]
		logger.debug('Found average azimuthal average array from all events')

	def calc_peak_and_intensity(self):
		"""Fit the gaussian with linear offset and get the peak 
		and integrated intensity around peak
		"""
		logger.debug('Fitting Guassian of average azimuthal binned array')
		azav_len = len(self._ave_azav)
		# Fit the line
		m, b = fit_line(self._ave_azav)
			
		# Estimate mean and std
		x = np.arange(azav_len)
		mean = sum(x * self.ave_azav) / sum(self.ave_azav)
		std = np.sqrt(sum((x - mean) ** 2 / azav_len))

		# Fit Gaussian and get center and integrated intensity
		popt, _ = curve_fit(gaussian, x, self.ave_azav, p0=[max(self.ave_azav), mean, std, m, b])
		self.peak = int(round(popt[1]))
		self.rad_start = self.peak - RADIAL_RANGE
		self.rad_end = self.peak + RADIAL_RANGE
		self.intensity = round(self.ave_azav[self.rad_start:self.rad_end].sum(axis=0) / self.i_dat.mean(), 2)
		logger.debug(f'Peak found at {self.peak}, with and intensity of {self.intensity}')

	def write_file(self):
		"""Write the data to the calib directory"""
		cal_dir = ''.join([self.exp_path, CAL_DIR])
		if not os.path.exists(cal_dir):
			logger.debug('calibration directory does not exist, creating')
			os.makedirs(cal_dir)
 
		logger.debug(f'Writing results to {cal_dir}')
		with open(''.join([cal_dir, CAL_FILE, '_', self.run]), 'w') as output:
			json.dump(self.results, output, indent=4, sort_keys=True)

	def intensity_fig(self):
		"""Generate a histogram of intensity monitor data with bounds"""
		fig = figure(
    		title=f'Used Intensity Distribution for {self.i_mon}',
    		x_axis_label='Intensity Values',
    		y_axis_label='Counts'
		)

		hist, edges = np.histogram(self.i_dat, bins=100)
		fig.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:])
		left_line = Span(location=edges[0], dimension='height', line_color='black')
		right_line = Span(location=edges[-1], dimension='height', line_color='black')
		fig.renderers.extend([left_line, right_line])

		return fig

	def idat_vs_peak_fig(self):
		"""Generate scatter plot of intensity vs peak of azav array"""
		intensity_data = self.i_dat
		# TODO: speed up with correct algorithm
		peak_vals = [azav[0][self.peak] for azav in self.azav_use]
		ratio = np.mean(peak_vals) / np.mean(self.i_dat)
		fig = figure(
			title=f'Intensity vs peak value.  Azav/Intensity ratio - {round(ratio, 2)}',
			x_axis_label='Intensity Values',
			y_axis_label='Azimuthal Average Peak'
		)
		fig.scatter(intensity_data, peak_vals)
		ratio_line = Span(location=ratio, dimension='width', line_color='red', line_width=3)
		fig.renderers.extend([ratio_line])

		return fig

	def azav_fig(self):
		"""Generate the azav fig for html file"""
		x_vals = np.arange(len(self.ave_azav))
		fig = figure(
			title=f'Average Azimuthal Binned Array: Center - {self.peak}, min/max - {self.rad_start}/{self.rad_end}, intensity - {round(self.intensity, 2)}',
			x_axis_label='Bins',
			y_axis_label='Intensity',
		)
		
		peak_line = Span(location=self.peak, dimension='height', line_color='green', line_width=2)
		lower_line = Span(location=self.rad_start, dimension='height', line_color='black')
		upper_line = Span(location=self.rad_end, dimension='height', line_color='black')
		ave_azav_curve = fig.scatter(x_vals, self.ave_azav)
		fig.renderers.extend([peak_line, lower_line, upper_line])

		azav_legend = Legend(items=[
			LegendItem(label='Azimuthal Average', renderers=[ave_azav_curve])
		])
		fig.add_layout(azav_legend)

		return fig

	def generate_html_file(self):
		"""Generate the html file for reporting"""
		gspec = pn.GridSpec(sizing_mode='stretch_both', max_height=1000)
		gspec[0:3, 0:3] = self.intensity_fig()
		gspec[4:6, 0:3] = self.azav_fig()
		gspec[7:9, 0:3] = self.idat_vs_peak_fig()		
		gspec.save('report.html')

if __name__ == '__main__':
	# This needs to be run as script
	exp = os.environ.get('EXPERIMENT', 'cxilr6716')
	run = os.environ.get('RUN_NUM', '139')
	parser = ArgumentParser()
	parser.add_argument('--exp', type=str, default=exp)
	parser.add_argument('--run', type=str, default=run)
	parser.add_argument('--imon', type=str, default='f_21_ENRC')
	parser.add_argument('--sig_low', type=float, default=1.5)
	parser.add_argument('--sig_hi', type=float, default=1.5)
	args = parser.parse_args()

	print('here is exp ', args)

	j = JTCal(args.exp, args.run)
	i_dat, low, hi = calc_intensity_bounds(j.f, args.imon, args.sig_low, args.sig_hi)
	j.calc_i_dat_events(i_dat, low, hi)
	j.det_azav()
	j.calc_peak_and_intensity()
	j.write_file()
	j.generate_html_file()
