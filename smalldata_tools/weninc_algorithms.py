import numpy as np
from scipy.signal import savgol_filter
import scipy.ndimage.measurements as smt
from scipy.ndimage.filters import maximum_filter


fp = np.array([[0, 1, 0],
               [1, 1, 1],
               [0, 1, 0]])

def find_droplets(img, seed_threshold, join_threshold):
    droplets, ndroplets = smt.label(img > seed_threshold)
    droplets = maximum_filter(droplets, footprint=fp)
    droplets[img < join_threshold] = 0
    droplets, ndroplets = smt.label(droplets)
    if ndroplets > 0:
        index = np.arange(1, ndroplets+1)
        adu = smt.sum(img, droplets, index)
        x, y = zip(*smt.center_of_mass(img, droplets, index))
        x = np.array(x, dtype=np.float32)
        y = np.array(y, dtype=np.float32)
    else:
        adu, x, y = None, None, None
    return ndroplets, x, y, adu

import scipy.ndimage.measurements as smt

def find_blobs(img, threshold, min_sum):
    blobs, nblobs = smt.label(img > threshold)
    if nblobs > 0:
        index = np.arange(1, nblobs+1)
        adu_sum = smt.sum(img, blobs, index)
        x, y = zip(*smt.center_of_mass(img, blobs, index))
        x = np.array(x, dtype=np.float32)
        y = np.array(y, dtype=np.float32)
        ind, = np.where(adu_sum > min_sum)
        nblobs = ind.size
        adu_sum = adu_sum[ind]
        x = x[ind]
        y = y[ind]
    else:
        x, y, adu_sum = None, None, None
    return nblobs, x, y, adu_sum

# Constant-Fraction Discriminator to find peaks in waveforms
def cfd(taxis, signal, fraction, delay, threshold, filter_length):
    bipolar = np.zeros(signal.shape)
    bipolar[delay:] = -fraction*signal[delay:] + signal[:-delay]
    bipolar = savgol_filter(bipolar, filter_length, 4)
    
    # find zero crossings from positive to negative
    pos = bipolar > 0
    ind, = (pos[:-1] & ~pos[1:]).nonzero()
    # only take zero crossing where signal is above threshold
    good_ind, = np.where(signal[ind] < -threshold)
    pos = ind[good_ind]
    
    # linear interpolation to find more precise zero crossing
    # y = y0 + (x - x0) * m  -->  x = x0 - y0/m
    dt = taxis[1] - taxis[0]
    m = (bipolar[pos+1] - bipolar[pos]) / dt
    return taxis[0] + pos*dt - bipolar[pos] / m