import numpy as np
from skimage import feature
from scipy import sparse, optimize
from skimage.measure import (CircleModel, ransac)
from scipy.signal import argrelextrema
from scipy.spatial import cKDTree
from numba import jit
from numba.typed import List as NTL
import itertools
from matplotlib import pyplot as plt
plt.ion()

def fitCircle(x,y,yerr=None, guess=None):
	"""
	fit a single circle. Transform input to lists to that fitCircles can be used.
	"""  
	x = [x]
	y = [y]
	rGuess = (np.nanmax(x)-np.nanmin(x)+np.nanmax(y)-np.nanmin(y))/4. #largest differences/2/2
	r = [rGuess]
	fitRes = _fit_circles(x,y,r,yerr=None, guess=None)
	#have only one circle.
	fitRes['R']=fitRes['R'][0]
	fitRes['residu']=fitRes['residu'][0]

	return fitRes

def find_edges(image, mask, sigma, hi_thresh, low_thresh):
    """Run the canny edge detection, probably doesn't need to live here

    Parameters
    ----------
    image: ndarray
        Scattering ring image to detect edges on
    mask: ndarray
        image mask with True/False values indicating whether to include in analysis
    sigma: float
        Standard deviation of the Gaussian filter
    hi_thresh: float
        Between 0 and 1. Upper bound for hysteresis thresholding (linking edges)
    low_thresh: float
        Between 0 and 1. Lower bound for hysteresis thresholding (linking edges)

    Returns
    -------
    edges: ndarray
        binary edges map of image
    sparse_edges: ndarray
        sparsified map of edges
    """
    edges = feature.canny(image, mask=mask.astype(bool), sigma=4 , low_threshold=0.92, \
                      high_threshold=0.98, use_quantiles=True)
    sparse_edges = sparse.coo_matrix(edges)
    
    return edges, sparse_edges

@jit(cache=True, nopython=True)
def _transform_hough_array(ar_hough, radii, zip_obj, center_x, center_y, dr, r_low, r_hi):
    """Transform Hough array on each iteration after narrowing search ranges
    Parameters
    ----------
    ar_hough: ndarray
        Initial Hough array to transform
    radii: np.array
        All the radii values from r ranges to check
    zip_obj: zipped object of np.arrays
        Used to iterate through the sparse_edges x, y, and data values
    center_x: np.array
        The x values of possible centers to iterate through
    center_y: np.array
        The x values of possible centers to iterate through
    dr: int
        The size of each radial step
    r_low: int
        The lower threshold for including in hough transform
    r_hi: int
        The upper threshold for including in hough transform
    """
    for row, col, data in zip_obj:
        for ix, cx in enumerate(center_x):
            dx = (row - cx) ** 2
            if dx < r_low or  dx > r_hi:
                # skip iteration if outside bounds
                continue
            for iy, cy in enumerate(center_y):
                dy = (col - cy) ** 2
                r = (dx + dy) ** 0.5
                ir = int((r - radii[0]) / dr)
                if ir >= 0 and ir < radii.shape[0]:
                    ar_hough[ir, ix, iy] += data

def _max_from_hough(ar_hough, radii, center_x, center_y):
    """Get maximum value from hough space for r, x, and y helper function
    
    Parameters
    ----------
    ar_hough: ndarray
        transformed hough array
    radii: np.array
        radii included in search
    center_x: np.array
        The x values of possible centers to iterate through
    center_y: np.array
        The x values of possible centers to iterate through

    Returns
    -------
    r_max: float
        max radius found in hough array radial plane
    x_max: float
        max x val found in hough array x plane
    y_max: float
        max y val found in hough array y plane
    maxdim_r: list(float)
        maximum values in radial plane along 0 axis of hough array
    """
    maxdim_r = [ar_hough[i,:,:].max() for i in range(ar_hough.shape[0])]
    maxdim_x = [ar_hough[:,i,:].max() for i in range(ar_hough.shape[1])]
    maxdim_y = [ar_hough[:,:,i].max() for i in range(ar_hough.shape[2])]

    r_max = radii[np.array(maxdim_r).argmax()]
    x_max = center_x[np.array(maxdim_x).argmax()]
    y_max = center_y[np.array(maxdim_y).argmax()]

    return r_max, x_max, y_max, maxdim_r

def iterate_center(sparse_edges, overfill, r_range, r_bin, c_bin, prec, red_factor, norm, min_dr):
    """Iterate through center finding until we are within defined precision, the main loop
    
    Parameters
    ----------
    sparse_edges: ndarray
        sparsified map of edges
    overfill: float
        mutliplier to extend range of search in x and y
    r_range: list(int)
        bounds for radial search
    r_bin: int
        divider used create number of bins for radial search
    c_bin: int
        divider used create number of bins for center search
    prec: float
        threshold for finishing center search
    red_factor: float
        factor we reduce search range by on each iteration
    norm: bool
        should we normalize extrema search for r values
    min_dr: int
        minimum delta r we need to include radius in final radius list

    Returns
    -------
    r_vals: list(int)
        accepted radii that pass threshold test
    x: float
        center x position
    y: float
        center y position
    """
    x_range = [sparse_edges.shape[0] * (1. - overfill), sparse_edges.shape[0] * overfill]
    y_range = [sparse_edges.shape[1] * (1. - overfill), sparse_edges.shape[1] * overfill]
    radii = np.arange(r_range[0], r_range[1], (r_range[1] - r_range[0]) / r_bin)
    dr = radii[1] - radii[0]
    r_low = radii[0] ** 2
    r_hi = radii[-1] ** 2
    center_x = np.arange(x_range[0], x_range[1], (x_range[1] - x_range[0]) / c_bin)
    center_y = np.arange(y_range[0], y_range[1], (y_range[1] - y_range[0]) / c_bin)
    ar_hough = np.zeros([radii.shape[0], center_x.shape[0], center_y.shape[0]])
    zip_obj_old = zip(sparse_edges.row, sparse_edges.col, sparse_edges.data)
    zip_obj = NTL()
    [zip_obj.append(x) for x in zip_obj_old]
    _transform_hough_array(ar_hough, radii, zip_obj, center_x, center_y, dr, r_low, r_hi)
    r, x, y, maxdim_r = _max_from_hough(ar_hough, radii, center_x, center_y)

    while (center_x[1] - center_x[0]) > prec:
        r_size_x = (x_range[1] - x_range[0]) / red_factor
        r_size_y = (y_range[1] - y_range[0]) / red_factor
        x_range = [x - r_size_x * 0.5, x + r_size_x * 0.5]
        y_range = [y - r_size_y * 0.5, y + r_size_y * 0.5]
        center_x = np.arange(x_range[0], x_range[1], (x_range[1] - x_range[0]) / c_bin)
        center_y = np.arange(y_range[0], y_range[1], (y_range[1] - y_range[0]) / c_bin)
        ar_hough = np.zeros([radii.shape[0], center_x.shape[0], center_y.shape[0]])
        _transform_hough_array(ar_hough, radii, zip_obj, center_x, center_y, dr, r_low, r_hi)
        r, x, y, maxdim_r = _max_from_hough(ar_hough, radii, center_x, center_y)

    x, y = y, x
    r_vals = _calc_r_vals(radii, maxdim_r, norm, min_dr)
    return r_vals, x, y

def _calc_r_vals(radii, maxdim_r, norm, min_dr):
    """Get the radii for the circles, helper function
    
    Parameters
    ----------
    radii: np.array
        possible radii values to check
    maxdim_r: list(float)
        maximum values in radial plane along 0 axis of hough array
    norm: bool
        should we normalize extrema search for r values
    min_dr: int
        minimum delta r we need to include radius in final radius list

    Returns
    -------
    nrad: list(int)
        All radii that pass threshold test
    """
    if norm:
        r_max = argrelextrema(np.array(maxdim_r) / radii, np.greater)
    else:
        r_max = argrelextrema(np.array(maxdim_r), np.greater)

    max_radii = radii[r_max[0]]

    if max_radii.shape[0] == 0:
        print('could not find a maxima in radii')
        return []

    res_at_radii = np.array(maxdim_r)[r_max[0]]
    res_at_radii, max_radii = (list(x) for x in zip(*sorted(zip(res_at_radii, max_radii))))

    nrad = []
    for rad in reversed(max_radii):
        cur_min_dr = 1e6
        for irrad in nrad:
            diff_val = np.abs(rad - irrad)
            if diff_val < cur_min_dr:
                cur_min_dr = diff_val
        if cur_min_dr > min_dr:
            nrad.append(rad)

    return nrad

def _fit_circles(x, y, r, yerr=False, guess=None):
    """
    simultanous least squares fitting of multiple concentric rings, helper function

    Parameters
    ----------
    x: float
        center in x
    y: float
        center in y
    r: list(int)
        accepted radii
    yerr: bool
        whether or not to get verbose output from scipy.optimize fit
    guess: list(float)
        provide a guess for the center

    Returns
    -------
    fit_res: dict
        all the parameters from least squares fit
    """  
    def calc_r(x, y, par):
        xc, yc = par
        return np.sqrt((x - xc) ** 2 + (y - yc) ** 2)

    def f(par, x, y):
        ri = calc_r(x, y, par)
        return ri - ri.mean()

    def f_global(par, mat_r, mat_x, mat_y):
        err = []
        for r, x, y in zip(mat_r, mat_x, mat_y):
            err_local = f(par, x, y)
            err = np.concatenate((err, err_local))          
        return err

    if guess is None or len(guess) != 2:
        x_m = np.mean(x[0])
        y_m = np.mean(y[0])
    else:
        x_m = guess[0]
        y_m = guess[1]
        
    center_estimate = x_m, y_m
    fit_res = {}  
    if yerr:      
        center, C, info, msg, success  = optimize.leastsq(f_global, \
            center_estimate, args=(r, x, y), full_output=True)
        fit_res['C'] = C
        fit_res['info'] = info
        fit_res['msg'] = msg
        fit_res['success'] = success
    else:
        center, ier = optimize.leastsq(f_global, center_estimate, args=(r, x, y))
        fit_res['ier'] = ier
    xc, yc = center
    fit_res['xCen'] = xc
    fit_res['yCen'] = yc
    rs=[]
    resids=[]
    for thisr, thisx, thisy in zip(r,x,y):
        ri     = calc_r(thisx, thisy, center)
        rs.append(ri.mean())
        resids.append(np.sum((ri - ri.mean()) ** 2))
    fit_res['residu'] = resids
    fit_res['R']      = rs

    return fit_res

def ransac_result(sparse_edges, center, r_vals, n_max, delta_r, min_points, min_samples, res_thresh, min_frac):
    """Find the number of edges inside proposed ring, if exceeds threshold, use ransac
    to fit the circle

    Parameters
    ----------
    sparse_edges: ndarray
        sparsified map of edges
    center: list(float)
        first element is center x value, second is y value
    r_vals: list(int)
        accepted radial guesses
    n_max: int
        maximum number of rings to check
    delta_r: int
        radial step size used in tree query
    min_points: int
        minimum absolute number of points in circle to be considered for final fit
    min_samples: int
        minimum number of samples to be drawn
    res_thresh: int
        allowed max residual to consider point being 'in'
    min_frac: float
        require fraction of points to pass RANSAC

    Returns
    -------
    comb_res: dict
        result from combined circle fitting
    ring_info: dict
        information about each ring that passed ransac fit
    sparse_edges: ndarray
        sparsified map of edges
    """
    tree = cKDTree(np.array([sparse_edges.col, sparse_edges.row]).T)
    ring_info = []
    for ir, r in enumerate(r_vals):
        if len(ring_info) >= n_max:
            break

        cur_ring_info = {}
        r_outer = tree.query_ball_point(center, r + delta_r)
        r_inner = tree.query_ball_point(center, r - delta_r)
        in_ring = set(r_outer).symmetric_difference(set(r_inner))
        cur_ring_info['r_input'] = r
        cur_ring_info['pointsInCircle'] = list(in_ring)

        if len(list(in_ring)) < min_points:
            continue

        model, inliers = ransac(np.array([sparse_edges.row[list(in_ring)], sparse_edges.col[list(in_ring)]]).T, \
            CircleModel, min_samples=min_samples, residual_threshold=res_thresh, max_trials=1000)
        
        cur_ring_info['in_frac'] = inliers.astype(int).sum() / float(len(list(in_ring)))
        cur_ring_info['inliers'] = inliers
        cur_ring_info['ransac_result'] = model

        if cur_ring_info['in_frac'] > min_frac and inliers.astype(int).sum() >= min_points:
            ring_info.append(cur_ring_info)

    all_x = []
    all_y = []
    all_r = []

    for ir, cur_ring_info in enumerate(ring_info):
        all_r.append(ir)
        all_x.append(sparse_edges.col[cur_ring_info['pointsInCircle']])
        all_y.append(sparse_edges.row[cur_ring_info['pointsInCircle']])

    try:
        comb_res = _fit_circles(all_x, all_y, all_r, yerr=True)
        return comb_res, ring_info, sparse_edges
    except Exception as e:
        print('Exception encountered during fitting: {0}'.format(e))
        return -1, ring_info, sparse_edges

def FindFitCenter(image, mask, **kwargs):
	"""arguments: image, mask (opt: inParams)
		image to use to find beam center (2-d)
		mask to be used (also 2-d)
		
		kwargs: dictionary of parameters used in finding the ebam center
			kwargs:
			#edge finding:
			sigma          #gaussian blurring for canny algorithm, def 1
			low_threshold  #threshold for only strongest features, using quantiles, def 0.92
			high_threshold #threshold for only strongest features, using quantiles, def 0.98
			#parameters for hough center finding
			precision #bin size in pixel for center in hough array, def 1
			nBinR     #number of bins for center, def 280
			overFac   #allow beam center to be out of the image by a factor of x, def 1.5 (50%)
			#parameters to find best rings
			norm      #normalize number of points/circle by its radius (point density), def True
			deltaR    #require new circle to be at least deltaR pixels away from last circle, def 5
			nMaxRing  #max number of rings to consider, def 6
			minInFrac #require 45% of points to pass RANSAC, def 0.45 
			minPoints #minimum absolute number of points in circle to be considered for final fit, def 40				
			RANSAC_residual_threshold #allowed max residual to consider point being 'in' (def 2)
			RANSAC_min_sample         #minimum number of samples to be drawn (def 10)
	"""
	sigma = kwargs.get('sigma', 4)
	hi_thresh = kwargs.get('high_threshold', 0.98)
	low_thresh = kwargs.get('low_threshold', 0.92)
	overfill = kwargs.get('overFac', 1.5)
	r_range = kwargs.get('r_range', [1, 1401])
	r_bin = kwargs.get('nBinR', 280)
	c_bin = kwargs.get('cBin', 100)
	prec = kwargs.get('precision', 1)
	red_factor = kwargs.get('deltaR', 5)
	norm = kwargs.get('norm', True)
	min_dr = kwargs.get('min_dr', -1)	
	n_max = kwargs.get('nMaxRing', 10)
	min_points = kwargs.get('minPoints', 40)
	min_samples = kwargs.get('RANSAC_min_sample', 10)
	res_thresh = kwargs.get('RANSAC_residual_threshold', 2)
	min_frac = kwargs.get('minInFrac', 0.45)

	edges, sparse_edges = find_edges(image, mask, sigma, hi_thresh, low_thresh)
	r_vals, x, y = iterate_center(sparse_edges, overfill, r_range, r_bin, c_bin, prec, red_factor, norm, min_dr)
	res, ring_info, edges = ransac_result(sparse_edges, [x, y], r_vals, n_max, red_factor, min_points, min_samples, res_thresh, min_frac)
	
	return res, ring_info, edges

def FindFitCenterPlot(image, mask, **kwargs):
        #plot image in grayscale
        combRes, ringInfo, arSparse = FindFitCenter(image, mask, **kwargs)

        #plot image in grayscale
        plt.figure(figsize=[12,12])
        plt.imshow(image, interpolation='none', cmap='gray',clim=[0,np.nanpercentile(image.flatten(),99.5)])
        #plot sparse images in blue.
        plt.plot(arSparse.col, arSparse.row,markersize=5,color='#ff9900',marker='.',linestyle='None')
        if combRes==-1:
            return -1
        greens=['#666600','#669900','#66cc00','#66cc99','#6699cc','#6633ff']
        print('center ',combRes['xCen'],combRes['yCen'])
        for ir,thisRingInfo,rFit in itertools.izip(itertools.count(),ringInfo,combRes['R']):
            #plot ransac selected data in green <-- consider different greens for first circles.
            plt.plot(arSparse.col[thisRingInfo['pointsInCircle']], arSparse.row[thisRingInfo['pointsInCircle']],marker='.',color=greens[ir],linestyle='None', markersize=4)
            circle1 = plt.Circle((combRes['xCen'],combRes['yCen']),rFit,color=greens[ir],fill=False,linestyle='dashed',linewidth=1)
            circle1b = plt.Circle((combRes['xCen'],combRes['yCen']),rFit,color='red',fill=False,linestyle='dotted',linewidth=2)
            plt.gca().add_artist(circle1)
            plt.gca().add_artist(circle1b)
        #circle1 = plt.Circle((combRes['xCen'],combRes['yCen']),ringInfo[0]['rInput'],color='r',fill=False,linestyle='dashed')
        #plt.gca().add_artist(circle1)
        plt.plot([combRes['xCen'],combRes['xCen']],[combRes['yCen']-20,combRes['yCen']+20],color='m')
        plt.plot([combRes['xCen']-20,combRes['xCen']+20],[combRes['yCen'],combRes['yCen']],color='m')
        limits=[[0,image.shape[0]],[0,image.shape[1]]]
        plt.xlim(limits[0])
        plt.ylim(limits[1])
        plt.show()
        for r in combRes['R']:
            print('aradii: ',r)
        print('center Final ',combRes['xCen'],combRes['yCen'])
        return combRes
 
