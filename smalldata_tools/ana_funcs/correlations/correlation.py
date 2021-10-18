import numpy as np
import matplotlib.pyplot as plt
from lmfit import models
from scipy import ndimage
from smalldata_tools.ana_funcs.correlations import utils


def _spatial_correlation_fourier(fim1, fim2_star, fmask, fmask_star):
    A_num1 = np.fft.irfft2(fim1*fim2_star)
    A_num2 = np.fft.irfft2(fmask*fmask_star)
    A_num2 *= A_num2>0.1 # remove Fourier components that are too small
    A = A_num1 * A_num2
    A_denom = np.fft.irfft2(fim1*fmask_star) * np.fft.irfft2(fim2_star*fmask) # symmetric normalization
    # make sure the normalization value isn't 0 otherwise the autocorr will 'explode'
    pos = np.where(np.abs(A_denom)!=0)
    A[pos] /= A_denom[pos]
    A = np.fft.fftshift(A)
    return A


def spatial_correlation_fourier(img, img2=None, mask=None):
    """ Compute spatial correlation between two images using Fourier transform.
    
    Args:
        img: first image
        img2: second image. If None, becomes the first image (autocorrelation)
        mask: mask spanning the region of interest
    
    Returns:
        A: 2d correlation matrix
    """
    if mask is None:
        mask = np.ones_like(img)
    if img2 is None:
        img2 = img
        
    # (i) restrain mask and imgs to a bounding box around the roi defined by the mask
    img, temp = utils.box_to_roi(img, mask)
    img2, mask = utils.box_to_roi(img2, mask)
    
    # (ii) compute the different terms
    fmask = np.fft.rfft2(mask)
    fmask_star = np.conjugate(fmask)
    fimg = np.fft.rfft2(img)
    fimg2_star = np.conjugate(np.fft.rfft2(img2))
    
    # (iii) compute correlation
    A = _spatial_correlation_fourier(fimg, fimg2_star, fmask, fmask_star)
    return A


def remove_central_corr(A, r=0):
    """ Remove the central part of the correlation, which is peaking too high. The range of data being removed
    can be adjusted with the parameter r.
    """
#     i,j = np.unravel_index(np.argmax(A), A.shape)
    center = utils.get_center(A)
    i,j = center[0], center[1]
    A[i-r:i+1+r, j-r:j+1+r] = np.nan
    return A


def correct_illumination(imgs, roi, kernel_size=5):
    """ Correct the detector images for non-uniform illumination.
    This implementation follows Part II in Duri et al. PHYS. REV. E 72, 051401 (2005).
    
    Args:
        imgs: stack of detector images
        roi: region of interest to consider. Important so that the normalization of the correction 
            is ~unity for the specific roi. Roi should be given as a boolean array (mask).
        kernel_size: size of the kernel for box-average. Can be None, in which case no kernel is 
            applied. The kernel is used to smooth out remaining speckly structure in the intensity correction.
        
    Returns:
        imgs: corrected images, cropped to an extended box aroung the roi
        roi: new roi for the cropped image
        bp: correction factor imgae
    """
    if kernel_size is None:
        extend = 0
    else:
        extend=2*kernel_size
    imgs, roi = utils.box_to_roi_extend(imgs, roi, extend=extend)
    if imgs.ndim==2:
        imgs = imgs[None, ...]
    bp = np.mean(imgs, axis=0)
    if kernel_size is not None:
        kernel = np.ones([kernel_size, kernel_size])/kernel_size/kernel_size
        bp = ndimage.convolve(bp, kernel)
    bp = bp / bp[roi].mean()
    zero = bp==0
    bp[zero] = 1e-6
    imgs_corr = imgs/bp
    return imgs_corr, roi, bp


def correct_illumination_gaussian(imgs, roi, kernel_size=5):
    """ Correct the detector images for non-uniform illumination by applying a
    Gaussina filter.
    
    Returns:
        imgs: corrected images, cropped to an extended box aroung the roi
        roi: new roi for the cropped image
        bp: correction factor image
    """
    extend=kernel_size
    imgs, roi = utils.box_to_roi_extend(imgs, roi, extend=extend)
    if imgs.ndim==2:
        imgs = imgs[None, ...]
    bp = np.mean(imgs, axis=0)
    if kernel_size is not None:
        bp = ndimage.gaussian_filter(bp, sigma=kernel_size/2, mode='mirror')
#     bp = bp / bp[roi].mean()
    zero = bp==0
    bp[zero] = 1e-6
    imgs_corr = imgs/bp
    return imgs_corr, roi, bp


def fit_correlation(A, r=None, ax=None, label='sigma'):
    xx = A.shape[0]//2
    yy = A.shape[1]//2
    if r is None:
        horiz = A[xx,:]
        vert = A[:,yy]
    else:
        horiz = A[xx, yy-r:yy+r]
        vert = A[xx-r:xx+r, yy]
    data = [vert, horiz]
    titles = ['vertical', 'horizontal']
    
    if ax is None:
        fig, ax = plt.subplots(nrows=2, figsize=(10,10))
    
    out = []
    for ii, dat in enumerate(data):
        center = np.where(np.isnan(dat))[0][0]
        
        gmodel = models.GaussianModel(nan_policy='omit')+models.ConstantModel(nan_policy='omit')
        params = gmodel.make_params()
        params['amplitude'].set(value=1.)
        params['sigma'].set(value=1.)
        params['center'].set(value=center)
        params['c'].set(value=1.)

        x = np.arange(dat.shape[0])
        res = gmodel.fit(dat, params, x=x, method='leastsq') #, fit_kws={'ftol':1e-10, 'xtol':1e-10})
        out.append(res)
        
        xfit = np.arange(0,2*center,0.1)
        if label=='sigma':
            ax[ii].plot(xfit, res.eval(x=xfit), color='purple', linewidth=2,
                label='sigma = {:.2f} +/- {:.2f}'.format(res.params['sigma'].value, res.params['sigma'].stderr))
        elif label=='contrast':
            ax[ii].plot(xfit, res.eval(x=xfit), color='purple', linewidth=2,
                label='contrast = {:.2f} +/- {:.2f}'.format(res.params['height'].value, res.params['height'].stderr))
        ax[ii].plot(dat, 'o', color='orange')
        ax[ii].legend(loc='upper right')
        ax[ii].set_title(titles[ii])
        ax[ii].set_xlabel('Pixel')
        ax[ii].set_ylabel('Correlation')
        print(res.fit_report())
        print('\n')

    plt.tight_layout()
    plt.show()
    return out
