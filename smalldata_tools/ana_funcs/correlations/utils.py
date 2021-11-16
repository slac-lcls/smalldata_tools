import numpy as np

def gauss(x, area, center, width):
    return abs(area)*np.exp(-np.power(x - center, 2.) / (2 * np.power(width, 2.))) / (width * np.sqrt(2 * np.pi))


def lorentzian(x, a, x0, gam):
    return a * gam**2 / ( gam**2 + (x-x0)**2)


def polarCoord(img_size, center):
    x, y = np.arange(img_size[1]), np.arange(img_size[0])
    xx, yy = np.meshgrid(x,y)
    cx = center[0]
    cy = center[1]
    
    rad = np.sqrt((xx-cx)**2+(yy-cy)**2)
    phi = np.rad2deg(-np.arctan2(yy-cy,xx-cx))
    return rad, phi


def polarROI(rad, phi, ra, rb, phia, phib):
    """ Polar ROI
    
    Args:
        rad, phi: output of polarCoord
        ra, rb: radial range
        phia, phib: azimuthal range
    """
    if phib<0 and phib<phia:
        ROI = np.logical_and((ra<=rad)*(rad<=rb), (phia<=phi)*(phi<=180))
        ROI = ROI + np.logical_and((ra<=rad)*(rad<=rb), (phi<=phib))
    else:
        ROI = np.logical_and((ra<=rad)*(rad<=rb), (phia<=phi)*(phi<=phib))
    return ROI


def box_to_roi(imgs, mask):
    """ Reduce the imgs to the smallest box around the roi defined by mask.
    
    Args:
        imgs: stack of images or single image (ndim= 2 or 3)
        mask: mask for imgs. Shape must correspond to the shape of the images
    Return:
        Cropped image and mask
    """
    if imgs.ndim==2:
        imgs = imgs[np.newaxis,...]
    nframe = imgs.shape[0]
    posx,posy = np.where(mask)
    shape = np.array([posx.max()-posx.min()+1,posy.max()-posy.min()+1])
    
    maskn = np.zeros(shape)
    maskn[posx-posx.min(),posy-posy.min()] = 1.
#     imgs_red = np.zeros([nframe, *maskn.shape])
    imgs_red = np.zeros([nframe, maskn.shape[0], maskn.shape[1]])
    imgs_red[:, posx-posx.min(),posy-posy.min()] = imgs[:,posx,posy].copy()
    return np.squeeze(imgs_red), maskn


def box_to_roi_extend(imgs, mask, extend=10):
    """ Reduce the imgs to an extended box around the roi defined by mask. Does not
    apply the mask.
    
    Args:
        imgs: stack of images or single image (ndim= 2 or 3)
        mask: mask for imgs. Shape must correspond to the shape of the images
        extend: number of pixel for the extension in each direction
    Return:
        Cropped image and mask
    """
    if imgs.ndim==2:
        imgs = imgs[np.newaxis,...]
    posx,posy = np.where(mask)
    posx_min = posx.min()-extend
    posx_max = posx.max()+extend
    posy_min = posy.min()-extend
    posy_max = posy.max()+extend
    
    maskn = mask[posx_min:posx_max, posy_min:posy_max]
    imgs_red = imgs[:,posx_min:posx_max, posy_min:posy_max]
    return np.squeeze(imgs_red), maskn


def get_center(img):
    """ Returns the center coordinates of img. """
    shape = img.shape
    center = (img.shape[0]//2,img.shape[1]//2)
    return center