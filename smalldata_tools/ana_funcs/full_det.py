import numpy as np
from skimage.measure import block_reduce
from smalldata_tools.DetObject import DetObjectFunc
from smalldata_tools.utilities import image_from_dxy

BIN_FCT = np.mean

class image_from_dat(DetObjectFunc):
    """ Smalldata analysis function to retrieve whole detector iamges
    Possibility to apply filters and bin image to make it smaller
    """
    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        name: str, optional
            Name of the function
        thresADU : None, float, optional
            Single pixel threshold value
        binned : None, int, float, 2d tuple, list or np.ndarray, optional
            How to bin the pixel.
        """
        self._name = kwargs.get('name','image_from_dat')
        super(image_from_dat, self).__init__(**kwargs)
        self.thresh = kwargs.get('thresADU', None)
        binned = kwargs.get('binned', None)
        if isinstance(binned,int) or isinstance(binned,float):
            self.binned = (int(binned), int(binned))
        elif isinstance(binned, tuple) or isinstance(binned, list), isinstance(binned, np.ndarray):
            if len(binned)==1:
                self.binned = (int(binned[0]), int(binned[0]))
            elif len(binned)==2:
                self.binned = (int(binned[0]), int(binned[1]))
        else:
            raise TypeError('Binned argument not valid. Int or tuple / array of length 2 only.')
        
    def setFromDet(self, det):
        # geometry
        self.ix = det.ix
        self.iy = det.iy
        
    def process(self, dat):
        image = image_from_dxy(dat, self.ix, self.iy)
        if self.thresh is not None:
            image[image<thresh] = 0
        if self.binned is not None:
            image = block_reduce(image, self.binned, BIN_FCT)
        output = {'image': np.asarray(image)}
        return output
    
    
    
