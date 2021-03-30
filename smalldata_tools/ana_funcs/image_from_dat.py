import numpy as np
from smalldata_tools.DetObject import DetObjectFunc
from smalldata_tools.utilities import image_from_dxy

class image_from_dat(DetObjectFunc):
    """ Return the whole image """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','image_from_dat')
        super(image_from_dat, self).__init__(**kwargs)
        
    def setFromDet(self, det):
        # geometry
        self.ix = det.ix
        self.iy = det.iy
        
    def process(self, dat):
        image = image_from_dxy(dat, self.ix, self.iy)
        output = {'image': np.asarray(image)}
        return output