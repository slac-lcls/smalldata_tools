import numpy as np
try:
    import pyFAI
except ModuleNotFoundError:
    print("pyFAI not available on LCLS-II env")

from smalldata_tools.DetObjectFunc import DetObjectFunc
from smalldata_tools.utilities import image_from_dxy

class azav_pyfai(DetObjectFunc):
    """
    Parameters
    ----------
    name: str
        Function name
    
    userMask: array, optional
        User defined mask. 1 for valid pixels.

    thres: int, float, optional
        Pixel intensity threshold.
    
    return2d: bool, optional
        Return a 2d (q,phi). Default: False
    
    poni_file: str, Path object, optional
        Path to a pyFAI calibration file
    
    ai_kwargs: dict, optional
        Arguments to pyFAI.AzimuthalIntegrator. Either this parameter or a calib file is necessary
        For arguments see: https://pyfai.readthedocs.io/en/master/api/pyFAI.html#module-pyFAI.azimuthalIntegrator
    
    pol_factor: float, optional
        Polarization factor. Default 1. Passed to integrate1d or integrate2d.
        
    npts_radial: int, optional
        Number of points for the radial binning. Default 256.
    
    npts_az: int, optional
        Number of points for the azimuthal binning. Default 360. Only used for the 2d integration.
        
    azav_kwargs: dict, optional
        Additonal arguments to pass to integrate1d or integrate2d.
        See https://pyfai.readthedocs.io/en/master/api/pyFAI.html#pyFAI.azimuthalIntegrator.AzimuthalIntegrator.integrate1d
    """
    
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','pyfai')
        super(azav_pyfai, self).__init__(**kwargs)
        
        self.mask = kwargs.pop('userMask',None)
        self.threshold = kwargs.pop('thres',None)
        self.return2d = kwargs.pop('return2d','False')
        
        # azimuthal integrator argument
        self.poni_file = kwargs.pop('poni_file',None)
        if self.poni_file is not None:
            self.ai = pyFAI.load(str(self.poni_file))
        else:
            self._ai_kwargs = kwargs.pop('ai_kwargs',None)
            assert self._ai_kwargs is not None, \
                "Need either a calibration file or a set of keywords arguments to instantiate the pyFAI integrator."
            # self.ai = pyFAI.azimuthalIntegrator.AzimuthalIntegrator(**self._ai_kwargs) # why does it not work?
            self.ai = pyFAI.AzimuthalIntegrator(**self._ai_kwargs)
        print(f'Created azimuthal integrator for {self._name}')
        
        # integration arguments
        self.pol_factor = kwargs.pop('polarization_factor',1)
        self.npts = kwargs.pop('npts',256)
        self.npts_az = kwargs.pop('npts_az',360)
        # self.radial_range = kwargs.pop('radial_range',None)
        # self.az_range = kwargs.pop('azimuthal_range',None)
        self._units = kwargs.pop('int_units','2th_deg') # q_A^-1, 2th_deg, r_mm
        self._azav_kwargs = kwargs.pop('azav_kwargs',{}) # additional arguments for integration
        return

    def setFromDet(self, det):
        super(azav_pyfai, self).setFromDet(det)
        # geometry
        self.ix = det.ix
        self.iy = det.iy
        
        if self.poni_file is None: # use calib value if calib file is given
            if len(det.pixelsize)==2:
                if det.pixelsize[0]>1:
                    self.ai.pixel1 = det.pixelsize[0]*1e-6
                    self.ai.pixel2 = det.pixelsize[1]*1e-6
                else:
                    self.ai.pixel1 = det.pixelsize[0]
                    self.ai.pixel2 = det.pixelsize[1]
            else:
                if det.pixelsize[0]>1:
                    self.ai.pixel1 = det.pixelsize[0]*1e-6
                    self.ai.pixel2 = det.pixelsize[0]*1e-6
                else:
                    self.ai.pixel1 = det.pixelsize[0]
                    self.ai.pixel2 = det.pixelsize[0]
        
        # use cmask by default. Note: pyFAI uses 0 for valid pixel.
        if self.mask is not None and det.cmask is not None and self.mask.shape==det.cmask.shape:
            print('Use user mask and cmask')
            self.mask = ~(self.mask.astype(bool)&det.cmask.astype(bool))
        else:
            print('Use default detector masks')
            try:
                self.mask = ~(det.cmask.astype(bool)&det.mask.astype(bool))
            except:
                if det.ped is None:
                    self.mask = None
                else:
                    try:
                        self.mask = ~(np.ones_like(det.ped).astype(bool))
                    except:
                        self.mask = None
        if self.mask.ndim==3:
            self.mask = image_from_dxy(self.mask, self.ix, self.iy)
            
        print(f'Azimuthal integrator:\n{self.ai}')
        return
            
    def setFromFunc(self, func=None):
        return
        
    def process(self, data):
        if self.threshold is not None:
            data[data<threshold] = 0
        if data.ndim==3:
            data = image_from_dxy(data, self.ix, self.iy)
        out = self._process(data)
        return out
    
    def _process(self, data):
        if self.return2d:
            I, q, az = self.ai.integrate2d(data, self.npts, self.npts_az, unit=self._units, 
                                           polarization_factor=self.pol_factor, mask=self.mask, method='cython', 
                                           **self._azav_kwargs)
            return {'azav':I, 'q': q, 'az':az}
        else:
            q, I = self.ai.integrate1d(data, self.npts, unit=self._units, 
                                       polarization_factor=self.pol_factor, mask=self.mask, method='cython', 
                                       **self._azav_kwargs)
            return {'azav':I, 'q': q}
