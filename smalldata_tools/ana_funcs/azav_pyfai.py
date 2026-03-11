import sys
import numpy as np

try:
    import pyFAI
except ModuleNotFoundError:
    print("pyFAI not available on LCLS-II env")

### Can we put LCLSGeom in environment?
sys.path.append("/sdf/home/l/lconreux/LCLSGeom")
from LCLSGeom.manager import get_geometry
from LCLSGeom.converter import PsanaToPyFAI

from smalldata_tools.common.detector_base import DetObjectFunc

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
        self._name = kwargs.get("name", "pyfai")
        super(azav_pyfai, self).__init__(**kwargs)
        self.mask = kwargs.pop("userMask", None)
        self.threshold = kwargs.pop("thres", None)
        self.return2d = kwargs.pop("return2d", "False")
        self.poni_file = kwargs.pop("poni_file", None)
        self._ai_kwargs = kwargs.pop("ai_kwargs", None)
        if self._ai_kwargs is None and self.poni_file is None:
            raise ValueError("Either a poni file or ai_kwargs must be provided to initialize the azimuthal integrator.")

        # integration arguments
        self.pol_factor = kwargs.pop("polarization_factor", -1)        
        self.npts = kwargs.pop("npts", 256)
        self.npts_az = kwargs.pop("npts_az", 360)
        # self.radial_range = kwargs.pop('radial_range',None)
        # self.az_range = kwargs.pop('azimuthal_range',None)
        self._units = kwargs.pop("int_units", "2th_deg")  # q_A^-1, 2th_deg, r_mm
        self._azav_kwargs = kwargs.pop(
            "azav_kwargs", {}
        )  # additional arguments for integration
        return

    def setFromDet(self, det):
        super(azav_pyfai, self).setFromDet(det)
        if (
            self.mask is not None
            and det.cmask is not None
            and self.mask.shape == det.cmask.shape
        ):
            print("Use user mask and cmask")
            self.mask = ~(self.mask.astype(bool) & det.cmask.astype(bool))
        else:
            print("Use default detector masks")
            try:
                self.mask = ~(det.cmask.astype(bool) & det.mask.astype(bool))
            except:
                if det.ped is None:
                    self.mask = None
                else:
                    try:
                        self.mask = ~(np.ones_like(det.ped).astype(bool))
                    except:
                        self.mask = None
        if self.mask is not None and self.mask.ndim == 3:
            self.mask = np.reshape(self.mask, (self.mask.shape[0]*self.mask.shape[1], self.mask.shape[2]))
      
        metrology = get_geometry(det.det.alias)
        # instantiate PyFAI detector first 
        # this adds the detector name mentioned in the poni file to pyFAI detector registry
        # reading the poni file will not throw an error if the detector name is not in the registry
        detector = PsanaToPyFAI.convert(metrology)
        # if detector is 3D (panels, slow, fast), then this scripts creates a fake 2D detector by
        # squeezing the panel and slow dimensions.
        # the 3D pixel layout is passed using the metrology data !
        if self._ai_kwargs is None:
            ai = pyFAI.load(self.poni_file)
            self._ai_kwargs = ai.getPyFAI()

        if self._units == "q_A^-1":
            if "wavelength" not in self._ai_kwargs:
                raise ValueError("Wavelength must be provided when using q_A^-1 units. Check poni file or ai_kwargs inputs.")
            self.ai = pyFAI.AzimuthalIntegrator(
                detector=detector,
                dist=self._ai_kwargs["dist"],
                poni1=self._ai_kwargs["poni1"],
                poni2=self._ai_kwargs["poni2"],
                rot1=self._ai_kwargs["rot1"],
                rot2=self._ai_kwargs["rot2"],
                rot3=self._ai_kwargs["rot3"],
                wavelength=self._ai_kwargs["wavelength"],
            )
        else:
            self.ai = pyFAI.AzimuthalIntegrator(
                detector=detector,
                dist=self._ai_kwargs["dist"],
                poni1=self._ai_kwargs["poni1"],
                poni2=self._ai_kwargs["poni2"],
                rot1=self._ai_kwargs["rot1"],
                rot2=self._ai_kwargs["rot2"],
                rot3=self._ai_kwargs["rot3"],
            )
        # we should have code that extract setup parameters
        #    in case people mess w/ the input
        # we should also extra q/theta/... bins for plotting later.
        azcfg = self.ai.get_config()
        dcfg = azcfg.pop("detector_config")
        ocfg = dcfg.pop("orientation")
        for k, v in azcfg.items():
            setattr(self, "azcfg_%s" % k, v)
        for k, v in dcfg.items():
            setattr(self, "azcfg_det_%s" % k, v)
        setattr(self, "azcfg_detorientation", ocfg.value)
        setattr(self, "azcfg_qvals", self.ai.get_qa())
        # setattr(self,'azcfg_qvals',self.ai.get_qa())
        print(f"Azimuthal integrator:\n{self.ai}")
        return

    def setFromFunc(self, func=None):
        return

    def process(self, data):
        if self.threshold is not None:
            data[data < self.threshold] = 0
        if data.ndim == 3:
            # instead of assembling the image, now we just squeeze the panel and slow dimensions
            # so that image and detector have same shape.
            data = np.reshape(data, (data.shape[0]*data.shape[1], data.shape[2]))
        out = self._process(data)
        return out

    def _process(self, data):
        if self.return2d:
            I, q, az = self.ai.integrate2d(
                data,
                self.npts,
                self.npts_az,
                unit=self._units,
                polarization_factor=self.pol_factor,
                mask=self.mask,
                method="cython",
                **self._azav_kwargs,
            )
            return {"azav": I, "q": q, "az": az}
        else:
            q, I = self.ai.integrate1d(
                data,
                self.npts,
                unit=self._units,
                polarization_factor=self.pol_factor,
                mask=self.mask,
                method="cython",
                **self._azav_kwargs,
            )
            return {"azav": I, "q": q}
