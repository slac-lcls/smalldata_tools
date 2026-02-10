import sys
import numpy as np

try:
    import pyFAI
except ModuleNotFoundError:
    print("pyFAI not available on LCLS-II env")

### Can LCLSGeom be passed through LUTE or put in environment?
sys.path.append("/sdf/group/lcls/ds/tools/LCLSGeom")
from LCLSGeom.psana2.converter import PsanaToPyFAI

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
        # Hardcoded metrology file for now, but I can get it from the jungfrau psana Detector object, just unsure on how to do it !
        # That allows us to map pixel indices to physical coordinates, without assembling and casting the image in a wrong 2D plane.
        self.geometry_path = "/sdf/home/l/lconreux/geom/mfx100852324/geometry-def-jungfrau16M.data"
        converter = PsanaToPyFAI(self.geometry_path)
        pyFAI_detector = converter.detector
        self.poni_file = kwargs.pop("poni_file", None)
        if self.poni_file is not None:
            ai = pyFAI.load(str(self.poni_file))
            dico = ai.getPyFAI()
            self.ai = pyFAI.AzimuthalIntegrator(
                detector=pyFAI_detector,
                dist=dico.pop("dist"),
                poni1=dico.pop("poni1"),
                poni2=dico.pop("poni2"),
                rot1=dico.pop("rot1"),
                rot2=dico.pop("rot2"),
                rot3=dico.pop("rot3"),
                wavelength=dico.pop("wavelength"),
            )
        else:
            self._ai_kwargs = kwargs.pop("ai_kwargs", None)
            if self._ai_kwargs is not None:
                self.ai = pyFAI.AzimuthalIntegrator(
                    detector=pyFAI_detector,
                    **self._ai_kwargs
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

        # integration arguments
        self.pol_factor = kwargs.pop("polarization_factor", 1)
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
        print(f"Azimuthal integrator:\n{self.ai}")
        return

    def setFromFunc(self, func=None):
        return

    def process(self, data):
        if self.threshold is not None:
            data[data < self.threshold] = 0
        if data.ndim == 3:
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