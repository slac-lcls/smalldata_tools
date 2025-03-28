# from pathlib import Path
import h5py as h5
import glob
import os
import numpy as np
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

from smalldata_tools.ana_funcs.svd_waveform import svd_waveform_processing as proc
from smalldata_tools.common.detector_base import DetObjectFunc
from smalldata_tools.utilities import printR

from mpi4py import MPI

rank = MPI.COMM_WORLD.Get_rank()


class SvdFit(DetObjectFunc):
    """Performs fit of waveforms using singular value decomposition. The main utility is to get the
    intensity of multiple pulse in a single waveform, from a basis of single-pulse reference waveforms.
    Check the file make_waveform_basis.py for information on how to create the basis.
    See svd_waveform_processing.py for the underlying processing algorithm.
    """

    def __init__(self, **kwargs):
        """Set delays and hyperparameters for the regressor
        Args:
            name (str): DetObjectName, default: svdFit
            n_components: number of component to build the basis (max. 25)
            n_pulse (int): number of pulses to fit. Default: 1
            delay (list or array): delay between each pulse
            sampling (float): sampling of the waveform. Default: 1
            basis_file (str or Path object): if not given will take the latest one in the detector calib folder
            mode (str): 'max', 'norm' or 'both', method to calculate the pulse amplitudes. Default: 'max'
            return_reconstructed (bool): return reconstructed waveforms or not
        """
        self._name = kwargs.get("name", "svdFit")
        super().__init__(**kwargs)
        self.n_components = kwargs.get("n_components", 2)
        self.n_pulse = kwargs.get("n_pulse", 1)
        self.delay = kwargs.get("delay", [0])
        if isinstance(self.delay, int) or isinstance(self.delay, float):
            self.delay = [self.delay]
        self.sampling = kwargs.get("sampling", 1)
        self.basis_file = kwargs.get("basis_file", None)
        self._mode = kwargs.get("mode", "max")
        self._return_reconstructed = kwargs.get("return_reconstructed", False)
        self.channel = None

    def setFromDet(self, det):
        """Load basis, projector and other useful variables from calib h5 file"""
        super().setFromDet(det)
        # calibDir = det.det.env.calibDir()
        if self.basis_file is None:
            try:
                # Automatically find the latest waveform basis file in calibDir
                if rank == 0:
                    logger.info("No basis file given, default to latest")
                exp = det.run.ds.exp

                BASE = os.environ.get("SIT_PSDM_DATA", "/sdf/data/lcls/ds/")
                basis_files = Path(BASE) / f"{exp[:3]}/{exp}/hdf5/smalldata/svd_basis"
                basis_files = basis_files.glob("wave_basis_" + det.det.alias + "*.h5")
                self.basis_file = max(basis_files, key=os.path.getctime)
            except:
                logger.error("Unable to find a basis file. Return now")
                return

        if not os.path.isfile(self.basis_file):
            if rank == 0:
                logger.error("No basis file found. Exit")
            return
        else:
            if rank == 0:
                logger.info(
                    "{}: Basis file found at {}".format(self._name, self.basis_file)
                )

        with h5.File(self.basis_file, "r") as f:
            components = f["components"][()]
            self.roi = f["roi"][()]
            if self.roi is None:
                self.roi = [0, 1e6]
            self.bkg_idx = f["background_roi"][()]
            if "channel" in f.keys():
                self.channel = f["channel"][()]

        if self.channel is not None:
            self._name += f"_ch{self.channel}"

        projector = components[
            : self.n_components
        ]  # see svd_waveform_processing for details
        A = projector.transpose()
        A, proj = proc.multiPulseProjector(
            A, n_pulse=self.n_pulse, delay=self.delay, sampling=self.sampling
        )
        self.regressor = proc.WaveformRegressor(
            A=A, projector=proj, n_pulse=self.n_pulse
        )

    def process(self, waveform):
        """
        Fit waveform and output dictionary with coefficient, intensity and score
        """
        if self.channel is not None:
            waveform = waveform[self.channel]
        if waveform.ndim == 1:
            waveform = waveform[None, :]
        if self.bkg_idx is not None:
            waveform = waveform - np.mean(waveform[:, : self.bkg_idx], axis=1)
        waveform = waveform[:, self.roi[0] : self.roi[1]]
        intensities = self.regressor.get_pulse_intensity(waveform, mode=self._mode)
        score = self.regressor.score(waveform)
        output = {
            "intensities": np.squeeze(intensities),
            "score": np.squeeze(score),
            "coefficients": np.squeeze(self.regressor.coeffs_),
        }
        if self._return_reconstructed:
            output["reconstructed"] = np.squeeze(self.regressor.reconstruct())
        return output

    def process_with_alignment(self, waveform):
        """To be implemented"""
        return None
