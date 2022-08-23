# from pathlib import Path
import h5py as h5
import glob
import os
import numpy as np

from smalldata_tools.ana_funcs import svd_waveform_processing as proc
from smalldata_tools.DetObjectFunc import DetObjectFunc


class svdFit(DetObjectFunc):
    """ Performs fit of waveforms using singular value decomposition. The main utility is to get the 
    intensity of multiple pulse in a single waveform, from a basis of single-pulse reference waveforms.
    Check the file make_waveform_basis.py for information on how to create the basis.
    See svd_waveform_processing.py for the underlying processing algorithm.
    """
    def __init__(self, **kwargs):
        """ Set delays and hyperparameters for the regressor
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
        self._name = kwargs.get('name','svdFit')
        super(svdFit, self).__init__(**kwargs)
        self.n_components = kwargs.get('n_components',2)
        self.n_pulse = kwargs.get('n_pulse',1)
        self.delay = kwargs.get('delay',[0])
        if isinstance(self.delay, int) or isinstance(self.delay, float):
            self.delay = [self.delay]
        self.sampling = kwargs.get('sampling', 1)
        self.basis_file = kwargs.get('basis_file', None)
        self._mode = kwargs.get('mode', 'max')
        self._return_reconstructed = kwargs.get('return_reconstructed', False)
            
    
    def setFromDet(self, det):
        """ Load basis, projector and other useful variables from calib h5 file """
        super(svdFit, self).setFromDet(det)
        calibDir = det.det.env.calibDir()
        if self.basis_file is None:
            try:
                # Automatically find the latest waveform basis file in calibDir
#                 basis_files = glob.glob('./wave_basis_'+det.det.alias+'*.h5') # dir for local test
#                 basis_files = (Path(det.det.env.calibDir()).parent / 'hdf5/basis').glob('wave_basis_'+det.det.alias+'*.h5')
                basis_files = glob.glob('/'.join(det.det.env.calibDir().split('/')[:-1]) 
                                        + '/hdf5/basis/' + 'wave_basis_'+det.det.alias+'*.h5') # hack if pathlib not available
                self.basis_file = max(basis_files, key=os.path.getctime)
                l = '/'.join(det.det.env.calibDir().split('/')[:-1])+'/hdf5'
            except:
                pass
#         if not Path(self.basis_file).is_file():
        if not os.path.isfile(self.basis_file):
            print("No basis file found. Exit")
            return
        else:
            print('{}: basis file found at {}'.format(self._name, self.basis_file))
        
        with h5.File(self.basis_file,'r') as f:
            components = f['components'][()]
            self.roi = f['roi'][()]
            if self.roi is None:
                self.roi=[0, 1e6]
            self.bkg_idx = f['background_index'][()]
            self.channel = f['channel'][()]
        
        projector = components[:self.n_components] # see svd_waveform_processing for details
        A = projector.transpose()
        A, proj = proc.multiPulseProjector(
                A,
                n_pulse=self.n_pulse,
                delay=self.delay,
                sampling=self.sampling
            )
        self.regressor = proc.WaveformRegressor(A=A, projector=proj, n_pulse=self.n_pulse)
    
    
    def process(self, waveform):
        """
        Fit waveform and output dictionary with coefficient, intensity and score
        """
        if self.channel is not None:
            waveform = waveform[self.channel]
        if waveform.ndim==1:
            waveform = waveform[None,:]
        if self.bkg_idx is not None:
            waveform = waveform - np.mean(waveform[:,:self.bkg_idx], axis=1)
        waveform = waveform[:,self.roi[0]:self.roi[1]]
        intensities = self.regressor.get_pulse_intensity(waveform, mode=self._mode)
        score = self.regressor.score(waveform)
        output = {
            'intensities': np.squeeze(intensities),
            'score': np.squeeze(score),
            'coefficients': np.squeeze(self.regressor.coeffs_)
        }
        if self._return_reconstructed:
            output['reconstructed'] = np.squeeze(self.regressor.reconstruct())
        return output
    
    
    def process_with_alignment(self, waveform):
        """ To be implemented
        """
        return None
