import numpy as np
# from pathlib import Path

from scipy.signal import savgol_filter
from scipy.stats import pearsonr

from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

"""
General linear algebra comments on subspace projection:
If A = [u_1, ... , u_k] is an orthonormal basis of a subspace U (basis column vector arranged in a matrix), then
the projector is given by P = A(A.T). If the basis A is not orthonormal (as in the case of two-pulse basis 
for example), then the projector is P = A(A.TA)^-1 A.T = A pinv(A).

In this code, the exact definition and nomenclature are not always respected, in order to facilitate computation.
The projector is defined as A.T (or pinv A).
The relation Ax = (x.TA.T).T is also used to flip the matrix multiplication for the same reason.
"""


def get_basis_and_projector(waveforms, n_components=1, n_iter=20):
    """
    Returns the basis vector A, subspace projector and svd of the waveforms.
    
    Remark: although in the single pulse case, A and the projector are simply the transpose of 
    each other, they are still assigned to two different variables, as their relationship is not as 
    straighforward for the multi-pulse case. The WaveformRegressor can thus handle both cases the 
    same way.
    """
    
    """ (i) Perform SVD"""
    svd = TruncatedSVD(n_components=25, n_iter=n_iter)
    svd.fit(waveforms)
    
    """ (ii) Construct projector """
#     The projector is defined as the pseudo inverse of the basis vectors A:
#     projector = np.linalg.pinv(A)
#     However, if the matrix A is orthonormal (which it is here!),then the pseudoinverse becomes:
    projector = svd.components_[:n_components]
    A = projector.transpose()
    return A, projector, svd
        

def multiPulseProjector(singlePulseBasis, n_pulse=1, delay=None, sampling=1, method='pinv', **kwargs):
    """ Build the basis for a multiple waveform analysis, based on single pulse reference basis.
    
    Args:
        singlePulseBasis: single pulse basis A
        n_pulse: number of pulse
        delay: delay between the pulses (if the number of delay is equal to n_pulse-1, then 
            the first pulse is assumed to have dl=0
        sampling: sampling of the waveform. Useful if the delay is given in time units instead of indices
        method: 'pinv', 'QR', 'Ridge'
        kwargs: only for Ridge regressor (see code)
    
    Returns:
        Basis matrix A and projector matrices
    
    The projection is given by:
        P=A.dot(projector).dot(data)
    The coefficients projector onto the subspace A are:
        coeffs=projector.dot(data)
        
    Note: if method 'Ridge' is used, then a RidgeRegressor object is returned instead of a projector (matrix).
    The function fitPulse will take care of handling this difference.
    """
    
    if delay is None:
        if n_pulse==1:
            delay = [0]
        else:
            raise ValueError('No delay given for multipulse basis.')
    if type(delay) is float or type(delay) is int:
        delay = [delay]
    if len(delay)==n_pulse-1:
        delay = np.insert(delay,0,0)
    elif len(delay)==n_pulse:
        delay = np.asarray(delay)
    else:
        raise ValueError('The number of delays given does not match the number of pulses.')
        
    
    """ (i) build the multiplulse basis matrix """
    A0 = singlePulseBasis    
    A = []
    for ii in range(n_pulse):
        A.append(np.roll(A0,int(delay[ii]/sampling), axis=0))
    A = np.concatenate(A, axis=1)
    
    """ (ii) Construct the projector """
    if method=='pinv':
        projector = np.linalg.pinv(A)
        return A, projector
    elif method=='QR':
        Q, R = np.linalg.qr(A)
        projector = np.transpose(np.linalg.inv(A.transpose().dot(Q))).dot(Q)
        return A, projector
    elif method=='Ridge':
        if 'alpha' in kwargs:
            alpha = kwargs.pop('alpha')
        else:
            alpha=0
        projector = Ridge(alpha=alpha, fit_intercept=False) # is a RidgeRegressor instance, not a matrix
        return A, projector
    else:
        raise NameError('Method not implemented')


def construct_waveformRegressor(X_ref, n_components=1, n_pulse=1, delay=None, **kwargs):
    """ Construct waveform regressor based on a set of reference waveforms.
    
    Args:
        X_ref: reference waveform
        n_components: nubmer of SVD components to use for the fit
        n_pulse: number of pulse to fit in the waveform
        **kwargs: see function multiPulseProjector. If n_pulse>1, a kwarg 'delay' is mandatory.
    """
    A, projector, svd = get_basis_and_projector(X_ref, n_components=n_components)
    A, projector = multiPulseProjector(A, n_pulse=n_pulse, delay=delay, **kwargs)
    return WaveformRegressor(A=A, projector=projector, n_pulse=n_pulse)



class WaveformRegressor(BaseEstimator, RegressorMixin):
    """ Regressor compatible with sk-learn package (or not really...) """
    def __init__(self, A=None, projector=None, n_pulse=1, roi=None):
        """
        Args:
            A: Basis vectors of the subspace in matrix form (column vectors)
                projector: projector on the subspace A
        
        Remarks:
            Because the basis A can be built artificially from non-orthogonal vectors, its projector is not necessarily 
            trivial (as in simply A.T), hence it is calculated separately.
            
            Also, in order to facilitate the fit of multiple waveforms at once, both matrices are transposed. The projection 
            and reconstruction are calculated thus as coeffs=X.dot(T) and coeffs.dot(A) respectively.
            
            The regressor is not fully compatible with sklearn unfortunately. This is because the basis and the projector 
            are external parameters that the package does not know how to handle properly.
        
        Construct basis A and projector using the function 'get_basis_projector' or 'construct_2PulseProjector'
        """
        self.A = A.T # transpose to allow fit of many waveforms at once, following sklearn convention
        self.projector = projector.T # transpose to allow fit of many data at once
        self.n_pulse_ = n_pulse
        if roi is None:
            self.roi=[0, 1e6]
    
    
    def fit(self, X, y=None):
        """ y=None for sklearn compatibility reason """
        # Check validity of input X:
        if X.shape[-1] != self.A.shape[1]:
            print('Data and projector shapes dont match.')
            self.coeffs_ = np.zeros(self.A.shape[1])
            return self
        
        if isinstance(self.projector, Ridge):
            ridge = self.projector.fit(self.A, X)
            coeffs = ridge.coef_
        else:
#             coeffs = self.projector.dot(X.T) # old way
            coeffs = X.dot(self.projector)
        
        if len(X.shape)==1:
            self.coeffs_ = coeffs[None,:]
        else:
            self.coeffs_ = coeffs
        return self
    
    
    def reconstruct(self):
        try:
            getattr(self,"coeffs_")
        except AttributeError:
            raise RuntimeError("You must fit the waveform before reconstructing it!")
            
#         reconstructed = self.A.dot(self.coeffs_) # old way
        reconstructed = self.coeffs_.dot(self.A)
        return reconstructed.real
    
    
    def predict(self, X=None):
        """ Just for sklearn compatibility reason """
        return self.reconstruct()
    
    
    def score(self, X):
        """ Returns the r2 score of the projected waveforms (one score value per waveform)
        Must have called fit(X) or fit_reconstruct(X) before.
        """
        return r2_score(X.T, self.reconstruct().T, multioutput='raw_values')
            # .T: hack so that the score of multiple waveforms is computed correctly
    
    
    def pearsonr_coeff(self, X):
        """ Pearson correlation between the fit and the waveform X
        Must have called fit(X) or fit_reconstruct(X) before.
        """
        if len(X.shape)==1:
            X = X[None,:]
        return np.asarray([pearsonr(xx, rr)[0] for xx, rr in zip(X, self.reconstruct())])
    
    
    def fit_reconstruct(self, X, return_score=False):
        self.fit(X)
        if return_score:
            return self.reconstruct(), self.score(X)
        return self.reconstruct()
    
    
    def get_pulse_intensity(self, X, mode='norm'):
        """
        Inputs:
            - waveform X
            - mode: 'norm' or 'max'. Return the norm of fitted coefficients or the max of each pulse
                For multipulse basis it is recommended to use the 'max' method, as the basis may not be orthogonal.
        Ouputs:
            - intensities: individual pulse intensities
        """
        self.fit(X)
        nCoeff = int(self.coeffs_.shape[1]/self.n_pulse_)
        intensities = np.zeros((self.coeffs_.shape[0],self.n_pulse_))
        if mode=='both':
            intensities_max = np.zeros((self.coeffs_.shape[0],self.n_pulse_))
        for ii in range(self.n_pulse_):
            coeffs = self.coeffs_[:,ii*nCoeff:(ii+1)*nCoeff]
            if mode=='norm':
                intensities[:,ii] = np.linalg.norm(coeffs, axis=1)
            elif mode=='max':
                reconstructed_single = coeffs.dot(self.A[ii*nCoeff:(ii+1)*nCoeff,:])
                intensities[:,ii] = np.max(reconstructed_single, axis=1)
            elif mode=='both':
                intensities[:,ii] = np.linalg.norm(coeffs, axis=1)
                reconstructed_single = coeffs.dot(self.A[ii*nCoeff:(ii+1)*nCoeff,:])
                intensities_max[:,ii] = np.max(reconstructed_single, axis=1)
        
        if mode=='both':
            return intensities, intensities_max
        return intensities

    
    
    
    
    
    
    
    
    
    
""" 
%%%%%%%%%%%%%%%%%%%%% OLD FUNCTIONS %%%%%%%%%%%%%%%%%%%%% 
May not work properly anymore.
Kept for potential backward compatibility on older scripts.
"""

def construct_2PulseProjector(singlePulseBasis, delay=None, sampling=.125, method='pinv', **kwargs):
    """
    Gives the projector onto the subspace mapped by the chosen single-pulse SVD components for a two-pulse waveform
    Inputs:
        singlePulseBasis: single pulse basis A
        delay: delay between the two pulses
        nCoeff: number of single pulse basis vectors to take
        method: 'pinv', 'QR', 'Ridge'
    Returns:
        Basis matrix A and projector function
        The projector is given by:
            P=A.dot(projector).dot(data)
        The coefficients projector onto the subspace A are:
            coeffs=projector.dot(data)
        
        Note: if method 'Ridge' is used, then a RidgeRegressor object is returned instead of a projector (matrix).
        The function fitPulse will take care of handling this difference.
    """
    
    if delay is None:
        raise ValueError('Delay is None, give it a value!')
        
    
    """ (i) build the basis matrix """
    A0 = singlePulseBasis
    A1 = A0
    A2 = np.roll(A0,int(delay/sampling),axis=0)
    A = np.append(A1,A2,axis=1)
    
    """ (ii) Construct the projector """
    if method=='pinv':
        projector = np.linalg.pinv(A)
        return A, projector
    elif method=='QR':
        Q, R = np.linalg.qr(A)
        projector = np.transpose(np.linalg.inv(A.transpose().dot(Q))).dot(Q)
        return A, projector
    elif method=='Ridge':
        if 'alpha' in kwargs:
            alpha = kwargs.pop('alpha')
        else:
            alpha=0
        projector = Ridge(alpha=alpha, fit_intercept=False) # is a RidgeRegressor instance, not a matrix
        return A, projector
    else:
        raise NameError('Method not implemented')


def construct_waveformRegressor_old(X_ref, n_components=1, n_pulse=1, **kwargs):
    """ 
    Construct waveform regressor based on a set of reference waveforms.
    """
    A, projector, svd = get_basis_and_projector(X_ref, n_components=n_components)
    if n_pulse==2:
        A, projector = construct_2PulseProjector(A, **kwargs)
    return WaveformRegressor(A=A, projector=projector, n_pulse=n_pulse)