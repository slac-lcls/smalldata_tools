import numpy as np
import time
import h5py
from scipy.ndimage import rotate
from scipy.ndimage import gaussian_filter1d
import scipy.optimize
from smalldata_tools.common.detector_base import DetObjectFunc

def Gaussian1D(amplitude, center,sigma,bkgC,bkgSlope):
    """ Returns a 1D Gaussian function with input parameters. """
    Xc=center  #Center
    #Now lets define the 1D gaussian function
    def Gauss1D(x) :
        """ Returns the values of the defined 1d gaussian at x """
        return amplitude*np.exp(-(((x-Xc)/sigma)**2)/2) +bkgC+bkgSlope*x

    return Gauss1D


class mecttFunc(DetObjectFunc):
    """
    Parameters
    ----------
    background Run : int (default = None)
         run to take the whitefield from
    roi_n : bounds for normalizing region
    """

    def __init__(self, **kwargs):
        self._name = kwargs.get("name", "mectt")
        super(mecttFunc, self).__init__(**kwargs)
        self.bkgfile = kwargs.get("bkgfile", None)
        self.bkgdet = kwargs.get("bkgdet", "alvium_tt_calib")
        self.roi_n = kwargs.get("roi_n", None)
        self.roi_s = kwargs.get("roi_s", None)
        self.tilt = kwargs.get("tilt")
        self.sigma_TT = kwargs.get("sigma_TT",18)
        self.sm = kwargs.get("sm")
        self.ttRaw = kwargs.get("ttRaw", False)
        self._debug = kwargs.get("debug", False)
                        
        if isinstance(self.bkgfile, str):
            try:
                fbkg = h5py.File(self.bkgfile, 'r')
                self.background = np.array(fbkg['Sums'][self.bkgdet])
            except:
                print(self._name,' could not get background data')
                return
        
        if self.roi_n.shape[0] == 2 and len(self.roi_n.shape) == 2:
            self.bkg_n = self.background[
                self.roi_n[0, 0] : self.roi_n[0, 1],
                self.roi_n[1, 0] : self.roi_n[1, 1]
            ]
        else:
            print(self._name,'background shape not as expected!')
        
    def tt_ana(self, data):
        ret_dict = {}
        if self.ttRaw:
            ret_dict['ttRaw']=data
        
        sig_n = data[
            self.roi_n[0, 0] : self.roi_n[0, 1],
            self.roi_n[1, 0] : self.roi_n[1, 1]
        ]
        if self._debug:
            ret_dict['sig_n']=sig_n
        
        norm_ratio = np.mean(self.bkg_n)/np.mean(sig_n)

        tt_im = data/self.background * norm_ratio
        #bckg = np.average(wf_map[roi_n[0]:roi_n[1],roi_n[2]:roi_n[3]])/np.average(run_map[roi_n[0]:roi_n[1],roi_n[2]:roi_n[3]])
        #tt_im = run_map/wf_map * bckg
        tt_im_rot = rotate(tt_im, angle=self.tilt, reshape ='false')
        
        sig = tt_im_rot[
            self.roi_s[0, 0] : self.roi_s[0, 1],
            self.roi_s[1, 0] : self.roi_s[1, 1]
        ]
        if self._debug:
            ret_dict['norm_ratio']=norm_ratio
            ret_dict['tt_im']=tt_im
            ret_dict['tt_im_rot']=tt_im_rot
            ret_dict['sig']=sig
            
        sig_pj = np.sum(sig, axis=1)
        ret_dict['sig_pj']=sig_pj

        if (np.max(sig_pj)-np.min(sig_pj) < 0.03):# | (np.average(sig) > 1):
            time_px = np.nan
            width_px = np.nan
            err =-1
            bkgC = np.nan
            bkgS = np.nan
        else:
            grad = np.gradient(gaussian_filter1d (sig_pj, self.sm))
            if self._debug:
                ret_dict['grad']=grad

            [sol, err] = self.FitGauss1D(-grad)

    # return amplitude*np.exp(-(((x-Xc)/sigma)**2)/2) +bkgC+bkgSlope*x
    # sol = [amplitude,center,sigma,bkgC,bkgSlope],[success]
    # err =1,2,3,4 if success fit
    #   print(sol)
            ampl_px = sol[0]# peak amplitude
            time_px = np.round(sol[1])# peak position gaussian fit
            width_px = np.round(sol[2])# sigma gaussian fit
            bkgC = np.round(sol[3])# background constant
            bkgS = np.round(sol[4])# background slopw
        ret_dict['time_px']=time_px
        ret_dict['ampl_px']=ampl_px
        ret_dict['width_px']=width_px
        ret_dict['err']=err
        ret_dict['bkg_const']=bkgC
        ret_dict['bkg_slope']=bkgS
        return ret_dict
            
    def process(self, data):
        ret_dict = self.tt_ana(data)
        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict[f"{k}_{kk}"] = subfuncResults[k][kk]
        return ret_dict

    def moments1D(self,inpData):
        """ Returns the (amplitude, center,sigma, bkgC, bkgSlope) estimated from moments in the 1d input array inpData """

        bkgC=inpData[0]  #Taking first point as intercept and fitting a straght line for slope
        bkgSlope=(inpData[-1]-inpData[0])/(len(inpData)-1.0)
    
        Data=np.ma.masked_less(inpData-(bkgC+bkgSlope*np.arange(len(inpData))),0)
        #Removing the background for calculating moments of pure gaussian
        #We also masked any negative values before measuring moments

        amplitude=Data.max()

        total=float(Data.sum())
        Xcoords=np.arange(Data.shape[0])

        center= np.argmax(Data)
        #(Xcoords*Data).sum()/total
    
        sigma = self.sigma_TT
        #sigma=np.sqrt(np.ma.sum((Data*(Xcoords-center)**2))/total)

        return amplitude,center,sigma,bkgC,bkgSlope

    def FitGauss1D(self,Data,ip=None):
        """ Fits 1D gaussian to Data with optional Initial conditions ip=(amplitude, center, sigma, bkgC, bkgSlope)
        Example: 
        >>> X=np.arange(40,dtype=np.float)
        >>> Data=np.exp(-(((X-25)/5)**2)/2) +1+X*0.5
        >>> FitGauss1D(Data)
        (array([  1. ,  25. ,   5. ,   1. ,   0.5]), 2)
        """
        if ip is None:   #Estimate the initial parameters from moments 
            ip=self.moments1D(Data)
        
        def errfun(ip):
            return np.ravel(Gaussian1D(*ip)(*np.indices(Data.shape)) - Data)

        p, success = scipy.optimize.leastsq(errfun, ip)

        return p,success

