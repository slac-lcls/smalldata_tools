import numpy as np
import scipy 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ---> SHOULD BE IN PYPSALGO
def gauss(x, mean, sigma, height=1.0, pedestal=0.0) :
    if sigma==0:
        print('sigma is 0: ',sigma,' ',mean)
        return 1e32
    return  (height * np.exp(-0.5 * ((x-mean)/sigma)**2)) + pedestal



def create_test_gaussian(x,noise_amplitude=1.0) :

    # create a gaussian for testing fiting routine
    mean     = np.random.random() * 511.0
    sigma    = np.random.random() * 5.0
    height   = np.random.random() * 100.0
    pedestal = np.random.random() * 10.0

    gauss_data = gauss(x,mean,sigma,height,pedestal) + \
        (np.random.normal(0.0,1.0,512)*noise_amplitude)

    return gauss_data,mean,sigma,height,pedestal


# ---> SHOULD BE IN PYPSALGO
def FWHM(data) :
    """
    Start from tallest bin and then go up in increasing x until bin
    content is half.  Then repeat in oppostite direction.  Take the
    average of the two as measure of FWHM.
    
    """

    # Find peak value and its index
    peakIndex = np.argmax(data)
    peakValue = data[peakIndex]
    
    # Use mean of lowest 10% of data to estimate pedestal
    lowestTenPercent = int(0.1 * data.size)
    lowestValue = np.sort(data)[:lowestTenPercent]
    pedestal = np.mean(lowestValue)

    # Value of height
    height = peakValue - pedestal
    
    # Split array to values before and after peak
    # --> lowerHalf array is reversed, so index is distance from peak
    lowerHalf = (data[:peakIndex])[::-1]
    upperHalf = data[peakIndex:]
    
    # Lower half FWHM
    #- Find largest index that is <= height/2.0 +pedestal
    lowerFWHMIndex = np.argwhere(lowerHalf <= (height/2.0 + pedestal)).flatten()
    lowerFWHM = lowerFWHMIndex[0] if lowerFWHMIndex.size>0 else None
    
    # Upper half FWHM
    #- Find smallest index that is <= height/2.0 +pedestal
    upperFWHMIndex = np.argwhere(upperHalf <= (height/2.0 + pedestal)).flatten()
    upperFWHM = upperFWHMIndex[0] if upperFWHMIndex.size>0 else None
        
    # Calculate FWHM, taking into account if lower or upper FWHM were
    # contained within data
    fwhm = None

    # - no upper or lower FWHM found ==> take entire width as estimate
    if (lowerFWHM is None) and (upperFWHM is None) :
        fwhm = data.size
    
    # - Only lower FWHM was found
    if (lowerFWHM is not None) and (upperFWHM is None) :
        fwhm = lowerFWHM

    # - Only upper FWHM was found
    if (lowerFWHM is None) and (upperFWHM is not None) :
        fwhm = upperFWHM

    # - Lower and Upper FWHM found - take average of the two
    if (lowerFWHM is not None) and (upperFWHM is not None) :
        fwhm = (lowerFWHM + upperFWHM) / 2.0

    return fwhm



def gauss_params_estimate(data) :
    """
    Estimate of gaussian parameters without fitting
    """
    # Use index of largest value in data to estimate mean
    mean = np.argmax(data)
    
    # Use mean of lowest 10% of data to estimate pedestal
    lowestTenPercent = int(0.1 * data.size)   
    lowestValue = np.sort(data)[:lowestTenPercent]
    pedestal = np.mean(lowestValue)
    
    # Use difference of largest value in data and pedestal to estimate
    # height
    height = data[mean] - pedestal
    
    #print 'DEBUG esimate: ',data[mean], pedestal, height
    # Estimate sigma from FWHM, which is 2.36*sigma
    sigma = FWHM(data) / 2.36

    return mean,sigma,height,pedestal

 

def GaussFit(data,xaxis=None, debug=False) :
    """
    General interface to fit data to gaussian
    """

    if xaxis is None:
        # use array index of data to define x-axis
        xaxis = np.arange(data.size)

    params_estimate = gauss_params_estimate(data)
    if (params_estimate[1]==0):
        print('params_estimate failed, will return nans.', params_estimate[0],' ', params_estimate[1])
        fit_results = {'mean':np.nan,
                       'sigma':np.nan,
                       'height':np.nan,
                       'pedestal':np.nan,
                       'mean_error':np.nan,
                       'sigma_error':np.nan,
                       'height_error':np.nan,
                       'pedestal_error':np.nan,
                       'mean_estimate':np.nan,
                       'sigma_estimate':np.nan,
                       'height_estimate':np.nan,
                       'pedestal_estimate':np.nan}
        return fit_results

    try:
        fit_params, fit_cov = curve_fit(gauss,xaxis,data,params_estimate)
        #print('ESTIMATE ', params_estimate)
        #print('fit_params ', fit_params)
        #print('fit_cov ', fit_cov)
        try:
            dCov = np.diag(fit_cov)
            fit_params_error = np.sqrt(dCov)
        except:
            if debug:
                print('np diag failed')
            fit_params_error = [-1.,-1.,-1., -1.]
    except:
        if debug:
            print('curve fit failed failed')
        fit_params = [-1.,-1.,-1., -1.]
        fit_params_error = [-1.,-1.,-1., -1.]        

    fit_results = {'mean':float(fit_params[0]),
                   'sigma':float(fit_params[1]),
                   'height':float(fit_params[2]),
                   'pedestal':float(fit_params[3]),
                   'mean_error':float(fit_params_error[0]),
                   'sigma_error':float(fit_params_error[1]),
                   'height_error':float(fit_params_error[2]),
                   'pedestal_error':float(fit_params_error[3]),
                   'mean_estimate':float(params_estimate[0]),
                   'sigma_estimate':float(params_estimate[1]),
                   'height_estimate':float(params_estimate[2]),
                   'pedestal_estimate':float(params_estimate[3])}

    return fit_results
    

def fitPeaks(peakDat, peakWindow, means, sigmas,heights, showPlot=False):
    debug=False
    lowEdge=0
    highEdge=peakWindow[0]

    while highEdge<peakDat.size and lowEdge<highEdge-peakWindow[0]*0.25:
        if debug:
            print('lowEdge ',lowEdge,' high ',highEdge)
        fitAr = GaussFit(peakDat[lowEdge:highEdge])
        if showPlot:
            plt.plot(peakDat[lowEdge:highEdge])
        if fitAr is not None:
           if fitAr['sigma']<0:
               if (debug):
                   print('Fit sigma is zero!', fitAr)
               lowEdge+=peakWindow[0]*0.5
               highEdge=min(peakDat.size,lowEdge+peakWindow[0])
           #before appending, check if peakis at or beyond range, if so, add 1/2 window and redo
           elif fitAr['mean']>(peakWindow[0]-2*fitAr['sigma']):
               lowEdge+=(peakWindow[0]*0.5)
               highEdge=min(peakDat.size,lowEdge+peakWindow[0])
           #if res is close to a previous one, add 1/2 window and redo
           elif len(means)>0 and np.abs((fitAr['mean']+lowEdge)-means[-1])<(2*fitAr['sigma']):
               lowEdge+=(peakWindow[0]*0.5)
               highEdge=min(peakDat.size,lowEdge+peakWindow[0])
           else: #we have a fit with a decent fitting window here.
               means.append(fitAr['mean']+lowEdge)
               sigmas.append(fitAr['sigma'])
               heights.append(fitAr['height'])
               if showPlot:
                   plt.plot([fitAr['mean'],fitAr['mean']],[fitAr['pedestal']*0.9,(fitAr['height']+fitAr['pedestal'])*1.05],'r')
               if debug:
                   print('mean: ',fitAr['mean'],' ',fitAr['height'])
               #check: if peak found, make sure that low edge of next fit is at least peak+2 sigma
               if (fitAr['mean']+2*fitAr['sigma']) < (peakWindow[0]-peakWindow[1]):
                   lowEdge+=peakWindow[0]-peakWindow[1]
               else:
                   lowEdge+=fitAr['mean']+2*fitAr['sigma']
                   highEdge=min(peakDat.size,lowEdge+peakWindow[0])
        else: #just advance the edges by half a window
            if (debug):
                print('Issues in GaussFit: returns None.')
            lowEdge+=peakWindow[0]*0.5
        highEdge=min(peakDat.size,lowEdge+peakWindow[0])
        if showPlot:
            plt.show()
    #print 'Done With Fits: ',len(means) 


if __name__ == "__main__" :

    xaxis = np.arange(512)
    noise_amplitude = 1.0

    params = np.empty(4)
    data,params[0],params[1],params[2],params[3] \
        = create_test_gaussian(xaxis,noise_amplitude)

    realgauss = gauss(xaxis, *params)
    
    params_estimate = gauss_params_estimate(data)
    crude_fit = gauss(xaxis,*params_estimate)
        
    print("Initial:",params)
    print("Estimate:",params_estimate)

    print("Fitting")
    fit_params, fit_cov = curve_fit(gauss,xaxis,data,params_estimate)
    fit_params_error = np.sqrt(np.diag(fit_cov))

    print("Final Fit(errors):",zip(fit_params,fit_params_error))
    print("Covariance Matrix:",fit_cov)
    gauss_fit = gauss(xaxis, *fit_params)    

    print("Parameter \t Fit \t Error \t Diff \t Status")
    for p,fit,fit_error in zip(params,fit_params,fit_params_error) :
        diff = abs(p-fit)
        fitStatus = "GOOD" if diff < fit_error else "BAD"
        print(p,"\t",fit,"\t",fit_error,"\t",diff,"\t",fitStatus)
        
    print("Plotting")
    plt.figure(1)
    plt.clf()
    plt.plot(data,'.', realgauss,'p',crude_fit,'.',gauss_fit,'-')
    plt.draw()
    plt.show()
    
