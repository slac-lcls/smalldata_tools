import numpy as np
import itertools
from skimage import feature
from numba import jit
from scipy import sparse
from scipy import optimize
from skimage.draw import circle_perimeter
from skimage.measure import (CircleModel, ransac)
from scipy.signal import argrelextrema
from scipy.spatial import cKDTree

def applyCanny(ar, mask, sigma=1, thres=0.9, thresH=0.995):
    """
    use canny to find edges in image to be used in hough algorithm 
    """      
    if thres==-999:
        arThres = feature.canny(ar, sigma=sigma,mask=mask)
    else:
        arThres = feature.canny(ar, sigma=sigma,mask=mask, low_threshold=thres, high_threshold=thresH, use_quantiles=True)
    return arThres,sparse.coo_matrix(arThres)

def fitCircle(x,y,yerr=None, guess=None):
    """
    fit a single circle. Transform input to lists to that fitCircles can be used.
    """  
    x = [x]
    y = [y]
    rGuess = (np.nanmax(x)-np.nanmin(x)+np.nanmax(y)-np.nanmin(y))/4. #largest differences/2/2
    r = [rGuess]
    fitRes = fitCircles(x,y,r,yerr=None, guess=None)
    #have only one circle.
    fitRes['R']=fitRes['R'][0]
    fitRes['residu']=fitRes['residu'][0]
    return fitRes

def fitCircles(x,y,r,yerr=None, guess=None):
  """
  simultanous least squares fitting of multiple concentric rings.
  """  
  def calc_R(x,y, par):
      xc,yc = par
      return np.sqrt((x-xc)**2 + (y-yc)**2)

  def f(par,x,y):
      Ri = calc_R(x,y,par)
      return Ri - Ri.mean()

  def f_global(par,Mat_r,Mat_x,Mat_y):
      err = []
      for r,x,y in zip( Mat_r, Mat_x, Mat_y):
          errLocal = f(par,x,y)
          err = np.concatenate((err, errLocal))          
      return err

  if guess is None or len(guess)!=2:
      x_m = np.mean(x[0])
      y_m = np.mean(y[0])
  else:
      x_m = guess[0]
      y_m = guess[1]
      
  center_estimate = x_m, y_m
  fitRes={}  
  if yerr is not None:      
      center, C, info, msg, success  = optimize.leastsq(f_global, center_estimate, args=(r,x,y), full_output=True)
      fitRes['C'] = C
      fitRes['info'] = info
      fitRes['msg'] = msg
      fitRes['success'] = success
  else:
      center, ier = optimize.leastsq(f_global, center_estimate, args=(r,x,y))
      #print(center)
      fitRes['ier'] = ier
  xc, yc = center
  fitRes['xCen'] = xc
  fitRes['yCen'] = yc
  ##now this will be a list.
  Rs=[]
  Resids=[]
  for thisr, thisx, thisy in zip(r,x,y):
      Ri     = calc_R(thisx, thisy, center)
      Rs.append(Ri.mean())
      Resids.append(np.sum((Ri - Ri.mean())**2))
  fitRes['residu'] = Resids
  fitRes['R']      = Rs
  return fitRes

# --------------------------------------------------------------            
#
# Hough Transform implementation.
#
# --------------------------------------------------------------            

@jit(forceobj=True)
def addToHough(x,y,arHough, hough_radii, center_x, center_y, wt=1):
    """
    add a single point in x-y space to hough space 
    """
    dr=hough_radii[1]-hough_radii[0]
    for icx,cx in enumerate(center_x):
        dx = (x-cx)*(x-cx)
        if dx < hough_radii[0]*hough_radii[0] or  dx > hough_radii[-1]*hough_radii[-1]:
            continue
        for icy,cy in enumerate(center_y):
            dy = (y-cy)*(y-cy)
            r = np.sqrt(dx+dy)
            ir = int((r-hough_radii[0])/dr)
            if ir>=0 and ir < hough_radii.shape[0]:
                arHough[ir, icx, icy]+=wt
#do this to force compilation
_ = addToHough(0,0,np.zeros([2,2,2]),np.arange(2), np.arange(2), np.arange(2))

@jit(forceobj=True)
def transformImage(arSparse,arHough, hough_radii, center_x, center_y):
    """
    transform a sparsified imaged to an array in hough space
    """
    assert arHough.shape[0]==hough_radii.shape[0]
    assert arHough.shape[1]==center_x.shape[0]
    assert arHough.shape[2]==center_y.shape[0]
    for trow,tcol,tdat in itertools.izip(arSparse.row, arSparse.col,arSparse.data):
        addToHough(trow,tcol, arHough,  hough_radii, center_x, center_y, tdat)
#do this to force compilation
_img = np.zeros((10, 10), dtype=np.uint8)
_arHough = np.zeros((10, 10, 10), dtype=np.uint8)
_rr, _cc = circle_perimeter(4, 4, 3)
_img[_rr, _cc] = 1
_imgSparse = sparse.coo_matrix(_img)
transformImage(_imgSparse,_arHough, np.arange(10), np.arange(10), np.arange(10))

def findCenter(arSparse, rBound, xcenBound, ycenBound, nbin=100, retHoughArray=False, nBinR=None):
    """
    find the beam center via hough transform, given a sparsified imaged 
    input: arSparse, rBound, xcenBound, ycenBound, nbin=100, retHoughArray=False, nBinR=None):
           arSparse  (sparsified image)
           rBound    ([r_min, r_max])
           xcenBound ([x_center_ min, x_center_max])
           ycenBound ([y_center_ min, y_center_max])
           nbin      (number of bins for center coordinates in hough space - default 100)
           nBinR     (number of bins for radius in hough space - defaults to nbin)
           retHoughArray (return the full filled array in hough space, other return projections, 
                          defaults to False)
    """
    if nBinR is None:
        nBinR = nbin
    radii = np.arange(rBound[0], rBound[1],(rBound[1]-rBound[0])/nBinR)
    centerx = np.arange(xcenBound[0],xcenBound[1],(xcenBound[1]-xcenBound[0])/nbin)
    centery = np.arange(ycenBound[0],ycenBound[1],(ycenBound[1]-ycenBound[0])/nbin)
    arHough = np.zeros([radii.shape[0], centerx.shape[0], centery.shape[0]])
    transformImage(arSparse, arHough, radii, centerx, centery)
    maxdim0 = [ arHough[i,:,:].max() for i in range(arHough.shape[0])]
    maxdim1 = [ arHough[:,i,:].max() for i in range(arHough.shape[1])]
    maxdim2 = [ arHough[:,:,i].max() for i in range(arHough.shape[2])]
    fitRes={}
    fitRes['R'] = radii[np.array(maxdim0).argmax()]
    fitRes['xCen'] = centerx[np.array(maxdim1).argmax()]
    fitRes['yCen'] = centery[np.array(maxdim2).argmax()]
    fitRes['radii'] = radii
    fitRes['centerx'] = centerx
    fitRes['centery'] = centery
    if retHoughArray:
        fitRes['houghArray'] = arHough
    else:
        fitRes['houghArray_projR'] = maxdim0
        fitRes['houghArray_projX'] = maxdim1
        fitRes['houghArray_projY'] = maxdim2
    return fitRes

def iterateCenter(arSparse, ar_shape, rRange, nbin=100, prec=1, redFac=5., retHoughArray=False, printProgress=False, overfillFactor=1, nBinR=None):
    """
    use hough method to find the center iteratively ro keep computing time dowb compared to a 
    very large 3-d hough space.
    input: arSparse, ar_shape, rRange, nbin=100, prec=1, redFac=5., retHoughArray=False, printProgress=False, overfillFactor=1, nBinR=Non
           arSparse (sparsified image)
           ar_shape (shape of the original image in #pixel)
           rRange   ()
           nbin     (number of bins for center coordinates in hough space - default 100)
    """
    xRange=[ar_shape[0]*(1.-overfillFactor),ar_shape[0]*overfillFactor]
    yRange=[ar_shape[1]*(1.-overfillFactor),ar_shape[1]*overfillFactor]
    fitRes = findCenter(arSparse, rRange, xRange, yRange, nbin=nbin, retHoughArray=retHoughArray, nBinR=nBinR)
    while (fitRes['centerx'][1]-fitRes['centerx'][0]) > prec:
        maxR = fitRes['R']
        maxX = fitRes['xCen']
        maxY = fitRes['yCen']
        rSizeX = (xRange[1]-xRange[0])/redFac
        rSizeY = (yRange[1]-yRange[0])/redFac
        xRange=[maxX-rSizeX*0.5, maxX+rSizeX*0.5]
        yRange=[maxY-rSizeY*0.5, maxY+rSizeY*0.5]
        if printProgress:
            print('finding the center in ranges: ',rRange, xRange, yRange)
        fitRes = findCenter(arSparse, rRange, xRange, yRange, nbin=nbin, retHoughArray=retHoughArray, nBinR=nBinR)
        if printProgress:
            print('found center: ',maxX, maxY, maxR)

    temp = fitRes['xCen']
    fitRes['xCen'] = fitRes['yCen']
    fitRes['yCen'] = temp
    return fitRes

# --------------------------------------------------------------            
#
# select maxima
#
# --------------------------------------------------------------            

def getMaxR(fitRes, norm=False, minDr=-1):
    """
    sort the radii of the circle found in hough transform step.
    arguments: fitRes, (opt: norm, minDr)
               fitRes: result of iterateCenter (as dictionary)
               norm: normalize the height of the peak in hough space by the radius
               minDr: request circles to be apart by at least minDr pixels
    """
    try:
        radii = fitRes['radii']
        rRes = fitRes['houghArray_projR']
    except:
        print('the passed fit Result does not have the necessary keys')
        return[]

    if norm:
        maxR = argrelextrema(np.array(rRes)/radii, np.greater)
    else:
        maxR = argrelextrema(np.array(rRes), np.greater)
    print('maxR: ',maxR)
    maxRadii = radii[maxR[0]]
    if maxRadii.shape[0]==0:
        print('no maxima could be found!')
        return []
    resAtRadii = np.array(rRes)[maxR[0]]
    resAtRadii,maxRadii = (list(x) for x in zip(*sorted(zip( resAtRadii, maxRadii))))

    #print('DEBUG: ',maxRadii)
    #print('DEBUG: ',resAtRadii)
    
    nrad=[]
    #loop inverse!
    for rad in reversed(maxRadii):
        thisminDr=1e6
        for irrad in nrad:
            if np.abs(rad-irrad)<thisminDr:
                thisminDr = np.fabs(rad-irrad)
        if thisminDr>minDr:
            nrad.append(rad)

    #return maxRadii
    return nrad

def findPointsInRing(arSparse, center, maxR1, ringInfo, inParams={}):
    """
    create list of points for each ring. 
    use cKDTree first to select points in 'donot' region arund input radius, 
         then use RANSAC to throw out points that do not fit with CircleModel
    parameters: arSparse, maxR1, (opt params={})
        arSparse: sparified image
        center: beam center found in Hough step  [xCenter, yCenter]
        maxR1: list of round maximu in Hough step
        ringInfo: list of rings
        params:
            nMaxRing: maximum number of rings to use in final fit, def 6
            deltaR: require new circle to be at least deltaR pixels away from last circle, def 5
            minPoints: minimum absolute number of points in circle to be considered for final fit, def 40
            minInFrac #require frav of points to pass RANSAC, def 0.45 (45%)
            RANSAC_residual_threshold: allowed max residual to consider point being 'in' (def 2)
            RANSAC_min_sample: minimu number of samples to be drawn (def 10)
    """
    #use cKDTree to select points in ring of width 2*deltaR around circle
    print('select points in circles using cKDTree and select good candidates for final fit using RANSAC')
    tree = cKDTree(np.array([arSparse.col, arSparse.row]).T)
    for ir, r in enumerate(maxR1):
        if len(ringInfo)>=inParams['nMaxRing']:
            break
        thisRingInfo={}
        rOuter = tree.query_ball_point(center, r+inParams['deltaR'])
        rInner = tree.query_ball_point(center, r-inParams['deltaR'])
        inRing = set(rOuter).symmetric_difference(set(rInner))
        thisRingInfo['rInput']=r
        thisRingInfo['pointsInCircle']=list(inRing)
        print('input radius, #points ',r,len(list(inRing)))
        if len(list(inRing))<inParams['minPoints']:
            continue
        #use ransac to pick points most likely to belong to circle
        model_robust, inliers = ransac(np.array([arSparse.row[list(inRing)], arSparse.col[list(inRing)]]).T, CircleModel, min_samples=inParams['RANSAC_min_samples'], residual_threshold=inParams['RANSAC_residual_threshold'], max_trials=1000)
        thisRingInfo['inFrac']=inliers.astype(int).sum()/float(len(list(inRing)))
        thisRingInfo['inliers']=inliers
        thisRingInfo['ransac_result']=model_robust
        print('ring %d, infrac %f, #points %d' %(ir, thisRingInfo['inFrac'], len(thisRingInfo['inliers'])))
        if thisRingInfo['inFrac']>inParams['minInFrac'] and thisRingInfo['inliers'].astype(int).sum()>=inParams['minPoints']:
            ringInfo.append(thisRingInfo)

# --------------------------------------------------------------            
#
# main function.
#
# --------------------------------------------------------------            

def FindFitCenter(image, mask, inParams={}):
    """
    main function to find the beam center
    arguments: image, mask (opt: inParams)
               image to use to find beam center (2-d)
               mask to be used (also 2-d)
               inParams: dictionary of parameters used in findcing the ebam center
    inParams:
             #edge finding:
             sigma          #gaussian blurring for canny algorithm, def 1
             low_threshold  #threshold for only strongest features, using quantiles, def 0.95
             high_threshold #threshold for only strongest features, using quantiles, def 0.0995
             #parameters for hough center finding
             precision #bin size in pixel for center in hough array, def 1
             nBinR     #number of bins for center, def 280
             overFac   #allow beam center to be out of the image by a factor of x, def 1.5 (50%)
             #parameters to find best rings
             norm      #normalize number of points/circle by its radius (point density), def True
             deltaR    #require new circle to be at least deltaR pixels away from last circle, def 5
             nMaxRing  #max number of rings to consider, def 6
             minInFrac #require 45% of points to pass RANSAC, def 0.45 
             minPoints #minimum absolute number of points in circle to be considered for final fit, def 40
             RANSAC_residual_threshold #allowed max residual to consider point being 'in' (def 2)
             RANSAC_min_sample         #minimum number of samples to be drawn (def 10)
    """
    ar = image
    if len(ar.shape)!=2 or len(mask.shape)!=2:
        print('image or mask are not 2-d arrays, cannot find centers')
        return

    #parameters found to be typcially right
    #parameters for canny
    params={'sigma':1}      
    params['low_threshold']=0.95  
    params['high_threshold']=0.995
    #parameters for hough center finding
    params['precision']=1  
    params['nBinR']=280    
    params['overFac']=1.5   
    #parameters to find best rings
    params['norm']=True    
    params['deltaR']=5     
    params['nMaxRing'] = 6 
    params['minPoints']=40 
    params['RANSAC_min_samples']=10 
    params['RANSAC_residual_threshold']=2 
    params['minInFrac']=0.45 

    for inkey in inParams:
        try:
            params[inkey]=inParams[inkey]
        except:
            print('fit parameters do not have key ',inkey,', available are: ',params.keys())
            pass

    print('now find edges')
    arThres,arSparse = applyCanny(ar, mask.astype(bool), sigma=params['sigma'], thres=params['low_threshold'], thresH=params['high_threshold'])

    print('use hough transform to find center & ring candidates')
    res = iterateCenter(arSparse, [arThres.shape[0], arThres.shape[1]], [1,1401], prec=params['precision'], overfillFactor=params['overFac'], nBinR=params['nBinR'])
    maxR1 =  getMaxR(res, norm=params['norm'], minDr=params['deltaR'])
    center = [ res['xCen'], res['yCen']]
    print('maxR1: ',maxR1)

    ringInfo=[]
    if len(maxR1) == 0:
        return -1, ringInfo, arSparse

    findPointsInRing(arSparse, center, maxR1, ringInfo, params)
 
    print('now prepare for the final fit')
    allX=[]
    allY=[]
    allR=[]
    for ir,thisRingInfo in enumerate(ringInfo):
        allR.append(ir)
        allX.append(arSparse.col[thisRingInfo['pointsInCircle']])
        allY.append(arSparse.row[thisRingInfo['pointsInCircle']])
    try:
        combRes = fitCircles(allX,allY,allR,yerr=True)
        return combRes, ringInfo, arSparse
    except:
        print('combined fit failed')
        return -1, ringInfo, arSparse
