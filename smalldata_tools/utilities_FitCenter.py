import numpy as np
from scipy import sparse
import itertools
from skimage import feature
from numba import jit
from scipy import sparse
from scipy import optimize
from skimage.draw import circle_perimeter
from skimage.measure import (CircleModel, ransac)
from scipy.signal import argrelextrema
from scipy.spatial import cKDTree

#
# stuff for new beam center fit
#

def applyCanny(ar, mask, sigma=1, thres=95):
    if thres==-999:
        arThres = feature.canny(ar, sigma=1,mask=mask)
    else:
        arThres = feature.canny(ar, sigma=1,mask=mask, low_threshold=thres)
    return arThres,sparse.coo_matrix(arThres)

def fitCircles(x,y,r,yerr=None, guess=None):
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
      print center
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
        
@jit
def addToHough(x,y,arHough, hough_radii, center_x, center_y, wt=1):
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

@jit
def transformImage(arSparse,arHough, hough_radii, center_x, center_y):
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
            print 'finding the center in ranges: ',rRange, xRange, yRange
        fitRes = findCenter(arSparse, rRange, xRange, yRange, nbin=nbin, retHoughArray=retHoughArray, nBinR=nBinR)
        if printProgress:
            print 'found center: ',maxX, maxY, maxR

    temp = fitRes['xCen']
    fitRes['xCen'] = fitRes['yCen']
    fitRes['yCen'] = temp
    return fitRes

def getMaxR(fitRes, norm=False, minDr=-1):
    try:
        radii = fitRes['radii']
        rRes = fitRes['houghArray_projR']
    except:
        print 'the passed fit Result does not have the necessary keys'
        return[]

    if norm:
        maxR = argrelextrema(np.array(rRes)/radii, np.greater)
    else:
        maxR = argrelextrema(np.array(rRes), np.greater)
    maxRadii = radii[maxR[0]]
    resAtRadii = np.array(rRes)[maxR[0]]
    resAtRadii,maxRadii = (list(x) for x in zip(*sorted(zip( resAtRadii, maxRadii))))

    #print 'DEBUG: ',maxRadii
    #print 'DEBUG: ',resAtRadii
    
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


def FindFitCenter(image, mask, inParams={}):
    ar = image
    if len(ar.shape)!=2 or len(mask.shape)!=2:
        print 'image of mask are not 2-d arrays, cannot find centers'
        return

    #parameters found to be typcially right
    #parameters for canny
    params={'sigma':1}       #gaussian blurring for canny algorithm
    params['threshold']=95   #threshold for only strongest features
    #parameters for hough center finding
    params['precision']=1   #bin size in pixel for center in hough array
    params['nBinR']=280     #number of bins for center
    params['overFac']=1.5   #allow beam center to be out of the image by a factor of 50%
    #parameters to find best rings
    params['norm']=True    #normalize number of points/circle by its radius (point density)
    params['deltaR']=5     #require new circle to be at least deltaR pixels away from last circle
    params['nMaxRing'] = 6 #max number of rings to consider.
    params['minInFrac']=0.45 #require 45% of points to pass RANSAC 
    params['minPoints']=40 #minimum absolute number of points in circle to be considered for final fit. 
                           #Should maybe be a fraction of highest number.
    for inkey in inParams:
        try:
            params[inkey]=inParams[inkey]
        except:
            print 'fit parameter do not have key ',inkey,', available are: ',params.keys()
            pass

    print 'now find edges'
    arThres,arSparse = applyCanny(ar, mask.astype(bool), sigma=params['sigma'], thres=params['threshold'])

    print 'use hough transform to find center & ring candidates'
    res = iterateCenter(arSparse, [arThres.shape[0], arThres.shape[1]], [1,1401], prec=params['precision'], overfillFactor=params['overFac'], nBinR=params['nBinR'])
    maxR1 =  getMaxR(res, norm=params['norm'], minDr=params['deltaR'])
    print 'maxR1: ',maxR1

    pointsInCircles=[]
    pointsInCirclesFitted=[]
    fitResultsRansac=[]
    residuals=[]
    ringInfo=[]
    plotDetail=False
    #use cKDTree to select points in ring of width 2*deltaR around circle
    print 'select points in circles using cKDTree and select good candidates for final fit using RANSAC'
    tree = cKDTree(np.array([arSparse.col, arSparse.row]).T)
    for ir, r in enumerate(maxR1):
        if len(ringInfo)>=params['nMaxRing']:
            break
        thisRingInfo={}
        rOuter = tree.query_ball_point([res['xCen'],res['yCen']], r+params['deltaR'])
        rInner = tree.query_ball_point([res['xCen'],res['yCen']], r-params['deltaR'])
        inRing = set(rOuter).symmetric_difference(set(rInner))
        thisRingInfo['rInput']=r
        thisRingInfo['pointsInCircle']=list(inRing)
        print 'input radius, #points ',r,len(list(inRing))
        if len(list(inRing))<params['minPoints']:
            continue
        #use ransac to pick points most likely to belong to circle
        model_robust, inliers = ransac(np.array([arSparse.row[list(inRing)], arSparse.col[list(inRing)]]).T, CircleModel, min_samples=10, residual_threshold=2, max_trials=1000)
        thisRingInfo['inFrac']=inliers.astype(int).sum()/float(len(list(inRing)))
        thisRingInfo['inliers']=inliers
        thisRingInfo['ransac_result']=model_robust
        print 'ring %d, infrac %f, #points %d' %(ir, thisRingInfo['inFrac'], len(thisRingInfo['inliers']))
        if thisRingInfo['inFrac']>params['minInFrac'] and thisRingInfo['inliers'].astype(int).sum()>=params['minPoints']:
            ringInfo.append(thisRingInfo)
 
    print 'now prepare for the final fit'
    allX=[]
    allY=[]
    allR=[]
    for ir,thisRingInfo in enumerate(ringInfo):
        allR.append(ir)
        allX.append(arSparse.col[thisRingInfo['pointsInCircle']])
        allY.append(arSparse.row[thisRingInfo['pointsInCircle']])
    combRes = fitCircles(allX,allY,allR,yerr=True)
    return combRes, ringInfo, arSparse


