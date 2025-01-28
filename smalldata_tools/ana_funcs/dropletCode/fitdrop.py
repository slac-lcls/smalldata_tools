# add a photon to an image
import numpy as np
from scipy.optimize import leastsq


def placephots(pkx, pky, sz, aduspphot):
    timg = np.zeros(sz)
    # wid=0.5  #tail size, note 1.0 is worse
    # a,b=np.indices(sz)
    # nn=aduspphot/10.0 #tail intensity
    # split center
    pkx = np.maximum(0.01, np.minimum(sz[0] - 1.01, pkx))
    pky = np.maximum(0, np.minimum(sz[1] - 1.01, pky))
    ii = np.int32(pkx)
    jj = np.int32(pky)

    pkx = pkx - ii
    pky = pky - jj
    for t in range(len(ii)):  # since positions could overlap
        timg[ii[t], jj[t]] += aduspphot * (1 - pkx[t]) * (1 - pky[t])
        timg[ii[t] + 1, jj[t]] += aduspphot * pkx[t] * (1 - pky[t])
        timg[ii[t], jj[t] + 1] += aduspphot * (1 - pkx[t]) * pky[t]
        timg[ii[t] + 1, jj[t] + 1] += aduspphot * pkx[t] * pky[t]
    # add tails as this helps convergence from farther away
    # timg += nn /(1.0+((a-pkx[t])/wid)**2+((b-pky[t])/wid)**2)
    return timg


# function for leastsq
def photres(p, y, err, aduspphot):
    l = len(p) // 2
    return np.ravel(y - placephots(p[0:l], p[l:], y.shape, aduspphot)) / err


# fit a droplet
def fitdrop(posi, posj, img, sigi, aduspphot, prflag=None):
    # fit to discrete model
    chimin = 1e6
    l = len(posi)
    reps = 50
    chigood = 0.5
    sz = 0.75
    p0 = np.append(posi, posj)
    if prflag != None:
        svv = np.zeros((reps, 4 * l + 2))
    for i in range(reps):  # try various random starting points
        d0 = 2 * sz * (np.random.rand(2 * len(posi)) - 0.5)
        plsq = leastsq(photres, p0 + d0, args=(img, sigi, aduspphot), full_output=1)
        chisq = np.sum(plsq[2]["fvec"] ** 2)
        if chisq < chimin:
            chimin = chisq
            pmin = plsq[0]
            # mini=plsq[2]['nfev'] #no function calls
            mini = i
        if prflag != None:
            svv[i, 0] = chisq
            svv[i, 1] = plsq[2]["nfev"]
            svv[i, 2 : 2 + 2 * l] = p0
            svv[i, 2 + 2 * l :] = plsq[0]
            print(i, svv[i])
        if chimin < chigood:
            break  # stop as is good enough
    if prflag != None:
        return (chimin, pmin[0:l], pmin[l:], mini, svv)
    else:
        return (chimin, pmin[0:l], pmin[l:], mini)
