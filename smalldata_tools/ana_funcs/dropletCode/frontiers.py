# return boundary sums for droplets.
import numpy as np


def frontiers(pp0, pp1, npixs, img):
    pos = np.append(0, np.cumsum(npixs))
    fadus = np.zeros(len(npixs))
    for k in xrange(len(npixs)):
        # make a sub-image of for droplet
        i = pp0[pos[k] : pos[k + 1]]
        j = pp1[pos[k] : pos[k + 1]]
        mi = np.min(i)
        mj = np.min(j)
        ii = 1 + i - mi
        jj = 1 + j - mj
        n, m = (np.max(i) - mi + 1, np.max(j) - mj + 1)
        img1 = np.zeros((n + 2, m + 2))
        img1[ii + 1, jj] = 1
        img1[ii, jj + 1] = 1
        img1[ii - 1, jj] = 1
        img1[ii, jj - 1] = 1
        img1[ii, jj] = 0
        fadus[k] = np.sum(img1 * img[mi - 1 : mi + n + 1, mj - 1 : mj + m + 1])
    return fadus
