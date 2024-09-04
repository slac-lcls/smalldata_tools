import numpy as np
import time as time
from smalldata_tools.ana_funcs.dropletCode.droplet import dropletfind, dropletanal


def convert_img(img, thres, photpts, mask=None):
    # setup arrays
    idhold = np.zeros(1 + np.size(img))
    # lg = 0
    bins = np.linspace(0, 2000, 1001)
    if mask is not None:
        img *= mask
    npeaks, dimg = dropletfind((img > thres).astype(int))
    dimghold = dimg.copy()
    if npeaks != 0:
        npix, xcen, ycen, xsig, ysig, xysig, adus, idlist = dropletanal(
            img.astype(int), dimg, npeaks
        )
        h, b = np.histogram(adus, bins=bins)
        # check idlist is sorted
        # if np.any(idlist[:-1]>=idlist[1:]): print("not sorted ",fr)

        # handle 1 photon droplets
        w = np.where((adus > photpts[1]) & (adus <= photpts[2]))[0]
        n = len(w)
        onephots = np.zeros((n, 5))
        # onephots[:,0] = fr
        onephots[:, 1] = xcen[w]
        onephots[:, 2] = ycen[w]
        onephots[:, 3] = adus[w]
        onephots[:, 4] = npix[w]
        nophots = n

        # more than single photon droplets
        # find ids of >1 photon drops

        # lg+=np.sum(adus>photpts[-1])
        w = np.where((adus > photpts[2]) & (adus <= photpts[-1]))[0]
        idhold[idlist[w]] = idlist[w]
        pp = np.where(idhold[dimghold])
        idhold[idlist[w]] = 0

        n = len(w)
        twos = np.zeros((n, 5))
        # twos[:,0] = fr
        twos[:, 1] = xcen[w]
        twos[:, 2] = ycen[w]
        twos[:, 3] = adus[w]
        twos[:, 4] = npix[w]

        n1 = len(pp[0])
        pixtwos = np.zeros((n1, 3), dtype=np.int16)
        ss = np.argsort(dimghold[pp[0], pp[1]])  # sort by droplet id
        pixtwos[:, 0] = pp[0][ss]
        pixtwos[:, 1] = pp[1][ss]
        pixtwos[:, 2] = img[pp[0][ss], pp[1][ss]]
        # no2pix = n1
        # notwos = n
    return onephots, twos, pixtwos, h, b
