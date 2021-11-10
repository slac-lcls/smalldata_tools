#finding photon positions without/with fit, code from Mark Sutton
#fitting is generally fast with intensity level 1e-3
from smalldata_tools.ana_funcs.dropletCode.greedyguess import greedyguess
from smalldata_tools.ana_funcs.dropletCode.fitdrop import *
import numpy as np


def loopdrops(onephots,twos,pixtwos,aduspphot,photpts):
    pos=np.append(np.array([0]),np.cumsum(twos[:,4].astype(np.int32)))
    nph = np.digitize(twos[:,3],photpts)-1
    photonlist = np.zeros((np.sum(nph),3))
    ppos = np.append(np.array([0]),np.cumsum(nph))
    aduerr = 20.0
    chisqs = np.zeros(twos.shape[0])
#     print "\n number of droplets = ", twos.shape[0]
    #loop over all the droplets and find photon position
    drops=range(twos.shape[0])
    for count,drop in enumerate(drops):
        i = pixtwos[pos[drop]:pos[drop+1],0]
        j = pixtwos[pos[drop]:pos[drop+1],1]
        adus = pixtwos[pos[drop]:pos[drop+1],2].copy()
        int0 = twos[drop,3]
        #nophots=np.int32((int0+aduspphot-offset)/aduspphot)
        nophots = np.digitize(int0, photpts)-1		
        npix = len(i)
        #trivial case: single pixel
        if(len(i))==1:
            photonlist[ppos[drop]:ppos[drop+1],0] = twos[drop,0]
            photonlist[ppos[drop]:ppos[drop+1],1] = i
            photonlist[ppos[drop]:ppos[drop+1],2] = j
            continue
        # make a sub-image of only the droplet
        #(i,j) in full image, (ii,jj) in sub-image
        # could probaly use: 
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
        mi = np.min(i)
        mj = np.min(j)
        ii = 1+i-mi
        jj = 1+j-mj
        n,m = (np.max(i)-np.min(i)+1,np.max(j)-np.min(j)+1)
        img = np.zeros((n+2,m+2))
        img[ii,jj] = adus
        posi,posj = greedyguess(img,nophots,aduspphot)
        fposi = posi
        fposj = posj
        chisq = np.sum(((img-placephots(posi,posj,img.shape,aduspphot))/aduerr)**2)
        chisqs[drop] = chisq
        photonlist[ppos[drop]:ppos[drop+1],0] = twos[drop,0]
        photonlist[ppos[drop]:ppos[drop+1],1] = fposi-1+mi
        photonlist[ppos[drop]:ppos[drop+1],2] = fposj-1+mj
        # photonlist = np.append(onephots[:,[0,2,1]], photonlist, axis=0)
    return (photonlist)
