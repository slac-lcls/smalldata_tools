#finding photon positions without/with fit, code from Mark Sutton
#fitting is generally fast with intensity level 1e-3
from smalldata_tools.ana_funcs.dropletCode.greedyguess import greedyguess
#from smalldata_tools.ana_funcs.dropletCode.fitdrop import *
import numpy as np
from numba import jit
from numba.typed import List as NTL

@jit(nopython=True, cache=True)
def loopdrops(twos,pixtwos,pixtwoadus,aduspphot,photpts):
    """
    place single photons into a multiphoton droplets
    arguments: arrays of multiphoton droplets, array of pixels assigned to photons,
        array of pixel ADU assigned to photons, expected energy/photon, photon energy boundaries
    
    creates the images of multiphoton droplets & uses greedyguess to place photons
    """
    pos=np.append(np.array([0]),np.cumsum(twos[:,4].astype(np.int32)))
    nph = np.digitize(twos[:,3],photpts)-1
    photonlist = np.zeros((np.sum(nph),3))
    ppos = np.append(np.array([0]),np.cumsum(nph))

    tiled = pixtwos.shape[1]>2
    #loop over all the droplets and find photon position
    for drop in range(twos.shape[0]):
        i = pixtwos[pos[drop]:pos[drop+1],pixtwos.shape[1]-2]
        j = pixtwos[pos[drop]:pos[drop+1],pixtwos.shape[1]-1]
        adus = pixtwoadus[pos[drop]:pos[drop+1]].copy()
        int0 = twos[drop,3]
        nophots = nph[drop]
        npix = len(i)
        #trivial case: single pixel
        if(npix)==1:
            photonlist[ppos[drop]:ppos[drop+1],0] = twos[drop,0]
            photonlist[ppos[drop]:ppos[drop+1],1] = i
            photonlist[ppos[drop]:ppos[drop+1],2] = j
            continue
        # make a sub-image of only the droplet
        img = np.zeros((np.max(i)-np.min(i)+3,np.max(j)-np.min(j)+3))
        zip_obj_old = zip(i+1,j+1,adus)
        zip_obj = NTL()
        [zip_obj.append(x) for x in zip_obj_old]

        mi = np.min(i)
        mj = np.min(j)
        for ti,tj,tadu in zip_obj:
            img[ti-mi,tj-mj]=tadu

        posi,posj = greedyguess(img,nophots,aduspphot)
        fposi = posi
        fposj = posj

        photonlist[ppos[drop]:ppos[drop+1],0] = twos[drop,0]
        photonlist[ppos[drop]:ppos[drop+1],1] = fposi-1+mi
        photonlist[ppos[drop]:ppos[drop+1],2] = fposj-1+mj

    return (photonlist)
