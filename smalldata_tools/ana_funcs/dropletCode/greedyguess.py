import numpy as np
#Greedy assignment of photons
#from fitdrop import placephots
from numba import jit

@jit(nopython=True, cache=True)
def greedyguess(img,nophots,aduspphot):
    pxs=np.zeros(nophots)
    pys=pxs.copy()
    
    fimg=img.copy()/aduspphot
    pimg=fimg//1
    timg=fimg-pimg

    nophots_assigned=0
    #pimg->sparse array, fill pxs&pys, figure out what to do with photons > 1.
    for ii in range(pimg.shape[0]):
        for jj in range(pimg.shape[1]):
            for photon in range(pimg[ii,jj]):
                pxs[nophots_assigned] = ii
                pys[nophots_assigned] = jj
                nophots_assigned=nophots_assigned+1

    for n in range(nophots_assigned,nophots):
        #find max intensity for next photon guess
        maximg=-1
        for ii in range(timg.shape[0]):
            for jj in range(timg.shape[1]):
                if timg[ii,jj]>maximg:
                    i=ii
                    j=jj
                    maximg = timg[ii,jj]
        #i,j=np.unravel_index(np.argmax(timg),timg.shape)
        #put a photon here if it fits
        t = 1.0-timg[i,j]  #fraction of photon not in i,j
        timg[i,j]=0.0 #use up central pixel

        #find the max neighboring pixel
        if timg[i-1,j]>max(timg[i+1,j],timg[i,j-1],timg[i,j+1]):
            pxs[n]=i-t
            pys[n]=j
            timg[i-1,j]=np.maximum(0.0,timg[i-1,j]-t)
        elif timg[i+1,j]>max(timg[i,j-1],timg[i,j+1]):
            pxs[n]=i+t
            pys[n]=j
            timg[i+1,j]=np.maximum(0.0,timg[i+1,j]-t)
        elif timg[i,j-1]>timg[i,j+1]:
            pxs[n]=i
            pys[n]=j-t
            timg[i,j-1]=np.maximum(0.0,timg[i,j-1]-t)
        else:
            pxs[n]=i
            pys[n]=j+t
            timg[i,j+1]=np.maximum(0.0,timg[i,j+1]-t)

    return(pxs,pys)
