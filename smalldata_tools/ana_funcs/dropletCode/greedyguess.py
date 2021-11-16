import numpy as np
#Greedy assignment of photons
#from fitdrop import placephots
#from numba import jit

#@jit
def greedyguess(img,nophots,aduspphot):
    #relative indices for the 4 tiles, order is carefully chosen
    xs=np.array([[0,0],[0,0],[-1,0],[1,0]])
    ys=np.array([[-1,0],[1,0],[0,0],[0,0]])

    timg=img.copy()
    pxs=np.zeros(nophots)
    pys=pxs.copy()
    for n in range(nophots):
        #find max intensity for next photon guess
        i,j=np.unravel_index(np.argmax(timg),timg.shape)
        #put a photon here if it fits
        if(timg[i,j]>aduspphot):
            timg[i,j] -= aduspphot
            pxs[n] = i
            pys[n] = j
            continue #go to next photon
        Is=np.sum(timg[i+xs,j+ys],axis=1)#sum the four pairs at (i,j)
        #find max pair
        mi=np.argmax(Is)
        #put all of pixel i,j into photon (not C.O.M.)
        t = 1.0-timg[i,j]/aduspphot  #fraction of photon not in i,j
        timg[i,j]=0.0 #use up central pixel
        if(mi<2): #vertical pair
            p1i=0.0
            p1j=t*ys[mi,0] #vertical shift
            timg[i,j+ys[mi,0]]=np.maximum(0.0,timg[i,j+ys[mi,0]]-t*aduspphot)
        else: #horizontal pair
            p1i=t*xs[mi,0] #horizontal shift
            p1j=0.0
            timg[i+xs[mi,0],j]=np.maximum(0.0,timg[i+xs[mi,0],j]-t*aduspphot)
        pxs[n]=i+p1i
        pys[n]=j+p1j
    return(pxs,pys)
