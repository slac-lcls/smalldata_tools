import numpy as np 
import time as time
from droplet import dropletfind,dropletanal
from frontiers import *
def convert(a,goodimgs,thres,photpts,avgflag=True,mask=None,prflag=None):
    #loop
    # setup arrays
    onephots=np.zeros((10000000,5))
    nophots=0
    if avgflag:
        avgimg=np.zeros_like(a[0,:,:])
        histimg=avgimg.copy()
    idhold=np.zeros(1+np.size(a[0,:,:]))
    twos=np.zeros((6000000,6))
    notwos=0
    pixtwos=np.zeros((30000000,3),dtype=np.int16)
    no2pix=0
    t0=time.time()
    lg = 0
    bins = np.linspace(0,2000,1001)
    hh = np.zeros((1000))
    for fr in range(np.size(goodimgs)):
    #for fr in range(10):
        if(goodimgs[fr]):
            img=a[fr,:,:]
            if mask is not None: img *= mask
            npeaks,dimg=dropletfind((img>thres).astype(int))
            dimghold=dimg.copy()
            if(npeaks != 0):
                    npix,xcen,ycen,xsig,ysig,xysig,adus,idlist=dropletanal(
                        img.astype(int),dimg,npeaks)
                    h,b = np.histogram(adus,bins = bins)
                    hh += h
                    #check idlist is sorted
                    if np.any(idlist[:-1]>=idlist[1:]): print("not sorted ",fr)
                    
                    #handle 1 photon droplets
                    w=np.where((adus>photpts[1])&(adus<=photpts[2]))[0]
                    n=len(w)
                    onephots[nophots:nophots+n,0]=fr
                    onephots[nophots:nophots+n,1]= xcen[w]
                    onephots[nophots:nophots+n,2]= ycen[w]
                    onephots[nophots:nophots+n,3]= adus[w]
                    onephots[nophots:nophots+n,4]= npix[w]
                    nophots += n
                    if avgflag:
                        idhold[idlist[w]]=idlist[w]
                        pp=np.where(idhold[dimghold])
                        avgimg[pp[0],pp[1]] += img[pp[0],pp[1]]
                        histimg[pp[0],pp[1]] += 1
                        idhold[idlist[w]]=0

                    # more than single photon droplets
                    #find ids of >1 photon drops
                    lg +=np.sum(adus>photpts[-1])
                    w=np.where((adus>photpts[2])&(adus<=photpts[-1]))[0]
                    idhold[idlist[w]]=idlist[w]
                    pp=np.where(idhold[dimghold])
                    if avgflag:
                        avgimg[pp[0],pp[1]] += img[pp[0],pp[1]]
                        histimg[pp[0],pp[1]] += 1
                    idhold[idlist[w]]=0

                    n=len(w)
                    twos[notwos:notwos+n,0]=fr
                    twos[notwos:notwos+n,1]=xcen[w]
                    twos[notwos:notwos+n,2]=ycen[w]
                    twos[notwos:notwos+n,3]=adus[w]
                    twos[notwos:notwos+n,4]=npix[w]
                    

                    n1=len(pp[0])
                    ss=np.argsort(dimghold[pp[0],pp[1]]) #sort by droplet id
                    pixtwos[no2pix:no2pix+n1,0]=pp[0][ss]
                    pixtwos[no2pix:no2pix+n1,1]=pp[1][ss]
                    pixtwos[no2pix:no2pix+n1,2]=img[pp[0][ss],pp[1][ss]]
                    no2pix += n1
#                     fadus=frontiers(pp[0][ss],pp[1][ss],npix[w],img)
#                     twos[notwos:notwos+n,5]=fadus
                    notwos += n


        if prflag is not None: 
            if(0==fr%500): 
                print('frame: ',fr,time.time()-t0)
    onephots=onephots[0:nophots,:]
    twos=twos[0:notwos,:]
    pixtwos=pixtwos[0:no2pix,:]
    if prflag is not None: print 'number of large droplet', lg
    if avgflag:
        return onephots,twos,pixtwos,avgimg,histimg, hh, b
    else:
        return onephots,twos,pixtwos
