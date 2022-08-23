import numpy as np
#use that??
from smalldata_tools.ana_funcs.dropletCode.convert_img import convert_img
from smalldata_tools.ana_funcs.dropletCode.loopdrops import *
from smalldata_tools.ana_funcs.dropletCode.getProb import *
import scipy.ndimage.measurements as measurements
#use that??
import skimage.measure as measure
import scipy.ndimage.filters as filters
from scipy import sparse
import time
from smalldata_tools.DetObject import DetObjectFunc
from numba import jit
from numba.typed import List as NTL

class droplet2Photons(DetObjectFunc):
    """
    Uses loopdrops to find the photons in the droplets (don't forget to append the ones)
    
    Counts the number of photons at each (rounded) coordinate.
    """
    
    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        aduspphot: float
            Expected single photon energy
        offset: float (default = 0.5*aduspphot)
            Offset to convert number of photons into energy boundaries.
        photpts: float
            energy boundaries for given number of photons (defaults to offset + N*aduspphot)
        one_photon_limits: list
            Energy range for a single photon ((defaults to [offset, offset+N*aduspphot])
        mask: np.ndarray
            Pass a mask in here (array-form). If None: use mask stored in DetObject.
        one_photon_info: bool (default = False)
            Only returns info for one-photon droplets. DOES NOT WORK FOR NOW.
        cputime: bool (default = False)
            Store CPU tiime of different steps in hdf5
        """
        self._name =  kwargs.get('name','droplet2phot')
        super(droplet2Photons, self).__init__(**kwargs)
        
        self.aduspphot = kwargs.get('aduspphot', 0)
        self.offset = kwargs.get('offset', self.aduspphot*0.5)
        self.photpts = np.arange(1000000)*self.aduspphot-self.aduspphot+self.offset
        
        # self.Np = kwargs.get('Np', None)
        self.cputime = kwargs.get('cputime', False)
        self.one_photon_info = kwargs.get('one_photon_info', False)
        
        self._footprintnbr = np.array([[0,1,0],[1,0,1],[0,1,0]])
        self.one_photon_limits = kwargs.get('one_photon_limits', None)
        
        # not sure if this should not just be photpts[:2]?
        if self.one_photon_limits is None:
            self.one_photon_limits = [self.offset, self.offset+self.aduspphot]
        return
    
    
    @jit(nopython=True, cache=True)
    def piximg(self, i, j, adu, pad=True):
        """
        Make a sub-image of only the droplet.
        Takes coordinates & values, equivalent to sparse matrix w/ optionial edge pixels.
        Parameters
        ----------
        i,j: int
            Coordinates of the pixel in the droplet
        adu: float
            Per pixel ADU
        pad: bool (default = True)
            Pad droplet with 3 zeros on each side.
        """
        if pad:
            img = np.zeros((np.max(i)-np.min(i)+3, np.max(j)-np.min(j)+3))
            zip_obj_old = zip(i+1, j+1, adus)
        else:
            img = np.zeros((np.max(i)-np.min(i)+1, np.max(j)-np.min(j)+1))
            zip_obj_old = zip(i, j, adus)
        zip_obj = NTL()
        [zip_obj.append(x) for x in zip_obj_old]

        mi = np.min(i)
        mj = np.min(j)
        for ti, tj, tadu in zip_obj:
            img[ti-mi, tj-mj] = tadu
        return

            
#    @jit(nopython=True, cache=True)
    def onephoton_max(self, pixones, npix_drop, ppos):
        """
        pixones: adus/droplet
        npix_drop: number of pixels/droplet
        ppos: array of pixel positons of droplets
        returns arrays of the energy of the maximum pixel & highest two pixels for 1-photon droplets
        """
        maxPixAdu=[]
        twoPixAdu=[]
        twobnrPixAdu=[]
        for drop in range(len(npix_drop)):
            i = pixones[ppos[drop]:ppos[drop+1], 0]
            j = pixones[ppos[drop]:ppos[drop+1], 1]
            adus = pixonesadu[ppos[drop]:ppos[drop+1]]
            maxPixAdu.append(np.nanmax(adus))
            ##they need to be neighbors, this is not required here....
            #twoPixAdu.append((np.sort(adus)[-1*(min(2, len(adus))):]).sum())
            piximg1=sparse.coo_matrix( (adus, (i, j)) ).todense()
            #find maximum, find maximum of footprint_nbradd neighbor, maximum.
            nbrpix = filters.maximum_filter(piximg1, mode='constant', footprint=self._footprintnbr).flatten()[np.argmax(piximg)]
            twonbrPixAdu.append(maxPixAdu[-1]+nbrpix)
        return maxPixAdu, twonbrPixAdu


    def onephoton(self, image, imgDrop, detail=True):
        """
        image: image
        imgDrop: droplet image (result of label)
        detail: return energy of (two) highest pixels
        returns dict with keys:
            adu: array of energy per droplet
            pos: array of positions (3d: [tile or zeros], 1st&2nd coord of image)
            npix: array of number of pixels in droplet
            n1: number of 1-photon droplets
            [maxpixadu]: highest pixel ADU
            [twopixadu]: two highest pixel ADU
        """
        drop_ind = np.arange(1, np.nanmax(imgDrop)+1)
        adu_drop_all = measurements.sum(image, imgDrop, drop_ind)
        vThres = np.where(
            (adu_drop_all>=self.one_photon_limits[0])
            & (adu_drop_all<self.one_photon_limits[1])
        )[0]
        
        adu_drop = np.array(measurements.sum(image, imgDrop, drop_ind[vThres])) #only to check!
        pos_drop = np.array(measurements.center_of_mass(image, imgDrop, drop_ind[vThres])) 
        npix_drop = np.array(measurements.sum(image.astype(bool).astype(int),imgDrop, drop_ind[vThres])).astype(int)
        ones_pos = np.zeros((len(npix_drop),3))
        
        if len(image.shape) == 2:
            ones_pos[:,1] = pos_drop[:,0]
            ones_pos[:,2] = pos_drop[:,1]
        else:
            ones_pos[:,1] = pos_drop[:,1]
            ones_pos[:,2] = pos_drop[:,2]
            ones_pos[:,0] = pos_drop[:,0] #tile.

        ones_dict = {'adu': adu_drop}
        ones_dict['pos'] = ones_pos
        ones_dict['npix'] = npix_drop
        ones_dict['n1'] = len(adu_drop)

        if not detail: 
            return ones_dict
        
        imgDrop1 = imgDrop.copy()
        vThres_1_veto = np.where((adu_drop>=self.one_photon_limits[0])|(adu_drop<self.one_photon_limits[1]))[0]
        vetoed = np.in1d(imgDrop.ravel(), (vThres_1_veto+1)).reshape(imgDrop.shape)
        imgDrop1[vetoed] = 0
        drop_ind_thres1 = np.delete(drop_ind, vThres_1_veto)

        pp = np.where(imgDrop1>0)
        ss = np.argsort(imgDrop1[pp])
        npixtot = len(pp[0])#number of pixels in 1-photon droplets
        pixones = np.zeros((npixtot,len(image.shape)),dtype=np.int16)
        pixonesadu = np.zeros(npixtot)
        for i in range(len(image.shape)):
            pixones[:,i] = pp[i][ss]
        
        # TODO: error below, fix if one_photon_info is needed
        pixonesadu = image[pp[0][ss], pp[1][ss], pp[2][ss]] # dim dependant!!!
        # pixonesadu = image[pp[0][ss], pp[1][ss]] # dim dependant!!! # Seems to fix for 2d detectors...
        ppos1 = np.append(np.array([0]), np.cumsum(npix_drop_1))
        maxPixAdu, twoPixAdu = self.onephoton_max(pixones, npix_drop_1, ppos1)
        ones_dict['maxpixadu'] = maxPixAdu
        ones_dict['twopixadu'] = twoPixAdu
        return ones_dict


    def multi_photon(self, image, imgDrop, drop_ind):
        """
        image: image
        imgDrop: droplet image (result of label)
        drop_ind: list of droplet indices
        returns:
            twos: droplet coordinates (5d: [tile or zeros], 1st&2nd coord of image, energy, npixel)
            pixtwos: pixel coordinates
            pixtwoadus: adus for pixels in droplets
        """
        adu_drop = measurements.sum(image, imgDrop, drop_ind)
        
        #veto below threshold
        vThres = np.where(
            (adu_drop<self.one_photon_limits[1])
            | (adu_drop>self.photpts[-1])
        )[0]
        vetoed = np.in1d(imgDrop.ravel(), (vThres+1)).reshape(imgDrop.shape)
        imgDrop[vetoed] = 0
        drop_ind_thres = np.delete(drop_ind,vThres)

        npix_drop = (measurements.sum(image.astype(bool).astype(int),imgDrop, drop_ind_thres)).astype(int) #need 
        adu_drop = np.array(measurements.sum(image,imgDrop, drop_ind_thres))
        pos_drop = np.array(measurements.center_of_mass(image,imgDrop, drop_ind_thres))

        twos = np.zeros((len(npix_drop),5))
        twos[:,3] = adu_drop
        twos[:,4] = npix_drop
        if len(image.shape)==2:
            twos[:,1] = pos_drop[:,0]
            twos[:,2] = pos_drop[:,1]
        else:
            twos[:,1] = pos_drop[:,1]
            twos[:,2] = pos_drop[:,2]
            twos[:,0] = pos_drop[:,0] #tile.

        pp = np.where(imgDrop>0)
        ss = np.argsort(imgDrop[pp])
        npixtot = len(pp[0])

        pixtwos = np.zeros((npixtot,len(image.shape)),dtype=np.int16)
        pixtwosadu = np.zeros(npixtot)
        for i in range(len(image.shape)):
            pixtwos[:,i] = pp[i][ss]
        if len(image.shape) == 3:
            pixtwosadu = image[pp[0][ss],pp[1][ss],pp[2][ss]]
        else:
            pixtwosadu = image[pp[0][ss],pp[1][ss]]

        #this not being returned. I assume the idea is to return this instead of the multiple 
        #    arguments in a future iteration
        multpixdict = {'pixtwos':pixtwos}
        multpixdict['pixtwosadu'] = pixtwosadu

        return twos, pixtwos, pixtwosadu

    
    def process(self, data):
        #this will be a dictionary.
        if (not isinstance(data, dict)) or (data.get('_imgDrop',None) is None): 
            print('droplet2photons expects a dictionary with imgDrop and image keys!')
            return
        
        time0 = time.time()
        img = data['_image']
        imgDrop = data['_imgDrop']
        drop_ind = np.arange(1,np.nanmax(imgDrop)+1)

        if self.one_photon_info:
            ones_dict = self.onephoton(img, imgDrop)
        else:
            ones_dict =  self.onephoton(img, imgDrop, detail=False)

        #alternative and likely better version for twos:
        #confirms its the same
        w = np.where((data['data']>self.photpts[2])&(data['data']<=self.photpts[-1]))[0]
        n = len(w)
        ntwos = np.zeros((n,5))
        if 'tile' in data:
            ntwos[:,0] = data['tile'][w]
        
        ntwos[:,1] = data['row'][w]
        ntwos[:,2] = data['col'][w]
        ntwos[:,3] = data['data'][w]
        ntwos[:,4] = data['npix'][w]

        twos, multpixs, multpixadus = self.multi_photon(img, imgDrop, drop_ind)

        #photon_list is array of [tiles, x, y]
        #photon_list = loopdrops(twos,pixtwos,aduspphot,photpts)
        photonlist = loopdrops(twos, multpixs, multpixadus, self.aduspphot, self.photpts)

        ###
        # get the ones
        #
        # figure out the indexing if needed!
        # -- from cython versionindexes are inverted for ones because of c vs python indexing
        ###
        #ones_dict =  self.onephoton(img, imgDrop, detail=False)
        timed = time.time()
        #ones_dict[pos] = [ tile? row col]
        #this is a list of 3 values
        if len(ones_dict['pos'])>0:
            photonlist = np.append(ones_dict['pos'], photonlist, axis=0) 
 
        timep = time.time()        
        p = getProb_img(photonlist, data['_mask'], 12)
        # output dictionary
        ret_dict = {'prob': np.squeeze(p)}
        
        ###
        # do I really want to make this the image? That will not work for tiles. Better to keep
        #   this a dictionay.....
        ###
        #nx, ny = img.shape
        #phot_img, xedges, yedges = np.histogram2d(photonlist[:,1]+0.5, photonlist[:,2]+0.5, bins=[np.arange(nx+1),np.arange(ny+1)])
        #self.dat = np.ma.masked_array(phot_img, mask=(~(self.mask.astype(bool))).astype(np.uint8))
            
        photon_dict={'tile': photonlist[:,0]}
        photon_dict['col'] = photonlist[:,1]
        photon_dict['row'] = photonlist[:,2]
        photon_dict['data'] = np.ones(photonlist[:,0].shape[0])

        if self.cputime: ret_dict['cputime'] = np.array([time.time()-time0, timed-time0, timep-timed, time.time()-timep]) 
        self.dat = photon_dict                                             

        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]

        return ret_dict
