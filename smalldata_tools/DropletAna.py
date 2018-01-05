import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse
from utilities import hist2d

class droplets(object):
    def __init__(self,h5file,detName='epix',dropName='droplet'):
        self._detName=detName
        self._dropName=dropName
        if self._dropName[-1]!='_':
            self._dropName+='_'
        self._h5 = h5file
        self._h5dir = self._h5.get_node('/'+self._detName)

    def getKeys(self):
        keyList=[]
        for data in self.__dict__.keys():
            if not data[0]=='_':
                keyList.append(data)
        return keyList

    def fillDropArrays(self, only_XYADU=False, Filter=None):
        for node in self._h5dir._f_list_nodes():
            if not (node.name.find(self._dropName)>=0):
                continue
            #likely a different dropName set of variables. Need own fillDropArray
            if node.name.replace(self._dropName,'').find('_')>=0:
                continue
            h5Name = self._dropName+node.name.replace(self._dropName,'')
            keyName=node.name.replace(self._dropName,'')
            if not only_XYADU:
                print 'fill drop ',h5Name
                if Filter is not None and Filter.shape[0] == self._h5.get_node('/'+self._detName, h5Name).shape[0]:
                    self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()[Filter]
                else:
                    self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()
            elif (keyName=='X' or keyName=='Y' or keyName=='adu' or keyName=='npix'):
                print 'fill drop ',h5Name
                if Filter is not None and Filter.shape[0] == self._h5.get_node('/'+self._detName, h5Name).shape[0]:
                    self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()[Filter]
                else:
                    self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()

    def flattenDropArray(self, filterArray=None):
        if filterArray is None:
            if 'adau' in self.getKeys():
                filterArray = (self.__dict__['adu']>0)
            else:
                print 'did not find adu, will not flatten'
                return

        for key in self.getKeys():
            if filterArray.shape == self.__dict__[key].shape:
                self.__dict__[key] = self.__dict__[key][filterArray].flatten()

    def getDropPixels(self, ievt, debug=False):
        print 'will need to be reimplemented'
        return
        if not 'Pix' in self.__dict__.keys():
            nDroplets = self._h5[self._detName][self._dropName+'nDroplets'][ievt]
            Pix = self._h5[self._detName][self._dropName+'Pix'][ievt]
            sizeX = self._h5[self._detName][self._dropName+'bbox_x1'][ievt] - self._h5[self._detName][self._dropName+'bbox_x0'][ievt]
            sizeY = self._h5[self._detName][self._dropName+'bbox_y1'][ievt] - self._h5[self._detName][self._dropName+'bbox_y0'][ievt]
            adu = self._h5[self._detName][self._dropName+'adu'][ievt]
            dX =  self._h5[self._detName][self._dropName+'X'][ievt]
            dY =  self._h5[self._detName][self._dropName+'Y'][ievt]
        else:
            nDroplets = self.nDroplets[ievt]
            Pix = self.Pix[ievt]
            sizeX = self.bbox_x1[ievt] - self.bbox_x0[ievt]
            sizeY = self.bbox_y1[ievt] - self.bbox_y0[ievt]
            adu = self.adu[ievt]
            dX =  self.X[ievt]
            dY =  self.Y[ievt]
        sizes = (sizeX*sizeY)
        imgs=[]
        idxPix=0
        for iDrop in range(0,nDroplets):
            img= np.array(Pix[idxPix:(idxPix+sizes[iDrop])]).reshape(sizeX[iDrop],sizeY[iDrop])
            if debug:
                print 'adu ',adu[iDrop],img.sum()
            imgs.append(img)
            idxPix+=sizes[iDrop]
        ret_dict = {'images' : imgs}
        ret_dict['adu']=adu[:len(imgs)]
        ret_dict['X']=dX[:len(imgs)]
        ret_dict['Y']=dY[:len(imgs)]
        return ret_dict

    def getDropPixelsRoi(self, ievt, mask, debug=False):
        dropInfo = self.getDropPixelsRoi(ievt, debug=debug)
        imgsROI = []
        for img,x,y in zip(dropInfo['images'],dropInfo['X'],dropInfo['Y']):
            if mask(int(x), int(y))==1:
                imgsROI.append(img)
        return imgsROI

    def plotSpectrum(self, plotLog=True, aduRange=[]):
        if len(aduRange)==0:
            aduRange=[0,700]
        elif len(aduRange)==1:
            aduRange=[0,aduRange[0]]
        elif len(aduRange)==2:
            aduRange=[aduRange[0], aduRange[1]]
            
        hst = np.histogram(self.__dict__['adu'], np.arange(aduRange[0],aduRange[1]))
        if plotLog:
            plt.semilogy(hst[1][1:],hst[0],'o')
        else:
            plt.plot(hst[1][1:],hst[0],'o')
            
    def plotAduX(self, ADUrange=[120,180],maxLim=99.5):
        if len(self.__dict__['X'].shape)>1:
            self.flattenDropArray()
        plt.figure()
        hist2d(self.__dict__['X'][self.__dict__['Y']<1000],self.__dict__['adu'][self.__dict__['Y']<1000], numBins=[702,180], histLims=[0,702,ADUrange[0], ADUrange[1]],limits=[1,maxLim])

    def plotAduY(self, ADUrange=[120,180],maxLim=99.5):    
        if len(self.__dict__['Y'].shape)>1:
            self.flattenDropArray()
        plt.figure()
        plt.subplot(211)
        hist2d(self.__dict__['Y'][self.__dict__['X']<351],self.__dict__['adu'][self.__dict__['X']<351], numBins=[766,180], histLims=[0,766,ADUrange[0], ADUrange[1]],limits=[1,maxLim])
        plt.subplot(212)
        hist2d(self.__dict__['Y'][self.__dict__['X']>353],self.__dict__['adu'][self.__dict__['X']>353], numBins=[766,180], histLims=[0,766,ADUrange[0], ADUrange[1]],limits=[1,maxLim])

    def plotXY(self, ADUrange=[120,180], npix=0):            
        allX = self.__dict__['X'][self.__dict__['adu']>ADUrange[0]]
        allY = self.__dict__['Y'][self.__dict__['adu']>ADUrange[0]]
        alladu = self.__dict__['adu'][self.__dict__['adu']>ADUrange[0]]
        if npix!=0:
            allNpix=self.__dict__['npix'][self.__dict__['adu']>ADUrange[0]]
            
        if ADUrange[1]>ADUrange[0]:
            allX = allX[alladu<ADUrange[1]]
            allY = allY[alladu<ADUrange[1]]
            if npix!=0:
                allNpix = allNpix[alladu<ADUrange[1]]
            alladu = alladu[alladu<ADUrange[1]]
        if npix!=0:
            if npix>0:
                allX = allX[allNpix==npix] 
                allY = allY[allNpix==npix] 
                alladu = alladu[allNpix==npix] 
            else:
                allX = allX[allNpix>=abs(npix)] 
                allY = allY[allNpix>=abs(npix)] 
                alladu = alladu[allNpix>=abs(npix)] 

        plt.figure()
        ndrop_int = max(1,490000./allX.shape[0])
        hist2d(allX,allY, numBins=[int(702/ndrop_int),int(766/ndrop_int)])

    def aduSlices(self,axis='y', ADUrange=[0,-1], npix=0):
        allX = self.__dict__['X'][self.__dict__['adu']>ADUrange[0]]
        allY = self.__dict__['Y'][self.__dict__['adu']>ADUrange[0]]
        alladu = self.__dict__['adu'][self.__dict__['adu']>ADUrange[0]]
        if npix!=0:
            allNpix=self.__dict__['npix'][self.__dict__['adu']>ADUrange[0]]
        if ADUrange[1]>ADUrange[0]:
            allX = allX[alladu<ADUrange[1]]
            allY = allY[alladu<ADUrange[1]]
            if npix!=0:
                allNpix = allNpix[alladu<ADUrange[1]]
            alladu = alladu[alladu<ADUrange[1]]
        if npix!=0:
            if npix>0:
                allX = allX[allNpix==npix] 
                allY = allY[allNpix==npix] 
                alladu = alladu[allNpix==npix] 
            else:
                allX = allX[allNpix>=abs(npix)] 
                allY = allY[allNpix>=abs(npix)] 
                alladu = alladu[allNpix>=abs(npix)] 
                
        if axis=='y':
            nSlice=16
            allY+=(allX>351).astype(int)*768
            sliceSize=768*2./nSlice
            binVar=allY
        elif axis=='x':
            nSlice=14
            sliceSize=704./nSlice
            binVar=allX
            
        aduS=[]
        for i in range(0,nSlice):
            aduS.append([])
        for adu,bv in zip(alladu,binVar):
            aduS[int(bv/sliceSize)].append(adu)
        return aduS
 
