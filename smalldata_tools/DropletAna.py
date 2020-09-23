import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse
from utilities import hist2d
import holoviews as hv

class droplets(object):
    def __init__(self,h5file,detName='epix',dropName='droplet', plotWith='matplotlib'):
        self._detName=detName
        self._dropName=dropName
        if self._dropName[-1]!='_':
            self._dropName+='_'
        self._h5 = h5file
        self._h5dir = self._h5.get_node('/'+self._detName)
        self._plotWith=plotWith

    def getKeys(self):
        keyList=[]
        for data in self.__dict__.keys():
            if not data[0]=='_':
                keyList.append(data)
        return keyList

    def fillDropArrays(self, only_XYADU=False, Filter=None):
        #get shape of detector, right now EPIX is hardcoded...
        try:
            dX = self._h5.get_node('/UserDataCfg/'+self._detName, 'iX').read().max() - self._h5.get_node('/UserDataCfg/'+self._detName, 'iX').read().min()
            dY = self._h5.get_node('/UserDataCfg/'+self._detName, 'iY').read().max() - self._h5.get_node('/UserDataCfg/'+self._detName, 'iY').read().min()
        except:
            try:
                dX = self._h5.get_node('/UserDataCfg/'+self._detName, 'ix').read().max() - self._h5.get_node('/UserDataCfg/'+self._detName, 'ix').read().min()
                dY = self._h5.get_node('/UserDataCfg/'+self._detName, 'iy').read().max() - self._h5.get_node('/UserDataCfg/'+self._detName, 'iy').read().min()
            except:
                print('could not find x/y arrays, quit')
                return
        self.__dict__['detSize'] = [dX, dY] #708, 772 or epix
        for node in self._h5dir._f_list_nodes():
            if not (node.name.find(self._dropName)>=0):
                continue
            ##likely a different dropName set of variables. Need own fillDropArray
            #if (node.name.replace(self._dropName,'').find('_')>=0):
            #    continue
            h5Name = self._dropName+node.name.replace(self._dropName,'')
            keyName=node.name.replace(self._dropName,'')
            if not only_XYADU:
                print('fill drop all %s'%h5Name)
                if Filter is not None and Filter.shape[0] == self._h5.get_node('/'+self._detName, h5Name).shape[0]:
                    self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()[Filter]
                else:
                    self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()
            elif (keyName=='X' or keyName=='Y' or keyName=='adu' or keyName=='npix'):
                print('fill drop (only x,y,adu) ',h5Name)
                #if Filter is not None:
                #    print 'DEBUG Shapes', h5Name, Filter.shape[0], self._h5.get_node('/'+self._detName, h5Name).shape[0]
                if Filter is not None and Filter.shape[0] == self._h5.get_node('/'+self._detName, h5Name).shape[0]:
                    print('DEBUG: ',h5Name, self._h5.get_node('/'+self._detName, h5Name).read().shape, Filter.shape)
                    self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()[Filter]
                else:
                    self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()

    def flattenDropArray(self, filterArray=None):
        self.dataname=None
        if 'adu' in self.getKeys():
            self.dataName='adu'
            self.xName='X'
            self.yName='Y'
        elif 'data' in self.getKeys():
            self.dataName='data'
            self.xName='row'
            self.yName='col'
        elif 'sparse_data' in self.getKeys():
            self.dataName='sparse_data'
            self.xName='sparse_row'
            self.yName='sparse_col'
        else:
            print('did not find adu, will not flatten')
            return

        if filterArray is None:
            print('filterarray....', filterArray)
            filterArray = (self.__dict__[self.dataName]>0)
            print('filterarray---....', filterArray)
        else:
            if len(filterArray.shape)==1 and filterArray.shape[0] == self.__dict__[self.dataName].shape[0]:
                filterArray = np.array([filterArray ,]*self.__dict__[self.dataName].shape[1]).transpose() & (self.__dict__[self.dataName]>0)

        for key in self.getKeys():
            if key == 'detSize':
                continue
            if isinstance(self.__dict__[key], np.ndarray) and (np.array(filterArray).shape == self.__dict__[key].shape):
                self.__dict__[key] = self.__dict__[key][filterArray].flatten()

    def getDropPixels(self, ievt, debug=False):
        print('will need to be reimplemented')
        return
        if not 'Pix' in self.__dict__.keys():
            nDroplets = self._h5[self._detName][self._dropName+'nDroplets'][ievt]
            Pix = self._h5[self._detName][self._dropName+'Pix'][ievt]
            sizeX = self._h5[self._detName][self._dropName+'bbox_x1'][ievt] - self._h5[self._detName][self._dropName+'bbox_x0'][ievt]
            sizeY = self._h5[self._detName][self._dropName+'bbox_y1'][ievt] - self._h5[self._detName][self._dropName+'bbox_y0'][ievt]
            adu = self._h5[self._detName][self._dropName+self.dataName][ievt]
            dX =  self._h5[self._detName][self._dropName+self.xName][ievt]
            dY =  self._h5[self._detName][self._dropName+self.yName][ievt]
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
                print('adu ',adu[iDrop],img.sum())
            imgs.append(img)
            idxPix+=sizes[iDrop]
        ret_dict = {'images' : imgs}
        ret_dict['adu']=adu[:len(imgs)]
        ret_dict[self.xName]=dX[:len(imgs)]
        ret_dict[self.yName]=dY[:len(imgs)]
        return ret_dict

    def getDropPixelsRoi(self, ievt, mask, debug=False):
        dropInfo = self.getDropPixelsRoi(ievt, debug=debug)
        imgsROI = []
        for img,x,y in zip(dropInfo['images'],dropInfo[self.xName],dropInfo[self.yName]):
            if mask(int(x), int(y))==1:
                imgsROI.append(img)
        return imgsROI

    def plotSpectrum(self, plotLog=True, aduRange=[], ROI=None):
        if len(aduRange)==0:
            aduRange=[0.1,700]
        elif len(aduRange)==1:
            aduRange=[0.1,aduRange[0]]
        elif len(aduRange)==2:
            aduRange=[aduRange[0], aduRange[1]]
            
        if ROI is None:
            dadu = self.__dict__[self.dataName]
        else:
            dadu = self.__dict__[self.dataName]
            if self.dataName == 'adu':
                dx = self.__dict__['x'].flatten()
                dy = self.__dict__['y'].flatten()
            elif self.dataName == 'adu':
                dx = self.__dict__['row'].flatten()
                dy = self.__dict__['col'].flatten()
            else:
                return
            inROIx = np.logical_and(dx > np.array(self.ROI).flatten()[0],dx < np.array(self.ROI).flatten()[1]) 
            inROIy = np.logical_and(dy > np.array(self.ROI).flatten()[2],dy < np.array(self.ROI).flatten()[3]) 
            inROI = np.logical_and(inROIx, inROIy)
            dadu = dadu.flatten()[inROI]

        hst = np.histogram(dadu, np.arange(aduRange[0],aduRange[1]))
        if self._plotWith=='matplotlib':
            if plotLog:
                plt.semilogy(hst[1][1:],hst[0],'o')
            else:
                plt.plot(hst[1][1:],hst[0],'o')
        elif self._plotWith=='holoviews':
            if plotLog:
                return hv.Scatter((hst[1][1:],np.log(hst[0])))#,'ADU','log_nDrops')
            else:
                return hv.Scatter((hst[1][1:],hst[0]))#,'ADU','nDrops')
        else:
            return hst
            
    def plotAduX(self, aduRange=[120,180],maxLim=99.5):
        #self.__dict__['detSize'] = [dX, dY] #708, 772 or epix
        dX = self.__dict__['detSize'][0]
        dY = self.__dict__['detSize'][1]
        if len(self.__dict__[self.xName].shape)>1:
            self.flattenDropArray()
        plt.figure()
        hist2d(self.__dict__[self.xName][self.__dict__[self.yName]<=dY],self.__dict__[self.dataName][self.__dict__[self.yName]<=dY], numBins=[dX,180], histLims=[0,dX,aduRange[0], aduRange[1]],limits=[1,maxLim])

    def plotAduY(self, aduRange=[120,180],maxLim=99.5):    
        dX = self.__dict__['detSize'][0]
        dY = self.__dict__['detSize'][1]
        if len(self.__dict__[self.yName].shape)>1:
            self.flattenDropArray()
        plt.figure()
        if self._detName.find('epix')>=0 or self._detName.find('Epix')>=0:
            plt.subplot(211)
            hist2d(self.__dict__[self.yName][self.__dict__[self.xName]<351],self.__dict__[self.dataName][self.__dict__[self.xName]<351], numBins=[dY,180], histLims=[0,dY,aduRange[0], aduRange[1]],limits=[1,maxLim])
            plt.subplot(212)
            hist2d(self.__dict__[self.yName][self.__dict__[self.xName]>353],self.__dict__[self.dataName][self.__dict__[self.xName]>353], numBins=[dY,180], histLims=[0,dY,aduRange[0], aduRange[1]],limits=[1,maxLim])
        else:
            hist2d(self.__dict__[self.yName][self.__dict__[self.xName]<=dX],self.__dict__[self.dataName][self.__dict__[self.xName]<=dX], numBins=[dY,180], histLims=[0,dY,aduRange[0], aduRange[1]],limits=[1,maxLim])

    def plotXY(self, aduRange=[120,180], npix=0):            
        dX = self.__dict__['detSize'][0]
        dY = self.__dict__['detSize'][1]
        allX = self.__dict__[self.xName][self.__dict__[self.dataName]>aduRange[0]]
        allY = self.__dict__[self.yName][self.__dict__[self.dataName]>aduRange[0]]
        alladu = self.__dict__[self.dataName][self.__dict__[self.dataName]>aduRange[0]]
        if npix!=0:
            allNpix=self.__dict__['npix'][self.__dict__[self.dataName]>aduRange[0]]
            
        if aduRange[1]>aduRange[0]:
            allX = allX[alladu<aduRange[1]]
            allY = allY[alladu<aduRange[1]]
            if npix!=0:
                allNpix = allNpix[alladu<aduRange[1]]
            alladu = alladu[alladu<aduRange[1]]
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
        hist2d(allX,allY, numBins=[int(dX/ndrop_int),int(dY/ndrop_int)])

    def aduSlices(self,axis='y', aduRange=[0,-1], npix=0):
        allX = self.__dict__[self.xName][self.__dict__[self.dataName]>aduRange[0]]
        allY = self.__dict__[self.yName][self.__dict__[self.dataName]>aduRange[0]]
        alladu = self.__dict__[self.dataName][self.__dict__[self.dataName]>aduRange[0]]
        if npix!=0:
            allNpix=self.__dict__['npix'][self.__dict__[self.dataName]>aduRange[0]]
        if aduRange[1]>aduRange[0]:
            allX = allX[alladu<aduRange[1]]
            allY = allY[alladu<aduRange[1]]
            if npix!=0:
                allNpix = allNpix[alladu<aduRange[1]]
            alladu = alladu[alladu<aduRange[1]]
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
 
