import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse
from utilities import hist2d

class photons(object):
    def __init__(self,h5file,detName='epix',photName='photons'):
        self._detName=detName
        self._photName=photName
        if self._photName[-1]!='_':
            self._photName+='_'
        self._h5 = h5file
        self._h5dir = self._h5.get_node('/'+self._detName)
        self.shape = h5file.get_node('/UserDataCfg/'+detName+'/cmask').shape

    def getKeys(self):
        keyList=[]
        for data in self.__dict__.keys():
            if not data[0]=='_':
                keyList.append(data)
        return keyList

    def fillPhotArrays(self):
        for node in self._h5dir._f_list_nodes():
            if not (node.name.find(self._photName)>=0):
                continue

            if node.name.replace(self._photName,'').find('_')>=0:
                continue

            h5Name = self._photName+node.name.replace(self._photName,'')
            #print('h5name ',h5Name)
            keyName=node.name.replace(self._photName,'')
            self.__dict__[keyName] = self._h5.get_node('/'+self._detName, h5Name).read()
        
    def flattenPhotArray(self, filterArray=None):
        for key in self.getKeys():
            if filterArray is None:
                filterArray=np.ones(self.__dict__[key].shape[0])
            if filterArray.shape == self.__dict__[key].shape[0]:
                self[key] = self.__dict__[key][filterArray].flatten()

    def photImgEvt(self, iEvt):
        data = self._h5.get_node('/'+self._detName+'/'+self._photName+'data')[iEvt]
        row = self._h5.get_node('/'+self._detName+'/'+self._photName+'row')[iEvt]
        col = self._h5.get_node('/'+self._detName+'/'+self._photName+'col')[iEvt]
        if max(row)>=self.shape[0] or max(col)>=self.shape[1]:
            print('inconsistent shape ',self.shape, max(row), max(col))
        return sparse.coo_matrix( (data, (row, col)),shape=self.shape ).todense()

    def photImg(self, filterArray=None):
        if 'pHist' not in self.__dict__.keys():
            self.fillPhotArrays()
        if filterArray is None:
            filterArray=np.ones(self.pHist.shape[0]).astype(bool)
        data = self.data[filterArray].flatten()
        data = data[data>0]
        row = self.row[filterArray].flatten()[data>0]
        col = self.col[filterArray].flatten()[data>0]
        img = sparse.coo_matrix( (data, (row, col)) ).todense()
        return img

    def photonHist(self, filterArray=None):
        if 'pHist' not in self.__dict__.keys():
            self.fillPhotArrays()
        if filterArray is None:
            filterArray=np.ones(self.pHist.shape[0]).astype(bool)
        if filterArray.shape[0] == self.pHist.shape[0]:
            pHist = self.pHist[filterArray].sum(axis=0)
        else:
            return
        return pHist

