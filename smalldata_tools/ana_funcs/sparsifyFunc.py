import os
import numpy as np
from scipy import sparse
import numpy.ma as ma
import itertools

import time
from smalldata_tools.utilities import unsparsify
from smalldata_tools.DetObjectFunc import DetObjectFunc

#
# for now, this works in "raw" coordinates.
# enable ROI in other coordinates: e.g. x, y: save ROI in x/y save across tiles.
# add asImg. Is this an option or a separate function?
#
#effectitely a projection onto a non-spatial coordinate.
#TEST ME WITH MASKED ARRAY  - WAS ALSO WRITTEN FOR MASKED ARRAYS
class sparsifyFunc(DetObjectFunc):
    """
    Function to sparsify a passed array (2 or 3-d input)
    nData: if passed, make output array rectangular (for storing in event based smlData)
    if a dictionary w/ data, row, col is passed, only make rectangular
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','sparse')
        super(sparsifyFunc, self).__init__(**kwargs)
        self.nData =  kwargs.get('nData',None)
        self._flagMasked = kwargs.get('flagMasked',False)
        self._needProps = kwargs.get('needProps',False)
        self._saveint = kwargs.get('saveInt',True)
        self._saveintadu = kwargs.get('saveIntADU',False)
        
        if type(self.nData) is float:
            self.nData = int(self.nData)

    def process(self, data):
        #apply mask - set to zero, so pixels will fall out in sparify step.
        if isinstance(data, np.ma.masked_array):            
            data = data.filled(fill_value=0)

        #already have dict w/ data
        if  isinstance(data, dict):
            if 'data' not in data or 'row' not in data or 'col' not in data:
                print('cannot make a make a sparse, rectangular array of', data)
                return
            ret_dict = data
        
        #sparsify image
        if  isinstance(data, np.ndarray):
            photonsImg = data.copy()
            if len(photonsImg.shape)>2: #tiled detector!
                data=[]
                row=[]
                col=[]
                tile=[]
                for iTile,photonTile in enumerate(photonsImg):
                    sImage = sparse.coo_matrix(photonTile)
                    data = list(itertools.chain.from_iterable([data, sImage.data.tolist()]))
                    row = list(itertools.chain.from_iterable([row, sImage.row.tolist()]))
                    col = list(itertools.chain.from_iterable([col, sImage.col.tolist()]))
                    tile = list(itertools.chain.from_iterable([tile, (np.ones_like(sImage.data)*iTile).tolist()]))
                data = np.array(data)
                row = np.array(row)
                col = np.array(col)
                tile= np.array(tile)
            else:
                sImage = sparse.coo_matrix(photonsImg)
                data = sImage.data
                row = sImage.row
                col = sImage.col
                tile = np.zeros_like(data)
            ret_dict = {'data':data}
            ret_dict['row'] = row
            ret_dict['col'] = col
            ret_dict['tile'] = tile
            
            
        #now fix shape of data in dict.
        data_dict = ret_dict
        ret_dict = {}
        if self.nData == 0:
            for key, d in data_dict.items():
                if key[0]=='_': continue
                ret_dict[f'ragged_{key}'] = d
        elif self.nData is not None:                
            for key in data_dict.keys():
                if key[0]=='_': continue
                if data_dict[key].shape[0] >= self.nData:
                    ret_dict[key]=data_dict[key][:self.nData]
                else:
                    ret_dict[key]=(np.append(data_dict[key], np.zeros(self.nData-len(data_dict[key]))))
                    if self._saveint:
                        if not self._saveintadu and key == 'data': continue
                        ret_dict[key] = ret_dict[key].astype(int)
        else:
            for key, d in data_dict.items():
                if key[0]=='_': continue
                ret_dict[f'var_{key}'] = d

        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]

        return ret_dict



class unsparsifyFunc(DetObjectFunc):
    """
    Function to create a 'dense'/raw data shaped array out of the sparsified dict
    """
    def __init__(self, **kwargs):
        self._name = kwargs.get('name','dense')
        super(unsparsifyFunc, self).__init__(**kwargs)

    def setFromDet(self, det):
        super(unsparsifyFunc, self).setFromDet(det)
        self._shape = det.mask.shape

    def process(self, data):
        #make sure we have dict w/ data
        if not isinstance(data, dict):
            if not (isinstance(data, np.ma.masked_array) or
                    isinstance(data, np.ndarray)):
                logger.warning('Try to unsparsify data that is neither a dict nor alraedy an array: {}'.format(data.dtype))
            return

        if 'data' not in data or 'row' not in data or 'col' not in data:
            print('cannot make a make a sparse, rectangular array of', data)
            return

        #now unsparsify the data
        self.dat = unsparsify(data, self._shape)
        ret_dict={}
        subfuncResults = self.processFuncs()
        for k in subfuncResults:
            for kk in subfuncResults[k]:
                ret_dict['%s_%s'%(k,kk)] = subfuncResults[k][kk]

        return ret_dict
