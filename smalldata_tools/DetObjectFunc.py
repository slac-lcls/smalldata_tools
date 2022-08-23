import os
import copy
import numpy as np

from future.utils import iteritems
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
try:
    basestring
except NameError:
    basestring = str

    
#epix10k: thermisotr value to temp in C
def getThermistorTemp(x):
    if x==0: return 0.
    u = x/16383.0 * 2.5
    i = u / 100000
    r = (2.5 - u)/i
    try:
        l = np.log(r/10000)
        t = 1.0 / (3.3538646E-03 + 2.5654090E-04 * l + 1.9243889E-06 * (l*l) + 1.0969244E-07 * (l*l*l))
        return t - 273.15
    except:
        return np.nan

#this class is a container which will hold the event based data. It will be created in the getData step.            
class event(object):
    pass

#
# implement DetObjectFunc as member of DetObjectFunc: e.g. apply a projection on a rotated ROI
# need to make sure which results should be stored temporary and which should go into the hdf5 file.
#
class DetObjectFunc(object):
    def __init__(self, **kwargs):
        self._debug = False
        self._proc = True # Will be set to False in cube only
        if '_name' not in kwargs and '_name' not in self.__dict__.keys():
            print('Function needs to have _name as parameter, will return None')
            return None
        for key in kwargs:
            self.__dict__[key] = kwargs[key]
    def setFromDet(self, det):
        for k, sfunc in iteritems(self.__dict__): 
            if isinstance(sfunc, DetObjectFunc):
                sfunc.setFromDet(det) #pass parameters from det (rms, geometry, .....)

    def setFromFunc(self, parentFunc=None):
        for k, sfunc in iteritems(self.__dict__): 
            if isinstance(sfunc, DetObjectFunc):
                sfunc.setFromFunc(self) #pass parameters from function (rms, boundaries, .....)

    def params_as_dict(self):
        """returns parameters as dictionary to be stored in the hdf5 file (once/file)"""
        #print('DEBUG get params for func -- params_as_dict : ',self._name)
        subFuncs = [ self.__dict__[key] for key in self.__dict__ if isinstance(self.__dict__[key], DetObjectFunc) ]
        funcPars={}
        for tfunc in subFuncs:
            sfuncDict = tfunc.params_as_dict()
            for key in sfuncDict:
                #print ('DEBUG get params for func - subfunc - key: ',self._name, tfunc._name, key)
                #funcPars['%s_%s_%s'%(self._name,tfunc._name,key)] = sfuncDict[key]
                funcPars['%s'%(key)] = sfuncDict[key]

        parList =  {key: self.__dict__[key] for key in self.__dict__ if (isinstance(getattr(self,key), (basestring, int, float, np.ndarray)) and key[0]!='_')}
        parList.update({key: np.array(self.__dict__[key]) for key in self.__dict__ if (isinstance(getattr(self,key), list) and key[0]!='_')})
        remKeys = [key for key in self.__dict__ if (key not in parList)]
        for key, value in iteritems(parList):
            funcPars['%s_%s'%(self._name, key)] = value
        if self._debug: print('DEBUG: keys which are not parameters:',self._name, remKeys)
        #print 'return for ',self._name, funcPars.keys()
        return funcPars

    def setDebug(self, debug):
        if isinstance(debug, bool):
            self._debug = debug
    def process(self, data):
        """returns results as dictionary to be stored in the hdf5 file (each event)"""
        return {}
    def addFunc(self, func):
        self.__dict__[func._name] = func
            
    def setKeyData(self, key, data):
        try:
            setattr(self, key, data)
        except:
            print('cound not set attribute %s of %s to:'%(key, self._name), data)
    def processFuncs(self):
        subFuncs = [ self.__dict__[key] for key in self.__dict__ if isinstance(self.__dict__[key], DetObjectFunc) ]
        subFuncResults={}
        if 'dat' not in self.__dict__.keys() and len(subFuncs)>0:
            print('cannot process subfunctions for %s as data is not being passed'%self.name)
            return
        for tfunc in subFuncs:
            subFuncResults[tfunc._name] = tfunc.process(self.dat)
        return subFuncResults

