import logging
import numpy as np
from dataclasses import dataclass
from future.utils import iteritems
from abc import ABCMeta, abstractmethod

logger = logging.getLogger(__name__)


class DetObject_base(metaclass=ABCMeta):
    pass


@dataclass
class Event:
    # This class is a container which will hold the event based data.
    # It will be created in the getData step.
    pass


def detData(detList, evt):
    """
    Retrieve all data for a given event and set of detectors
    Parameters
    ----------
        detList (list): List of detector. Each of them must be a subclass of
                        DefaultDetector_base
        evt: psana.Event
    """
    data = {}
    for det in detList:
        try:
            data[det.name] = det.data(evt)
        except:
            # print('could not get data in this event for detector ',det.name)
            pass
    return data


class DefaultDetector_base(metaclass=ABCMeta):
    def __init__(self, detname, name):
        self.name = name
        self.detname = detname
        self._debug = False
        self.det = None  # Move definition to relevant sub-classes?

        if self.in_run():
            self.det = self.get_psana_det()

    @abstractmethod
    def in_run(self):
        """Returns whether the detector is available for this data source."""
        pass
        # dNames=[]
        # try:
        #     detnames = psana.DetNames()
        #     for dn in detnames:
        #         for dnn in dn:
        #             if dnn!='':
        #                 dNames.append(dnn)
        # except:
        #     detnames = self._run.detinfo
        #     for dn in detnames:
        #         dNames.append(dn[0])
        # if self.detname in dNames:
        #     return True
        # return False

    @abstractmethod
    def get_psana_det(self):
        """Method to return the underlying detector object"""
        pass

    def _setDebug(self, debug):
        self._debug = debug

    def params_as_dict(self):
        """Returns parameters as dictionary to be stored in the hdf5 file (once/file)"""
        parList = {
            key: self.__dict__[key]
            for key in self.__dict__
            if (
                key[0] != "_"
                and isinstance(getattr(self, key), (str, int, float, np.ndarray, tuple))
            )
        }
        parList.update(
            {
                key: np.array(self.__dict__[key])
                for key in self.__dict__
                if (
                    key[0] != "_"
                    and isinstance(getattr(self, key), list)
                    and len(getattr(self, key)) > 0
                    and isinstance(getattr(self, key)[0], (str, int, float, np.ndarray))
                )
            }
        )
        return parList

    @abstractmethod
    def data(self, evt):
        """Method that should return a dict of values from event"""
        pass


class DetObjectFunc(object):
    def __init__(self, **kwargs):
        self._debug = False
        self._proc = True  # Will be set to False in cube only
        if "_name" not in kwargs and "_name" not in self.__dict__.keys():
            print("Function needs to have _name as parameter, will return None")
            return None
        for key in kwargs:
            self.__dict__[key] = kwargs[key]

    def setFromDet(self, det):
        for k, sfunc in iteritems(self.__dict__):
            if isinstance(sfunc, DetObjectFunc):
                sfunc.setFromDet(det)  # pass parameters from det (rms, geometry, .....)

    def setFromFunc(self, parentFunc=None):
        for k, sfunc in iteritems(self.__dict__):
            if isinstance(sfunc, DetObjectFunc):
                sfunc.setFromFunc(
                    self
                )  # pass parameters from function (rms, boundaries, .....)

    def params_as_dict(self):
        """returns parameters as dictionary to be stored in the hdf5 file (once/file)"""
        # print('DEBUG get params for func -- params_as_dict : ',self._name)
        subFuncs = [
            self.__dict__[key]
            for key in self.__dict__
            if isinstance(self.__dict__[key], DetObjectFunc)
        ]
        funcPars = {}
        for tfunc in subFuncs:
            sfuncDict = tfunc.params_as_dict()
            for key in sfuncDict:
                # print ('DEBUG get params for func - subfunc - key: ',self._name, tfunc._name, key)
                # funcPars['%s_%s_%s'%(self._name,tfunc._name,key)] = sfuncDict[key]
                funcPars["%s" % (key)] = sfuncDict[key]

        parList = {
            key: self.__dict__[key]
            for key in self.__dict__
            if (
                isinstance(getattr(self, key), (str, int, float, np.ndarray))
                and key[0] != "_"
            )
        }
        parList.update(
            {
                key: np.array(self.__dict__[key])
                for key in self.__dict__
                if (isinstance(getattr(self, key), list) and key[0] != "_")
            }
        )
        remKeys = [key for key in self.__dict__ if (key not in parList)]
        for key, value in iteritems(parList):
            funcPars["%s_%s" % (self._name, key)] = value
        if self._debug:
            print("DEBUG: keys which are not parameters:", self._name, remKeys)
        # print 'return for ',self._name, funcPars.keys()
        return funcPars

    def setDebug(self, debug):
        if isinstance(debug, bool):
            self._debug = debug

    def process(self, data):
        """
        Returns results as dictionary to be stored in the hdf5 file (each event).
        To be overwritten in sub-classes
        """
        return {}

    def addFunc(self, func):
        self.__dict__[func._name] = func

    def setKeyData(self, key, data):
        try:
            setattr(self, key, data)
        except:
            print("cound not set attribute %s of %s to:" % (key, self._name), data)

    def processFuncs(self):
        subFuncs = [
            self.__dict__[key]
            for key in self.__dict__
            if isinstance(self.__dict__[key], DetObjectFunc)
        ]
        subFuncResults = {}
        if "dat" not in self.__dict__.keys() and len(subFuncs) > 0:
            print(
                "Cannot process subfunctions for %s as data is not being passed"
                % self.name
            )
            return
        for tfunc in subFuncs:
            subFuncResults[tfunc._name] = tfunc.process(self.dat)
        return subFuncResults


# Utility functions to retrive data etc
def getUserData(det):
    """Return dictionary with event-based user data from input detector, apply mask here."""
    det_dict = {}
    try:
        userData_keys = [
            key for key in det.evt.__dict__.keys() if key.find("_write_") >= 0
        ]
        for key in userData_keys:
            if isinstance(getattr(det.evt, key), dict):
                for kkey in getattr(det.evt, key).keys():
                    fieldName = ("%s_%s" % (key, kkey)).replace("_write_", "")
                    if isinstance(getattr(det.evt, key)[kkey], tuple):
                        det_dict[fieldName] = np.array(getattr(det.evt, key)[kkey])
                    else:
                        det_dict[fieldName] = getattr(det.evt, key)[kkey]
            elif isinstance(det.evt.__dict__[key], np.ndarray):
                if isinstance(det.evt.__dict__[key], np.ma.masked_array):
                    data = det.evt.__dict__[key].data
                    if np.issubdtype(data.dtype, np.integer):
                        data[det.evt.__dict__[key].mask] = 0
                    else:
                        data[det.evt.__dict__[key].mask] = np.nan
                    det_dict[key.replace("_write_", "")] = data
                else:
                    det_dict[key.replace("_write_", "")] = det.evt.__dict__[key]
            else:
                det_dict[key.replace("_write_", "")] = np.array(det.evt.__dict__[key])
    except:
        logger.info(f"Could not retrieve data for det {det}")
        pass
    newDict = {}
    for key in det_dict:
        if key.find("ragged") >= 0:
            newDict["ragged_%s" % (key.replace("ragged_", ""))] = det_dict[key]
        elif key.find("var") >= 0:
            fname = key.split("_")[-1]
            dname = key.replace(f"_{fname}", "").replace("_var", "")
            if "var_%s" % dname in newDict:
                newDict[f"var_{dname}"][fname] = det_dict[key]
            else:
                newDict[f"var_{dname}"] = {fname: det_dict[key]}
        else:
            newDict["%s" % (key)] = det_dict[key]
    det_dict = newDict
    return det_dict


def getUserEnvData(det):
    """Return environment data for a detectors as a dictionary (temperatures,....)"""
    env_dict = {}
    try:
        envData_keys = [key for key in det.evt.__dict__.keys() if key.find("env_") >= 0]
        for key in envData_keys:
            if isinstance(det.evt.__dict__[key], np.ndarray):
                if isinstance(det.evt.__dict__[key], np.ma.masked_array):
                    data = det.evt.__dict__[key].data
                    data[det.evt.__dict__[key].mask] = np.nan
                    env_dict[key.replace("env_", "")] = data
                else:
                    env_dict[key.replace("env_", "")] = det.evt.__dict__[key]
            else:
                env_dict[key.replace("env_", "")] = np.array(det.evt.__dict__[key])
    except:
        pass
    return env_dict
