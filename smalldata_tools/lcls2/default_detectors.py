# importing generic python modules
import psana
import numpy as np
from abc import ABCMeta, abstractmethod

from smalldata_tools.common.detector_base import DefaultDetector_base


class DefaultDetector(DefaultDetector_base, metaclass=ABCMeta):
    def __init__(self, detname, name, run):
        self._run = run
        super().__init__(detname, name)
        # self.name = name
        # self.detname = detname
        # self._debug = False

        # if self.in_run():
        # self.det = self._run.Detector(detname)
        if not hasattr(self, "_veto_fields"):
            self._veto_fields = ["TypeId", "Version", "config"]

    def in_run(self):
        detnames = self._run.detnames
        for dn in detnames:
            if dn == self.detname:
                return True
        return False

    def get_psana_det(self):
        return self._run.Detector(self.detname)

    def get_fields(self, top_field):
        """
        Parameters
        ----------
        Top field: str
            Usually 'fex' or 'raw'
        """
        top_field_obj = getattr(self.det, top_field)
        if top_field_obj is None:
            return top_field_obj, None
        fields = [
            field
            for field in dir(top_field_obj)
            if (field[0] != "_" and field not in self._veto_fields)
        ]
        return top_field_obj, fields

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

    # @abstractmethod
    # def data(self, evt):
    #     """ Method that should return a dict of values from event """
    #     pass


def detOnceData(det, data, ts, noarch):
    if not noarch:
        # If we've missed something, look for it in the archiver.
        # ts is our timestamp.
        data = addArchiverData(det, data, ts)
    return data


class genericDetector(DefaultDetector):
    def __init__(self, detname, run=None, name=None):
        if name is None:
            name = detname
        super().__init__(detname, name, run)
        if self.in_run():
            self._data_field, self._sub_fields = self.get_fields("raw")

    def data(self, evt):
        dl = {}
        # raw = getattr( self.det, 'raw')
        # vetolist = ['TypeId', 'Version', 'config']
        if self._data_field is not None:
            #    fields = [ field for field in dir(raw) if (field[0]!='_' and field not in vetolist) ]
            for field in self._sub_fields:
                dat = getattr(self._data_field, field)(evt)
                if dat is None:
                    continue
                dl[field] = dat
                if isinstance(dl[field], list):
                    dl[field] = np.array(dl[field])
        return dl


class ttDetector(DefaultDetector):
    def __init__(self, detname, run=None, saveTraces=False):
        self._veto_fields = ["TypeId", "Version", "calib", "image", "raw", "config"]
        super().__init__(detname, "tt", run)
        self.saveTraces = saveTraces
        if self.in_run():
            self._data_field, self._sub_fields = self.get_fields("ttfex")
            if self.saveTraces:
                self.proj_field, self.proj__sub_fields = self.get_fields("ttproj")

    def data(self, evt):
        dl = {}
        # fex = getattr(self.det, 'ttfex')
        # _veto_fields = ['TypeId', 'Version', 'calib', 'image', 'raw', 'config' ]
        if self._data_field is not None:
            # fields = self.get_fields('ttfex')
            # fields = [ field for field in dir(fex) if (field[0]!='_' and field not in _veto_fields) ]
            for field in self._sub_fields:
                dat = getattr(self._data_field, field)(evt)
                if dat is None:
                    continue
                dl[field] = dat
                if isinstance(dl[field], list):
                    dl[field] = np.array(dl[field])

        if self.saveTraces:
            # fex = getattr( self.det, 'ttproj')
            # _veto_fields = ['TypeId', 'Version', 'calib', 'image', 'raw', 'config' ]
            if self.proj_field is not None:
                #    fields = [ field for field in dir(fex) if (field[0]!='_' and field not in _veto_fields) ]
                for field in self.proj__sub_fields:
                    dat = getattr(self.proj_field, field)(evt)
                    if dat is None:
                        continue
                    dl[field] = dat
                    if isinstance(dl[field], list):
                        dl[field] = np.array(dl[field])
        return dl


class fimfexDetector(DefaultDetector):
    def __init__(self, detname, run=None, name=None):
        if name is None:
            name = detname
        self._veto_fields = ["TypeId", "Version", "calib", "image", "raw", "config"]
        super().__init__(detname, name, run)
        if self.in_run():
            self._data_field, self._sub_fields = self.get_fields("fex")

    def data(self, evt):
        dl = {}
        if self._data_field is not None:
            for field in self._sub_fields:
                dat = getattr(self._data_field, field)(evt)
                if dat is None:
                    continue
                dl[field] = dat
                if isinstance(dl[field], list):
                    dl[field] = np.array(dl[field])
        return dl


class lightStatus(DefaultDetector):
    def __init__(self, destination, laser_codes, run):
        """
        Parameters
        ----------
        destination : <enum 'BeamDestination'>
            X-ray beam destination. See enum definition in hutch_default.py
        laser_codes : list of int
            Positive values are event codes marking dropped shots
            Negative values are event codes marking requested shots
        """
        super().__init__("timing", "lightStatus", run)
        self.destination = destination
        print(f"Beam destination target: {destination}")
        self.laserCodes_drop = [c for c in laser_codes if c > 0]
        self.laserCodes_req = [-c for c in laser_codes if c > 0]

    def data(self, evt):
        xfel_status, laser_status = (0, 1)  # Can we make laser default to 0 as well?
        # x-ray
        if self.det.raw.destination(evt) == self.destination.value:
            xfel_status = 1

        # laser
        event_codes = self.det.raw.eventcodes(evt)

        for lOff in self.laserCodes_drop:
            if evtCodes[lOff]:
                laser_status = 0
        if len(self.laserCodes_req) > 0:
            laser_status = 0
            for lOff in self.laserCodes_req and laser_status == 1:
                if evtCodes[lOff]:
                    laser_status = 1

        dl = {"xray": xfel_status, "laser": laser_status}
        return dl


class epicsDetector(DefaultDetector):
    ### TO REVISE / REVIEW for (det)name, super, etc. Does it inherit DefaultDetector?
    def __init__(self, PVlist=[], run=None, name="epics"):
        # super().__init__('epics', name, run)
        self.name = name
        self.detname = "epics"
        self.PVlist = []
        self.missing = []
        self.missingPV = []
        self.pvs = []
        # run.epicsinfo is an alias or (alias, pvname) --> pvname dict.
        # Let's rearrange this somewhat to make it more useful.
        self.alias_2_pv = {}
        self.pv_to_alias = {}

        for k, pv in run.epicsinfo:
            if type(k) == tuple:
                al = k[0]
            else:
                al = k
            self.alias_2_pv[al] = pv
            self.pv_to_alias[pv] = al

        for p in PVlist:
            try:
                if type(p) == tuple:
                    # User specified PV and alias.
                    pv = p[0]
                    al = p[1]
                    self.pv_to_alias[pv] = al
                    self.alias_2_pv[al] = pv
                elif p in self.alias_2_pv.keys():
                    # Known Alias
                    al = p
                    pv = self.alias_2_pv[al]
                elif p in self.pv_to_alias.keys():
                    # Known PV
                    pv = p
                    al = self.pv_to_alias[pv]
                else:
                    # We don't know.  Assume it's a PV we'll find later.
                    pv = p
                    al = p
                    self.pv_to_alias[pv] = al
                    self.alias_2_pv[al] = pv
                self.pvs.append(run.Detector(pv))
                self.addPV(al, pv)
            except:
                print("could not find LCLS2 EPICS PV %s in data" % pv)
                self.missing.append(al)
                self.missingPV.append(pv)
        # Add these now so they don't interfere in the above loop.
        for al, pv in run.epicsinfo:
            self.alias_2_pv[pv] = pv

    def addPV(self, al, pv):
        self.PVlist.append(al)

    def in_run(self):
        if len(self.pvs) > 0:
            return True
        return False

    def data(self, evt):
        dl = {}
        for pvname, pv in zip(self.PVlist, self.pvs):
            try:
                if pv(evt) is not None:
                    dl[pvname] = pv(evt)
                    if isinstance(dl[pvname], str):
                        dl[pvname] = np.nan
            except:
                # print('we have issues with %s in this event'%pvname)
                pass
        return dl

    def params_as_dict(self):
        """returns parameters as dictionary to be stored in the hdf5 file (once/file)"""
        parList = {
            key: self.__dict__[key]
            for key in self.__dict__
            if (
                key[0] != "_"
                and isinstance(getattr(self, key), (str, int, float, np.ndarray, tuple))
            )
        }
        PVlist = getattr(self, "PVlist")
        parList.update(
            {"PV_%d" % ipv: pv for ipv, pv in enumerate(PVlist) if pv is not None}
        )
        parList.update(
            {
                "PVname_%d" % ipv: self.alias_2_pv[pv]
                for ipv, pv in enumerate(PVlist)
                if pv is not None
            }
        )
        return parList


class scanDetector(DefaultDetector):  # ##### NEED FIX
    def __init__(self, detname="scan", run=None, name=None):
        name = name or detname
        super().__init__(detname, name, run)
        self.scans = []
        self.scanlist = []
        vetolist = ["step_docstring"]
        try:
            scanlist = [k[0] for k in run.scaninfo if k[0] not in vetolist]
            for scan in scanlist:
                try:
                    self.scans.append(run.Detector(scan))
                    self.scanlist.append(scan)
                except:
                    print("could not find LCLS2 EPICS PV %s in data" % pv)
        except:
            pass

    def in_run(self):
        if len(self.scans) > 0:
            return True
        return False

    def data(self, evt):
        dl = {}
        for scanname, scan in zip(self.scanlist, self.scans):
            try:
                if scan(evt) is not None:
                    dl[scanname] = scan(evt)
                    if isinstance(dl[scanname], str):
                        dl[scanname] = np.nan
            except:
                # print('we have issues with %s in this event'%scanname)
                pass
        return dl
