# importing generic python modules
import psana
import numpy as np
from abc import ABCMeta, abstractmethod
import logging

logger = logging.getLogger(__name__)

from smalldata_tools.common.detector_base import DefaultDetector_base
from smalldata_tools.utilities import printR


from mpi4py import MPI

rank = MPI.COMM_WORLD.Get_rank()


class DefaultDetector(DefaultDetector_base, metaclass=ABCMeta):
    def __init__(self, detname, name, run):
        self._run = run
        super().__init__(detname, name)

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


class damageDetector(DefaultDetector):
    def __init__(self, detname="damage", run=None):
        self.name = detname
        self.detname = detname
        self.detNames = []
        self.damageDets = []
        if run is None:
            return
        # run.detnames: set[str]
        for dn in run.detnames:
            if dn in ("epicsinfo", "triginfo") or "pvdetinfo" in dn:
                continue
            self.detNames.append(dn)
            det = run.Detector(dn)
            self.damageDets.append(psana.detector.damage.Damage(det.raw))
        self.detAlias = [det for det in self.detNames]

    def in_run(self):
        return True

    def data(self, evt):
        dl = {}
        for idx, detName in enumerate(self.detNames):
            damage = self.damageDets[idx].count(evt)
            val = 0
            if damage is None:
                val = 0
            elif damage == {}:
                val = 1
            else:
                val = 0
            dl[detName] = val
        return dl


class genericDetector(DefaultDetector):
    """
    Detector to get whatever is under the 'raw' field, minus the default
    veto fields.
    """

    def __init__(self, detname, run=None, name=None, var_fields=None):
        if name is None:
            name = detname
        super().__init__(detname, name, run)
        self._var_fields = var_fields
        if self.in_run():
            self._data_field, self._sub_fields = self.get_fields("raw")

    def data(self, evt):
        dl = {}
        if self._data_field is not None:
            for field in self._sub_fields:
                if self._var_fields is not None and field in self._var_fields:
                    field_name = f"var_{field}"
                else:
                    field_name = field
                dat = getattr(self._data_field, field)(evt)
                if dat is None:
                    continue
                dl[field_name] = dat
                if isinstance(dl[field_name], list):
                    dl[field_name] = np.array(dl[field_name])
        return dl


class ttDetector(DefaultDetector):
    def __init__(self, detname, run=None, saveTraces=False, iocTimetool=False):
        self._veto_fields = ["TypeId", "Version", "calib", "image", "raw", "config"]
        super().__init__(detname, "tt", run)
        self.saveTraces = saveTraces
        self.iocTimetool = iocTimetool
        if self.in_run():
            if iocTimetool:
                fields_name = "raw"
            else:
                fields_name = "ttfex"
            self._data_field, self._sub_fields = self.get_fields(fields_name)
            if self.saveTraces and not self.iocTimetool:
                self.proj_field, self.proj_sub_fields = self.get_fields("ttproj")

    def data(self, evt):
        dl = {}
        if self._data_field is not None:
            for field in self._sub_fields:
                dat = getattr(self._data_field, field)(evt)
                if dat is None:
                    continue
                if self.iocTimetool:
                    if field == "value":
                        dl["fltpos"] = dat[0]
                        dl["fltpos_ps"] = dat[1]
                        dl["fltposfwhm"] = dat[2]
                        dl["ampl"] = dat[3]
                        dl["amplnxt"] = dat[4]
                        dl["refampl"] = dat[5]
                        return dl
                    else:
                        continue
                dl[field] = dat
                if isinstance(dl[field], list):
                    dl[field] = np.array(dl[field])

        if self.saveTraces:
            if self.proj_field is not None:
                for field in self.proj_sub_fields:
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


class usbEncoder(DefaultDetector):
    """
    Detector for US Digital USB encoder bld.
    """

    def __init__(self, detname, run=None, name="enc"):
        super().__init__(detname, name, run)
        if self.in_run():
            self._data_field, self._sub_fields = self.get_fields("raw")

    def data(self, evt):
        dl = {}
        if self._data_field is not None and self._sub_fields is not None:
            channel_descriptions = []
            data = []
            for field in self._sub_fields:
                if field == "descriptions":
                    channel_descriptions = getattr(self._data_field, field)(evt)
                else:
                    data = getattr(self._data_field, field)(evt)

            if data is not None:
                for i in range(4):
                    if (desc := channel_descriptions[i]) != "":
                        dl[desc] = data[desc]
        return dl


class interpolatedEncoder(DefaultDetector):
    """
    Detector for interpolated encoder.
    """

    def __init__(self, detname, run=None, name=None):
        if name is None:
            name = detname
        super().__init__(detname, name, run)
        if self.in_run():
            self._raw_field, self._raw_sub_fields = self.get_fields("raw")
            self._interp_field, self._interp_sub_fields = self.get_fields(
                "interpolated"
            )

    def data(self, evt):
        dl = {}
        if self._raw_field is not None:
            for field in self._raw_sub_fields:
                dat = getattr(self._raw_field, field)(evt)
                if dat is None:
                    continue
                dl[f"raw_{field}"] = dat
                if isinstance(dl[f"raw_{field}"], list):
                    dl[field] = np.array(dl[f"raw_{field}"])

        if self._interp_field is not None:
            for field in self._interp_sub_fields:
                dat = getattr(self._interp_field, field)(evt)
                if dat is None:
                    continue
                dl[f"interpolated_{field}"] = dat
                if isinstance(dl[f"interpolated_{field}"], list):
                    dl[field] = np.array(dl[f"interpolated_{field}"])
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
        if rank == 0:
            logger.info(f"Beam destination target: {destination}")
        self.laserCodes_drop = [c for c in laser_codes if c > 0]
        self.laserCodes_req = [-c for c in laser_codes if c < 0]

    def data(self, evt):
        xfel_status, laser_status = (0, 1)  # Can we make laser default to 0 as well?
        # x-ray
        if self.det.raw.destination(evt) == self.destination.value:
            xfel_status = 1

        # laser
        event_codes = self.det.raw.eventcodes(evt)

        for lOff in self.laserCodes_drop:
            laser_status = 1
            if event_codes[lOff]:
                laser_status = 0
        if len(self.laserCodes_req) > 0:
            laser_status = 0
            for lOff in self.laserCodes_req:
                if event_codes[lOff]:
                    laser_status = 1

        dl = {"xray": xfel_status, "laser": laser_status}
        return dl


class lightStatusLcls1Timing(DefaultDetector):
    def __init__(self, beam_las_codes, run):
        """
        Parameters
        ----------
        destination : <enum 'BeamDestination'>
            X-ray beam destination. See enum definition in hutch_default.py
        beam_las_codes : list of list of int
            First list is for x-ray codes, second list is for laser codes.
            Positive values are event codes marking dropped shots
            Negative values are event codes marking requested shots
        """
        super().__init__("timing", "lightStatus", run)

        beam_codes = beam_las_codes[0]
        laser_codes = beam_las_codes[1]

        self.laserCodes_drop = [c for c in laser_codes if c > 0]
        self.laserCodes_req = [-c for c in laser_codes if c < 0]

        self.beamCodes_drop = [c for c in beam_codes if c > 0]
        self.beamCodes_req = [-c for c in beam_codes if c < 0]

    def data(self, evt):
        xfel_status, laser_status = (0, 1)  # Can we make laser default to 0 as well?
        # laser
        event_codes = self.det.raw.eventcodes(evt)

        for lOff in self.laserCodes_drop:
            laser_status = 1
            if event_codes[lOff]:
                laser_status = 0
        if len(self.laserCodes_req) > 0:
            laser_status = 0
            for lOn in self.laserCodes_req:
                if event_codes[lOn]:
                    laser_status = 1

        for bOff in self.beamCodes_drop:
            xfel_status = 1
            if event_codes[bOff]:
                xfel_status = 0
        if len(self.beamCodes_req) > 0:
            xfel_status = 0
            for bOn in self.beamCodes_req:
                if event_codes[bOn]:
                    xfel_status = 1

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


class scanDetector(DefaultDetector):
    def __init__(self, detname="scan", run=None, name=None):
        self._veto_fields = ["step_docstring"]
        self.scans = []
        self.scanlist = []
        name = name or detname
        super().__init__(detname, name, run)

    def in_run(self):
        if self._run.scaninfo == {}:
            return False
        return True

    def get_psana_det(self):
        scanlist = [k[0] for k in self._run.scaninfo if k[0] not in self._veto_fields]
        for scan in scanlist:
            self.scans.append(self._run.Detector(scan))
            self.scanlist.append(scan)

    def data(self, evt):
        dl = {}
        for scanname, scan in zip(self.scanlist, self.scans):
            try:
                if scan(evt) is not None:
                    dl[scanname] = scan(evt)
                    if isinstance(dl[scanname], str):
                        dl[scanname] = np.nan
            except:
                pass
        return dl
