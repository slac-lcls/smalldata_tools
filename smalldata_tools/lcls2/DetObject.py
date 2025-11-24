import os
import copy
import numpy as np
import logging
from psana.pscalib.calib.MDBWebUtils import calib_constants
from smalldata_tools.common.detector_base import DetObjectFunc, Event
from smalldata_tools.utilities import cm_epix
from mpi4py import MPI

rank = MPI.COMM_WORLD.Get_rank()

import psana

logger = logging.getLogger(__name__)


def DetObject(srcName, run, **kwargs):
    if rank == 0:
        logger.info(f"Getting the detector for: {srcName}")
    det = None
    try:
        det = run.Detector(srcName)
    except:
        if rank == 0:
            logger.warning(f"failed to make detector for {srcName}: {run.detnames}")
        return NullDetObject(name=srcName)
    det.alias = srcName
    detector_classes = {
        "epix100": Epix100Object,
        "epix10ka": Epix10kObject,
        "opal": SimpleCameraObject,
        "hsd": HsdObject,
        "wave8": Wave8Object,
        "pv": PVObject,
        "piranha4": PiranhaObject,
        "archon": ArchonObject,
        "jungfrau": JungfrauObject,
        "axis": SimpleCameraObject,
        "generic_container": GenericContainerObject,
    }

    cls = detector_classes[det._dettype]
    return cls(det, run, **kwargs)


class DetObjectClass(object):
    def __init__(
        self, det, run, **kwargs
    ):  # name=None, common_mode=None, applyMask=0):
        self.det = det
        self._detid = det._detid
        self._name = kwargs.get("name", self.det._det_name)  # srcName)

        self.run = run
        self._storeSum = {}
        self.applyMask = kwargs.get("applyMask", 0)

        self.dataAccessTime = 0.0

    def params_as_dict(self):
        """Returns parameters as dictionary to be stored in the hdf5 file (once/file)"""
        parList = {
            key: self.__dict__[key]
            for key in self.__dict__
            if (
                key[0] != "_"
                and isinstance(getattr(self, key), (str, int, float, np.ndarray))
            )
        }
        parList.update(
            {
                key: np.array(self.__dict__[key])
                for key in self.__dict__
                if (
                    key[0] != "_"
                    and isinstance(getattr(self, key), tuple)
                    and isinstance(getattr(self, key)[0], (str, int, float, np.ndarray))
                )
            }
        )
        parList.update(
            {
                key: np.array(self.__dict__[key])
                for key in self.__dict__
                if (
                    key[0] != "_"
                    and isinstance(getattr(self, key), list)
                    and getattr(self, key)  # Check for empty list
                    and isinstance(getattr(self, key)[0], (str, int, float, np.ndarray))
                )
            }
        )

        # Add parameters of function to dict with composite keyt(sf._name, key)
        subFuncs = [
            self.__dict__[key]
            for key in self.__dict__
            if isinstance(self.__dict__[key], DetObjectFunc)
        ]
        # We cannot save arrays of chars/strings. There should be a proper test instead of this explicit rejection of coords.
        # I can't see how to do that in a single line, but that is not a great reason...
        for sf in subFuncs:
            sfPars = sf.params_as_dict()
            parList.update(
                {
                    "%s__%s" % (sf._name, key): value
                    for key, value in sfPars.items()
                    if (
                        key[0] != "_"
                        and isinstance(value, (str, int, float, np.ndarray, tuple))
                        and key.find("coords") < 0
                    )
                }
            )
        return parList

    # here we should get masks from the calibstore once we learn how to.
    def _getMasks(self):
        self.mask = None
        self.cmask = None

    def _applyMask(self):
        try:
            if self.applyMask == 1:
                self.evt.dat[self.mask == 0] = 0
            if self.applyMask == 2:
                self.evt.dat[self.cmask == 0] = 0
        except:
            print("Could not apply mask to data for detector ", self.name)

    def storeSum(self, sumAlgo=None):
        if sumAlgo is not None:
            self._storeSum[sumAlgo] = None
        else:
            return self._storeSum

    def setPed(self, ped):
        self.ped = ped

    def setMask(self, mask):
        self.mask = mask

    def setcMask(self, mask):
        self.cmask = np.amin(np.array([self.mask, mask]), axis=0)

    def setGain(self, gain):
        """
        Set a local gain.
        This file is supposed to be applied on top of whatever corrections DetObject will apply, given the common mode
        """
        self.local_gain = gain

    def getData(self, evt):
        try:
            getattr(self, "evt")
        except:
            self.evt = Event()
        self.evt.dat = None

    def addFunc(self, func):
        func.setFromDet(self)  #  Pass parameters from det (rms, geometry, .....)
        try:
            func.setFromFunc()  # Pass parameters from itself to children (rms, bounds, .....)
        except:
            print("Failed to pass parameters to children of ", func._name)
        self.__dict__[func._name] = func

    def processFuncs(self):
        if self.evt.dat is None:
            logger.debug("This event has no data to be processed for %s" % self._name)
            return
        for func in [
            self.__dict__[k]
            for k in self.__dict__
            if isinstance(self.__dict__[k], DetObjectFunc)
        ]:
            # try:
            if 1:
                retData = func.process(self.evt.dat)
                self.evt.__dict__["_write_%s" % func._name] = retData
            # except Exception as e:
            #     print('Could not run function %s on data of detector %s of shape'%(func._name, self._name), self.evt.dat.shape)
            #     print(e)

    def processSums(self):
        for key in self._storeSum.keys():
            asImg = False
            thres = -1.0e9
            for skey in key.split("_"):
                if skey.find("thresADU") >= 0:
                    thres = float(skey.replace("thresADU", ""))

            if self.evt.dat is None:
                return
            dat_to_be_summed = self.evt.dat
            if thres > 1e-9:
                dat_to_be_summed[self.evt.dat < thres] = 0.0

            if key.find("nhits") >= 0:
                dat_to_be_summed[dat_to_be_summed > 0] = 1

            if key.find("square") >= 0:
                dat_to_be_summed = np.square(dat_to_be_summed)

            if self._storeSum[key] is None:
                if dat_to_be_summed is not None:
                    self._storeSum[key] = dat_to_be_summed.copy()
            else:
                try:
                    if key.find("max") < 0:
                        self._storeSum[key] += dat_to_be_summed
                    else:
                        # Care about nan behaviour?
                        self._storeSum[key] = np.maximum(
                            self._storeSum[key], dat_to_be_summed
                        )
                except:
                    print("could not add ", dat_to_be_summed)
                    print("could not to ", self._storeSum[key])
            # print('%s'%key, self._storeSum[key] )


class NullDetObject:
    """
    A dummy detector object that does nothing.
    Useful to handle instantiation in case the detector is not present in the data.
    """

    def __init__(self, *args, **kwargs):
        self._name = kwargs.get("name", "NullDetObject")
        self.det = None
        self.run = None
        self.evt = Event()
        self._storeSum = {}
        self.applyMask = 0
        self.dataAccessTime = 0.0

    def addFunc(self, func):
        """
        Do nothing, as this is a null object.
        """
        pass


class GenericContainerObject(DetObjectClass):
    def __init__(self, det, run, **kwargs):
        super().__init__(det, run, **kwargs)
        self._epicsinfo = run.epicsinfo
        self._stashed_run = run

        self._is_tiled_camera = False
        self._veto_fields = ["TypeId", "Version", "config"]
        self._searched_for_calib = False

        if hasattr(self.det, "xtc1dump"):
            # Container is a generic coming from an XTC1 to XTC2 conversion
            self._data_field, self._sub_fields = self.get_fields("xtc1dump")
        else:
            raise RuntimeError(
                "Currently only xtc1-xtc2 GenericContainers are supported."
            )

        if "calib" in self._sub_fields:
            # We'll say that this is a camera-type generic
            self._is_tiled_camera = True
            self.ped = None
            self.gain = None
            self.mask = None
            self.cmask = None
            self.rms = None

            self.ix = None
            self.iy = None
            self.x = None
            self.y = None

    def getData(self, evt):
        super().getData(evt)
        if not self._searched_for_calib:
            for name, name0 in self._epicsinfo:
                det = self._stashed_run.Detector(name)
                if self._name in name:
                    dgrams = det._env_store.env_managers[0].dgrams
                    for dgram in dgrams:
                        raw = dgram.epics[0].raw
                        attrs = vars(raw)
                        for attr_name in attrs:
                            if "pedestals" in attr_name:
                                self.ped = getattr(raw, attr_name)
                            if "gain" in attr_name:
                                self.gain = getattr(raw, attr_name)
                            if "mask" in attr_name:
                                self.mask = getattr(raw, attr_name)
                                self.cmask = self.mask
                            if "pixel_index_map" in attr_name:
                                index_map = getattr(raw, attr_name)
                                self.ix = index_map[..., 0]
                                self.iy = index_map[..., 1]
                            if "pixel_position" in attr_name:
                                pixel_pos = getattr(raw, attr_name)
                                self.x = pixel_pos[..., 0]
                                self.y = pixel_pos[..., 1]
                    self._searched_for_calib = True
                    self.reprocessFuncs()
                    break

        if self._is_tiled_camera:
            # We save calibrated data using the conversion process
            # No need to do anything else
            self.evt.dat = self.det.xtc1dump.calib(evt)

    def addFunc(self, func):
        try:
            # Some of these may fail for these containers because mask and such
            # is not available until after first event. We therefor will reprocess
            # afterwards
            super().addFunc(func)
        except:
            ...

    def reprocessFuncs(self):
        """Rerun setup for DetObjectFuncs.

        DetObjectFuncs may need reprocessing for these containers because certain
        information that is available for all other detectors at Configure/BeginRun
        is no available until after the first event. Rather than rework everything for
        this edge-case, the edge-case will just catch errors and rerun the config.
        """
        for func in [
            self.__dict__[k]
            for k in self.__dict__
            if isinstance(self.__dict__[k], DetObjectFunc)
        ]:
            func.setFromDet(self)
            try:
                func.setFromFunc()
            except:
                print("Failed to pass parameters to children of ", func._name)

    def get_fields(self, top_field):
        """
        Parameters
        ----------
        Top field: str
            Usually 'fex' or 'raw'
        """
        # Not ideal duplicating this here.... But don't want to inherit
        # from default detector directly
        top_field_obj = getattr(self.det, top_field)
        if top_field_obj is None:
            return top_field_obj, None
        fields = [
            field
            for field in dir(top_field_obj)
            if (field[0] != "_" and field not in self._veto_fields)
        ]
        return top_field_obj, fields


class CameraObject(DetObjectClass):
    def __init__(self, det, run, **kwargs):
        super(CameraObject, self).__init__(det, run, **kwargs)
        self._common_mode_list = [0, -1, 30]  # none, raw, calib
        self.common_mode = kwargs.get("common_mode", self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if (
            self.common_mode not in self._common_mode_list
            and type(self) is CameraObject
        ):
            print(
                "Common mode %d is not an option for a CameraObject, please choose from: "
                % self.common_mode,
                self._common_mode_list,
            )
        self.pixelsize = None
        self.isGainswitching = False

        # try calibconst...
        # detrawid = det.raw._uniqueid
        # self.peds = calib_constants(detrawid, exp=run.expt, ctype='pedestals', run=run.runnum)[0]
        # self.rms = calib_constants(detrawid, exp=run.expt, ctype='pixel_rms', run=run.runnum)[0]

        try:
            self.ped = det.raw._pedestals()
        except:
            self.ped = None
        try:
            self.rms = det.raw._rms()
        except:
            self.rms = None
        try:
            self.gain = det.raw._gain()
        except:
            self.gain = None
        try:
            self.mask = det.raw._mask(calib=False, status=True, edges=True)
            self.cmask = det.raw._mask(calib=True, status=True, edges=True)
        except:
            self.mask = None
            self.cmask = None
        # self.det.calibconst['pop_rbfs'][1] return meta data for the calib data.
        # self.rms = self.det.rms(run)
        # self.gain_mask = self.det.gain_mask(run)
        # self.gain = self.det.gain(run)
        # self.common_mode_pars=self.det.common_mode(run)
        self.local_gain = None
        self._getImgShape()
        self._gainSwitching = False
        try:
            self.x, self.y, self.z = det.raw._pixel_coords(do_tilt=True, cframe=0)
            self.x = self.x.squeeze()
            self.y = self.y.squeeze()
            self.z = self.z.squeeze()
        except:
            self.x, self.y, self.z = None, None, None

    def getData(self, evt):
        super(CameraObject, self).getData(evt)

    # this is not important until we get to tiled detectors w/ geometry.
    def _getImgShape(self):
        self.imgShape = None


class SimpleCameraObject(CameraObject):
    def __init__(self, det, run, **kwargs):
        super().__init__(det, run, **kwargs)

        # Try to be smart about common mode and use calib if the detector has pedestals.
        if self.ped is not None and self.common_mode == 0:
            self.common_mode = 30

    def getData(self, evt):
        super().getData(evt)
        if self.common_mode < 0:
            self.evt.dat = self.det.raw.raw(evt)
        elif (
            self.common_mode == 0
        ):  # we need to figure out how to do this. Don't implement for, return raw
            self.evt.dat = self.det.raw.raw(evt)
        elif self.common_mode % 100 == 30:
            self.evt.dat = self.det.raw.calib(evt)

        # Override gain if desired
        if (
            self.local_gain is not None
            and self.common_mode in [0, 30]
            and self._gainSwitching is False
            and self.local_gain.shape == self.evt.dat.shape
        ):
            self.evt.dat *= self.local_gain  # apply own gain


class TiledCameraObject(CameraObject):
    def __init__(self, det, run, **kwargs):
        # super().__init__(det,env,run, **kwargs)
        super(TiledCameraObject, self).__init__(det, run, **kwargs)
        try:
            self.ix, self.iy = det.raw._pixel_coord_indexes()
        except:
            if rank == 0:
                print(
                    "failed to get geometry info, likely because we do not have a geometry file"
                )
            self.ix = self.x  # need to change this so ix & iy are integers!
            self.iy = self.y
        self.ix = self.ix.squeeze()
        self.iy = self.iy.squeeze()
        self._needsGeo = True  # FIX ME: not sure it should be here.

    def getData(self, evt):
        super(TiledCameraObject, self).getData(evt)


class Epix100Object(TiledCameraObject):
    def __init__(self, det, run, **kwargs):
        super().__init__(det, run, **kwargs)
        self._common_mode_list = [
            6,
            36,
            4,
            34,
            45,
            46,
            47,
            0,
            -1,
            30,
        ]  # Jacob (norm), Jacob, def, def(ami-like), mine, mine (norm), mine (norm-bank), none, raw, calib
        self.common_mode = kwargs.get("common_mode", self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print(
                "Common mode %d is not an option for as Epix detector, please choose from: "
                % self.common_mode,
                self._common_mode_list,
            )
        self.pixelsize = [50e-6]
        self.areas = None

        if self.rms is None or self.rms.shape != self.ped.shape:
            self.rms = np.ones_like(self.ped)
        try:
            # ok, this call does NOT work (yet) for LCLS2. Needs an event...
            # self.imgShape=det.raw.image(run.runnum, self.ped[0]).shape
            self.imgShape = (self.ix.max(), self.iy.max())
        except:
            if len(self.ped[0].squeeze().shape) == 2:
                self.imgShape = self.ped[0].squeeze().shape
            else:
                self.imgShape = None

        self.bankMasks = []
        if self.common_mode == 47:
            for i in range(0, 16):
                bmask = np.zeros_like(self.rms)
                bmask[
                    (i % 2) * 352 : (i % 2 + 1) * 352,
                    768 / 8 * (i / 2) : 768 / 8 * (i / 2 + 1),
                ] = 1
                self.bankMasks.append(bmask.astype(bool))

    def getData(self, evt):
        super().getData(evt)
        mbits = 0  # do not apply mask (would set pixels to zero)
        # mbits=1 #set bad pixels to 0

        # Common Mode currently NOT working for ePix100
        # Skip all options that pass cmpars
        if self.common_mode % 100 == 6:
            # Was using cmpars=[6] -> This doesn't work for psana2
            # Some of these other kwargs may not do anything
            # self.evt.dat = self.det.raw.calib(
            #    evt, cmpars=(0,6,100,10), mbits=mbits, rms=self.rms, normAll=True
            # )
            self.evt.dat = self.det.raw.calib(
                evt, mbits=mbits, rms=self.rms, normAll=True
            )
        elif self.common_mode % 100 == 36:
            # self.evt.dat = self.det.raw.calib(evt, cmpars=[6], mbits=mbits, rms=self.rms)
            self.evt.dat = self.det.raw.calib(evt, mbits=mbits, rms=self.rms)
        elif self.common_mode % 100 == 34:
            # self.evt.dat = self.det.raw.calib(evt, cmpars=(4, 6, 100, 100), mbits=mbits)
            self.evt.dat = self.det.raw.calib(evt, mbits=mbits)
        elif self.common_mode % 100 == 4:
            # self.evt.dat = self.det.raw.calib(evt, cmpars=(4, 6, 30, 10), mbits=mbits)
            self.evt.dat = self.det.raw.calib(evt, mbits=mbits)
        elif self.common_mode % 100 == 45:
            self.evt.dat = self.det.raw.raw(evt) - self.ped
            self.evt.dat = cm_epix(self.evt.dat, self.rms, mask=self.statusMask)
        elif self.common_mode % 100 == 46:
            self.evt.dat = self.det.raw.raw(evt) - self.ped
            self.evt.dat = cm_epix(
                self.evt.dat, self.rms, normAll=True, mask=self.statusMask
            )
        elif self.common_mode % 100 == 47:
            self.evt.dat = self.det.raw_data(evt) - self.ped
            for _, bMask in enumerate(self.bankMasks):
                self.evt.dat[bMask] -= np.median(self.evt.dat[bMask])
            self.evt.dat = cm_epix(self.evt.dat, self.rms, mask=self.statusMask)

        # override gain if desired
        if (
            self.local_gain is not None
            and self.local_gain.shape == self.evt.dat.shape
            and self.common_mode in [6, 36, 34, 3, 4, 45, 46, 47]
        ):
            self.evt.dat *= self.local_gain  # apply own gain
        elif (
            self.local_gain is None
            and self.gain is not None
            and self.gain.shape == self.evt.dat.shape
            and self.common_mode in [45, 46, 47]
        ):
            self.evt.dat *= self.gain  # apply gain after own common mode correction

        # correct for area of pixels.
        if self.areas is not None:
            self.evt.dat /= self.areas


class Epix10kObject(TiledCameraObject):
    def __init__(self, det, run, **kwargs):
        # super().__init__(det,env,run, **kwargs)
        super(Epix10kObject, self).__init__(det, run, **kwargs)
        self._common_mode_list = [
            80,
            0,
            -1,
            -2,
            30,
            83,
            84,
            85,
        ]  # calib-noCM, ped sub, raw, raw_gain, calib, calib-CMrow, calib-CMrowcol
        self.common_mode = kwargs.get("common_mode", self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print(
                "Common mode %d is not an option for as Epix detector, please choose from: "
                % self.common_mode,
                self._common_mode_list,
            )
        self.pixelsize = [100e-6]
        self.isGainswitching = True

        if self.rms is None or self.rms.shape != self.ped.shape:
            self.rms = np.ones_like(self.ped)
        try:
            # ok, this call does NOT work (yet) for LCLS2. Needs an event...
            # self.imgShape=det.raw.image(run.runnum, self.ped[0]).shape
            self.imgShape = (self.ix.max(), self.iy.max())
        except:
            if len(self.ped[0].squeeze().shape) == 2:
                self.imgShape = self.ped[0].squeeze().shape
            else:
                self.imgShape = None
        self._gainSwitching = True

    def getData(self, evt):
        super(Epix10kObject, self).getData(evt)
        mbits = 0  # do not apply mask (would set pixels to zero)
        # mbits=1 #set bad pixels to 0
        if self.common_mode < 0:
            if self.common_mode == -2:
                self.evt.dat = self.evt.dat
            else:
                self.evt.dat = self.evt.dat & 0x3FFF
            # self.evt.gainbit = (self.evt.dat&0xc000>0)
        elif self.common_mode == 0:
            ##########
            ### FIX ME epix10ka
            # will need to read gain bit from data and use right pedestal.
            # will hopefully get calib function for this.
            ##########
            if len(self.ped.shape) > 3:
                self.evt.dat = (self.det.raw.raw(evt) & 0x3FFF) - self.ped[0]
            else:
                self.evt.dat = (self.det.raw.raw(evt) & 0x3FFF) - self.ped
        elif self.common_mode == 30:
            self.evt.dat = self.det.raw.calib(evt)
        elif self.common_mode == 80:
            self.evt.dat = self.det.raw.calib(evt, cmpars=(7, 0, 100))
        elif self.common_mode == 83:
            self.evt.dat = self.det.raw.calib(evt, cmpars=(7, 1, 100, 10))
        elif self.common_mode == 84:
            self.evt.dat = self.det.raw.calib(evt, cmpars=(7, 2, 100, 10))
        elif self.common_mode == 85:
            self.evt.dat = self.det.raw.calib(evt, cmpars=(7, 3, 100, 10))

        # override gain if desired -- this looks like CsPad.
        if (
            self.local_gain is not None
            and self.local_gain.shape == self.evt.dat.shape
            and self.common_mode in [1, 5, 55, 10]
        ):
            self.evt.dat *= self.local_gain  # apply own gain


class JungfrauObject(TiledCameraObject):
    def __init__(self, det, run, **kwargs):
        super().__init__(det, run, **kwargs)
        self._common_mode_list = [
            0,
            7,
            71,
            72,
            -1,
            30,
        ]  # none, epix-style corr on row*col, row, col, raw, calib
        self.common_mode = kwargs.get("common_mode", self._common_mode_list[0])
        if self.common_mode is None:
            self.common_mode = self._common_mode_list[0]
        if self.common_mode not in self._common_mode_list:
            print(
                "Common mode %d is not an option for Jungfrau, please choose from: "
                % self.common_mode,
                self._common_mode_list,
            )
        self.pixelsize = [75e-6]
        self.isGainswitching = True

        try:
            self.imgShape = self.det.raw.image(run, self.ped[0]).shape
        except:
            pass
        self._gainSwitching = True
        # self._common_mode_list.append()

    def getData(self, evt):
        super(JungfrauObject, self).getData(evt)
        mbits = 0  # do not apply mask (would set pixels to zero)
        # mbits=1 #set bad pixels to 0
        if self.common_mode == 0:
            self.evt.dat = self.det.raw.calib(evt, cmpars=(7, 0, 10), mbits=mbits)
        elif self.common_mode % 100 == 71:
            self.evt.dat = self.det.raw.calib(
                evt, cmpars=(7, 1, 10, 10), mbits=mbits
            )  # correction in rows
        elif self.common_mode % 100 == 72:
            self.evt.dat = self.det.raw.calib(
                evt, cmpars=(7, 2, 10, 10), mbits=mbits
            )  # correction in columns
        elif self.common_mode % 100 == 7:
            self.evt.dat = self.det.raw.calib(
                evt, cmpars=(7, 3, 10, 10), mbits=mbits
            )  # correction in rows&columns

        # override gain if desired
        if (
            self.local_gain is not None
            and self.local_gain.shape == self.evt.dat.shape
            and self.common_mode in [7, 71, 72, 0]
        ):
            self.evt.dat *= self.local_gain  # apply own gain

        # correct for area of pixels.
        # if self.areas is not None:
        #    self.evt.dat /= self.areas
        # self._getRawShape()


class PVObject(CameraObject):
    def __init__(self, det, run, **kwargs):
        super(PVObject, self).__init__(det, run, **kwargs)

    def getData(self, evt):
        super(PVObject, self).getData(evt)
        if self.common_mode < 0:
            self.evt.dat = self.det.raw.value(evt)
        elif (
            self.common_mode == 0
        ):  # we need to figure out how to do this. Don't implement for, return raw
            self.evt.dat = self.det.raw.value(evt)
        elif self.common_mode % 100 == 30:
            self.evt.dat = self.det.raw.value(evt)


class WaveformObject(DetObjectClass):
    def __init__(self, det, run, **kwargs):
        super(WaveformObject, self).__init__(det, run, **kwargs)
        self.common_mode = kwargs.get("common_mode", -1)
        self.ped = None
        self.rms = None
        self.mask = None
        self.wfx = None
        self.gain = None

    def getData(self, evt):
        super(WaveformObject, self).getData(evt)
        # self.evt.dat = self.det.raw.waveform(evt)
        # if self.evt.dat is None:
        #    self.wfx = self.det.wftime(evt)


class HsdObject(WaveformObject):
    def __init__(self, det, run, **kwargs):
        super(HsdObject, self).__init__(det, run, **kwargs)
        self.cidx = [k for k in self.det.raw._seg_chans()]

    def getData(self, evt):
        """
        HSD data handling a bit particular.
        Because the different channels can have different sampling and lengths, the data are
        put in a contiguous array (i.e. 1 dimension vector) with all channels appended to
        one-another.
        The wfxlen variable is used to then isolated individual channel when needed.
        The information on the x axis is gathered on the first event, and cached for the rest
        of the run.
        """
        super(HsdObject, self).getData(evt)
        datadict = self.det.raw.waveforms(evt)
        if datadict is None:
            print("HSD data is None")
            return  # ensure that self.evt.dat is None

        if self.wfx is None:
            self.wfxlen = np.array([datadict[k]["times"].shape[0] for k in self.cidx])
            self.wfx = np.zeros(self.wfxlen.sum())
            startidx = 0
            for ik, k in enumerate(self.cidx):
                self.wfx[startidx : (startidx + self.wfxlen[ik])] = datadict[k]["times"]
                startidx += self.wfxlen[ik]

        if self.evt.dat is None:
            startidx = 0
            self.evt.dat = np.zeros(self.wfxlen.sum())
            for ik, k in enumerate(self.cidx):
                self.evt.dat[startidx : (startidx + self.wfxlen[ik])] = datadict[k][0]
                startidx += self.wfxlen[ik]
            try:
                self.evt.dat = np.array(self.evt.dat)
            except:
                print(
                    "HsdObject: could not cast waveform times to array ", self.evt.dat
                )


class Wave8Object(WaveformObject):
    def __init__(self, det, run, **kwargs):
        super(Wave8Object, self).__init__(det, run, **kwargs)
        self._chan_names = [name for name in dir(det.raw) if name[0] != "_"]
        if "raw_all" in self._chan_names:
            self._chan_names = ["raw_all"]
        self.chan_names = " ".join(self._chan_names)

    def getData(self, evt):
        super(Wave8Object, self).getData(evt)
        vetolist = ["config"]
        self.evt.dat = [
            getattr(self.det.raw, name)(evt)
            for name in self._chan_names
            if name not in vetolist
        ]
        try:
            self.evt.dat = np.squeeze(np.array(self.evt.dat)).astype(int)
        except:
            logger.debug(f"Wave8: Could not cast waveform times to array, set to None ")
            self.evt.dat = None


class PiranhaObject(WaveformObject):
    def __init__(self, det, run, **kwargs):
        super(PiranhaObject, self).__init__(det, run, **kwargs)
        return

    def getData(self, evt):
        super(PiranhaObject, self).getData(evt)
        try:
            self.evt.dat = self.det.raw.raw(evt).astype(int)
        except:
            print("piranha: Could not cast piranha to array, set to None")
            self.evt.dat = None
        return


class ImageObject(CameraObject):
    def __init__(self, det, run, **kwargs):
        """
        Extension of CameraObject for detector that should use the psana image method instead of the
        usual calib method.
        For now (1/22/2025), this is only applicable to the Archon detector.
        """
        super().__init__(det, run, **kwargs)
        self.params_to_img()
        return

    def params_to_img(self):
        """
        Convert all parameters from the calib shape to the image shape:
        """
        # pedestals
        self.ped = self._to_img_shape(self.ped)

        # rms
        self.rms = self._to_img_shape(self.rms)

        # gain
        self.gain = self._to_img_shape(self.gain)

        # masks
        self.mask = self._to_img_shape(self.mask)
        self.cmask = self._to_img_shape(self.cmask)
        return

    def _to_img_shape(self, arr):
        if arr is None:
            return None
        else:
            return self.det.raw.image(None, arr)

    def getData(self, evt):
        super().getData(evt)
        if self.common_mode < 0:
            dat = self.det.raw.raw(evt)
        elif self.common_mode in [0, 30]:
            dat = self.det.raw.calib(evt)
        if dat is None:
            return
        self.evt.dat = self.det.raw.image(None, dat)
        return


class ArchonObject(ImageObject):
    def __init__(self, det, run, **kwargs):
        super(ArchonObject, self).__init__(det, run, **kwargs)
        self.common_mode = kwargs.get("common_mode", 30)  # default to calib
        if self.common_mode is None:
            self.common_mode = 30
        return
