"""
droplet2photons_gpu.py — drop-in for smalldata_tools' droplet2Photons.

Subclasses DetObjectFunc and matches its process()/self.dat contract, but performs the
photon assignment with droplets_gpu (vectorized NumPy on CPU, or CuPy on GPU, including
the giant-blob path the reference lacks). Output self.dat = {tile,row,col,data} is photon-
count-identical to the reference.

Usage (exactly like droplet2Photons in a smalldata pipeline):
    from droplet2photons_gpu import droplet2photons_gpu
    f = droplet2photons_gpu(aduspphot=172, name='droplet2phot')      # CPU (default)
    f = droplet2photons_gpu(aduspphot=172, use_gpu=True)             # GPU path
It consumes the same data dict {_image, _imgDrop, _mask} produced upstream by dropletFunc.
"""

import warnings

import numpy as np
import scipy.ndimage as ndi

try:
    import cupy as cp

    _HAS_CUPY = True
except ImportError:
    cp = None
    _HAS_CUPY = False

from smalldata_tools.common.detector_base import DetObjectFunc

try:  # vendored inside smalldata_tools/ana_funcs/
    from smalldata_tools.ana_funcs.droplets_gpu import find_photons, _labeled_sum_gpu
except ImportError:  # or droplets_gpu on PYTHONPATH (dev)
    from droplets_gpu import find_photons, _labeled_sum_gpu


class droplet2photons_gpu(DetObjectFunc):
    def __init__(self, **kwargs):
        self._name = kwargs.get("name", "droplet2phot")
        super(droplet2photons_gpu, self).__init__(**kwargs)
        # aduspphot / offset / photpts are resolved ONCE, here, and the photon-finding
        # call receives explicit inputs -- the three supported cases (issue #289):
        #   1. aduspphot only            -> offset = aduspphot/2, uniform thresholds
        #   2. aduspphot + offset        -> uniform thresholds n*aduspphot - offset
        #   3. photpts (a list of edges) -> used VERBATIM (non-uniform thresholds, e.g. to
        #      cut a low-energy background); offset is not meaningful and is ignored with a
        #      warning if also given; aduspphot is still needed for multi-photon splitting,
        #      so it is derived from median(diff(photpts)) unless supplied (a supplied value
        #      that disagrees with the derived one gets a warning, since that combination is
        #      legitimate but easy to do by accident).
        self.aduspphot = kwargs.get("aduspphot", 0)
        self.photpts = kwargs.get("photpts", None)
        self._user_photpts = self.photpts is not None
        self._user_offset = "offset" in kwargs
        if self._user_photpts:
            self.photpts = np.asarray(self.photpts, dtype=np.float64)
            derived = (
                float(np.median(np.diff(self.photpts)))
                if self.photpts.size >= 2
                else 0.0
            )
            if self.aduspphot <= 0.0:
                self.aduspphot = derived
            elif derived > 0 and abs(self.aduspphot - derived) > 0.2 * derived:
                warnings.warn(
                    "droplet2photons_gpu: aduspphot=%g vs median(diff(photpts))=%g -- using "
                    "your aduspphot for photon splitting and your photpts for the "
                    "thresholds (fine if intentional, e.g. a background cutoff)"
                    % (self.aduspphot, derived),
                    RuntimeWarning,
                )
            if self._user_offset:
                warnings.warn(
                    "droplet2photons_gpu: offset= is ignored when photpts= is given",
                    RuntimeWarning,
                )
            self.offset = self.aduspphot * 0.5
        else:
            self.offset = kwargs.get("offset", self.aduspphot * 0.5)
        self.use_gpu = bool(kwargs.get("use_gpu", False)) and _HAS_CUPY
        self._photpts_dev = None  # device copy of a user photpts, uploaded once
        return

    def process(self, data):
        if (not isinstance(data, dict)) or (data.get("_imgDrop", None) is None):
            print("droplet2photons_gpu expects a dict with _image and _imgDrop keys!")
            return {}

        img = data["_image"]
        imgDrop = data["_imgDrop"]

        # build the droplet dict our find_photons expects, from the pre-labeled image
        if self.use_gpu:
            img = cp.asarray(img, dtype=cp.float32)
            img_drop = cp.asarray(imgDrop).astype(cp.int32)
            n_drop = int(img_drop.max()) if img_drop.size else 0
            drop_ind = cp.arange(1, n_drop + 1, dtype=cp.int32)
            adu_drop = (
                _labeled_sum_gpu(img, img_drop, n_drop)
                if n_drop
                else cp.zeros(0, dtype=cp.float32)
            )
        else:
            img = np.asarray(img, dtype=np.float32)
            img_drop = np.asarray(imgDrop).astype(np.int32)
            n_drop = int(img_drop.max()) if img_drop.size else 0
            drop_ind = np.arange(1, n_drop + 1, dtype=np.int32)
            adu_drop = (
                ndi.sum_labels(img, labels=img_drop, index=drop_ind).astype(np.float32)
                if n_drop
                else np.zeros(0, dtype=np.float32)
            )

        if n_drop == 0:
            self.dat = {k: np.array([]) for k in ("tile", "row", "col", "data")}
            return self._collect()

        droplet_dict = {
            "img": img,
            "img_drop": img_drop,
            "drop_ind": drop_ind,
            "adu_drop": adu_drop,
        }
        # photon_pts, resolved per the cases above (any argument passed IS used):
        #   user photpts -> ship the user's list verbatim (typically short; uploaded once);
        #   user offset  -> uniform edges n*aduspphot - offset, sized to THIS event's
        #                   droplets (a handful of floats, not the reference's 1e6 array);
        #   neither      -> photon_pts=None: find_photons' built-in default is the same
        #                   classification as n*aduspphot - aduspphot/2, with no array
        #                   shipped to the GPU at all.
        if self._user_photpts:
            if self.use_gpu:
                if self._photpts_dev is None:
                    self._photpts_dev = cp.asarray(self.photpts)
                photon_pts = self._photpts_dev
            else:
                photon_pts = self.photpts
        elif self._user_offset:
            xp = cp if self.use_gpu else np
            n_max = int(float(adu_drop.max()) / self.aduspphot) + 3
            photon_pts = (
                xp.arange(n_max, dtype=xp.float64) * self.aduspphot - self.offset
            )
        else:
            photon_pts = None
        photons = find_photons(
            droplet_dict, float(self.aduspphot), photon_pts=photon_pts
        )
        if self.use_gpu:
            photons = photons.get()

        n = int(photons.shape[0])
        self.dat = {
            "tile": np.zeros(n),
            "row": photons[:, 0],
            "col": photons[:, 1],
            "data": np.ones(n),
        }
        return self._collect()

    def _collect(self):
        ret_dict = {}
        sub = self.processFuncs()
        for k in sub:
            for kk in sub[k]:
                ret_dict["%s_%s" % (k, kk)] = sub[k][kk]
        return ret_dict
