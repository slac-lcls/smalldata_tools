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
import numpy as np
import scipy.ndimage as ndi

try:
    import cupy as cp
    _HAS_CUPY = True
except ImportError:
    cp = None
    _HAS_CUPY = False

from smalldata_tools.common.detector_base import DetObjectFunc
try:                                          # vendored inside smalldata_tools/ana_funcs/
    from smalldata_tools.ana_funcs.droplets_gpu import find_photons, _labeled_sum_gpu
except ImportError:                           # or droplets_gpu on PYTHONPATH (dev)
    from droplets_gpu import find_photons, _labeled_sum_gpu


class droplet2photons_gpu(DetObjectFunc):
    def __init__(self, **kwargs):
        self._name = kwargs.get("name", "droplet2phot")
        super(droplet2photons_gpu, self).__init__(**kwargs)
        self.aduspphot = kwargs.get("aduspphot", 0)
        self.offset = self.aduspphot * 0.5
        # photpts[n] = n*aduspphot - offset  (identical convention to droplet2Photons)
        photpts = np.arange(1000000) * self.aduspphot - self.offset
        self.photpts = kwargs.get("photpts", photpts)
        self.use_gpu = bool(kwargs.get("use_gpu", False)) and _HAS_CUPY
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
            adu_drop = (_labeled_sum_gpu(img, img_drop, n_drop) if n_drop
                        else cp.zeros(0, dtype=cp.float32))
        else:
            img = np.asarray(img, dtype=np.float32)
            img_drop = np.asarray(imgDrop).astype(np.int32)
            n_drop = int(img_drop.max()) if img_drop.size else 0
            drop_ind = np.arange(1, n_drop + 1, dtype=np.int32)
            adu_drop = (ndi.sum_labels(img, labels=img_drop, index=drop_ind).astype(np.float32)
                        if n_drop else np.zeros(0, dtype=np.float32))

        if n_drop == 0:
            self.dat = {k: np.array([]) for k in ("tile", "row", "col", "data")}
            return self._collect()

        droplet_dict = {"img": img, "img_drop": img_drop,
                        "drop_ind": drop_ind, "adu_drop": adu_drop}
        # photon_pts=None → find_photons sizes the (uniform) edges to the data; identical
        # classification to self.photpts, without shipping a 1e6 array to the GPU.
        photons = find_photons(droplet_dict, float(self.aduspphot), photon_pts=None)
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
