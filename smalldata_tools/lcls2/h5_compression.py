"""Optional lossless compression for the smalldata HDF5 files.

Why: psana2's smalldata writer (psana/smalldata.py, ``Server.new_dset``) creates every
dataset with ``h5py.Group.create_dataset(..., chunks=...)`` and has no compression
option. Integer and boolean per-shot data (counters, codes, flags, raw waveforms) usually
compress well with gzip level 1 and the shuffle filter, at little CPU cost. This module is
an opt-in stopgap until lcls2 supports compression natively.

``enable()`` wraps ``h5py.Group.create_dataset`` (``h5py.File`` inherits it) in the
current process. A call gets ``compression="gzip", compression_opts=1`` and
``shuffle=True`` (unless the caller set ``shuffle``) only when all of these hold:

- the dataset is chunked (``chunks`` is given; ``Server.new_dset`` always passes it);
- its dtype (``dtype=``, or the dtype of ``data``) is a plain signed/unsigned integer or
  bool; float, complex, string, object, compound and sub-array dtypes are left alone;
- the caller gave no ``compression``, ``compression_opts``, ``scaleoffset`` or
  ``external``.

Every other call is passed through unchanged. dtype, shape, maxshape and chunk shape are
kept, so the stored values are bit-identical and h5py readers decompress transparently.
Virtual datasets are not affected: ``Group.create_virtual_dataset`` builds them through
``VirtualLayout.make_dataset``, not through ``create_dataset``.

``enable()`` is idempotent; ``disable()`` restores the original method.
"""

import logging

import h5py
import numpy as np

logger = logging.getLogger(__name__)

GZIP_LEVEL = 1
INT_KINDS = "iub"  # signed integer, unsigned integer, bool

_orig_create_dataset = None


def _resolve_dtype(dtype, data):
    try:
        if dtype is not None:
            return np.dtype(dtype)
        if data is not None:
            dt = getattr(data, "dtype", None)
            return dt if dt is not None else np.asarray(data).dtype
    except (TypeError, ValueError):
        return None
    return None


def wants_compression(dtype, data, kwds):
    """Return True if a create_dataset call with these arguments gets gzip-1 + shuffle."""
    if kwds.get("chunks") in (None, False):
        return False
    for k in ("compression", "compression_opts", "scaleoffset", "external"):
        if kwds.get(k) is not None:
            return False
    dt = _resolve_dtype(dtype, data)
    if dt is None or dt.fields is not None or dt.subdtype is not None:
        return False
    return dt.kind in INT_KINDS


def _create_dataset(self, name, shape=None, dtype=None, data=None, **kwds):
    if wants_compression(dtype, data, kwds):
        kwds["compression"] = "gzip"
        kwds["compression_opts"] = GZIP_LEVEL
        kwds.setdefault("shuffle", True)
    return _orig_create_dataset(self, name, shape=shape, dtype=dtype, data=data, **kwds)


def enable():
    """Install the create_dataset wrapper (idempotent). Return True if it is active."""
    global _orig_create_dataset
    if _orig_create_dataset is None:
        _orig_create_dataset = h5py.Group.create_dataset
        h5py.Group.create_dataset = _create_dataset
        logger.debug(
            "h5_compression enabled: gzip-%d+shuffle for chunked integer/bool datasets "
            "(h5py %s, HDF5 %s)",
            GZIP_LEVEL,
            h5py.__version__,
            h5py.version.hdf5_version,
        )
    return is_enabled()


def disable():
    """Restore the original h5py.Group.create_dataset."""
    global _orig_create_dataset
    if _orig_create_dataset is not None:
        h5py.Group.create_dataset = _orig_create_dataset
        _orig_create_dataset = None


def is_enabled():
    return h5py.Group.create_dataset is _create_dataset
