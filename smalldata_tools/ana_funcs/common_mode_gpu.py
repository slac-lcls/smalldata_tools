"""
common_mode_gpu.py
Common-mode (baseline) subtraction on GPU or CPU, for detectors whose readout
leaves a per-ASIC / per-bank / per-row offset on top of the calibrated frame.

Standalone functions -- no DetObjectFunc dependency.

GPU path : arrays stay on GPU as CuPy arrays.
CPU path : arrays stay as NumPy arrays.
The dispatcher picks the path from the input array type, exactly like
droplets_gpu: pass a NumPy array to force CPU even when a GPU is present.

Why this module exists
----------------------
psana's ePixUHR calib does NOT apply a common-mode correction: the `cmpars`
argument is accepted and stored by `UtilsEpixUHR.calib_v02` but the call to
`common_mode_epix_multigain_apply` is commented out, so
`det.raw.calib(evt, cmpars=...)` is a silent no-op for that detector.  Users
therefore do the correction by hand after calib, and on a (2, 336, 576) frame a
per-ASIC median loop costs several ms per frame -- comparable to the read
itself.  These functions do the same arithmetic in one vectorized pass, which
the GPU path then runs ~4x faster than the loop.

The one algorithm, and the family it covers
-------------------------------------------
Every member of the common-mode family in LCLS detector code is the same
operation on a different tiling of the frame:

    partition the frame into tiles; from each tile subtract the median of
    its own pixels below `cut`.

`cut` keeps photons out of the baseline estimate -- set it well below one
photon.  So:

    tile = (asic_rows, asic_cols)   per-ASIC / "region" common mode
    tile = (1, cols_per_bank)       per-row-per-bank ("rowCommonMode")
    tile = (rows_per_bank, 1)       per-column-per-bank ("colCommonMode")
    tile = whole frame              one global offset

`common_mode_tiles` is that operation; the named wrappers below just pick a
tile.  Which tiling is right for a given detector is a calibration question,
not a software one -- see `common_mode_gpu_README.md`.

Typical usage
-------------
    import cupy as cp
    from smalldata_tools.ana_funcs.common_mode_gpu import common_mode_asic

    frame = det.raw.calib(evt)                 # (2, 336, 576) float32, NumPy
    corr  = common_mode_asic(frame, cut=2.0)   # CPU, vectorized

    corr  = common_mode_asic(cp.asarray(frame), cut=2.0)   # GPU, same answer

    # the per-tile offsets are themselves the drift signal:
    corr, cm = common_mode_asic(frame, cut=2.0, return_cm=True)
    # cm.shape == (2, 2, 3): one baseline per ASIC
"""

import warnings

import numpy as np

try:
    import cupy as cp

    _HAS_CUPY = True
except ImportError:
    cp = None
    _HAS_CUPY = False


# ePixUHR 3x2: each (336, 576) panel is a 2x3 grid of 168x192 ASICs.
UHR3X2_ASIC = (168, 192)


def _is_gpu(arr):
    return _HAS_CUPY and isinstance(arr, cp.ndarray)


def common_mode_tiles(frame, tile, cut, return_cm=False):
    """
    Subtract, from every tile of `frame`, the median of that tile's pixels
    below `cut`.  Works on CuPy or NumPy input; the output type matches.

    Parameters
    ----------
    frame     : array, shape (..., R, C).  Leading axes (panels, events, ...)
                are independent; the last two are the detector plane.
                Integer input is promoted to float32.
    tile      : (th, tw).  Tile height and width in pixels; both must divide
                R and C.  Use 0 or None for "the whole extent of that axis",
                so (0, 0) is a single global median and (1, 576) is per-row.
    cut       : float.  Only pixels strictly below `cut` enter the median.
                Units are whatever `frame` is in -- for a calib ePixUHR frame
                at ~1 ADU/keV, cut=2.0 is ~0.2 photon at 9.9 keV.
    return_cm : bool.  Also return the per-tile medians.

    Returns
    -------
    corrected                    if return_cm is False
    (corrected, cm)              if return_cm is True
        corrected : same shape and type as `frame`
        cm        : shape (..., R//th, C//tw), the subtracted baselines

    Notes
    -----
    A tile with no pixels below `cut` gets no correction (baseline 0), which
    is what `calibrationSuite.basicSuiteScript.rowCommonModeCorrection` does
    via its `if not np.isnan(rowCM)` guard.
    """
    xp = cp if _is_gpu(frame) else np

    a = xp.asarray(frame)
    if a.ndim < 2:
        raise ValueError("frame must have at least 2 dimensions, got %d" % a.ndim)
    if not xp.issubdtype(a.dtype, xp.floating):
        a = a.astype(xp.float32)

    lead = a.shape[:-2]
    nrow, ncol = a.shape[-2:]
    th, tw = tile
    th = nrow if not th else int(th)
    tw = ncol if not tw else int(tw)
    if nrow % th or ncol % tw:
        raise ValueError(
            "tile %r does not divide the frame plane (%d, %d)" % ((th, tw), nrow, ncol)
        )
    ntr, ntc = nrow // th, ncol // tw

    # (..., ntr, th, ntc, tw) -> (..., ntr, ntc, th*tw): one axis per tile
    blocked = a.reshape(*lead, ntr, th, ntc, tw)
    flat = xp.moveaxis(blocked, -2, -3).reshape(*lead, ntr, ntc, th * tw)

    sel = xp.where(flat < cut, flat, xp.asarray(float("nan"), dtype=flat.dtype))
    with warnings.catch_warnings():
        # an all-above-cut tile is a legitimate state, not a problem
        warnings.simplefilter("ignore", RuntimeWarning)
        cm = xp.nanmedian(sel, axis=-1)  # (..., ntr, ntc)
    cm = xp.nan_to_num(cm, nan=0.0)

    # broadcast back over the (ntr, th, ntc, tw) layout, then restore the plane
    corrected = (blocked - cm[..., :, None, :, None]).reshape(*lead, nrow, ncol)
    return (corrected, cm) if return_cm else corrected


def common_mode_asic(frame, asic=UHR3X2_ASIC, cut=2.0, return_cm=False):
    """
    Per-ASIC common mode: one median per ASIC, subtracted from that ASIC.

    This is the correction in use on the ePixUHR 3x2 at XPP (November 2026):
    `asic=(168, 192)` splits each (336, 576) panel into its 2x3 ASIC grid.
    It is the `regionCommonModeCorrection` of the beamtime calibration suite
    applied to the ASIC grid rather than to a single ROI.

    See common_mode_tiles for the parameters; `asic` is its `tile`.
    """
    return common_mode_tiles(frame, asic, cut, return_cm=return_cm)


def common_mode_rows(frame, cols_per_bank=None, cut=1000.0, return_cm=False):
    """
    Per-row-per-bank common mode ("rowCommonMode"): one median per row segment.

    `cols_per_bank=None` means the full row width (one bank per row).  This is
    the algorithm of `calibrationSuite.basicSuiteScript.rowCommonModeCorrection`,
    whose default `arbitraryCut` is 1000 because it runs on raw ADU, not calib.
    """
    return common_mode_tiles(frame, (1, cols_per_bank), cut, return_cm=return_cm)


def common_mode_cols(frame, rows_per_bank=None, cut=1000.0, return_cm=False):
    """
    Per-column-per-bank common mode ("colCommonMode"): one median per column
    segment.  `rows_per_bank=None` means the full column height.
    """
    return common_mode_tiles(frame, (rows_per_bank, 1), cut, return_cm=return_cm)


def common_mode_global(frame, cut=1000.0, return_cm=False):
    """
    One median for the whole detector plane, subtracted from all of it.
    This is `regionCommonModeCorrection` with the region being the full frame.
    """
    return common_mode_tiles(frame, (None, None), cut, return_cm=return_cm)
