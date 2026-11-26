"""
Self-contained equivalence tests for common_mode_gpu (no psana / no data files).

  - per-ASIC tiling is bit-identical to the per-ASIC median loop in use at XPP
  - per-row-per-bank tiling is bit-identical to
    calibrationSuite.basicSuiteScript.rowCommonModeCorrection
  - a global tile is one median per plane, with leading axes independent
  - a tile with no pixels below the cut is left alone (the suite's isnan guard)
  - the GPU path matches the CPU path bit for bit (skipped if CuPy/GPU absent)

Run:  pytest tests/test_common_mode_gpu.py -v
"""

import numpy as np
import pytest

from smalldata_tools.ana_funcs.common_mode_gpu import (
    UHR3X2_ASIC,
    common_mode_asic,
    common_mode_cols,
    common_mode_global,
    common_mode_rows,
    common_mode_tiles,
)

try:
    import cupy

    _HAS_GPU = cupy.cuda.runtime.getDeviceCount() > 0
except Exception:
    _HAS_GPU = False

CUT = 2.0
ADU_PER_PHOTON = 990.0  # ePixUHR at XPP, 9.9 keV, Yanwen's x100 scaling


def _frame(seed=0, shape=(2, 336, 576), nphot=4000):
    """A calib-like frame: noise baseline, a per-ASIC offset, sparse photons."""
    rng = np.random.default_rng(seed)
    f = rng.normal(0.0, 0.8, shape).astype(np.float32)
    ar, ac = UHR3X2_ASIC
    for p in range(shape[0]):  # a different baseline on every ASIC
        for i in range(shape[1] // ar):
            for j in range(shape[2] // ac):
                f[p, i * ar : (i + 1) * ar, j * ac : (j + 1) * ac] += rng.normal(0, 3)
    idx = rng.integers(0, shape[1] * shape[2], size=nphot)
    for p in range(shape[0]):
        f[p].reshape(-1)[idx] += ADU_PER_PHOTON
    return f


def _ref_asic(frame, cut=CUT, asic=UHR3X2_ASIC):
    """The per-ASIC median loop in use at XPP, November 2026."""
    ar, ac = asic
    rows = [slice(i * ar, (i + 1) * ar) for i in range(frame.shape[-2] // ar)]
    cols = [slice(j * ac, (j + 1) * ac) for j in range(frame.shape[-1] // ac)]
    out = frame.copy()
    for p in range(out.shape[0]):
        for rs in rows:
            for cs in cols:
                sub = out[p, rs, cs]
                sel = sub[sub < cut]
                if sel.size:
                    sub -= np.median(sel)
    return out


def _ref_rows(frame, cut=1000.0, cols_per_bank=288):
    """calibrationSuite.basicSuiteScript.rowCommonModeCorrection, per panel."""
    out = frame.copy()
    for p in range(out.shape[0]):
        for r in range(out.shape[1]):
            off = 0
            while off < out.shape[2]:
                seg = out[p, r, off : off + cols_per_bank]
                sel = seg[seg < cut]
                cm = np.median(sel) if sel.size else np.nan
                if not np.isnan(cm):
                    out[p, r, off : off + cols_per_bank] -= cm
                off += cols_per_bank
    return out


def test_asic_matches_reference_loop():
    f = _frame()
    assert np.array_equal(_ref_asic(f), common_mode_asic(f, cut=CUT))


def test_rows_matches_suite_loop():
    f = _frame(seed=1)
    got = common_mode_rows(f, cols_per_bank=288, cut=1000.0)
    assert np.array_equal(_ref_rows(f, cut=1000.0, cols_per_bank=288), got)


def test_cm_shape_and_values():
    f = _frame(seed=2)
    corr, cm = common_mode_asic(f, cut=CUT, return_cm=True)
    assert cm.shape == (2, 336 // 168, 576 // 192)
    # subtracting the reported baseline is exactly what was done
    assert np.array_equal(corr, common_mode_tiles(f, UHR3X2_ASIC, CUT))


def test_global_is_one_median_per_plane():
    f = _frame(seed=3)
    got = common_mode_global(f, cut=1000.0)
    for p in range(f.shape[0]):
        want = f[p] - np.median(f[p][f[p] < 1000.0])
        assert np.array_equal(got[p], want)


def test_tile_with_nothing_below_cut_is_untouched():
    z = np.full((1, 4, 4), 5.0, dtype=np.float32)
    assert np.array_equal(common_mode_tiles(z, (2, 2), cut=1.0), z)


def test_non_dividing_tile_raises():
    with pytest.raises(ValueError):
        common_mode_tiles(_frame(), (100, 100), CUT)


def test_integer_input_is_promoted():
    out = common_mode_tiles(np.arange(16).reshape(4, 4), (2, 2), 10.0)
    assert out.dtype == np.float32


def test_cols_tiling_runs_and_is_a_transpose_of_rows():
    f = _frame(seed=4)
    by_col = common_mode_cols(f, rows_per_bank=168, cut=1000.0)
    by_row_of_T = np.stack(
        [
            common_mode_rows(f[p].T[None], cols_per_bank=168, cut=1000.0)[0].T
            for p in range(f.shape[0])
        ]
    )
    assert np.array_equal(by_col, by_row_of_T)


@pytest.mark.skipif(not _HAS_GPU, reason="no CuPy/GPU")
@pytest.mark.parametrize(
    "fn,kw",
    [
        (common_mode_asic, {"cut": CUT}),
        (common_mode_rows, {"cols_per_bank": 288, "cut": 1000.0}),
        (common_mode_cols, {"rows_per_bank": 168, "cut": 1000.0}),
        (common_mode_global, {"cut": 1000.0}),
    ],
)
def test_gpu_matches_cpu(fn, kw):
    f = _frame(seed=5)
    cpu = fn(f, **kw)
    gpu = cupy.asnumpy(fn(cupy.asarray(f), **kw))
    assert np.array_equal(cpu, gpu)


@pytest.mark.skipif(not _HAS_GPU, reason="no CuPy/GPU")
def test_gpu_returns_cupy():
    out = common_mode_asic(cupy.asarray(_frame(seed=6)), cut=CUT)
    assert isinstance(out, cupy.ndarray)
