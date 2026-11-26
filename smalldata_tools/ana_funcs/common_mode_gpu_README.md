# `common_mode_gpu` — common-mode subtraction on GPU or CPU

Subtracts a per-tile baseline from a calibrated detector frame: for each tile,
the median of that tile's pixels below `cut`. One vectorized pass, CuPy or
NumPy, output type matching the input, exactly like `droplets_gpu`.

```python
from smalldata_tools.ana_funcs.common_mode_gpu import common_mode_asic

frame = det.raw.calib(evt)                              # (2, 336, 576) float32
corr  = common_mode_asic(frame, cut=2.0)                # CPU
corr  = common_mode_asic(cp.asarray(frame), cut=2.0)    # GPU, same answer
```

## Why this is here and not in psana

psana's ePixUHR `calib` does **not** apply a common-mode correction, and does
not say so. `epixuhr3x2.py::calib()` delegates to
`UtilsEpixUHR.calib_v02()`, which documents `cmpars` in its docstring and
stores it on the `Storage` object, but has both the parameter read and the
application commented out:

```python
#cmpars  = kwa.get('cmpars', None)
...
#if store.cmpars is not None:
#    common_mode_epix_multigain_apply(arrf, gmaps, store)
```

`common_mode_epix_multigain_apply` exists and is in production use for
epix10ka (`UtilsEpix10ka.py:854`), and `epix100.py::calib` takes and applies
`cmpars` normally. So `det.raw.calib(evt, cmpars=...)` on an ePixUHR is a
**silent no-op**: no warning, no error, no correction. That is why UHR users do
the correction by hand downstream, and why this module exists.

Uncommenting it is not a one-liner: `gmaps` is not constructed in
`calib_v02` (it computes the gain-bit index `igr` instead), so the gain maps
for the UHR's two gain ranges have to be built first. The proper fix belongs
upstream in psana; this module is the near-term path, and is also useful on its
own for anyone who wants the correction on a GPU or wants the baselines back as
a diagnostic.

## The family is one operation

Every common-mode variant in the LCLS detector code is the same arithmetic on a
different tiling:

| tiling | `tile` | name elsewhere |
|---|---|---|
| per-ASIC / region | `(168, 192)` for UHR 3x2 | `regionCommonModeCorrection` |
| per-row, per bank | `(1, cols_per_bank)` | `rowCommonModeCorrection` |
| per-column, per bank | `(rows_per_bank, 1)` | `colCommonModeCorrection` |
| one global offset | whole plane | `regionCommonModeCorrection`, region = frame |

`common_mode_tiles(frame, tile, cut)` is the operation; `common_mode_asic`,
`common_mode_rows`, `common_mode_cols` and `common_mode_global` are one-line
wrappers that pick a tile. Leading axes (panel, event, ...) are independent, so
a whole stack can go through in one call.

## Which tiling is right for the ePixUHR 3x2 is not settled

This matters more than the speed, so it is stated plainly:

- The correction in use at XPP in November 2026 is **per-ASIC**, `cut=2.0` on
  calib ADU — the 2x3 ASIC grid inside each `(336, 576)` panel.
- The beamtime calibration suite's own `detectorInfo.setup_epixUHR3x2()` sets
  `nBanksRow = 1` and `nBanksCol = 1` with the comments *"common mode seems
  global"*, and sets `preferredCommonMode = "clusterCommonMode"  ## not yet
  defined`. So for this detector the suite has no settled choice, and its row
  and column variants would run with one bank — per-row across the full
  576-pixel width.
- The per-ASIC choice is finer than anything the suite currently configures for
  the 3x2, and is well motivated: the ASIC grid is where independent readout
  offsets would actually arise.

The `cut` defaults differ for the same reason: `1000.0` on the row/column/global
wrappers follows the suite, which runs on **raw ADU**; `2.0` on
`common_mode_asic` follows the XPP usage on **calib** output at roughly
1 ADU/keV. Check which your frames are in before trusting either default.

## The baselines are the drift signal

`return_cm=True` also returns the subtracted medians, shaped
`(..., R//th, C//tw)`:

```python
corr, cm = common_mode_asic(stack, cut=2.0, return_cm=True)
# stack (N, 2, 336, 576) -> cm (N, 2, 2, 3): one baseline per ASIC per event
```

Plotting `cm` against event number is a direct look at the per-ASIC baseline
drift — the thing the correction is there to remove — at no extra cost.

## Measured

One A100 and one core of the same node (`sdfampere027`), 60 calib frames of `xpp102335627` run 45,
`epixuhr1M_0`, both panels, per-ASIC, `cut=2.0`:

| | ms/frame | vs the reference loop |
|---|---|---|
| reference per-ASIC median loop, CPU | 6.40 | 1.00x |
| `common_mode_asic`, CPU | 6.84 | 0.94x |
| `common_mode_asic`, GPU, one frame per call | 0.92 | 7.0x |
| same, including the 1.55 MB H2D copy | 0.94 | 6.8x |
| same, all 60 frames in one call | 0.29 | 22x |

**The CPU path is not a speedup** — it is slightly slower than the loop it replaces, because 12
medians of 32k pixels is already NumPy-bound and the blocking copy costs more than it saves. It is
here so that the same call works in the psana environment, where CuPy is not available, and so the
two paths are one implementation. The GPU is the speedup, and batching the leading axis is worth
another 3x on top.

## Provenance

The tiling family is Philip Hart's, from `calibrationSuite/basicSuiteScript.py` in
`beamtime-calibration-suite` — `regionCommonModeCorrection`, `rowCommonModeCorrection`,
`colCommonModeCorrection`. No code was copied: the implementation here is a single tiled-median
pass, and the suite's functions are per-row and per-region Python loops. The suite is the reference
the tests check against, and the names in the table above are the suite's names so the
correspondence is traceable.

The per-ASIC application is Yanwen Sun's, written for this beamtime.

## Correctness

`tests/test_common_mode_gpu.py` (no psana, no data files) checks that

- per-ASIC tiling is **bit-identical** to the per-ASIC median loop in use at XPP;
- per-row-per-bank tiling is **bit-identical** to
  `calibrationSuite.basicSuiteScript.rowCommonModeCorrection`, including its
  `if not np.isnan(rowCM)` guard — a tile with nothing below `cut` is left
  alone rather than being filled with NaN;
- a global tile is one median per plane, leading axes independent;
- integer input is promoted to float32, and a tile that does not divide the
  frame raises rather than silently truncating;
- the GPU path is bit-identical to the CPU path for all four tilings.

The same checks were run on real data — 60 calib frames of `xpp102335627`
run 45, `epixuhr1M_0` — with the per-ASIC arm bit-identical to the reference
loop on every frame. Numbers are in the module's verification record.

Equality is bit-for-bit, not `allclose`, because `np.median` and
`cp.nanmedian` average the two middle values on an even count and both do the
same thing on the same multiset. That equality is what makes the GPU path a
drop-in; **if the correction changes, re-run the test before trusting the GPU
path again.** The check is one line:
`assert np.array_equal(reference(frame), common_mode_asic(frame, cut=cut))`.
