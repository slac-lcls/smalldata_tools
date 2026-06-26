# `droplet2photons_gpu` — a GPU/CPU drop-in for `droplet2Photons`

Handoff note for the smalldata_tools maintainers. This adds a contract-identical, much
faster implementation of the droplet → photon step, with an optional GPU path.

## Why

`droplet2Photons` does the multi-photon split with `loopdrops`/`greedyguess`, a per-droplet
**Python** loop. It's fine for sparse frames but collapses on dense ones: on a real XPP
`epix_alc5` frame (the direct-beam blob is one ~160 k-photon droplet) it takes **~2.8 s/frame**
(worse on the brightest shots). This drop-in (a) adds a dedicated path for saturated/giant
droplets, and (b) **compiles the greedyguess per-droplet loop** (numba-jitted driver, no
Python-level loop) — so the same frame is **~18 ms on CPU / ~5–8 ms on a GPU** with **identical
photon counts**. Compiling the driver alone is ~4× over the interpreted per-droplet loop
(single-threaded); an optional `parallel=True` adds numba `prange` threading, but that's only
worth it on dense frames *and* when the producer isn't already parallel across events, so the
serial compiled path is the default.

## What it adds (2 files, both in `smalldata_tools/ana_funcs/`)

| file | role |
|------|------|
| `droplet2photons_gpu.py` | `DetObjectFunc` subclass — the drop-in. Same `process(data)` / `self.dat` contract as `droplet2Photons`. |
| `droplets_gpu.py` | core library: `dropletize` + `find_photons`, NumPy on CPU **or** CuPy on GPU (auto-dispatch by array type). Vendored here like `dropletCode/`. |

CuPy is an **optional** import — if it isn't present the module still imports and the CPU path
works unchanged (`use_gpu=True` silently falls back). So the CPU drop-in is safe to land on the
existing psana CPU/MPI farm with no new hard dependency.

## How to use (same as `droplet2Photons`)

```python
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.droplet2photons_gpu import droplet2photons_gpu

drop = dropletFunc(threshold=..., thresholdLow=..., thresADU=..., useRms=True)
drop.addFunc(droplet2photons_gpu(aduspphot=172, name='photons'))   # CPU
# drop.addFunc(droplet2photons_gpu(aduspphot=172, use_gpu=True))   # GPU path
```

It consumes the upstream `_imgDrop` from `dropletFunc` and writes
`self.dat = {tile, row, col, data}` — same keys/format as `droplet2Photons`, so anything
downstream (HDF5 writer, `getProb`/XPCS sub-funcs) is unaffected.

## How to test

A self-contained pytest (synthetic frames — no psana, no data files):

```bash
pytest tests/test_droplet2photons_gpu.py -v
```

It checks, against the reference `droplet2Photons`:
- identical photon counts and `self.dat` contract, both sparse and with a saturated giant blob;
- identical photon positions when nothing is saturated;
- the GPU path matches the CPU path (auto-skipped if CuPy/GPU is absent).

(If the repo's `tests/conftest.py` fails to import in your environment, add `--noconftest` —
this test uses no fixtures.)

Optional live-data integration check:

```bash
python tests/example_pipeline_psana1.py [exp] [run] [detector]   # defaults: xpplw3920 130 epix_alc4
```

runs a psana1 `DataSource` through `dropletFunc` → `droplet2photons_gpu` and reports
droplets/event and photons/event (e.g. ~33 droplets → ~19 photons/event on the default run).

## Benchmarks (real XPP `xpplw3920` r130, 704×768, identical photon counts across all rows; end-to-end, incl. host↔device transfers)

| method | sparse speckle (`epix_alc4`) | dense · beam **masked** (`epix_alc5`) | dense · beam **unmasked** (`epix_alc5`) |
|---|---|---|---|
| `droplet2Photons` (reference, jitted) | 24 ms | 60 ms | ~2.8 s |
| this drop-in, **CPU** | 5.0 ms | 15 ms | 18 ms |
| this drop-in, **GPU** (A100 / H100) | 3.3 / 2.3 ms | 7.2 / 4.3 ms | 8.0 / 4.8 ms |

**Sparse speckle is the operating point** — that's where XPCS actually runs. The beam blob is
pile-up (≥2 photons/pixel, where photon counting is undefined) and is **masked** out of the
analysis, so the *unmasked* column is the worst case: the reference stalls at seconds/frame
(a DAQ hazard on a saturated shot), while the drop-in stays in the tens of ms regardless — the
giant-blob path is robustness insurance, not the headline. The everyday CPU speedup (~4–5×) is
the faster reductions, not the numba JIT (the reference's master jits its driver too); the GPU
adds the placement kernels.

Note: on small frames the GPU is launch/transfer-bound and can be **slower than the CPU**
(~2.4 vs ~1.1 ms/frame on a 256² frame) — it wins on the full detector and dense frames. Pick the
path per workload (`use_gpu=True` is opt-in). On dense frames an opt-in `parallel=True` shaves a
further ~1.4× on CPU.

## Known differences (by design)

- **Counts are identical** to `droplet2Photons` in every test.
- **Positions differ only inside saturated/giant droplets** (> 64 photons, e.g. the beam blob).
  There the drop-in uses integer placement + a centroid for the residual instead of running
  full greedyguess (which is the bulk of the ~2.8 s cost and where sub-pixel position is
  meaningless anyway).
  Outside such blobs, positions match 100% (sparse) / pixel-exact (dense). Mask the beam
  (`max_drop_adu=`) and the difference vanishes entirely.
- Multi-photon positions are sub-pixel floats (from greedyguess); round if integer pixels are
  required downstream.

## Notes for maintainers

- A genuine bug surfaced while validating: `droplet2Photons`'s default `photpts` is fine, but a
  fixed-length count-edge array silently caps a droplet at its last entry. This drop-in sizes the
  edges to the data (and uses float64 — float32 loses integer precision above ~1.7e7 ADU, which a
  saturated blob exceeds). Worth a look in the reference too.
- The greedyguess driver (`_greedyguess_driver_serial`) is numba-jitted and writes into
  pre-offset output slices, so it has no Python-level per-droplet overhead. The same routine
  is the reference's `greedyguess`, so the math is unchanged — only the loop is compiled. The
  threaded twin (`_greedyguess_driver_par`, `prange`) is byte-identical and opt-in.
- Vendoring `droplets_gpu.py` mirrors how `dropletCode/` is vendored. Alternatively make it a small
  installable package and import it.
- Registration: add `droplet2photons_gpu` wherever the available funcs are listed for producer
  configs (same place `droplet2Photons` is referenced).
