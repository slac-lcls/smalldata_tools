# azimuthalBinning_gpu

A drop-in for `azimuthalBinning` that moves the per-frame reduction to the GPU. It subclasses
`azimuthalBinning` and overrides **only** `doCake()`; the geometry, the q/phi/r binning, the
polarization and solid-angle corrections, the masking and the thresholds all remain the
parent's, unchanged.

## Why this is a small change

The parent's hot path is a single line:

```python
I = np.bincount(self.Cake_idxs, weights=img / self.correction, minlength=nradial * nphi)
```

A `bincount` over a fixed index array *is* a sparse matrix-vector product. Build `M` once with
`M[Cake_idxs[i], i] = 1` and `M @ img` is the same sum; folding `1/correction` into the matrix
values removes the per-frame divide too. `M` depends only on geometry, so it is built once and
reused for every event. On the GPU that matvec is cuSPARSE `csrmv`.

This is the same formulation as the LCLS DRP radial-integration primitive
(`slac-lcls/drp-benchmarks`, `radial_integration/radial.py`), which is itself the GPU form of the
pattern `azimuthalBinning` already uses on the CPU.

## Usage

```python
from smalldata_tools.ana_funcs.azimuthalBinning_gpu import azimuthalBinning_gpu

f = azimuthalBinning_gpu(center=..., dis_to_sam=..., eBeam=...)                 # CPU (default)
f = azimuthalBinning_gpu(center=..., dis_to_sam=..., eBeam=..., use_gpu=True)   # GPU path
```

Every other keyword is the parent's. `use_gpu=True` needs CuPy; without it the class warns once
and falls back to the parent's CPU path, so a pipeline written for the GPU still runs on a
CPU-only node.

For many frames at once:

```python
Icake = f.doCake_batch(imgs)     # (B, ...) -> (B, nphi, nradial) in one SpMM
```

## When it is worth using — and when it is not

Measured on one A100 against the parent, 3000x3000 (9 Mpix) frame:

| | time | vs parent |
| --- | --- | --- |
| `azimuthalBinning`, 1 CPU core | 32.6 ms | — |
| this, GPU-resident frame | 0.34 ms | ~95x |
| this, including a fresh host->device copy | ~2.7 ms | ~12x |

Two honest caveats on those numbers:

- **The frame should already be on the GPU.** A host->device copy of a 9 Mpix frame is ~2.4 ms,
  about seven times the compute. It is still a net win against a 32.6 ms CPU bincount, but the
  ~95x figure is only available to a pipeline that already has the data on the device.
- **Small frames do not saturate a GPU.** On a per-panel-sized frame (~0.19 Mpix) a single matvec
  is dominated by kernel-launch and per-call cuSPARSE setup rather than by pixels, so a per-event
  loop is a poor fit; the crossover for the batch path sits around B~32. Use `doCake_batch`, or
  stay on the CPU.

On the CPU, numpy's `bincount` is roughly 2x *faster* than a scipy CSR matvec, which is why the
CPU default here simply defers to the parent. This class is a GPU accelerator, not a CPU rewrite,
and `azimuthalBinning` remains the right CPU/MPI tool.

## Numerics

`M @ img` and `np.bincount` compute the same sum in a different order, so agreement is to
floating-point tolerance rather than bit-exact. Measured on an A100, against `bincount`:

```
scipy CSR   vs bincount : max rel 3.9e-16
cupyx CSR   vs bincount : max rel 4.9e-16
cupyx SpMM  vs bincount : max rel 5.2e-16
```

## Testing

`tests/test_azimuthalBinning_gpu.py` — 9 passing, no psana and no data files needed. It asserts
the CPU path is *identical* to the parent (not merely close), that the sparse operator reproduces
the parent's `bincount` to 1e-12 with masked pixels and multiple phi bins, that masked pixels
cannot contribute however bright, and that the no-CuPy fallback still returns the parent's answer.

The GPU-marked tests skip unless CuPy and a device are present. They could not be run in either
S3DF environment — `ana-4.0.68-py3` has smalldata_tools' `skimage` dependency but no CuPy, and
`xpp_drp_gpu_311` has CuPy but not `skimage` — so the GPU arithmetic was validated separately by
the standalone equivalence check whose numbers are quoted above.

## Running the GPU path in practice

There is currently no single S3DF environment with both smalldata_tools' dependency chain and
CuPy: `ana-4.0.68-py3` has `skimage`/psana but no CuPy, and `xpp_drp_gpu_311` has CuPy but not
`skimage`. So `use_gpu=True` only helps where CuPy happens to be importable alongside
smalldata_tools.

The proven way around this in a real producer is a shared-memory env bridge, prototyped for GPU
droplet finding (`py2py/shm_bridge/` on `slac-lcls/drp-benchmarks`, branch `add-droplet-gpu`): a
`DetObjectFunc` keeps running in the psana env and offloads the GPU work to a child process in
the CuPy env over a `/dev/shm` zero-copy ring, so the psana env never needs CuPy. (The related
`envbridge` package on that repo's `main` solves cross-env calls by copying arrays through a
socket; it does not provide the shared-memory frame path.) That pattern applies unchanged here —
the operator `M` is built once and only the matvec would cross — and is the natural route for
MPI production. This PR does not add it; it matches the direct-CuPy shape of
`droplet2photons_gpu`, which is what was upstreamed previously.

## One pre-existing bug this surfaced (not fixed here)

`azimuthalBinning.doCake()` selects `nradial = self.nr` when `rbin` is set and bincounts into
`nr * nphi` bins, then reshapes unconditionally to `(nphi, self.nq)`
(`azimuthalBinning.py:385`). Whenever `nr != nq` that raises `ValueError`, so the `rbin` path
does not currently work. `test_parent_rbin_path_is_broken_upstream` documents it.

This is the single place the GPU path deliberately differs from the parent: it reshapes to
`(nphi, nradial)`, i.e. it does the consistent thing and does not crash. The one-line fix to the
parent is deliberately left out of this PR to keep it focused — happy to fold it in.
