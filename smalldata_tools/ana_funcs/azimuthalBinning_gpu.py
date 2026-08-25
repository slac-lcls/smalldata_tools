"""
azimuthalBinning_gpu.py — drop-in for smalldata_tools' azimuthalBinning.

Subclasses azimuthalBinning and overrides ONLY the per-frame reduction. All of the geometry,
the q/phi/r binning, the polarization and solid-angle corrections, the masking and the
thresholds are the parent's — unchanged, and still the single source of truth.

The parent's hot path is one line::

    I = np.bincount(self.Cake_idxs, weights=img / self.correction, minlength=nradial * nphi)

A bincount over a fixed index array IS a sparse matrix-vector product: build ``M`` once with
``M[Cake_idxs[i], i] = 1`` and ``M @ img`` is the same sum. Folding ``1/correction`` into the
matrix values makes the per-frame divide disappear as well, so a frame costs one matvec —
cuSPARSE ``csrmv`` on the GPU, ``scipy.sparse`` on the CPU. ``M`` depends only on geometry,
so it is built once and reused for every event.

Usage (exactly like azimuthalBinning in a smalldata pipeline)::

    from smalldata_tools.ana_funcs.azimuthalBinning_gpu import azimuthalBinning_gpu
    f = azimuthalBinning_gpu(center=..., dis_to_sam=..., eBeam=...)                 # CPU (default)
    f = azimuthalBinning_gpu(center=..., dis_to_sam=..., eBeam=..., use_gpu=True)   # GPU path

``use_gpu=True`` requires CuPy; without it the class falls back to the CPU path and warns once,
so a pipeline written for the GPU still runs on a CPU-only node.

WHEN THIS IS WORTH IT
    The win is on large frames whose data is already on the GPU. Measured on one A100 against
    this module's parent, on a 3000x3000 (9 Mpix) frame: azimuthalBinning 32.6 ms on one CPU
    core, this 0.34 ms GPU-resident (~95x), or ~2.7 ms including a fresh host->device copy of
    the frame (~12x). On a per-panel frame (~0.19 Mpix) a single matvec is dominated by launch
    and occupancy overhead rather than by pixels — there the SpMM batch path (``doCake_batch``)
    is what pays, and a per-event loop is not. On the CPU, numpy's ``bincount`` is roughly 2x
    faster than a scipy CSR matvec, which is why the CPU default here simply defers to the
    parent: this class is a GPU accelerator, not a CPU rewrite.

NUMERICS
    ``M @ img`` and ``np.bincount`` compute the same sum in a different order, so results agree
    to floating-point tolerance rather than bit-exactly (typically <1e-12 relative on float64
    frames). ``tests/test_azimuthalBinning_gpu.py`` asserts that against the parent on random
    and structured frames, including masked pixels and multiple phi bins.

Origin: the "sparse matrix as accumarray" formulation is the same one used by the LCLS DRP
radial-integration primitive (slac-lcls/drp-benchmarks, ``radial_integration/radial.py``),
which is in turn the GPU form of the pattern this module's parent already uses on the CPU.
"""

import warnings

import numpy as np
import scipy.sparse

try:
    import cupy as cp
    import cupyx.scipy.sparse as cpsparse

    _HAS_CUPY = True
except ImportError:  # CPU-only node: the GPU path degrades to the parent
    cp = None
    cpsparse = None
    _HAS_CUPY = False

from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning


class azimuthalBinning_gpu(azimuthalBinning):
    """azimuthalBinning with the per-frame bincount replaced by a sparse matvec.

    Extra keyword
    -------------
    use_gpu : bool (default False)
        Run the per-frame reduction on the GPU via CuPy/cuSPARSE. Falls back to the parent's
        CPU path (with a one-time warning) if CuPy is unavailable.
    """

    def __init__(self, **kwargs):
        self.use_gpu = bool(kwargs.pop("use_gpu", False))
        self._name = kwargs.get("name", "azav")
        super(azimuthalBinning_gpu, self).__init__(**kwargs)
        # cache: applyCorrection -> operator. Built lazily, after the parent's _setup() has
        # produced Cake_idxs/correction (which happens in setFromDet/setFromFunc, not __init__).
        self._M = {}
        self._warned_no_cupy = False

    # ------------------------------------------------------------------ operator
    def _nradial(self):
        """Rows per phi bin, matching the parent's own choice in doCake()."""
        return self.nr if self.rbin is not None else self.nq

    def _on_gpu(self):
        if not self.use_gpu:
            return False
        if not _HAS_CUPY:
            if not self._warned_no_cupy:
                warnings.warn(
                    "azimuthalBinning_gpu: use_gpu=True but CuPy is not available; "
                    "falling back to the CPU bincount path.",
                    RuntimeWarning,
                )
                self._warned_no_cupy = True
            return False
        return True

    def _operator(self, applyCorrection):
        """CSR M with M[bin, pixel] = 1 (or 1/correction[pixel]), built once per geometry.

        Rows = nradial * nphi, columns = the unmasked pixels, i.e. exactly the length of the
        vector the parent's doCake() feeds to bincount.
        """
        gpu = self._on_gpu()
        key = (bool(applyCorrection), gpu)
        if key in self._M:
            return self._M[key]

        idx = np.asarray(self.Cake_idxs).ravel()
        npix = idx.size
        # The parent uses minlength=nradial*nphi for the corrected branch and nq*nphi for the
        # uncorrected one; reproduce that exactly rather than quietly normalising the two.
        nrows = (self._nradial() if applyCorrection else self.nq) * self.nphi
        nrows = max(nrows, int(idx.max()) + 1)

        if applyCorrection:
            vals = 1.0 / np.asarray(self.correction).ravel().astype(np.float64)
        else:
            vals = np.ones(npix, dtype=np.float64)

        M = scipy.sparse.csr_matrix(
            (vals, (idx, np.arange(npix))), shape=(nrows, npix), dtype=np.float64
        )
        if gpu:
            M = cpsparse.csr_matrix(M)
        self._M[key] = M
        return M

    # ------------------------------------------------------------------ per frame
    def doCake(self, img, applyCorrection=True):
        """One frame -> Icake, shape (nphi, nradial). Same contract as the parent."""
        if not self._on_gpu():
            return super(azimuthalBinning_gpu, self).doCake(
                img, applyCorrection=applyCorrection
            )

        # Mirror the parent's pre-processing exactly, in the same order.
        if self.darkImg is not None:
            img = img - self.darkImg
        if self.gainImg is not None:
            img = img / self.gainImg

        good = np.asarray(self._mask).ravel() == 0
        flat = (
            img.ravel()[good]
            if not isinstance(img, cp.ndarray)
            else img.ravel()[cp.asarray(good)]
        )
        x = cp.asarray(flat, dtype=cp.float64)

        M = self._operator(applyCorrection)
        I = M @ x

        nradial = self._nradial() if applyCorrection else self.nq
        I = I[: nradial * self.nphi].reshape(self.nphi, nradial)
        self.Icake = I / cp.asarray(self.Cake_norm)
        return self.Icake

    def doCake_batch(self, imgs, applyCorrection=True):
        """(B, ...) frames -> (B, nphi, nradial) in one SpMM.

        A single small frame does not saturate a GPU: the matvec is dominated by launch and
        per-call cuSPARSE setup rather than by pixels, so on per-panel-sized data the batch
        path is the one that pays (the crossover sits around B~32 for a ~0.19 Mpix panel).
        Frames must already be dark/gain corrected the same way doCake() would.
        """
        if not self._on_gpu():
            return np.stack(
                [
                    super(azimuthalBinning_gpu, self).doCake(
                        np.asarray(im), applyCorrection=applyCorrection
                    )
                    for im in imgs
                ]
            )

        good = np.asarray(self._mask).ravel() == 0
        X = cp.asarray(imgs, dtype=cp.float64).reshape(len(imgs), -1)[
            :, cp.asarray(good)
        ]
        M = self._operator(applyCorrection)
        # csrmm prefers an F-contiguous dense RHS; X.T of a C-contiguous X is exactly that.
        I = (M @ X.T).T

        nradial = self._nradial() if applyCorrection else self.nq
        I = I[:, : nradial * self.nphi].reshape(len(imgs), self.nphi, nradial)
        return I / cp.asarray(self.Cake_norm)
