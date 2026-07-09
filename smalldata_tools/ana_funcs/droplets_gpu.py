"""
droplets_gpu.py
GPU-accelerated droplet finding and photon extraction using CuPy,
with automatic CPU fallback via scipy/numpy when CuPy is unavailable.

Standalone functions — no DetObjectFunc dependency.

GPU path  : arrays stay on GPU as CuPy arrays.
CPU path  : arrays stay as NumPy arrays.
The dispatcher (dropletize / find_photons) picks the path based on
whether CuPy is importable AND the input array is a CuPy array.
Pass a NumPy array to force CPU even when a GPU is present.

Typical GPU usage
-----------------
    import cupy as cp
    from droplets_gpu import dropletize, find_photons

    img_gpu  = cp.asarray(raw_image)
    mask_gpu = cp.asarray(mask)
    comp_gpu = cp.asarray(rms / gain)

    drops   = dropletize(img_gpu, mask=mask_gpu, comp_data=comp_gpu,
                         threshold=5.0, thres_adu=0.5 * adu_per_photon)
    photons = find_photons(drops, adu_per_photon=adu_per_photon)
    # photons: cp.ndarray [N, 2]  →  (row, col) float32

Typical CPU usage (or when CuPy is unavailable)
------------------------------------------------
    import numpy as np
    from droplets_gpu import dropletize, find_photons

    drops   = dropletize(raw_image, mask=mask, comp_data=rms / gain,
                         threshold=5.0, thres_adu=0.5 * adu_per_photon)
    photons = find_photons(drops, adu_per_photon=adu_per_photon)
    # photons: np.ndarray [N, 2]  →  (row, col) float32
"""

import numpy as np
import scipy.ndimage as sndi

try:
    import cupy as cp
    import cupyx.scipy.ndimage as cundi

    _HAS_CUPY = True
except ImportError:
    cp = None
    cundi = None
    _HAS_CUPY = False

try:
    from numba import njit as _njit, prange as _prange

    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False

    def _njit(**kwargs):  # passthrough when numba unavailable
        return lambda fn: fn

    _prange = range  # plain range fallback


def _is_gpu(arr):
    return _HAS_CUPY and isinstance(arr, cp.ndarray)


# Cross-shaped 4-connectivity kernel (same as the original footprint)
_FOOTPRINT_NP = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.uint8)
_FOOTPRINT_GPU = (
    None  # created lazily on first GPU call to avoid import-time device init
)


def _get_footprint_gpu():
    global _FOOTPRINT_GPU
    if _FOOTPRINT_GPU is None:
        _FOOTPRINT_GPU = cp.asarray(_FOOTPRINT_NP)
    return _FOOTPRINT_GPU


# ── greedyguess CUDA kernel ────────────────────────────────────────────────────
# One thread per multi-photon droplet. Each thread:
#   1. Unpacks its pixels from a CSR-packed input (rows, cols, adus).
#   2. Builds a local 2D sub-image (with 1-pixel zero border).
#   3. Places photons from integer-part ADU first, then fractionally from residuals.
#      The fractional position splits the photon between the peak pixel and its
#      brightest 4-connected neighbor — matching Sutton/Lurio greedyguess.py.
_GREEDYGUESS_KERNEL_SRC = r"""
#define MAX_SUB_DIM 36
#define MAX_SUB     (MAX_SUB_DIM * MAX_SUB_DIM)
#define MAX_PHOT    64

extern "C" __global__ void greedyguess_kernel(
    const int*   d_rows,     // packed pixel rows   (global coords)
    const int*   d_cols,     // packed pixel cols
    const float* d_adus,     // packed pixel ADUs
    const int*   d_pstart,   // CSR start index per droplet
    const int*   d_npix,     // pixel count per droplet
    const int*   d_nphot,    // photon count per droplet
    const int*   d_ostart,   // output start index per droplet
    float*       out_rows,   // output photon rows (sub-pixel)
    float*       out_cols,   // output photon cols (sub-pixel)
    float        aduspphot,
    int          n_drops
) {
    int did = blockIdx.x * blockDim.x + threadIdx.x;
    if (did >= n_drops) return;

    int pstart = d_pstart[did];
    int npix   = d_npix[did];
    int nphot  = d_nphot[did];
    int ostart = d_ostart[did];

    if (nphot <= 0 || npix == 0) return;
    if (nphot > MAX_PHOT) nphot = MAX_PHOT;

    // Bounding box of droplet pixels
    int rmin = d_rows[pstart], rmax = rmin;
    int cmin = d_cols[pstart], cmax = cmin;
    for (int k = 1; k < npix; k++) {
        int r = d_rows[pstart + k], c = d_cols[pstart + k];
        if (r < rmin) rmin = r; if (r > rmax) rmax = r;
        if (c < cmin) cmin = c; if (c > cmax) cmax = c;
    }

    // Sub-image dimensions: bounding box + 1-pixel zero border on each side
    int snrows = rmax - rmin + 3;
    int sncols = cmax - cmin + 3;

    // Oversized droplet: fall back to center-of-mass for all photon positions
    if (snrows > MAX_SUB_DIM || sncols > MAX_SUB_DIM) {
        float sw = 0.f, sr = 0.f, sc = 0.f;
        for (int k = 0; k < npix; k++) {
            float a = d_adus[pstart + k];
            sw += a;
            sr += a * d_rows[pstart + k];
            sc += a * d_cols[pstart + k];
        }
        float cr = (sw > 0.f) ? sr / sw : (float)d_rows[pstart];
        float cc = (sw > 0.f) ? sc / sw : (float)d_cols[pstart];
        for (int n = 0; n < nphot; n++) {
            out_rows[ostart + n] = cr;
            out_cols[ostart + n] = cc;
        }
        return;
    }

    // Local fractional sub-image (zero-initialized)
    float timg[MAX_SUB];
    for (int k = 0; k < snrows * sncols; k++) timg[k] = 0.f;

    float pxs[MAX_PHOT], pys[MAX_PHOT];
    int nassigned = 0;

    // Fill sub-image; place photons for any pixel carrying >= 1 photon worth of ADU
    for (int k = 0; k < npix; k++) {
        int   ii  = d_rows[pstart + k] - rmin + 1;
        int   jj  = d_cols[pstart + k] - cmin + 1;
        float val = d_adus[pstart + k] / aduspphot;
        int   ip  = (int)val;                 // integer part
        timg[ii * sncols + jj] = val - (float)ip;  // fractional residual
        for (int p = 0; p < ip && nassigned < nphot; p++) {
            pxs[nassigned] = (float)ii;
            pys[nassigned] = (float)jj;
            nassigned++;
        }
    }

    // Greedy fractional assignment for remaining photons
    for (int n = nassigned; n < nphot; n++) {
        // Find max fractional pixel
        int   bi = 0, bj = 0;
        float bv = -1.f;
        for (int ii = 0; ii < snrows; ii++) {
            for (int jj = 0; jj < sncols; jj++) {
                float v = timg[ii * sncols + jj];
                if (v > bv) { bv = v; bi = ii; bj = jj; }
            }
        }

        int i = bi, j = bj;
        float t = 1.f - bv;       // fraction of photon NOT at central pixel
        timg[i * sncols + j] = 0.f;

        // Brightest 4-connected neighbor
        float vup    = (i > 0)        ? timg[(i-1)*sncols + j] : -1.f;
        float vdown  = (i < snrows-1) ? timg[(i+1)*sncols + j] : -1.f;
        float vleft  = (j > 0)        ? timg[i*sncols + (j-1)] : -1.f;
        float vright = (j < sncols-1) ? timg[i*sncols + (j+1)] : -1.f;

        if (vup >= vdown && vup >= vleft && vup >= vright) {
            pxs[n] = i - t; pys[n] = (float)j;
            if (i > 0) timg[(i-1)*sncols + j] = fmaxf(0.f, timg[(i-1)*sncols + j] - t);
        } else if (vdown >= vleft && vdown >= vright) {
            pxs[n] = i + t; pys[n] = (float)j;
            if (i < snrows-1) timg[(i+1)*sncols + j] = fmaxf(0.f, timg[(i+1)*sncols + j] - t);
        } else if (vleft >= vright) {
            pxs[n] = (float)i; pys[n] = j - t;
            if (j > 0) timg[i*sncols + (j-1)] = fmaxf(0.f, timg[i*sncols + (j-1)] - t);
        } else {
            pxs[n] = (float)i; pys[n] = j + t;
            if (j < sncols-1) timg[i*sncols + (j+1)] = fmaxf(0.f, timg[i*sncols + (j+1)] - t);
        }
    }

    // Write sub-image coords back to global coords
    for (int n = 0; n < nphot; n++) {
        out_rows[ostart + n] = pxs[n] - 1.f + (float)rmin;
        out_cols[ostart + n] = pys[n] - 1.f + (float)cmin;
    }
}
"""

_greedyguess_kernel = None  # compiled lazily on first GPU call


def _get_greedyguess_kernel():
    global _greedyguess_kernel
    if _greedyguess_kernel is None:
        _greedyguess_kernel = cp.RawKernel(
            _GREEDYGUESS_KERNEL_SRC, "greedyguess_kernel"
        )
    return _greedyguess_kernel


# ── image preparation ──────────────────────────────────────────────────────────


def prepare_image(img, mask, comp_data, threshold, pixmax=None):
    """
    Apply mask, optional pixel-value clip, and threshold.

    Parameters
    ----------
    img       : cp.ndarray, 2D detector frame on GPU.
    mask      : cp.ndarray or None, uint8/bool — 1 = valid, 0 = bad pixel.
    comp_data : cp.ndarray, per-pixel threshold scale (e.g. rms/gain).
    threshold : float, threshold in units of comp_data.
    pixmax    : float or None, zero pixels above this value before thresholding.

    Returns
    -------
    cp.ndarray — copy of img with sub-threshold and masked pixels set to 0.
    """
    out = img.copy()
    if mask is not None:
        out *= mask.astype(out.dtype)
    if pixmax is not None:
        out[out > pixmax] = 0.0
    out[out < comp_data * threshold] = 0.0
    return out


# ── labeled reductions via scatter-add (bincount) ───────────────────────────────
# NOTE: cupyx.scipy.ndimage.sum(labels=, index=) is very slow — it dispatches to a
# generic sort-based labeled reduction. Profiling showed it dominated GPU runtime
# (>70% at low occupancy). cp.bincount(weights=) is a direct scatter-add and is
# 1–2 orders of magnitude faster, so all labeled sums go through it.


def _labeled_sum_gpu(img_thr, img_drop, n_drop):
    """
    Total intensity per label (1..n_drop) over the thresholded image.
    Only bright (labeled) pixels are scattered, avoiding bin-0 contention.

    Returns
    -------
    cp.ndarray of length n_drop, dtype float32 — sum for labels 1..n_drop.
    """
    flat_lab = img_drop.ravel()
    nz = flat_lab > 0
    labs = flat_lab[nz]
    vals = img_thr.ravel()[nz].astype(cp.float64)
    sums = cp.bincount(labs, weights=vals, minlength=n_drop + 1)
    return sums[1 : n_drop + 1].astype(cp.float32)


def _center_of_mass_batch(img, img_drop, drop_ind):
    """
    Intensity-weighted center of mass for a batch of labeled regions.

    Sparse scatter-add: extract only labeled pixels, map their label to a group
    index, then three cp.bincount reductions (total, row·w, col·w). Avoids the
    slow cundi.sum path entirely.

    Returns
    -------
    cp.ndarray of shape [N, 2], dtype float32 — (row, col) per droplet.
    """
    n = int(drop_ind.size)
    if n == 0:
        return cp.zeros((0, 2), dtype=cp.float32)

    r, c = cp.where(img_drop > 0)
    if r.size == 0:
        return cp.zeros((n, 2), dtype=cp.float32)

    labels_all = img_drop[r, c]
    adus_all = img[r, c].astype(cp.float64)

    # map label value → group index in drop_ind (-1 if that droplet was filtered out)
    max_label = int(img_drop.max())
    lut = cp.full(max_label + 1, -1, dtype=cp.int32)
    lut[drop_ind] = cp.arange(n, dtype=cp.int32)
    group = lut[labels_all]

    valid = group >= 0
    g = group[valid]
    w = adus_all[valid]
    rv = r[valid].astype(cp.float64)
    cv = c[valid].astype(cp.float64)

    total = cp.bincount(g, weights=w, minlength=n)
    row_m = cp.bincount(g, weights=w * rv, minlength=n)
    col_m = cp.bincount(g, weights=w * cv, minlength=n)

    safe_total = cp.where(total > 0, total, 1.0)
    return cp.stack([row_m / safe_total, col_m / safe_total], axis=1).astype(cp.float32)


# ── droplet finding (GPU) ─────────────────────────────────────────────────────


def _dropletize_gpu(
    img,
    mask=None,
    comp_data=None,
    threshold=5.0,
    threshold_low=None,
    thres_adu=None,
    relabel=True,
    pixmax=None,
    max_drop_adu=None,
):
    """GPU implementation of dropletize(). Input/output are CuPy arrays."""
    if threshold_low is None:
        threshold_low = threshold
    if comp_data is None:
        comp_data = cp.ones_like(img)

    img_thr = prepare_image(img, mask, comp_data, threshold, pixmax)
    # _FOOTPRINT_NP passed as numpy so label() doesn't have to convert internally
    img_drop, n_drop = cundi.label(img_thr, structure=_FOOTPRINT_NP)

    if threshold != threshold_low:
        # Expand each droplet to neighboring pixels that pass the lower threshold
        neighbor = cundi.maximum_filter(img_drop, footprint=_get_footprint_gpu())
        img_low = prepare_image(img, mask, comp_data, threshold_low, pixmax)
        if relabel:
            neighbor[img_low == 0] = 0
            img_drop, n_drop = cundi.label(neighbor, structure=_FOOTPRINT_NP)
        else:
            img_drop = neighbor

    n_all = int(n_drop)
    if n_all == 0:
        return {
            "nDroplets_all": 0,
            "nDroplets": 0,
            "img": img_thr,
            "img_drop": img_drop,
            "drop_ind": cp.array([], dtype=cp.int32),
            "adu_drop": cp.array([], dtype=cp.float32),
        }

    drop_ind = cp.arange(1, n_drop + 1, dtype=cp.int32)
    adu_drop = _labeled_sum_gpu(img_thr, img_drop, n_drop)

    if thres_adu is not None or max_drop_adu is not None:
        keep = cp.ones(adu_drop.shape, dtype=cp.bool_)
        if thres_adu is not None:
            keep &= adu_drop >= thres_adu
        if max_drop_adu is not None:  # veto saturated/bright blobs
            keep &= adu_drop <= max_drop_adu
        # Zero out vetoed droplets via a label lookup table — O(n_pixels), no Python loop
        lut = cp.zeros(n_drop + 1, dtype=img_drop.dtype)
        kept_labels = drop_ind[keep]
        lut[kept_labels] = kept_labels
        img_drop = lut[img_drop]
        drop_ind = kept_labels
        adu_drop = adu_drop[keep]

    return {
        "nDroplets_all": n_all,
        "nDroplets": int(drop_ind.size),
        "img": img_thr,
        "img_drop": img_drop,
        "drop_ind": drop_ind,
        "adu_drop": adu_drop,
    }


# ── photon extraction (GPU) ───────────────────────────────────────────────────


def _photon_count_per_drop(adu_drop, photon_pts):
    """
    Return integer photon count for each element of adu_drop.

    photon_pts[n] is the lower ADU boundary for n photons, so:
        n_photons = searchsorted(photon_pts, adu, side='right') - 1
    """
    return (cp.searchsorted(photon_pts, adu_drop, side="right") - 1).astype(cp.int32)


def _multi_photon_positions(
    img, img_drop, drop_ind, n_photons_per_drop, adu_per_photon
):
    """
    Assign sub-pixel photon positions for multi-photon droplets via greedyguess.

    One CUDA thread per droplet. Pixels are CSR-packed (sorted by label),
    then the kernel builds a local sub-image and runs the greedy fractional
    assignment algorithm (Sutton/Lurio) matching the reference greedyguess.py.

    Returns
    -------
    cp.ndarray of shape [total_photons, 2], dtype float32 — (row, col).
    """
    n_drops = int(drop_ind.size)
    if n_drops == 0:
        return cp.zeros((0, 2), dtype=cp.float32)

    # --- pixels belonging to multi-photon droplets ---
    max_label = int(img_drop.max())
    is_multi_lut = cp.zeros(max_label + 1, dtype=cp.bool_)
    is_multi_lut[drop_ind] = True
    rows, cols = cp.where(is_multi_lut[img_drop])
    if rows.size == 0:
        return cp.zeros((0, 2), dtype=cp.float32)

    labels = img_drop[rows, cols]
    adus = img[rows, cols].astype(cp.float32)

    # --- CSR packing: sort pixels by label so each droplet is contiguous ---
    # Stability not required: the kernel rebuilds each droplet's sub-image from
    # all its pixels regardless of intra-group order, so a plain sort suffices.
    sort_idx = cp.argsort(labels)
    labels_s = labels[sort_idx].astype(cp.int32)
    rows_s = rows[sort_idx].astype(cp.int32)
    cols_s = cols[sort_idx].astype(cp.int32)
    adus_s = adus[sort_idx]

    # start/end index in the packed arrays for each droplet (searchsorted on sorted labels)
    di32 = drop_ind.astype(cp.int32)
    pstart = cp.searchsorted(labels_s, di32, side="left").astype(cp.int32)
    pend = cp.searchsorted(labels_s, di32, side="right").astype(cp.int32)
    npix_arr = (pend - pstart).astype(cp.int32)
    nphot_arr = n_photons_per_drop.astype(cp.int32)

    # --- output offsets (cumulative sum of per-droplet photon counts) ---
    out_start = cp.zeros(n_drops + 1, dtype=cp.int32)
    cp.cumsum(nphot_arr, out=out_start[1:])
    total_phot = int(out_start[-1])
    if total_phot == 0:
        return cp.zeros((0, 2), dtype=cp.float32)

    out_rows = cp.empty(total_phot, dtype=cp.float32)
    out_cols = cp.empty(total_phot, dtype=cp.float32)

    # --- launch: one thread per droplet ---
    kernel = _get_greedyguess_kernel()
    threads = 128
    blocks = (n_drops + threads - 1) // threads
    kernel(
        (blocks,),
        (threads,),
        (
            rows_s,
            cols_s,
            adus_s,
            pstart,
            npix_arr,
            nphot_arr,
            out_start[:n_drops],
            out_rows,
            out_cols,
            np.float32(adu_per_photon),
            np.int32(n_drops),
        ),
    )

    return cp.stack([out_rows, out_cols], axis=1)


# Droplets with more photons than the greedyguess kernel can hold (MAX_PHOT) are
# routed to the vectorized giant-droplet path instead.
GIANT_NPHOT = 64


def _run_length_index(counts):
    """
    GPU-safe equivalent of the index that np.repeat(values, counts) would gather:
    returns src such that values[src] == repeat(values, counts).
    (cp.repeat does not accept an array of repeats, so we build the index ourselves.)
    """
    ends = cp.cumsum(counts)
    total = int(ends[-1]) if counts.size else 0
    if total == 0:
        return cp.zeros(0, dtype=cp.int64)
    return cp.searchsorted(ends, cp.arange(total), side="right")


def _giant_droplet_positions(img, img_drop, giant_ind, giant_nphot, adu_per_photon):
    """
    Photon placement for droplets too large for the per-thread greedyguess kernel
    (saturated/beam blobs with up to ~1e5 photons). Fully vectorized, no per-droplet
    loop:
      • integer photons  → floor(adu_px / adu_per_photon) placed AT each pixel
      • per-droplet residual (nphot − Σ integer) → placed at the droplet centroid
    Photon COUNT is exact (matches the digitize count); positions inside the blob are
    approximate, which is appropriate for saturated regions where sub-pixel assignment
    is not meaningful.
    """
    n = int(giant_ind.size)
    if n == 0:
        return cp.zeros((0, 2), dtype=cp.float32)

    max_label = int(img_drop.max())
    is_g = cp.zeros(max_label + 1, dtype=cp.bool_)
    is_g[giant_ind] = True
    r, c = cp.where(is_g[img_drop])

    labels = img_drop[r, c]
    adus = img[r, c].astype(cp.float64)
    rf = r.astype(cp.float64)
    cf = c.astype(cp.float64)

    lut = cp.full(max_label + 1, -1, dtype=cp.int32)
    lut[giant_ind] = cp.arange(n, dtype=cp.int32)
    grp = lut[labels]

    # integer photons per pixel, expanded via run-length index
    k = cp.maximum(cp.floor(adus / adu_per_photon).astype(cp.int64), 0)
    idx_int = _run_length_index(k)
    int_rows = r.astype(cp.float32)[idx_int]
    int_cols = c.astype(cp.float32)[idx_int]

    # per-droplet integer total + intensity-weighted centroid
    sum_k = cp.bincount(grp, weights=k.astype(cp.float64), minlength=n)
    total = cp.bincount(grp, weights=adus, minlength=n)
    row_m = cp.bincount(grp, weights=adus * rf, minlength=n)
    col_m = cp.bincount(grp, weights=adus * cf, minlength=n)
    safe = cp.where(total > 0, total, 1.0)
    com_r = (row_m / safe).astype(cp.float32)
    com_c = (col_m / safe).astype(cp.float32)

    # residual photons placed at the centroid (exact count restored)
    resid = cp.maximum(giant_nphot.astype(cp.int64) - sum_k.astype(cp.int64), 0)
    idx_res = _run_length_index(resid)
    res_rows = com_r[idx_res]
    res_cols = com_c[idx_res]

    return cp.stack(
        [cp.concatenate([int_rows, res_rows]), cp.concatenate([int_cols, res_cols])],
        axis=1,
    )


def _find_photons_gpu(droplet_dict, adu_per_photon, photon_pts=None):
    """GPU implementation of find_photons(). Input/output are CuPy arrays."""
    img = droplet_dict["img"]
    img_drop = droplet_dict["img_drop"]
    drop_ind = droplet_dict["drop_ind"]
    adu_drop = droplet_dict["adu_drop"]

    if drop_ind.size == 0:
        return cp.zeros((0, 2), dtype=cp.float32)

    if photon_pts is None:
        # Size the count edges to the data (float64 — float32 loses integer precision
        # above ~1.7e7 ADU, which a saturated blob can exceed). photon_pts[k]=k·APH−APH/2.
        offset = adu_per_photon * 0.5
        amax = float(adu_drop.max()) if adu_drop.size else 0.0
        n_max = int(amax / adu_per_photon) + 3
        photon_pts = cp.arange(n_max, dtype=cp.float64) * adu_per_photon - offset

    n_photons = _photon_count_per_drop(adu_drop, photon_pts)  # [nDroplets]

    is_single = n_photons == 1
    is_multi = (n_photons > 1) & (n_photons <= GIANT_NPHOT)
    is_giant = n_photons > GIANT_NPHOT

    chunks = []

    # --- single-photon droplets: center of mass
    single_ind = drop_ind[is_single]
    if single_ind.size > 0:
        pos = _center_of_mass_batch(img, img_drop, single_ind)  # [N1, 2]
        chunks.append(pos)

    # --- multi-photon droplets (≤ MAX_PHOT): greedyguess kernel
    multi_ind = drop_ind[is_multi]
    multi_n = n_photons[is_multi]
    if multi_ind.size > 0:
        pos = _multi_photon_positions(img, img_drop, multi_ind, multi_n, adu_per_photon)
        if pos.size > 0:
            chunks.append(pos)

    # --- giant droplets (saturated blobs): vectorized integer + centroid placement
    giant_ind = drop_ind[is_giant]
    giant_n = n_photons[is_giant]
    if giant_ind.size > 0:
        pos = _giant_droplet_positions(
            img, img_drop, giant_ind, giant_n, adu_per_photon
        )
        if pos.size > 0:
            chunks.append(pos)

    if not chunks:
        return cp.zeros((0, 2), dtype=cp.float32)

    return cp.concatenate(chunks, axis=0)


# ══════════════════════════════════════════════════════════════════════════════
# CPU implementations (NumPy + SciPy)
# ══════════════════════════════════════════════════════════════════════════════


def _prepare_image_cpu(img, mask, comp_data, threshold, pixmax=None):
    out = img.copy()
    if mask is not None:
        out *= mask.astype(out.dtype)
    if pixmax is not None:
        out[out > pixmax] = 0.0
    out[out < comp_data * threshold] = 0.0
    return out


def _center_of_mass_batch_cpu(img, img_drop, drop_ind):
    """
    Sparse CoM: extract only labeled pixels then use np.bincount for groupby.
    Avoids 3 × full-image scipy.ndimage.sum passes — O(n_labeled + n_drops)
    instead of O(3 × n_pixels).
    """
    # one pass to find all labeled pixels (any droplet label > 0)
    r, c = np.where(img_drop > 0)
    if r.size == 0:
        return np.zeros((len(drop_ind), 2), dtype=np.float32)

    labels_all = img_drop[r, c]
    adus_all = img[r, c].astype(np.float64)

    # map label value → group index in drop_ind (drop_ind is sorted 1..N after relabel)
    max_label = int(img_drop.max())
    lut = np.full(max_label + 1, -1, dtype=np.int32)
    lut[drop_ind] = np.arange(len(drop_ind), dtype=np.int32)
    group = lut[labels_all]

    # discard pixels whose droplet was filtered out
    valid = group >= 0
    if not valid.any():
        return np.zeros((len(drop_ind), 2), dtype=np.float32)
    g = group[valid]
    w = adus_all[valid]
    rv = r[valid].astype(np.float64)
    cv = c[valid].astype(np.float64)

    n = len(drop_ind)
    total = np.bincount(g, weights=w, minlength=n)
    row_m = np.bincount(g, weights=w * rv, minlength=n)
    col_m = np.bincount(g, weights=w * cv, minlength=n)

    safe_total = np.where(total > 0, total, 1.0)
    return np.stack([row_m / safe_total, col_m / safe_total], axis=1).astype(np.float32)


@_njit(cache=True)
def _greedyguess_single_cpu(prows, pcols, padus, nphot, adu_per_photon):
    """
    Greedyguess for one droplet. Explicit loops → compatible with numba njit.
    When numba is unavailable the decorator is a no-op and this runs as plain Python.
    """
    npx = len(prows)
    rmin = int(prows[0])
    rmax = rmin
    cmin = int(pcols[0])
    cmax = cmin
    for k in range(1, npx):
        r = int(prows[k])
        c = int(pcols[k])
        if r < rmin:
            rmin = r
        if r > rmax:
            rmax = r
        if c < cmin:
            cmin = c
        if c > cmax:
            cmax = c

    snrows = rmax - rmin + 3
    sncols = cmax - cmin + 3

    timg = np.zeros((snrows, sncols))
    pxs = np.zeros(nphot)
    pys = np.zeros(nphot)
    nassigned = 0

    for k in range(npx):
        ii = int(prows[k]) - rmin + 1
        jj = int(pcols[k]) - cmin + 1
        val = float(padus[k]) / adu_per_photon
        ip = int(val)
        timg[ii, jj] = val - ip
        for _ in range(ip):
            if nassigned < nphot:
                pxs[nassigned] = float(ii)
                pys[nassigned] = float(jj)
                nassigned += 1

    for n in range(nassigned, nphot):
        # Explicit max-scan (np.argmax not supported by numba in nopython mode)
        bi = 0
        bj = 0
        bv = -1.0
        for ii in range(snrows):
            for jj in range(sncols):
                if timg[ii, jj] > bv:
                    bv = timg[ii, jj]
                    bi = ii
                    bj = jj

        i = bi
        j = bj
        t = 1.0 - bv
        timg[i, j] = 0.0

        vup = timg[i - 1, j] if i > 0 else -1.0
        vdown = timg[i + 1, j] if i < snrows - 1 else -1.0
        vleft = timg[i, j - 1] if j > 0 else -1.0
        vrite = timg[i, j + 1] if j < sncols - 1 else -1.0
        best = max(vup, vdown, vleft, vrite)

        if best == vup:
            pxs[n] = i - t
            pys[n] = float(j)
            if i > 0:
                timg[i - 1, j] = max(0.0, timg[i - 1, j] - t)
        elif best == vdown:
            pxs[n] = i + t
            pys[n] = float(j)
            if i < snrows - 1:
                timg[i + 1, j] = max(0.0, timg[i + 1, j] - t)
        elif best == vleft:
            pxs[n] = float(i)
            pys[n] = j - t
            if j > 0:
                timg[i, j - 1] = max(0.0, timg[i, j - 1] - t)
        else:
            pxs[n] = float(i)
            pys[n] = j + t
            if j < sncols - 1:
                timg[i, j + 1] = max(0.0, timg[i, j + 1] - t)

    return pxs - 1.0 + rmin, pys - 1.0 + cmin


@_njit(cache=True)
def _greedyguess_driver_serial(
    rows_s,
    cols_s,
    adus_s,
    pstart,
    pend,
    nphot_arr,
    offsets,
    adu_per_photon,
    out_r,
    out_c,
):
    """Fused serial driver: the whole per-droplet loop compiled, writing into
    pre-offset output slices. Removes the Python per-droplet call + list-append
    overhead of the old interpreted loop (~4× faster single-threaded)."""
    n = len(nphot_arr)
    for k in range(n):
        s = pstart[k]
        e = pend[k]
        np_ = nphot_arr[k]
        if e <= s or np_ == 0:
            continue
        pr, pc = _greedyguess_single_cpu(
            rows_s[s:e], cols_s[s:e], adus_s[s:e], np_, adu_per_photon
        )
        off = offsets[k]
        for j in range(np_):
            out_r[off + j] = pr[j]
            out_c[off + j] = pc[j]


@_njit(parallel=True, cache=True)
def _greedyguess_driver_par(
    rows_s,
    cols_s,
    adus_s,
    pstart,
    pend,
    nphot_arr,
    offsets,
    adu_per_photon,
    out_r,
    out_c,
):
    """OpenMP-style (numba prange) driver: one droplet per thread, writing into
    disjoint pre-offset output slices so there are no races. Mirrors the GPU
    one-thread-per-droplet RawKernel decomposition. Opt-in (parallel=True): only
    worth it for dense frames AND when not already parallel across events."""
    n = len(nphot_arr)
    for k in _prange(n):
        s = pstart[k]
        e = pend[k]
        np_ = nphot_arr[k]
        if e <= s or np_ == 0:
            continue
        pr, pc = _greedyguess_single_cpu(
            rows_s[s:e], cols_s[s:e], adus_s[s:e], np_, adu_per_photon
        )
        off = offsets[k]
        for j in range(np_):
            out_r[off + j] = pr[j]
            out_c[off + j] = pc[j]


def _multi_photon_positions_cpu(
    img, img_drop, drop_ind, n_photons_per_drop, adu_per_photon, parallel=False
):
    """Greedyguess photon placement for all multi-photon droplets.

    The per-droplet loop is fully compiled (numba) and writes into disjoint
    pre-offset output slices. `parallel=True` runs it over threads (numba
    prange); the serial fused path is the default since most of the speedup
    comes from the fusion, not the threads, and the producer is usually already
    parallel across events. Assumes every kept droplet is a valid multi (np_>=2,
    e>s) — guaranteed by the nph>1 selection — so the cumsum offsets are exact.
    """
    max_label = int(img_drop.max())
    is_multi_lut = np.zeros(max_label + 1, dtype=bool)
    is_multi_lut[drop_ind] = True

    rows, cols = np.where(is_multi_lut[img_drop])
    if rows.size == 0:
        return np.zeros((0, 2), dtype=np.float32)

    labels = img_drop[rows, cols]
    adus = img[rows, cols].astype(np.float32)
    sort_idx = np.argsort(labels, kind="stable")
    labels_s = labels[sort_idx]
    rows_s = rows[sort_idx].astype(np.int64)
    cols_s = cols[sort_idx].astype(np.int64)
    adus_s = adus[sort_idx]

    nphot_arr = np.asarray(n_photons_per_drop, dtype=np.int64)
    di64 = np.asarray(drop_ind, dtype=np.int64)
    pstart = np.searchsorted(labels_s, di64, side="left").astype(np.int64)
    pend = np.searchsorted(labels_s, di64, side="right").astype(np.int64)

    offsets = np.zeros(len(di64), dtype=np.int64)
    if len(di64) > 1:
        offsets[1:] = np.cumsum(nphot_arr[:-1])
    total = int(nphot_arr.sum())
    out_r = np.empty(total, dtype=np.float64)
    out_c = np.empty(total, dtype=np.float64)

    driver = _greedyguess_driver_par if parallel else _greedyguess_driver_serial
    driver(
        rows_s,
        cols_s,
        adus_s,
        pstart,
        pend,
        nphot_arr,
        offsets,
        float(adu_per_photon),
        out_r,
        out_c,
    )
    return np.stack([out_r, out_c], axis=1).astype(np.float32)


def _multi_photon_positions_cpu_par(
    img, img_drop, drop_ind, n_photons_per_drop, adu_per_photon
):
    """Back-compat alias for the threaded path (= parallel=True)."""
    return _multi_photon_positions_cpu(
        img, img_drop, drop_ind, n_photons_per_drop, adu_per_photon, parallel=True
    )


def _giant_droplet_positions_cpu(img, img_drop, giant_ind, giant_nphot, adu_per_photon):
    """CPU twin of _giant_droplet_positions (integer placement + centroid residual)."""
    n = int(giant_ind.size)
    if n == 0:
        return np.zeros((0, 2), dtype=np.float32)

    max_label = int(img_drop.max())
    is_g = np.zeros(max_label + 1, dtype=bool)
    is_g[giant_ind] = True
    r, c = np.where(is_g[img_drop])

    labels = img_drop[r, c]
    adus = img[r, c].astype(np.float64)
    rf, cf = r.astype(np.float64), c.astype(np.float64)

    lut = np.full(max_label + 1, -1, dtype=np.int32)
    lut[giant_ind] = np.arange(n, dtype=np.int32)
    grp = lut[labels]

    k = np.maximum(np.floor(adus / adu_per_photon).astype(np.int64), 0)
    int_rows = np.repeat(r.astype(np.float32), k)
    int_cols = np.repeat(c.astype(np.float32), k)

    sum_k = np.bincount(grp, weights=k.astype(np.float64), minlength=n)
    total = np.bincount(grp, weights=adus, minlength=n)
    row_m = np.bincount(grp, weights=adus * rf, minlength=n)
    col_m = np.bincount(grp, weights=adus * cf, minlength=n)
    safe = np.where(total > 0, total, 1.0)
    com_r = (row_m / safe).astype(np.float32)
    com_c = (col_m / safe).astype(np.float32)

    resid = np.maximum(
        np.asarray(giant_nphot, dtype=np.int64) - sum_k.astype(np.int64), 0
    )
    res_rows = np.repeat(com_r, resid)
    res_cols = np.repeat(com_c, resid)

    return np.stack(
        [np.concatenate([int_rows, res_rows]), np.concatenate([int_cols, res_cols])],
        axis=1,
    ).astype(np.float32)


def _dropletize_cpu(
    img,
    mask,
    comp_data,
    threshold,
    threshold_low,
    thres_adu,
    relabel,
    pixmax,
    max_drop_adu=None,
):
    if comp_data is None:
        comp_data = np.ones_like(img)

    img_thr = _prepare_image_cpu(img, mask, comp_data, threshold, pixmax)
    img_drop, n_drop = sndi.label(img_thr, structure=_FOOTPRINT_NP)

    if threshold != threshold_low:
        neighbor = sndi.maximum_filter(img_drop, footprint=_FOOTPRINT_NP)
        img_low = _prepare_image_cpu(img, mask, comp_data, threshold_low, pixmax)
        if relabel:
            neighbor[img_low == 0] = 0
            img_drop, n_drop = sndi.label(neighbor, structure=_FOOTPRINT_NP)
        else:
            img_drop = neighbor

    n_all = int(n_drop)
    if n_all == 0:
        return {
            "nDroplets_all": 0,
            "nDroplets": 0,
            "img": img_thr,
            "img_drop": img_drop,
            "drop_ind": np.array([], dtype=np.int32),
            "adu_drop": np.array([], dtype=np.float32),
        }

    drop_ind = np.arange(1, n_drop + 1, dtype=np.int32)
    adu_drop = sndi.sum(img_thr, labels=img_drop, index=drop_ind)

    if thres_adu is not None or max_drop_adu is not None:
        keep = np.ones(adu_drop.shape, dtype=bool)
        if thres_adu is not None:
            keep &= adu_drop >= thres_adu
        if max_drop_adu is not None:  # veto saturated/bright blobs
            keep &= adu_drop <= max_drop_adu
        lut = np.zeros(n_drop + 1, dtype=img_drop.dtype)
        kept_labels = drop_ind[keep]
        lut[kept_labels] = kept_labels
        img_drop = lut[img_drop]
        drop_ind = kept_labels
        adu_drop = adu_drop[keep]

    return {
        "nDroplets_all": n_all,
        "nDroplets": int(drop_ind.size),
        "img": img_thr,
        "img_drop": img_drop,
        "drop_ind": drop_ind,
        "adu_drop": np.asarray(adu_drop, dtype=np.float32),
    }


def _find_photons_cpu(droplet_dict, adu_per_photon, photon_pts):
    img = droplet_dict["img"]
    img_drop = droplet_dict["img_drop"]
    drop_ind = droplet_dict["drop_ind"]
    adu_drop = droplet_dict["adu_drop"]

    if drop_ind.size == 0:
        return np.zeros((0, 2), dtype=np.float32)

    if photon_pts is None:
        offset = adu_per_photon * 0.5
        amax = float(adu_drop.max()) if adu_drop.size else 0.0
        n_max = int(amax / adu_per_photon) + 3
        photon_pts = np.arange(n_max, dtype=np.float64) * adu_per_photon - offset

    n_photons = (np.searchsorted(photon_pts, adu_drop, side="right") - 1).astype(
        np.int32
    )

    is_single = n_photons == 1
    is_multi = (n_photons > 1) & (n_photons <= GIANT_NPHOT)
    is_giant = n_photons > GIANT_NPHOT
    chunks = []

    single_ind = drop_ind[is_single]
    if single_ind.size > 0:
        chunks.append(_center_of_mass_batch_cpu(img, img_drop, single_ind))

    multi_ind = drop_ind[is_multi]
    multi_n = n_photons[is_multi]
    if multi_ind.size > 0:
        pos = _multi_photon_positions_cpu(
            img, img_drop, multi_ind, multi_n, adu_per_photon
        )
        if pos.size > 0:
            chunks.append(pos)

    giant_ind = drop_ind[is_giant]
    giant_n = n_photons[is_giant]
    if giant_ind.size > 0:
        pos = _giant_droplet_positions_cpu(
            img, img_drop, giant_ind, giant_n, adu_per_photon
        )
        if pos.size > 0:
            chunks.append(pos)

    if not chunks:
        return np.zeros((0, 2), dtype=np.float32)

    return np.concatenate(chunks, axis=0)


# ══════════════════════════════════════════════════════════════════════════════
# Public dispatchers — automatically select GPU or CPU path
# ══════════════════════════════════════════════════════════════════════════════


def dropletize(
    img,
    mask=None,
    comp_data=None,
    threshold=5.0,
    threshold_low=None,
    thres_adu=None,
    relabel=True,
    pixmax=None,
    max_drop_adu=None,
):
    """
    Find droplets in a 2D detector image.
    Automatically uses the GPU path if img is a CuPy array, otherwise CPU.

    Parameters
    ----------
    img           : cp.ndarray or np.ndarray, 2D detector frame.
    mask          : array or None, 1=valid 0=bad pixel. Use to exclude bright/beam
                    ROIs so they never form giant droplets (fast, exact path).
    comp_data     : array or None, per-pixel threshold scale (e.g. rms/gain).
                    Defaults to all-ones (threshold is in raw ADU).
    threshold     : float, high threshold in units of comp_data.
    threshold_low : float or None, low threshold for droplet expansion.
                    Defaults to threshold (no expansion).
    thres_adu     : float or None, minimum total ADU per droplet (lower veto).
    relabel       : bool, re-run label after expanding with threshold_low.
    pixmax        : float or None, clip pixels above this before thresholding
                    (per-pixel hot-pixel/saturation ceiling).
    max_drop_adu  : float or None, drop droplets whose total ADU exceeds this
                    (per-droplet upper veto — excludes saturated/beam blobs from
                    photon counting).

    Returns
    -------
    dict with keys: nDroplets_all (int), nDroplets (int),
                    img, img_drop, drop_ind, adu_drop.
    Array types match the input (CuPy on GPU path, NumPy on CPU path).
    """
    if threshold_low is None:
        threshold_low = threshold

    if _is_gpu(img):
        return _dropletize_gpu(
            img,
            mask,
            comp_data,
            threshold,
            threshold_low,
            thres_adu,
            relabel,
            pixmax,
            max_drop_adu,
        )
    return _dropletize_cpu(
        img,
        mask,
        comp_data,
        threshold,
        threshold_low,
        thres_adu,
        relabel,
        pixmax,
        max_drop_adu,
    )


def find_photons(droplet_dict, adu_per_photon, photon_pts=None):
    """
    Extract photon positions from the output of dropletize().

    Single-photon droplets  → intensity-weighted center of mass (sub-pixel).
    Multi-photon droplets   → top-N pixels by ADU, where N = photon count.

    Parameters
    ----------
    droplet_dict   : dict, output of dropletize().
    adu_per_photon : float, expected ADU for one photon.
    photon_pts     : array or None
                     ADU bin edges: photon_pts[n] = lower bound for n photons.
                     Defaults to n * adu_per_photon - 0.5 * adu_per_photon.

    Returns
    -------
    array of shape [N_photons, 2], dtype float32 — (row, col).
    Shape (0, 2) if no photons found.
    Array type matches the arrays in droplet_dict.
    """
    if _is_gpu(droplet_dict["img"]):
        return _find_photons_gpu(droplet_dict, adu_per_photon, photon_pts)
    return _find_photons_cpu(droplet_dict, adu_per_photon, photon_pts)
