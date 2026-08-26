"""
Self-contained equivalence tests for azimuthalBinning_gpu (no psana / no data files).

  - the CPU default is the parent, byte-for-byte (the subclass must not perturb it)
  - the sparse operator reproduces the parent's bincount to floating-point tolerance,
    with masked pixels, several phi bins and applyCorrection=False
  - the parent's rbin path is documented as broken upstream (pre-existing)
  - the GPU path matches the CPU path (skipped if CuPy/GPU absent)
  - doCake_batch matches a loop of doCake

Run:  pytest tests/test_azimuthalBinning_gpu.py -v
"""

import numpy as np
import pytest

from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning
from smalldata_tools.ana_funcs.azimuthalBinning_gpu import azimuthalBinning_gpu

try:
    import cupy  # noqa: F401

    _HAS_GPU = cupy.cuda.runtime.getDeviceCount() > 0
except Exception:
    _HAS_GPU = False

SHAPE = (128, 160)
PIXSIZE = 100e-6  # m


def _geometry(shape=SHAPE, masked=True, seed=0):
    """A centred detector plus an optional mask, in the form the parent's setFromFunc wants."""
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(np.float64)
    x = (xx - nx / 2.0) * PIXSIZE * 1e6  # micron, as the parent expects
    y = (yy - ny / 2.0) * PIXSIZE * 1e6
    mask = np.zeros(shape, dtype=bool)
    if masked:
        rng = np.random.default_rng(seed)
        mask[rng.random(shape) < 0.05] = True  # 5% bad pixels
        mask[:, 0] = True  # a whole dead column
    return x, y, mask


def _make(cls, phiBins=1, rbin=None, masked=True, **kw):
    x, y, mask = _geometry(masked=masked)
    f = cls(
        center=(0.0, 0.0),
        dis_to_sam=100.0,
        eBeam=9.5,
        phiBins=phiBins,
        qbin=5e-3,
        rbin=rbin,
        userMask=~mask,  # parent inverts: userMask is GOOD pixels
        **kw,
    )
    # _setup() consumes x/y/z_off/_mask off the instance, so set them before calling it. We
    # call _setup() directly rather than setFromFunc(), which would run it before x/y exist.
    # setFromDet() would normally supply z/z_off; for a flat detector normal to the beam
    # z_off = nanmean(z) - z is identically zero.
    f.x, f.y = x.ravel(), y.ravel()
    f.z = np.zeros_like(f.x)
    f.z_off = np.zeros_like(f.x)
    f._mask = mask.ravel()
    f._setup()
    return f


def _frames(n=4, seed=1):
    rng = np.random.default_rng(seed)
    base = rng.poisson(30, (n,) + SHAPE).astype(np.float64)
    yy, xx = np.mgrid[0 : SHAPE[0], 0 : SHAPE[1]]
    r = np.hypot(yy - SHAPE[0] / 2, xx - SHAPE[1] / 2)
    base += 200.0 * np.exp(-((r - 30.0) ** 2) / 8.0)  # a ring, so q structure is real
    return base


@pytest.mark.parametrize("phiBins,rbin", [(1, None), (8, None)])
@pytest.mark.parametrize("applyCorrection", [True, False])
def test_cpu_matches_parent(phiBins, rbin, applyCorrection):
    """The CPU default defers to the parent, so it must be identical, not merely close."""
    ref = _make(azimuthalBinning, phiBins=phiBins, rbin=rbin)
    new = _make(azimuthalBinning_gpu, phiBins=phiBins, rbin=rbin)
    img = _frames(1)[0]
    a = ref.doCake(img.copy(), applyCorrection=applyCorrection)
    b = new.doCake(img.copy(), applyCorrection=applyCorrection)
    np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("phiBins", [1, 8])
def test_operator_reproduces_bincount(phiBins):
    """M @ img must equal the parent's bincount to float tolerance (different summation order)."""
    ref = _make(azimuthalBinning, phiBins=phiBins)
    new = _make(azimuthalBinning_gpu, phiBins=phiBins)
    img = _frames(1)[0]

    expected = ref.doCake(img.copy())
    M = new._operator(True)  # CPU scipy operator
    good = new._mask.ravel() == 0
    got = (M @ img.ravel()[good])[: new.nq * new.nphi].reshape(new.nphi, new.nq)
    got = got / new.Cake_norm
    np.testing.assert_allclose(got, expected, rtol=1e-12, atol=0)


def test_masked_pixels_excluded():
    """A masked pixel must not contribute, however bright."""
    new = _make(azimuthalBinning_gpu, masked=True)
    img = _frames(1)[0]
    a = new.doCake(img.copy())
    bad = np.where(new._mask.ravel())[0][0]
    img2 = img.copy()
    img2.ravel()[bad] += 1e9
    b = new.doCake(img2)
    np.testing.assert_array_equal(a, b)


@pytest.mark.skipif(not _HAS_GPU, reason="CuPy/GPU not available")
@pytest.mark.parametrize("phiBins", [1, 8])
def test_gpu_matches_cpu(phiBins):
    cpu = _make(azimuthalBinning_gpu, phiBins=phiBins, use_gpu=False)
    gpu = _make(azimuthalBinning_gpu, phiBins=phiBins, use_gpu=True)
    img = _frames(1)[0]
    a = cpu.doCake(img.copy())
    b = gpu.doCake(img.copy())
    assert isinstance(b, np.ndarray)  # host contract: small_data rejects device arrays
    np.testing.assert_allclose(b, a, rtol=1e-11, atol=0)


@pytest.mark.skipif(not _HAS_GPU, reason="CuPy/GPU not available")
def test_batch_matches_loop():
    gpu = _make(azimuthalBinning_gpu, use_gpu=True)
    imgs = _frames(6)
    loop = np.stack([gpu.doCake(im.copy()) for im in imgs])
    batch = gpu.doCake_batch(imgs)
    assert isinstance(batch, np.ndarray)
    np.testing.assert_allclose(batch, loop, rtol=1e-11, atol=0)


@pytest.mark.skipif(not _HAS_GPU, reason="CuPy/GPU not available")
def test_gpu_applies_calibration_images_to_device_frame():
    """darkImg/gainImg configured + a GPU-resident frame: the advertised fast path. This used
    to raise (device frame minus host darkImg); now the calibration images are device-cached
    and applied on the GPU, matching the parent's CPU result."""
    rng = np.random.default_rng(7)
    dark = rng.random(SHAPE) * 5.0
    gainim = rng.random(SHAPE) * 0.5 + 0.75
    cpu = _make(azimuthalBinning_gpu, use_gpu=False, darkImg=dark, gainImg=gainim)
    gpu = _make(azimuthalBinning_gpu, use_gpu=True, darkImg=dark, gainImg=gainim)
    img = _frames(1)[0]
    a = cpu.doCake(img.copy())
    b = gpu.doCake(cupy.asarray(img))  # device-resident input
    assert isinstance(b, np.ndarray)
    np.testing.assert_allclose(b, a, rtol=1e-11, atol=0)


@pytest.mark.skipif(not _HAS_GPU, reason="CuPy/GPU not available")
def test_batch_applies_calibration_like_loop():
    """doCake_batch must produce what a loop of doCake produces when darkImg/gainImg are set
    (the GPU batch path used to skip both while the CPU fallback applied them)."""
    rng = np.random.default_rng(8)
    dark = rng.random(SHAPE) * 5.0
    gainim = rng.random(SHAPE) * 0.5 + 0.75
    gpu = _make(azimuthalBinning_gpu, use_gpu=True, darkImg=dark, gainImg=gainim)
    imgs = _frames(5)
    loop = np.stack([gpu.doCake(im.copy()) for im in imgs])
    batch = gpu.doCake_batch(imgs)
    np.testing.assert_allclose(batch, loop, rtol=1e-11, atol=0)


def test_batch_cpu_fallback_matches_loop_and_preserves_inputs():
    """The CPU fallback of doCake_batch, with darkImg/gainImg set: must equal a loop of
    doCake AND must not mutate the caller's batch (the parent's doCake corrects in place;
    the fallback therefore has to hand it copies)."""
    rng = np.random.default_rng(9)
    dark = rng.random(SHAPE) * 5.0
    gainim = rng.random(SHAPE) * 0.5 + 0.75
    f = _make(azimuthalBinning_gpu, use_gpu=False, darkImg=dark, gainImg=gainim)
    imgs = _frames(4)
    keep = imgs.copy()
    loop = np.stack([f.doCake(im.copy()) for im in imgs])
    batch = f.doCake_batch(imgs)
    np.testing.assert_array_equal(imgs, keep)  # caller's data untouched
    np.testing.assert_allclose(batch, loop, rtol=0)


def test_falls_back_without_cupy(monkeypatch):
    """use_gpu=True on a CPU-only node must warn once and still produce the parent's answer."""
    import smalldata_tools.ana_funcs.azimuthalBinning_gpu as mod

    monkeypatch.setattr(mod, "_HAS_CUPY", False)
    ref = _make(azimuthalBinning)
    new = _make(azimuthalBinning_gpu, use_gpu=True)
    img = _frames(1)[0]
    with pytest.warns(RuntimeWarning):
        b = new.doCake(img.copy())
    np.testing.assert_array_equal(b, ref.doCake(img.copy()))


def test_parent_rbin_path_is_broken_upstream():
    """Documents a PRE-EXISTING bug in the parent, unrelated to this change.

    azimuthalBinning.doCake() picks nradial = self.nr when rbin is set and bincounts into
    nr * nphi bins, but then reshapes unconditionally to (nphi, self.nq)
    (azimuthalBinning.py:385). Whenever nr != nq that raises. The GPU path here reshapes to
    (nphi, nradial), i.e. it does the consistent thing -- so this is the one place the two
    deliberately differ, and it differs by not crashing. Left as a separate one-line fix
    rather than folded into this PR.
    """
    ref = _make(azimuthalBinning, rbin=20.0)
    assert ref.nr != ref.nq, "test geometry no longer exercises the mismatch"
    with pytest.raises(ValueError, match="reshape"):
        ref.doCake(_frames(1)[0])
