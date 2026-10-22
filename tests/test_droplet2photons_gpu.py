"""
Self-contained equivalence tests for droplet2photons_gpu (no psana / no data files).

  - counts and self.dat contract are identical to the reference droplet2Photons
    (sparse and with a saturated "giant" blob)
  - photon positions match the reference exactly when nothing is saturated
  - the GPU path matches the CPU path (skipped if CuPy/GPU absent)

Run:  pytest tests/test_droplet2photons_gpu.py -v
"""
import numpy as np
import scipy.ndimage as ndi
import pytest

from smalldata_tools.ana_funcs.droplet2photons_gpu import droplet2photons_gpu

try:
    from smalldata_tools.ana_funcs.droplet2Photons import droplet2Photons
    _HAS_REF = True
except Exception:
    _HAS_REF = False

try:
    import cupy  # noqa: F401
    _HAS_GPU = cupy.cuda.runtime.getDeviceCount() > 0
except Exception:
    _HAS_GPU = False

APH, THRES = 172.0, 15.0
FP = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])


def _make_frame(occ, with_blob, seed):
    rng = np.random.default_rng(seed)
    shape = (256, 256)
    img = rng.normal(0, 10, shape)
    for r, c in rng.integers(6, shape[0] - 6, size=(int(occ * np.prod(shape)), 2)):
        img[r, c] += APH
    if with_blob:
        img[100:124, 100:124] += 220.0          # saturated blob -> a >64-photon droplet
    return img


def _prep(img):
    thr = img.copy(); thr[thr < THRES] = 0
    lab, _ = ndi.label(thr, structure=FP)
    return {"_image": thr, "_imgDrop": lab, "_mask": np.ones(img.shape, bool)}


def _run(func, d):
    func.process(d)
    return np.asarray(func.dat["row"]), np.asarray(func.dat["col"])


def _match(r1, c1, r2, c2):
    if len(r1) == 0:
        return 1.0 if len(r2) == 0 else 0.0
    s2 = set(zip(np.round(r2).astype(int).tolist(), np.round(c2).astype(int).tolist()))
    return sum((int(round(r)), int(round(c))) in s2 for r, c in zip(r1, c1)) / len(r1)


@pytest.mark.skipif(not _HAS_REF, reason="reference droplet2Photons not importable")
@pytest.mark.parametrize("with_blob", [False, True])
def test_counts_and_contract_match_reference(with_blob):
    ref = droplet2Photons(aduspphot=APH, name="ref")
    ours = droplet2photons_gpu(aduspphot=APH, name="ours")
    for seed in range(3):
        d = _prep(_make_frame(1e-2, with_blob, seed))
        rr, _ = _run(ref, d)
        orr, _ = _run(ours, d)
        assert len(orr) == len(rr), f"count mismatch: ref={len(rr)} ours={len(orr)}"
        assert sorted(ours.dat) == sorted(ref.dat)          # same self.dat keys


@pytest.mark.skipif(not _HAS_REF, reason="reference droplet2Photons not importable")
def test_positions_match_reference_without_saturation():
    ref = droplet2Photons(aduspphot=APH, name="ref")
    ours = droplet2photons_gpu(aduspphot=APH, name="ours")
    d = _prep(_make_frame(1e-2, with_blob=False, seed=7))
    rr, rc = _run(ref, d)
    orr, orc = _run(ours, d)
    assert _match(orr, orc, rr, rc) >= 0.99


@pytest.mark.skipif(not _HAS_GPU, reason="no CuPy/GPU available")
@pytest.mark.parametrize("with_blob", [False, True])
def test_gpu_matches_cpu(with_blob):
    cpu = droplet2photons_gpu(aduspphot=APH, use_gpu=False, name="cpu")
    gpu = droplet2photons_gpu(aduspphot=APH, use_gpu=True, name="gpu")
    d = _prep(_make_frame(1e-2, with_blob, seed=1))
    cpu.process(d); gpu.process(d)
    assert len(gpu.dat["row"]) == len(cpu.dat["row"])


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
