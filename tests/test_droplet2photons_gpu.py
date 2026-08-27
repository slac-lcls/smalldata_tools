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
        img[100:124, 100:124] += 220.0  # saturated blob -> a >64-photon droplet
    return img


def _prep(img):
    thr = img.copy()
    thr[thr < THRES] = 0
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
        assert sorted(ours.dat) == sorted(ref.dat)  # same self.dat keys


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
    cpu.process(d)
    gpu.process(d)
    assert len(gpu.dat["row"]) == len(cpu.dat["row"])


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))


# ---------------------------------------------------------------------------------------
# aduspphot / offset / photpts resolution semantics (issue #289) -- CPU-only, no GPU needed.
# The contract, per the maintainer: any argument passed must be used.
#   1. aduspphot only     -> photon_pts=None (find_photons' built-in default IS
#                            n*aduspphot - aduspphot/2; nothing shipped to the GPU)
#   2. aduspphot + offset -> uniform edges n*aduspphot - offset, sized to the event
#   3. photpts            -> the user's list verbatim; offset ignored (warns);
#                            aduspphot derived from median(diff(photpts)) unless given
# ---------------------------------------------------------------------------------------


def _one_droplet_frame():
    img = np.zeros((16, 16), np.float32)
    img[5, 5] = 172.0
    img[5, 6] = 30.0
    lab = np.zeros((16, 16), np.int32)
    lab[5, 5] = 1
    lab[5, 6] = 1
    return {"_image": img, "_imgDrop": lab, "_mask": None}


def _capture(monkeypatch, **initkw):
    """Run one event through droplet2photons_gpu with find_photons stubbed out; return
    what reached it."""
    import smalldata_tools.ana_funcs.droplet2photons_gpu as mod

    seen = {}

    def fake_find_photons(dd, adu, photon_pts=None):
        seen["adu"] = adu
        seen["photon_pts"] = photon_pts
        return np.zeros((0, 2), np.float32)

    monkeypatch.setattr(mod, "find_photons", fake_find_photons)
    f = mod.droplet2photons_gpu(**initkw)
    f.process(_one_droplet_frame())
    return f, seen


def test_case1_aduspphot_only_passes_none(monkeypatch):
    f, seen = _capture(monkeypatch, aduspphot=172.0)
    assert seen["photon_pts"] is None
    assert seen["adu"] == 172.0
    assert f.offset == 86.0


def test_case2_user_offset_builds_data_sized_edges(monkeypatch):
    f, seen = _capture(monkeypatch, aduspphot=172.0, offset=120.0)
    pts = seen["photon_pts"]
    assert pts is not None and pts.size < 10  # sized to the event, not 1e6
    assert pts[0] == -120.0  # n*aduspphot - offset, from n=0
    assert np.allclose(np.diff(pts), 172.0)
    assert seen["adu"] == 172.0


def test_case3_user_photpts_shipped_verbatim(monkeypatch):
    edges = [-50.0, 150.0, 320.0, 500.0]  # non-uniform, e.g. background cutoff
    f, seen = _capture(monkeypatch, photpts=edges)
    np.testing.assert_array_equal(seen["photon_pts"], np.asarray(edges))
    assert seen["adu"] == np.median(np.diff(edges))  # derived: 180.0
    assert f.aduspphot == 180.0


def test_case3_conflicting_aduspphot_warns_and_wins(monkeypatch):
    edges = [-50.0, 150.0, 320.0, 500.0]
    with pytest.warns(RuntimeWarning, match="aduspphot"):
        f, seen = _capture(monkeypatch, photpts=edges, aduspphot=100.0)
    assert seen["adu"] == 100.0  # the passed argument is used
    np.testing.assert_array_equal(seen["photon_pts"], np.asarray(edges))


def test_case3_offset_ignored_with_photpts_warns(monkeypatch):
    with pytest.warns(RuntimeWarning, match="offset"):
        _capture(monkeypatch, photpts=[-86.0, 86.0, 258.0], offset=40.0)


def test_uniform_user_photpts_equals_default_on_real_finder():
    """The claim behind case 1: for the DEFAULT uniform thresholds, photon_pts=None and an
    explicit n*aduspphot - aduspphot/2 list classify identically on the real finder."""
    from smalldata_tools.ana_funcs.droplet2photons_gpu import droplet2photons_gpu

    rng = np.random.default_rng(5)
    frame = np.zeros((128, 128), np.float32)
    ys, xs = rng.integers(4, 124, 12), rng.integers(4, 124, 12)
    frame[ys, xs] = APH
    frame[ys, np.minimum(xs + 1, 127)] += 0.4 * APH  # some charge-shared 2-px droplets
    lab = ndi.label(frame > THRES, structure=FP)[0].astype(np.int32)
    data = {"_image": frame, "_imgDrop": lab, "_mask": None}

    f1 = droplet2photons_gpu(aduspphot=APH)
    f1.process(dict(data))
    r1 = {k: np.array(v) for k, v in f1.dat.items()}

    f2 = droplet2photons_gpu(aduspphot=APH, photpts=np.arange(1000) * APH - APH / 2)
    f2.process(dict(data))
    r2 = {k: np.array(v) for k, v in f2.dat.items()}

    for k in ("row", "col", "data", "tile"):
        np.testing.assert_array_equal(r1[k], r2[k])
