"""
Example: run droplet2photons_gpu in a live smalldata pipeline.
    psana1 DataSource → dropletFunc (live droplet finding) → droplet2photons_gpu (drop-in)
The drop-in receives the LIVE _imgDrop produced by dropletFunc each event — no synthetic labeling.

Usage (in a psana1 / "ana" environment):
    python example_pipeline_psana1.py [exp] [run] [detector]
    # defaults: xpplw3920 130 epix_alc4

Note: this drives dropletFunc standalone (no DetObject), so _compData and _is_tiled —
normally set by DetObject.setFromFunc in a producer run — are set by hand below.
"""
import sys, numpy as np, psana

from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.droplet2photons_gpu import droplet2photons_gpu

EXP = sys.argv[1] if len(sys.argv) > 1 else 'xpplw3920'
RUN = int(sys.argv[2]) if len(sys.argv) > 2 else 130
DET = sys.argv[3] if len(sys.argv) > 3 else 'epix_alc4'
APH, THRES, N = 172, 15, 50

ds  = psana.DataSource(f'exp={EXP}:run={RUN}')
det = psana.Detector(DET, ds.env())

df = None
nphot, ndrop, n = [], [], 0
for evt in ds.events():
    calib = det.calib(evt)
    if calib is None:
        continue
    img = (calib * 17).astype(np.float64)          # keV → ADU (detector-specific scale)
    if df is None:
        mask = np.ones(img.shape, dtype=int)
        df = dropletFunc(threshold=THRES, thresholdLow=THRES, useRms=False, name='droplet', mask=mask)
        df._compData = np.ones_like(img)            # normally set by DetObject.setFromFunc
        df._is_tiled = (img.ndim == 3)
        df.addFunc(droplet2photons_gpu(aduspphot=APH, name='photons'))  # CPU; use_gpu=True for GPU
        subs = [k for k, v in df.__dict__.items() if hasattr(v, 'process') and k != 'dat']
        print(f"pipeline wired: dropletFunc → sub-funcs {subs}")
    df.process(img)                                 # find (live) → processFuncs → drop-in
    nphot.append(len(df.photons.dat['row']))
    nd = df.dat['_imgDrop']
    ndrop.append(int(np.nanmax(nd)) if nd.size else 0)
    n += 1
    if n >= N:
        break

print(f"\n{EXP} r{RUN} {DET}:  {n} events")
print(f"  droplets/event (live dropletFunc): {np.mean(ndrop):.1f}")
print(f"  photons/event  (drop-in)         : {np.mean(nphot):.1f}")
print(f"  self.dat keys                    : {sorted(df.photons.dat)}")
