import psana
import numpy as np
from weninc_algorithms import cfd, find_droplets

def process_waveform(evt, config):
    det = config['det']
    wf = det.waveform(evt)
    if wf is None:
        return [np.zeros(0, dtype=np.float32) for c in config['channels']]
    taxis = det.wftime(evt)
    peaks = [cfd(taxis[c], wf[c], config['fraction'], config['delay'], config['threshold'], 21).astype(np.float32) for c in config['channels']]
    return peaks

class QuadAnodeDLD:
    def __init__(self, dld_config, mcp_config):
        self.dld_config = dld_config
        self.mcp_config = mcp_config
        
    def process_event(self, evt, batch):
        peaks = process_waveform(evt, self.dld_config)
        keys = ['x1_peaks', 'x2_peaks', 'y1_peaks', 'y2_peaks']
        for i, k in enumerate(keys):
            batch.add_data(k, peaks[i])
        
        keys = ['x1_counts', 'x2_counts', 'y1_counts', 'y2_counts']
        for i, k in enumerate(keys):
            batch.add_data(k, peaks[i].size)
        
        peaks, = process_waveform(evt, self.mcp_config)
        batch.add_data('mcp_counts', peaks.size)
        batch.add_data('mcp_peaks', peaks)
        
class VonHamos:
    def __init__(self, config, run):
        self.config = config
        self.mask = self.config['det'].mask(run, status=True, edges=True, central=True)
        
    def missing_values(self, batch):
        batch.add_data('ndroplets', 0)
        batch.add_data('x', np.zeros(0, dtype=np.float32))
        batch.add_data('y', np.zeros(0, dtype=np.float32))
        batch.add_data('adu', np.zeros(0, dtype=np.float32))
        
    def process_event(self, evt, batch):
        img = self.config['det'].calib(evt)
        if img is None:
            print 'Missing detector in VonHamos'
            self.missing_values(batch)
            return
        img *= self.mask
        img = np.squeeze(img) # remove leading dimension with size 1 of jungfrau
        ndroplets, x, y, adu = find_droplets(img, self.config['seed_threshold'], self.config['join_threshold'])
        if ndroplets > 0:
            batch.add_data('ndroplets', ndroplets)
            batch.add_data('x', x)
            batch.add_data('y', y)
            batch.add_data('adu', adu)
        else:
            self.missing_values(batch)


class VonHamosNoDroplet:
    def __init__(self, config, run):
        self.config = config
        self.mask = self.config['det'].mask(run, status=True, edges=True, central=True)
        
    def missing_values(self, batch):
        batch.add_data('ndroplets', 0)
        batch.add_data('x', np.zeros(0, dtype=np.float32))
        batch.add_data('y', np.zeros(0, dtype=np.float32))
        batch.add_data('adu', np.zeros(0, dtype=np.float32))
        
    def process_event(self, evt, batch):
        img = self.config['det'].calib(evt)
        if img is None:
            print 'Missing detector in VonHamos'
            self.missing_values(batch)
            return
        img *= self.mask
        img = np.squeeze(img)
        
        img[:, :250] = 0.0
        ind = np.where(img > self.config['threshold'])
        
        if ind:
            x, y = ind
            npixels = x.shape[0]
            adu = img[ind]
            batch.add_data('ndroplets', npixels)
            batch.add_data('x', x)
            batch.add_data('y', y)
            batch.add_data('adu', adu)
        else:
            self.missing_values(batch)
            
