import numpy as np
import argparse
import sys
import psana as ps
from pathlib import Path
import datetime
import h5py as h5

import smalldata_tools.DetObject as dobj
import smalldata_tools.ana_funcs.svd_waveform_processing as proc


def make_basis(exp_name, run, det_name, nWaveforms=500, bkg_idx=None, n_c=2, roi=None, channel=None):
    """ Create the SVD basis for the waveform fitting. The output is saved in a h5 file stored
    in the calib directory of the relevant experiment.
    
    Args:
        exp_name (str): experiment name
        run (int): run number
        det_name (str): psana name for the detector (alias)
        nWaveforms (int): number of waveform to use to build the basis
        bkg_idx (int): background index. np.mean(waveform[:bkg_idx]) will be subtracted
        roi (list, array): two index specifying the region of interest of the waveform
        channel (int): detector channel index. Necessary if the digitizer has more than 1 channel
    """
    hutch = exp_name[:3]
    savePath = Path('/reg/d/psdm/{}/{}/hdf5/basis/'.format(hutch, exp_name))
    if not savePath.exists():
        savePath.mkdir()

    dstr = 'exp={}:run={}'.format(exp_name, run)
    print('\n'+dstr+'\n')
    ds = ps.MPIDataSource(dstr)
    det = dobj.DetObject(det_name, ds.env(), int(run))

    if roi is None:
        roi = [0,int(1e6)]

    """ ---------------------------- GET BASIS ---------------------------- """
    waveforms = []
    ii = 0
    for nevt,evt in enumerate(ds.events()):
        if ii>nWaveforms:
            break

        try:
            det.getData(evt)
            wave = np.squeeze(det.evt.dat)
            if channel is not None:
                wave = wave[channel]
#             print(wave.shape)
        except:
            continue
        
        if np.asarray(wave).ndim>1:
            raise ValueError('The dimension of the waveform is larger than 1. Perhaps a channel argument is missing?')
            
        if bkg_idx is not None:
            wave = wave - np.mean(wave[:bkg_idx])
        wave = wave[roi[0]:roi[1]]
        waveforms.append(wave)
        ii+=1

    waveforms = np.asarray(waveforms)
#     print(waveforms.shape)
    A, proj, svd = proc.get_basis_and_projector(waveforms, n_components=n_c)


    """ ---------------------------- SAVE DATA ---------------------------- """
    fname = 'wave_basis_' + det_name + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")+ '.h5'
    print(fname)
    fname = savePath / fname
#     fname = Path('./') / fname
    with h5.File(fname,'w') as f:
        dset = f.create_dataset("projector", data=proj)
        dset = f.create_dataset("components", data=svd.components_)
        dset = f.create_dataset("singular_values", data=svd.singular_values_)
        dset = f.create_dataset("ref_waveforms", data=waveforms)
        dset = f.create_dataset("roi", data=roi)
        dset = f.create_dataset("background_index", data=bkg_idx)
        dset = f.create_dataset("channel", data=channel)
    print('Basis file saved as {}.'.format(fname))
    return fname



""" ---------------------------- PARSER ---------------------------- """
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Create waveform basis file for given detector.')
    parser.add_argument('-e','--exp', type=str, help='Experiment name')
    parser.add_argument('-r','--run', type=int, help='Run to be analyzed')
    parser.add_argument('-d', '--detector', type=str, help='psana detector name')
    parser.add_argument('-n', '--nWaveform', type=int, nargs='?', default=500, 
                        help='Number of waveforms taken to build the basis')
    parser.add_argument('-b', '--baseline', type=int, nargs='?', default=None, 
                        help='Subtract b-average baseline of the waveform')
    parser.add_argument('-nc', '--nComponents', type=int, nargs='?', default=2, 
                        help='Number of SVD components in the basis')
    parser.add_argument('--roi', type=int, nargs='*', default=None, help='ROI: indx1 idx2')
    parser.add_argument('--channel', type=int, nargs=1, default=None, 
                        help='Channel of the digitizer (need if there are several channels in the detector)')

    args = parser.parse_args()

    exp_name = args.exp
    run = args.run
    det_name = args.detector
    nWaveform = args.nWaveform
    b = args.baseline
    n_c = args.nComponents
    roi = args.roi
    channel = args.channel
    
    make_basis(exp_name, run, det_name, nWaveform=nWaveform, bkg_idx=b, n_c=n_c, roi=roi, channel=channel)