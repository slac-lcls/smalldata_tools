import logging
import os

import sys
import psana
import pytest 

import smalldata_tools
from smalldata_tools.DetObject import DetObject
from tests.conftest import datasource, detector, getPhotons, getDroplets, getD2P, getSparsify
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.droplet2Photons import droplet2Photons
from smalldata_tools.ana_funcs.sparsifyFunc import sparsifyFunc
from smalldata_tools.SmallDataUtils import getUserData

import tables

logger = logging.getLogger(__name__)
#logger.setLevel(logging.DEBUG)

logger.info('Loading detector: Epix630')

@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 631}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'epix_alc3'}], indirect=True)
def test_detector_type(datasource, detector):
    logger.info('Running detector type test')
    det = detector
    assert(isinstance(det, smalldata_tools.DetObject.EpixObject))
    logger.info('Pass the test')
    
@pytest.mark.parametrize('getPhotons', [{'ADU_per_photon': 9.5, 'thresADU': 6.8}], indirect=True) 
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 631}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'epix_alc3'}], indirect=True)
def test_Photons(datasource, detector, getPhotons):
    logger.info('Running test for Photons function to check output format.')
    try:
        func_kwargs = getPhotons
    except Exception as e:
        print(f'Can\'t instantiate Photon args: {e}')
        func_kwargs = []
    # Analysis the function
    photon = photonFunc(**func_kwargs)
    det = detector
    ds, run = datasource
    det.addFunc(photon)
                
    userDict = {}
    small_data = ds.small_data('./test_photon.h5', gather_interval=5) # file to save data to
    ds.break_after(5) # stop event loop after 5 events
    
    # dats = []

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        # dats.append(det.evt.dat)
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_photon.h5"
    h5explorer = tables.File(f).root
    
    #checking for the output
    assert(list(h5explorer.epix_alc3._v_children)==['photon_nPhot'])
    logger.info('correct name and number of outputs')
    
    #checking for the shape
    
    assert(h5explorer.epix_alc3.photon_nPhot.shape == (5,))
    logger.info('Correct data shape for photon_nPhot')
    
    #checking for the type
    for i in range(5):        
        assert(h5explorer.epix_alc3.photon_nPhot[i].dtype == 'uint64')
    logger.info('Correct data type for photon_nPhot')
     
    logger.info('Pass the photon_function test')
    
    tables.file._open_files.close_all()
    
    os.remove('test_photon.h5')

    
@pytest.mark.parametrize('getDroplets', [{'theshold': 5,'thresholdLow':5, 'thresADU': 60, 'useRms':True, 'nData': 1e5}], indirect=True) 
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 631}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'epix_alc3'}], indirect=True)
def test_Droplets(datasource, detector, getDroplets):
    logger.info('Running test for Photons function to check output format.')
    try:
        func_kwargs = getDroplets
    except Exception as e:
        print(f'Can\'t instantiate Photon args: {e}')
        func_kwargs = []
    # Analysis the function
    print(func_kwargs)
    droplet = dropletFunc(**func_kwargs)
    det = detector
    ds, run = datasource
    det.addFunc(droplet)
                
    userDict = {}
    small_data = ds.small_data('./test_droplet.h5', gather_interval=5) # file to save data to
    ds.break_after(5) # stop event loop after 5 events

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        # dats.append(det.evt.dat)
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_droplet.h5"
    h5explorer = tables.File(f).root
    
    #checking for the output
    assert(list(h5explorer.epix_alc3._v_children)==['droplet_nDroplets', 'droplet_nDroplets_all'])
    logger.info('correct name and number of outputs')
    
    #checking for the shape
    
    assert(h5explorer.epix_alc3.droplet_nDroplets.shape == (5,))
    assert(h5explorer.epix_alc3.droplet_nDroplets_all.shape == (5,))
    logger.info('Correct data shape for droplet_nDroplets_all and droplet_nDroplets')
    
    #checking for the type
    for i in range(5):        
        assert(h5explorer.epix_alc3.droplet_nDroplets[i].dtype == 'int64')
        assert(h5explorer.epix_alc3.droplet_nDroplets_all[i].dtype == 'int64')
        
    logger.info('Correct data type for droplet_nDroplets, droplet_nDroplets_all')
     
    logger.info('Pass the Droplet_function test')
    
    tables.file._open_files.close_all()
    
    os.remove('test_droplet.h5')
    

@pytest.mark.parametrize('getD2P', [{'aduspphot': 162, 'cputime': True}], indirect=True) 
@pytest.mark.parametrize('getDroplets', [{'threshold': 10, 'thresholdLow':3, 'thresADU':10, 'useRms':True}], indirect=True) 
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 631}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'epix_alc1'}], indirect=True)
@pytest.mark.parametrize('getSparsify',[{'nData': 30000}], indirect=True)
def test_D2P(datasource, detector, getDroplets, getD2P, getSparsify):
    logger.info('Running test for D2P function to check output format.')
    try:
        func_kwargs = getD2P
    except Exception as e:
        print(f'Can\'t instantiate D2P args: {e}')
        func_kwargs = []
    d2p = droplet2Photons(**func_kwargs)
    
    try:
        func_kwargs = getDroplets
    except Exception as e:
        print(f'Can\'t instantiate Droplet args: {e}')
        func_kwargs = []
    droplet = dropletFunc(**func_kwargs)
    
    try:
        func_kwargs = getSparsify
    except Exception as e:
        print(f'Can\'t instantiate Droplet args: {e}')
        func_kwargs = []
    sparsify = sparsifyFunc(**func_kwargs)
    nData = func_kwargs['nData']

    d2p.addFunc(sparsify)
    droplet.addFunc(d2p)
    
    det = detector
    ds, run = datasource
    det.addFunc(droplet)
                
    userDict = {}
    small_data = ds.small_data('./test_d2p.h5', gather_interval=5) # file to save data to
    ds.break_after(5) # stop event loop after 5 events

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        # dats.append(det.evt.dat)
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_d2p.h5"
    h5explorer = tables.File(f).root
    
    #checking for the output
    if nData == None:
        assert(list(h5explorer.epix_alc1._v_children)==['var_droplet_droplet2phot_sparse', 
                                                        'droplet_droplet2phot_cputime', 
                                                        'droplet_droplet2phot_prob', 
                                                        'droplet_nDroplets', 
                                                        'droplet_nDroplets_all', 
                                                        'var_droplet_droplet2phot_sparse_len'])
        assert(list(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse._v_children)==['col', 'data', 'row', 'tile'])
    else:
        assert(list(h5explorer.epix_alc1._v_children)==['droplet_droplet2phot_cputime',
                                                        'droplet_droplet2phot_prob',
                                                        'droplet_droplet2phot_sparse_col',
                                                        'droplet_droplet2phot_sparse_data',
                                                        'droplet_droplet2phot_sparse_row',
                                                        'droplet_droplet2phot_sparse_tile',
                                                        'droplet_nDroplets',
                                                        'droplet_nDroplets_all'])
    logger.info('correct name and number of outputs')
    
    #checking for the shape
    if nData == None:
        assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse.col.shape == (99506,))
        assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse.data.shape == (99506,))
        assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse.row.shape == (99506,))
        assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse.tile.shape == (99506,))
        assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse_len.shape == (5,))
    else:
        assert(h5explorer.epix_alc1.droplet_droplet2phot_sparse_col.shape == (5, nData))
        assert(h5explorer.epix_alc1.droplet_droplet2phot_sparse_data.shape == (5, nData))
        assert(h5explorer.epix_alc1.droplet_droplet2phot_sparse_row.shape == (5, nData))
        assert(h5explorer.epix_alc1.droplet_droplet2phot_sparse_tile.shape == (5, nData))
    assert(h5explorer.epix_alc1.droplet_droplet2phot_cputime.shape == (5,4))
    assert(h5explorer.epix_alc1.droplet_droplet2phot_prob.shape == (5,12))
    assert(h5explorer.epix_alc1.droplet_nDroplets.shape == (5,))
    assert(h5explorer.epix_alc1.droplet_nDroplets_all.shape == (5,))
    
    logger.info('Correct data shape for all outputs')
    
    #checking for the type
    if nData==None:
        assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse.col.dtype == 'float64')
        assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse.data.dtype == 'float64')
        assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse.row.dtype == 'float64')
        assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse.tile.dtype == 'float64')
        for i in range(5):
            assert(h5explorer.epix_alc1.var_droplet_droplet2phot_sparse_len[i].dtype == 'int64') 
    else:
        for i in range(5):
            assert(h5explorer.epix_alc1.droplet_droplet2phot_sparse_col[i].dtype == 'int64')
            assert(h5explorer.epix_alc1.droplet_droplet2phot_sparse_data[i].dtype == 'float64')
            assert(h5explorer.epix_alc1.droplet_droplet2phot_sparse_row[i].dtype == 'int64')
            assert(h5explorer.epix_alc1.droplet_droplet2phot_sparse_tile[i].dtype == 'int64')
    for i in range(5): 
        assert(h5explorer.epix_alc1.droplet_droplet2phot_cputime[i].dtype == 'float64')
        assert(h5explorer.epix_alc1.droplet_droplet2phot_prob[i].dtype == 'float64')
        assert(h5explorer.epix_alc1.droplet_nDroplets[i].dtype == 'int64')
        assert(h5explorer.epix_alc1.droplet_nDroplets_all[i].dtype == 'int64')
    logger.info('Correct data type for all outputs')
     
    logger.info('Pass the d2p_function test')
    
    tables.file._open_files.close_all()
    
    os.remove('test_d2p.h5')