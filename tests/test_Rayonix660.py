import logging
import os

import sys
import psana
import pytest 

import smalldata_tools
from smalldata_tools.DetObject import DetObject
from tests.conftest import datasource, detector, getAzavPyfais
from smalldata_tools.ana_funcs.azav_pyfai import azav_pyfai
from smalldata_tools.SmallDataUtils import getUserData

import tables


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

logger.info('Loading detector: Rayonix660')

@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 660}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'Rayonix'}], indirect=True)
def test_detector_type(datasource, detector):
    logger.info('Running detector type test')
    det = detector
    assert(isinstance(det, smalldata_tools.DetObject.RayonixObject))
    logger.info('Pass the test')
    

pix_size = 1.760e-04
@pytest.mark.parametrize('getAzavPyfais', [{'userMask': None, 'thres': None, 'return2d': True, 'ai_kwargs': {'dist':1, 'poni1':960*pix_size, 'poni2':960*pix_size}, 'pol_factor': 1, 'npts_radial': 256,'npts_az': 360 }], indirect=True)
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 660}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'Rayonix'}], indirect=True)
def test_Azav_pyfai(getAzavPyfais, detector, datasource):
    logger.info('Running test for azimuthal integration to check output format.')
    try:
        func_kwargs = getAzavPyfais
    except Exception as e:
        print(f'Can\'t instantiate azimuthal integration args: {e}')
        func_kwargs = []
    return2d = func_kwargs['return2d']
    azav = azav_pyfai(**func_kwargs)
    det = detector
    ds, run = datasource
    det.addFunc(azav)
    
    userDict = {}
    small_data = ds.small_data('./test_azav.h5', gather_interval=5) # file to save data to
    ds.break_after(5) # stop event loop after 5 events

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_azav.h5"
    h5explorer = tables.File(f).root
    
    
    #checking for the output
    if return2d:
        assert(list(h5explorer.Rayonix._v_children)==['pyfai_az', 'pyfai_azav', 'pyfai_q'])
    else:
        assert(list(h5explorer.Rayonix._v_children)==['pyfai_azav', 'pyfai_q'])
    logger.info('Correct name and number of outputs')
    
    #checking for the shape
    assert(h5explorer.Rayonix.pyfai_azav.shape == (5, 360, 256))
    assert(h5explorer.Rayonix.pyfai_q.shape == (5, 256))
    if return2d:
        assert(h5explorer.Rayonix.pyfai_az.shape == (5, 360))
    logger.info('Correct data shape for pyfai_azav, pyfai_q, pyfai_az')
    
    #checking for the type
    for i in range(5):        
        assert(h5explorer.Rayonix.pyfai_azav[i].dtype == 'float32')
        assert(h5explorer.Rayonix.pyfai_q[i].dtype == 'float64')
        if return2d:
            assert(h5explorer.Rayonix.pyfai_az[i].dtype == 'float64')
    logger.info('Correct data type for pyfai_azav, pyfai_q, pyfai_az')
    
    logger.info('Pass the azimuthal integration test')
    
    tables.file._open_files.close_all()
    
    os.remove('test_azav.h5')