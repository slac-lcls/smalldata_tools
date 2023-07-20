import logging
import os

import sys
import psana
import pytest 

import smalldata_tools
from smalldata_tools.DetObject import DetObject
from tests.conftest import datasource, detector, getROIs
import tests.conftest
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, projectionFunc, spectrumFunc
from smalldata_tools.SmallDataUtils import getUserData

import tables

logger = logging.getLogger(__name__)
#logger.setLevel(logging.DEBUG)

logger.info('Loading detector: jungfrau1M')


@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 650}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'jungfrau1M'}], indirect=True)
def test_detector_type(datasource, detector):
    logger.info('Running detector type test')
    det = detector
    assert(isinstance(det, smalldata_tools.DetObject.JungfrauObject))
    logger.info('Pass the Detector_type test')

    
@pytest.mark.parametrize('getROIs', [{'ROI': [[1,2], [157,487], [294,598]], 'writeArea': True, 'thresADU': None, }], indirect=True)
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 650}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'jungfrau1M'}], indirect=True)
def test_ROI(getROIs, detector, datasource):
    logger.info('Running test for ROI function to check output format.')
    try:
        func_kwargs = getROIs
    except Exception as e:
        print(f'Can\'t instantiate ROI args: {e}')
        func_kwargs = []
    # Analysis the function
    [[x1, x2], [y1, y2], [z1, z2]] = func_kwargs['ROI']
    writeArea = func_kwargs['writeArea']
    roi = ROIFunc(**func_kwargs)
    det = detector
    ds, run = datasource
    det.addFunc(roi)
    
                
    userDict = {}
    small_data = ds.small_data('./test_roi.h5', gather_interval=5) # file to save data to
    ds.break_after(5) # stop event loop after 5 events
    
    # dats = []

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        # dats.append(det.evt.dat)
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_roi.h5"
    h5explorer = tables.File(f).root

    for i in range(5):
        if writeArea:
            assert(list(h5explorer.jungfrau1M._v_children)==['ROI_area', 'ROI_com', 'ROI_max', 'ROI_mean', 'ROI_sum'])
        else:
            assert(list(h5explorer.jungfrau1M._v_children)==['ROI_com', 'ROI_max', 'ROI_mean', 'ROI_sum'])
        logger.info('Correct name and number of outputs')
        
        # dat = dats[i]    
        area = h5explorer.jungfrau1M.ROI_area[i]
        com = h5explorer.jungfrau1M.ROI_com[i]
        max = h5explorer.jungfrau1M.ROI_max[i]
        mean = h5explorer.jungfrau1M.ROI_mean[i]
        sum = h5explorer.jungfrau1M.ROI_sum[i]
        
        #checking for the type
        assert(com.dtype == 'float64')
        assert(max.dtype == 'float32')
        assert(mean.dtype == 'float32')
        assert(sum.dtype == 'float32')
        if writeArea:
            assert(area.dtype == 'float32')
        logger.info('Correct data type for ROI_area, ROI_com, ROI_max, ROI_mean, ROI_sum')
        
        #checking for the size/shape
        assert(com.shape == (2,))
        assert(max.size == 1)
        assert(mean.size == 1)
        assert(sum.size == 1)
        #assert(area == (330,304), True)
        if writeArea:
            assert(area.shape == (y2-y1,z2-z1))
        logger.info('Correct shape for ROI_area, ROI_com, ROI_max, ROI_mean, ROI_sum')

    logger.info('Pass the ROI_function test')
    
    tables.file._open_files.close_all()
    
    os.remove('test_roi.h5')



@pytest.mark.parametrize('getProjections', [{'axis': 2, 'mean': True, 'thresADU': 1e-6, 'thresRms': 1e-6, 'singlePhoton': False}], indirect=True)
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 650}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'jungfrau1M'}], indirect=True)
def test_Projection(getProjections, detector, datasource):
    logger.info('Running test for Projection function to check output format.')
    try:
        func_kwargs = getProjections
    except Exception as e:
        print(f'Can\'t instantiate Projection args: {e}')
        func_kwargs = []
    axis = func_kwargs['axis']
    proj = projectionFunc(**func_kwargs)
    det = detector
    ds, run = datasource
    det.addFunc(proj)
    
    userDict = {}
    small_data = ds.small_data('./test_proj.h5', gather_interval=5) # file to save data to
    ds.break_after(5) # stop event loop after 5 events

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_proj.h5"
    h5explorer = tables.File(f).root
    
    #checking for the output
    assert(list(h5explorer.jungfrau1M._v_children)==['proj_data'])
    logger.info('Correct name and number of outputs')
    
    #checking for the shape
    if axis==0:
        assert(h5explorer.jungfrau1M.proj_data.shape == (5, 512, 1024))
    elif axis==1:
        assert(h5explorer.jungfrau1M.proj_data.shape == (5, 2, 1024))
    elif axis==2:
        assert(h5explorer.jungfrau1M.proj_data.shape == (5, 2, 512))
    else:
        assert(h5explorer.jungfrau1M.proj_data.shape == (5, ))
    logger.info('Correct data shape for proj_data')
    
    #checking for the type
    for i in range(5):        
        assert(h5explorer.jungfrau1M.proj_data[i].dtype == 'float32')
    logger.info('Correct data type for proj_data')
    
    logger.info('Pass the Projection_function test')
    
    tables.file._open_files.close_all()
    
    os.remove('test_proj.h5')


@pytest.mark.parametrize('getROIs', [{'ROI': [[1,2], [157,487], [294,598]], 'writeArea': True, 'thresADU': None, }], indirect=True)
@pytest.mark.parametrize('getProjections', [{'axis': 1, 'mean': True, 'thresADU': 1e-6, 'thresRms': 1e-6, 'singlePhoton': False}], indirect=True)
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 650}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'jungfrau1M'}], indirect=True)
def test_ROI_Projection(getROIs, getProjections, detector, datasource):
    logger.info('Running test for ROI_Projection to check output format.')
    #roi function
    try:
        func_kwargs = getROIs
    except Exception as e:
        print(f'Can\'t instantiate ROI args: {e}')
        func_kwargs = []
    writeArea = func_kwargs['writeArea']
    [[x1, x2], [y1, y2], [z1, z2]] = func_kwargs['ROI']
    func = ROIFunc(**func_kwargs)
    #project function
    try:
        func_kwargs = getProjections
    except Exception as e:
        print(f'Can\'t instantiate Projection args: {e}')
        func_kwargs = []
    axis = func_kwargs['axis']
    proj = projectionFunc(**func_kwargs)
    
    if axis >= 2:
        print('axis out of bound')
        logger.info('Test failed because of wrong axis')
        return
    
    func.addFunc(proj)
    
    det = detector
    ds, run = datasource
    det.addFunc(func)
    
    # Process a couple events
    userDict = {} # dictionary to store detector data
    small_data = ds.small_data('./test_roi_proj.h5', gather_interval=5) # file to save data to
    ds.break_after(5) # stop event loop after 5 events

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_roi_proj.h5"
    h5explorer = tables.File(f).root
    
    #checking for the output
    if writeArea:
        assert(list(h5explorer.jungfrau1M._v_children)==['ROI_area', 'ROI_com', 'ROI_max', 'ROI_mean', 'ROI_proj_data', 'ROI_sum'])
    else:
        assert(list(h5explorer.jungfrau1M._v_children)==['ROI_com', 'ROI_max', 'ROI_mean', 'ROI_proj_data', 'ROI_sum'])
    logger.info('correct name and number of outputs')
    
    #checking for the shape
    if axis==0:
        assert(h5explorer.jungfrau1M.ROI_proj_data.shape == (5, x2-x1))
    elif axis==1:
        assert(h5explorer.jungfrau1M.ROI_proj_data.shape == (5, y2-y1))
    elif axis==2:
        assert(h5explorer.jungfrau1M.ROI_proj_data.shape == (5, z2-z1))
    else:
        assert(h5explorer.jungfrau1M.ROI_proj_data.shape == (5, ))
    assert(h5explorer.jungfrau1M.ROI_com.shape == (5,2))
    if writeArea:
        assert(h5explorer.jungfrau1M.ROI_area.shape == (5,y2-y1,z2-z1))
    assert(h5explorer.jungfrau1M.ROI_max.shape == (5,))
    assert(h5explorer.jungfrau1M.ROI_mean.shape == (5,))
    assert(h5explorer.jungfrau1M.ROI_sum.shape == (5,))
    logger.info('Correct data shape for ROI_area, ROI_com, ROI_max, ROI_mean, ROI_sum and ROI_proj_data')
    
    #checking for the type
    for i in range(5):        
        assert(h5explorer.jungfrau1M.ROI_proj_data[i].dtype == 'float32')
        assert(h5explorer.jungfrau1M.ROI_com[i].dtype == 'float64')
        assert(h5explorer.jungfrau1M.ROI_max[i].dtype == 'float32')
        assert(h5explorer.jungfrau1M.ROI_mean[i].dtype == 'float32')
        assert(h5explorer.jungfrau1M.ROI_sum[i].dtype == 'float32')
        if writeArea:
            assert(h5explorer.jungfrau1M.ROI_area[i].dtype == 'float32')
    logger.info('Correct data type for ROI_area, ROI_com, ROI_max, ROI_mean, ROI_sum and ROI_proj_data')
     
    logger.info('Pass the ROI_Projection_function test')
    
    tables.file._open_files.close_all()
    
    os.remove('test_roi_proj.h5')
    

@pytest.mark.parametrize('getSpectrums', [{'bins': [1,2,3,4,5,6,7,8]}], indirect=True)
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 650}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'jungfrau1M'}], indirect=True)
def test_Spectrum(getSpectrums, detector, datasource):
    logger.info('Running test for Spectrum function to check output format.')
    try:
        func_kwargs = getSpectrums
    except Exception as e:
        print(f'Can\'t instantiate Projection args: {e}')
        func_kwargs = []
    bins = len(func_kwargs['bins'])
    spec = spectrumFunc(**func_kwargs)
    det = detector
    ds, run = datasource
    det.addFunc(spec)
    
    userDict = {}
    small_data = ds.small_data('./test_spec.h5', gather_interval=5) # file to save data to
    ds.break_after(5) # stop event loop after 5 events

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_spec.h5"
    h5explorer = tables.File(f).root
    
    #checking for the output
    assert(list(h5explorer.jungfrau1M._v_children)==['spectrum_histogram'])
    logger.info('name and number of outputs')
    
    #checking for the shape
    assert(h5explorer.jungfrau1M.spectrum_histogram.shape==(5, bins-1))
    logger.info('Correct shape for spectrum_histogram')
    
    #checking for the type
    for i in range(5):        
        assert(h5explorer.jungfrau1M.spectrum_histogram[i].dtype == 'int64')
    logger.info('Correct data type for spectrum_histogram')
    
    logger.info('Pass the Spectrum_function test')
    
    tables.file._open_files.close_all()
    
    os.remove('test_spec.h5')
    
@pytest.mark.parametrize('getROIs', [{'ROI': [[1,2], [157,487], [294,598]], 'writeArea': True, 'thresADU': None, }], indirect=True)
@pytest.mark.parametrize('getSpectrums', [{'bins': [1,2,3,4,5,6,7,8]}], indirect=True)
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 650}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'jungfrau1M'}], indirect=True)
def test_ROI_Spectrum(getROIs, getSpectrums, detector, datasource):
    logger.info('Running test for ROI_Spectrum function to check output format.')
    #roi function
    try:
        func_kwargs = getROIs
    except Exception as e:
        print(f'Can\'t instantiate ROI args: {e}')
        func_kwargs = []
    writeArea = func_kwargs['writeArea']
    [[x1, x2], [y1, y2], [z1, z2]] = func_kwargs['ROI']
    func = ROIFunc(**func_kwargs)
    #spectrum function
    try:
        func_kwargs = getSpectrums
    except Exception as e:
        print(f'Can\'t instantiate Projection args: {e}')
        func_kwargs = []
    bins = len(func_kwargs['bins'])
    spec = spectrumFunc(**func_kwargs)
    
    func.addFunc(spec)
    
    det = detector
    ds, run = datasource
    det.addFunc(func)
    
    # Process a couple events
    userDict = {} # dictionary to store detector data
    small_data = ds.small_data('./test_roi_spec.h5', gather_interval=5) # file to save data to
    ds.break_after(5) # stop event loop after 5 events

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_roi_spec.h5"
    h5explorer = tables.File(f).root
    
    #checking for the output
    if writeArea:
        assert(list(h5explorer.jungfrau1M._v_children)==['ROI_area', 'ROI_com', 'ROI_max', 'ROI_mean', 'ROI_spectrum_histogram', 'ROI_sum'])
    else:
        assert(list(h5explorer.jungfrau1M._v_children)==['ROI_com', 'ROI_max', 'ROI_mean', 'ROI_sum', 'ROI_spectrum_histogram'])
    logger.info('correct name and number of outputs')
    
    #checking for the shape
    assert(h5explorer.jungfrau1M.ROI_com.shape == (5,2))
    if writeArea:
        assert(h5explorer.jungfrau1M.ROI_area.shape == (5,y2-y1,z2-z1))
    assert(h5explorer.jungfrau1M.ROI_max.shape == (5,))
    assert(h5explorer.jungfrau1M.ROI_mean.shape == (5,))
    assert(h5explorer.jungfrau1M.ROI_sum.shape == (5,))
    assert(h5explorer.jungfrau1M.ROI_spectrum_histogram.shape==(5, bins-1))
    logger.info('Correct data shape for ROI_area, ROI_com, ROI_max, ROI_mean, ROI_sum and ROI_spectrum_histogram_data')
    
    #checking for the type
    for i in range(5):        
        assert(h5explorer.jungfrau1M.ROI_spectrum_histogram[i].dtype == 'int64')
        assert(h5explorer.jungfrau1M.ROI_com[i].dtype == 'float64')
        assert(h5explorer.jungfrau1M.ROI_max[i].dtype == 'float32')
        assert(h5explorer.jungfrau1M.ROI_mean[i].dtype == 'float32')
        assert(h5explorer.jungfrau1M.ROI_sum[i].dtype == 'float32')
        if writeArea:
            assert(h5explorer.jungfrau1M.ROI_area[i].dtype == 'float32')
    logger.info('Correct data type for ROI_area, ROI_com, ROI_max, ROI_mean, ROI_sum and ROI_spectrum_histogram_data')
     
    logger.info('Pass the ROI_Spectrum_function test')
    
    tables.file._open_files.close_all()
    
    os.remove('test_roi_spec.h5')
    