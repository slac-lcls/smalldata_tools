import logging

import sys
import psana
import pytest 

import smalldata_tools
from smalldata_tools.DetObject import DetObject
from tests.conftest import datasource, detector, getROIs
import tests.conftest
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, projectionFunc
from smalldata_tools.SmallDataUtils import getUserData

import tables

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

logger.info('Loading detector: jungfrau1M')


@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 650}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'jungfrau1M'}], indirect=True)
def test_detector_type(datasource, detector):
    logger.debug('Running detector type test')
    det = detector
    assert(isinstance(det, smalldata_tools.DetObject.JungfrauObject))
    logger.debug('Pass the Detector_type test')

    
# x1 = 0
# x2 = 1
# y1 = 100
# y2 = 120
# z1 = 50
# z2 = 60
# writeArea = True
    
@pytest.mark.parametrize('getROIs', [{'ROI': [[1,2], [157,487], [294,598]], 'writeArea': True, 'thresADU': None, }], indirect=True)
@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 650}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'jungfrau1M'}], indirect=True)
#@pytest.mark.parametrize('getROIs', [{'ROIs': [ [[x1,x2], [y1,y2], [z1,z2]] ], 'writeArea': writeArea, 'thresADU': None, }], indirect=True)
def test_ROI(getROIs, detector, datasource):
    logger.debug('Running test for ROI function')
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
    
    dats = []

    for nevt,evt in enumerate(ds.events()): # usual psana event loop
        det.getData(evt) # get the detector data
        det.processFuncs() # process the attached functions
        userDict[det._name]=getUserData(det) # get the function results
        dats.append(det.evt.dat)
        small_data.event(userDict) # write data to h5
    
    # Test the function
    f = "test_roi.h5"
    h5explorer = tables.File(f).root
    
    print(list(h5explorer.jungfrau1M._v_children))

    for i in range(5):
        if writeArea:
            assert(list(h5explorer.jungfrau1M._v_children)==['ROI_area', 'ROI_com', 'ROI_max', 'ROI_mean', 'ROI_sum'], True)
        else:
            assert(list(h5explorer.jungfrau1M._v_children)==['ROI_com', 'ROI_max', 'ROI_mean', 'ROI_sum'], True)
        logger.debug('Correct output')
        
        dat = dats[i]    
        area = h5explorer.jungfrau1M.ROI_area[i]
        com = h5explorer.jungfrau1M.ROI_com[i]
        max = h5explorer.jungfrau1M.ROI_max[i]
        mean = h5explorer.jungfrau1M.ROI_mean[i]
        sum = h5explorer.jungfrau1M.ROI_sum[i]
        
        #checking for the type
        assert(com.dtype == 'float32', True)
        assert(max.dtype == 'float32', True)
        assert(mean.dtype == 'float32', True)
        assert(sum.dtype == 'float32', True)
        logger.debug('Correct type')
        
        #checking for the size/shape
        assert(com.shape == (2,), True)
        assert(max.size == 1, True)
        assert(mean.size == 1, True)
        assert(sum.size == 1, True)
        #assert(area == (330,304), True)
        assert(area.shape == (y2-y1,z2-z1), True)
        logger.debug('Correct shape')
        
    logger.debug('Pass the ROI_function test')
    
    

        
        




        
        
        
        
        
        
        
        

        # assert(area[0][0] == dat[1][157][294], True)
        # assert(area[329][303] == dat[1][486][597], True)
#         assert(area[0][0] == dat[x1][y1][z1], True)
#         assert(area[0][z2-z1-1] == dat[x2-1][y1][z2-1], True)
#         assert(area[y2-y1-1][0] == dat[x2-1][y2-1][z1], True)
#         assert(area[y2-y1-1][z2-z1-1] == dat[x2-1][y2-1][z2-1], True)
#         logger.debug('Correct boundary')

#         assert(area.sum() == h5explorer.jungfrau1M.ROI_sum[4], True)
    
    
    
    
    
    
    
