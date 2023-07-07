import logging

import sys
import psana
import pytest 

import smalldata_tools
from smalldata_tools.DetObject import DetObject
from tests.conftest import datasource, detector

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
    
