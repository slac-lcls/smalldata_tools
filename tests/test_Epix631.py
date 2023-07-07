import logging

import sys
import psana
import pytest 

import smalldata_tools
from smalldata_tools.DetObject import DetObject
from tests.conftest import datasource, detector

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

logger.info('Loading detector: Epix630')

@pytest.mark.parametrize('datasource', [{'exp': 'xpptut15', 'run': 631}], indirect=True)
@pytest.mark.parametrize('detector', [{'name': 'epix_alc3'}], indirect=True)
def test_detector_type(datasource, detector):
    logger.info('Running detector type test')
    det = detector
    assert(isinstance(det, smalldata_tools.DetObject.EpixObject))
    logger.info('Pass the test')