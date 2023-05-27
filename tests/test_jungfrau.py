import logging
import sys
import psana
import pytest 
import smalldata_tools
from smalldata_tools.DetObject import DetObject

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

logger.info('Loading detector: jungfrau1M')

exp = 'xpptut15'
run = 650
ds_str = f"exp={exp}:run={run}"
ds = psana.MPIDataSource(ds_str)

det = DetObject('jungfrau1M', ds.env(), run)

def test_detector_type():
    logger.debug('Running detector type test')
    assert(isinstance(det, smalldata_tools.DetObject.JungfrauObject))
    logger.debug('Pass the test')
