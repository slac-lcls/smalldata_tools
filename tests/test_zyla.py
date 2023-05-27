import logging
import sys
import psana
import pytest 
import smalldata_tools
from smalldata_tools.DetObject import DetObject

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

logger.info('Loading detector: zyla_0')

exp = 'xpptut15'
run = 630
ds_str = f"exp={exp}:run={run}"
ds = psana.MPIDataSource(ds_str)

det = DetObject('zyla_0', ds.env(), run)

def test_detector_type():
    logger.debug('Running detector type test')
    assert(isinstance(det, smalldata_tools.DetObject.ZylaObject))
    logger.debug('Pass the test')
