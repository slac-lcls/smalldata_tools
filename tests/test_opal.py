import logging
import sys
import psana
import pytest 
import smalldata_tools
from smalldata_tools.DetObject import DetObject

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

logger.info('Loading detector: Opal_2')

exp = 'xpptut15'
run = 650
ds_str = f"exp={exp}:run={run}"
ds = psana.MPIDataSource(ds_str)

det = DetObject('Opal_2', ds.env(), run)

def test_detector_type():
    logger.debug('Running detector type test')
    assert(isinstance(det, smalldata_tools.DetObject.OpalObject))
    logger.debug('Pass the test')
