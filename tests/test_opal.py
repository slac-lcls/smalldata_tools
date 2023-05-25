import logging

import sys
import psana
import pytest

from ..smalldata_tools import DetObject

logger = logging.getLogger(__name__)

print('Loading detector: Opal650')
exp = 'xpptut15'
run = 650
ds_str = f"exp={exp}:run={run}"
ds = psana.MPIDataSource(ds_str)

det = DetObject('Opal_2', ds.env(), run)

@pytest
def test_detector_type():
    print('Running detector type test')
    assert(isinstance(det, DetObject.OpalObject))
    print('Pass the test')
