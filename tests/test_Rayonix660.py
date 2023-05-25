import logging

import sys
import psana
import pytest

from ..smalldata_tools import DetObject

logger = logging.getLogger(__name__)

print('Loading detector: Rayonix660')
exp = 'xpptut15'
run = 660
ds_str = f"exp={exp}:run={run}"
ds = psana.MPIDataSource(ds_str)
sys.path.append('/sdf/home/h/hzhang12/smalldata_tools/')
import smalldata_tools 
from smalldata_tools.DetObject import DetObjects
det = DetObject('Rayonix', ds.env(), run)

@pytest
def test_detector_type():
    print('Running detector type test')
    assert(isinstance(det, DetObject.RayonixObject))
    print('Pass the test')
