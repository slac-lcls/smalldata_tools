import sys
import psana
import pytest

sys.path.append('/sdf/home/h/hzhang12/smalldata_tools/')
import smalldata_tools 

@pytest.fixture
def test_opal_object(det):
    assert(isinstance(det, smalldata_tools.DetObject.OpalObject))
    
def test_opal_object(det):
    assert(isinstance(det, smalldata_tools.DetObject.EpixObject))
    
def test_jungfrau_object(det):
    assert(isinstance(det, smalldata_tools.DetObject.JungfrauObject))
    
def test_rayonix_object(det):
    assert(isinstance(det, smalldata_tools.DetObject.RayonixObject))

def test_zyla_object(det):
    assert(isinstance(det, smalldata_tools.DetObject.ZylaObject))