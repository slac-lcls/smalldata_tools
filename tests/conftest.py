import pytest
import psana

from smalldata_tools.DetObject import DetObject

@pytest.fixture(scope="module")
def datasource(request):
    exp = request.param.get('exp')
    run = request.param.get('run')
    ds_str = f"exp={exp}:run={run}"
    ds = psana.MPIDataSource(ds_str)
    yield ds

    ds.close()

@pytest.fixture(scope="module")
def detector(request):
    
    detector_name = request.param.get('name')
    detector = DetObject(detector_name)
    yield detector

    detector.cleanup()
