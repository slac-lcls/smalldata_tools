import pytest
import psana

from smalldata_tools.DetObject import DetObject

@pytest.fixture(scope="module")
def datasource(request):
    exp = request.param.get('exp')
    run = request.param.get('run')
    ds_str = f"exp={exp}:run={run}"
    ds = psana.MPIDataSource(ds_str)
    return (ds, run)

@pytest.fixture(scope="module")
def detector(request, datasource):
    ds, run = datasource
    detector_name = request.param.get('name')
    detector = DetObject(detector_name, ds.env(), run)
    return detector


   