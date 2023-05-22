import sys
import psana
from smalldata_tools.DetObject import DetObject
import pytest
exp = 'xpptut15'
run = 650
ds_str = f"exp={exp}:run={run}"
ds = psana.MPIDataSource(ds_str)
sys.path.append('/sdf/home/h/hzhang12/smalldata_tools/')

det = DetObject('jungfrau1M', ds.env(), run)

assert(isinstance(det, smalldata_tools.DetObject.JungfrauObject))