import sys
import psana
import pytest

exp = 'xpptut15'
run = 650
ds_str = f"exp={exp}:run={run}"
ds = psana.MPIDataSource(ds_str)
sys.path.append('/sdf/home/h/hzhang12/smalldata_tools/')
import smalldata_tools 
from smalldata_tools.DetObject import DetObjects

det = DetObject('Opal_2', ds.env(), run)

assert(isinstance(det, smalldata_tools.DetObject.OpalObject))