import pytest
import psana

from smalldata_tools.DetObject import DetObject


@pytest.fixture(scope='module')
def datasource(request):
    exp = request.param.get('exp')
    run = request.param.get('run')
    ds_str = f"exp={exp}:run={run}"
    ds = psana.MPIDataSource(ds_str)
    yield (ds, run)
    

@pytest.fixture(scope='module')
def detector(request, datasource):
    ds, run = datasource
    detector_name = request.param.get('name')
    detector = DetObject(detector_name, ds.env(), run)
    yield detector
    

@pytest.fixture(scope="function")
def getROIs(request):
    """ Set parameter for ROI analysis. Set writeArea to True to write the full ROI in the h5 file.
    See roi_rebin.py for more info
    """
    func_kwargs = {}
    func_kwargs['ROI'] = request.param.get('ROI')
    func_kwargs['writeArea'] = request.param.get('writeArea')
    func_kwargs['thresADU'] = request.param.get('thresADU')
    yield func_kwargs
    
@pytest.fixture(scope="function")  
def getProjections(request):
    """ Set parameter for Projection analysis.
    See roi_rebin.py for more info
    """
    roi_dict = {}
    roi_dict['name'] = 'proj'
    roi_dict['axis'] = request.param.get('axis')
    roi_dict['mean'] = request.param.get('mean')
    roi_dict['thresADU'] = request.param.get('thresADU')
    roi_dict['thresRms'] = request.param.get('thresRms')
    roi_dict['singlePhoton'] = request.param.get('singlePhoton')
    yield roi_dict

@pytest.fixture(scope="function")  
def getSpectrums(request):
    """ Set parameter for Spectrum analysis.
    See roi_rebin.py for more info
    """
    roi_dict = {}
    roi_dict['bins'] = request.param.get('bins')
    
    yield roi_dict

   
@pytest.fixture(scope="function")  
def getAzavPyfais(request):
    """ Set parameter for Spectrum analysis.
    See roi_rebin.py for more info
    """
    roi_dict = {}
    roi_dict['name'] = 'pyfai'
    roi_dict['userMask'] = request.param.get('userMask')
    roi_dict['return2d'] = request.param.get('return2d')
    roi_dict['ai_kwargs'] = request.param.get('ai_kwargs')
    roi_dict['pol_factor'] = request.param.get('pol_factor')
    roi_dict['npts_radial'] = request.param.get('npts_radial')
    roi_dict['npts_az'] = request.param.get('npts_az')
    
    yield roi_dict

