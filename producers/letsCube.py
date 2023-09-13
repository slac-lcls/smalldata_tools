#!/usr/bin/env python

import sys
import os
import time
from datetime import datetime
import numpy as np
import argparse
import socket
import logging 
import re
import requests
from requests.auth import HTTPBasicAuth
import tables
from pathlib import Path
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

fpath=os.path.dirname(os.path.abspath(__file__))
fpathup = '/'.join(fpath.split('/')[:-1])
sys.path.append(fpathup)
print(fpathup)
import smalldata_tools.cube.cube_mpi_fun as mpi_fun

##########################################################
# Custom exception handler to make job abort if a single rank fails.
# Avoid jobs hanging forever and report actual error message to the log file.
import traceback as tb

def global_except_hook(exctype, value, exc_traceback):
    tb.print_exception(exctype, value, exc_traceback)
    sys.stderr.write("except_hook. Calling MPI_Abort().\n")
    sys.stdout.flush() # Command to flush the output - stdout
    sys.stderr.flush() # Command to flush the output - stderr
    # Note: mpi4py must be imported inside exception handler, not globally.
    import mpi4py.MPI
    mpi4py.MPI.COMM_WORLD.Abort(1)
    sys.__excepthook__(exctype, value, exc_traceback)
    return

sys.excepthook = global_except_hook
##########################################################

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run', type=str, default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument("--indirectory", help="directory w/ smallData file if not default")
parser.add_argument("--outdirectory", help="directory w/ smallData for cube if not same as smallData", default='')
parser.add_argument("--nevents", help="number of events/bin", default=-1)
parser.add_argument("--postRuntable", help="postTrigger for seconday jobs", action='store_true')
parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
parser.add_argument('--config', help='Name of the config file module to use (without .py extension)', default=None, type=str)
parser.add_argument("--optimize_cores", help="split processing over more cores than bins", action='store_true')
args = parser.parse_args()
    
exp = args.experiment
run = args.run

if rank==0:
    from smalldata_tools.SmallDataAna import SmallDataAna
    from smalldata_tools.SmallDataAna_psana import SmallDataAna_psana
    import smalldata_tools.cube.cube_rank0_utils as utils
    
    # Constants
    HUTCHES = [
        'XPP',
        'XCS',
        'MFX',
        'CXI',
        'MEC'
    ]

    S3DF_BASE = Path('/sdf/data/lcls/ds/')
    FFB_BASE = Path('/cds/data/drpsrcf/')
    PSANA_BASE = Path('/cds/data/psdm/')
    PSDM_BASE = Path(os.environ.get('SIT_PSDM_DATA', S3DF_BASE))
    SD_EXT = Path('./hdf5/smalldata/')
    logger.info(f"PSDM_BASE={PSDM_BASE}")

    begin_prod_time = datetime.now().strftime('%m/%d/%Y %H:%M:%S')

    ##### START SCRIPT ########

    # Define hostname
    hostname = socket.gethostname()

    # Parse hutch name from experiment and check it's a valid hutch
    exp = args.experiment
    run = args.run
    logger.debug(f'Analyzing data for EXP:{args.experiment} - RUN:{args.run}')

    hutch = exp[:3].upper()
    if hutch not in HUTCHES:
        print(f'Could not find {hutch} in list of available hutches')
        comm.Abort()
    # load config file
    if args.config is None:
        if hutch=='XPP':
            import cube_config_xpp as config
        elif hutch=='XCS':
            import cube_config_xcs as config
        elif hutch=='MFX':
            import cube_config_mfx as config
        elif hutch=='CXI':
            import cube_config_cxi as config
    else:
        print(f'Importing custom config {args.config}')
        config = importlib.import_module(args.config)


    dirname=''
    if args.indirectory:
        dirname = args.indirectory
        if dirname.find('/')<=0:
            dirname+='/'
        
    ana = None
    anaps = SmallDataAna_psana(exp,run,dirname)
    if anaps and anaps.sda is not None and 'fh5' in anaps.sda.__dict__.keys():
        print('/nCreate ana module from anaps')
        ana = anaps.sda
    else:
        print('We will now try to open the smallData file directly')
        ana = SmallDataAna(exp, run, dirname, fname)
        if 'fh5' not in ana.__dict__.keys():
            ana = None
    if ana is None:
        print('No ana instance found. Abort')
        comm.Abort()

    ana.printRunInfo()

    for filt in config.filters:
        if rank==0: print('Filter: ',filt)
        ana.addCut(*filt)

    varList = config.varList

    ####
    # if you are interested in laser-xray delay, please select the delay of choice here!
    ####
    ana.setDelay(use_ttCorr=config.use_tt, addEnc=False, addLxt=False, reset=True)
    
    cubeName='cube' #initial name
    scanName, scanValues = ana.getScanValues()
    binSteps=[]
    binName=''
    
    if isinstance(scanName, list):
       scanName = scanName[0]
       scanValues = scanValues[0]

    if scanName!='':
        if scanName.find('delay')<0:
            scanSteps = np.unique(scanValues)
            scanSteps = np.append(scanSteps, scanSteps[-1]+abs(scanSteps[1]-scanSteps[0])) #catch value at right edge?
            scanSteps = np.append(scanSteps[0]-abs(scanSteps[1]-scanSteps[0]),scanSteps) #catch values at left edge
            binSteps = scanSteps
            cubeName = scanName
            if 'lxt' in scanName:
                bins = config.binBoundaries(run)
                #only use special binning if we use the timetool correction!
                if bins is not None and config.use_tt:
                    binSteps = bins
                print('bin data using ',scanName,' and bins: ',binSteps)
                binName='delay'
            else:
                print('bin data using ',scanName,' and bins: ',scanSteps)
                binName='scan/%s'%scanName
        else:
            binSteps = config.binBoundaries(run)
            cubeName = 'delay'
            if binSteps is None:
                #assume a fast delay scan here.
                cubeName = 'lasDelay'
                enc = ana.getVar('enc/lasDelay')
                binSteps = np.arange( (int(np.nanmin(enc*10)))/10., 
                                     (int(np.nanmax(enc*10)))/10., 0.1)
            binName = 'delay'
    else:
        cubeName='randomTest'
        binName = [key for key in ana.Keys() if re.search("ipm.*/sum",key)][0]
        # binName='ipm2/sum'
        binVar=ana.getVar(binName)
        binSteps=np.nanpercentile(binVar,[0,25,50,75,100])

    print(f'Bin name: {binName}, bins: {binSteps}')

    if int(args.nevents)>0: # not really implemented
        cubeName+='_%sEvents'%args.nevents
    
    cubeName_base = cubeName
    for filterName in ana.Sels:
        if filterName!='filter1':
            cubeName = cubeName_base+'_'+filterName
        else:
            cubeName = cubeName_base
        ana.addCube(cubeName, binName, binSteps, filterName)
        ana.addToCube(cubeName, varList)
            
    nBins = binSteps.shape[0]
    try:
        addBinVars = config.get_addBinVars(run)
    except Exception as e:
        print('Not setting additional binVar. Will assume that 1D cube is requested. Fix if this is not the case.')
        logger.info(f'addBinVars error info: {e}')
        addBinVars = None
    if addBinVars is not None and isinstance(addBinVars, dict):
        for cubeName, cube in ana.cubes.items():
            cube.add_BinVar(addBinVars)
        nBins *= len(addBinVars[1])
            
    #figure out how many bins & cores we have to maybe add a fake variable.
    nCores = size
    nSplit = int(nCores/nBins)
    #ok, we have more cores than bins. add a random variable to the data & add bins for that to the cube.
    if nSplit > 1 and args.optimize_cores:
        my_randomvar = np.random.random(ana.getVar('fiducials').shape[0])
        ana.addVar('random_remove',my_randomvar)
        randomBins = np.linspace(0,1.,nSplit+1)
        for cubeName, cube in ana.cubes.items():
            cube.add_BinVar({'random_remove': randomBins.tolist()})

    anaps._broadcast_xtc_dets(cubeName) # send detectors dict to workers. All cubes MUST use the same det list.
    
    cube_infos = []
    for ii, cubeName in enumerate(ana.cubes):
        print('Cubing {}'.format(cubeName))
        comm.bcast('Work!', root=0)
        time.sleep(1) # is this necessary? Just putting it here in case...
        if config.laser:
            #request 'on' events base on input filter (add optical laser filter)
            cubeName, bins, nEntries = anaps.makeCubeData(
                cubeName,
                onoff=1,
                nEvtsPerBin=args.nevents,
                dirname=args.outdirectory
            )
            cube_infos.append([f'{cubeName}_on', bins, nEntries])
            comm.bcast('Work!', root=0)
            time.sleep(1) # is this necessary? Just putting it here in case...
            
            #request 'off' events base on input filter (switch optical laser filter, drop tt
            cubeName, bins, nEntries = anaps.makeCubeData(
                cubeName,
                onoff=0, 
                nEvtsPerBin=args.nevents,
                dirname=args.outdirectory
            )
            cube_infos.append([f'{cubeName}_off', bins, nEntries])
        
        else:
            # no laser filters
            cubeName, bins, nEntries = anaps.makeCubeData(
                cubeName,
                onoff=2,
                nEvtsPerBin=args.nevents,
                dirname=args.outdirectory
            )
            cube_infos.append([cubeName, bins, nEntries])
    comm.bcast('Go home!', root=0)
            
    if nSplit > 1 and args.optimize_cores:
        cubedirectory= args.outdirectory
        cubeFNames = []
        for cubeName in ana.cubes:
            cubeFName = f'Cube_{args.experiment}_Run{int(run).zfill(4)}_{cubeName}.h5'
            if config.laser:
                cubeFNames.append(cubeFName.replace('.h5','_on.h5'))
                cubeFNames.append(cubeFName.replace('.h5','_off.h5'))
            else:
                cubeFNames.append(cubeFName)
        print('.....remove random_remove....')
        for cubeFName in cubeFNames:
            dat = tables.open_file('%s/%s'%(cubedirectory, cubeFName ),'r+')
            bins_remove = getattr(dat.root,'bins_random_remove', None)
            if bins_remove is None:
                print('Did not find the axis to be removed, quit.')
                continue
            nBins_remove = bins_remove.shape[0]
            print('Remove extra dimension: ', nSplit, nBins_remove)
            dat.remove_node('/bins_random_remove')
            #get all keys (not directories) & find binning variable random_remove
            h5_variables = [k for k in dir(dat.root) if k[0]!='_']
            #get data            
            for var in h5_variables:
                try:
                    org_data = getattr(dat.root,var).read()
                    if len(org_data.shape)<=1:
                        #print('Dim not available')
                        continue
                    if org_data.shape[1]!=nBins_remove:
                        #print('Did not find the right variable to sum over ', org_data.shape, nBins_remove)
                        continue
                    #sum to remove dim
                    new_data = org_data.sum(axis=1)
                    #remove dataset
                    dat.remove_node('/%s'%var)
                    #add dataset
                    dat.create_array('/',var, obj=new_data)
                except:
                    pass

    if config.save_tiff is True:
        tiffdir='/cds/data/drpsrcf/mec/meclu9418/scratch/run%s'%args.run
        cubeFName = 'Cube_%s_Run%04d_%s'%(args.experiment, int(run), cubeName )
        cubedirectory= args.outdirectory
        if cubedirectory=='':
            cubedirectory = '/cds/data/drpsrcf/mec/meclu9418/scratch/hdf5/smalldata'
        dat = tables.open_file('%s/%s.h5'%(cubedirectory, cubeFName )).root
        detnames = ['Epix10kaQuad2','epix100a_1_IXS','epix100a_2_XRTS']
        for detname in detnames:
            cubedata = np.squeeze(getattr(dat, detname+'_data'))
            if cubedata.ndim==2:
                im = Image.fromarray(cubedata)
                tiff_file = '%s/Run_%s_%s_filter_%s.tiff'%(tiffdir, run, detname, filterName)
                im.save(tiff_file)
 
    if os.environ.get('ARP_JOB_ID', None) is not None:
        requests.post(os.environ["JID_UPDATE_COUNTERS"], json=[{"key": "<b>Cube </b>", "value": "Done"}])
    
    # Make summary plots
    logger.info('###### MAKE REPORT #####')
    tabs = utils.make_report(anaps, cube_infos, config.hist_list, config.filters, config.varList, exp, run)
        
else:
    work = 1
    binWorker = mpi_fun.BinWorker(run, exp)
    while(work):
        time.sleep(1) # is this necessary? Just putting it here in case...
        amIStillWorking = comm.bcast(None, root=0)
        if amIStillWorking=='Go home!':
            work = 0
            logger.debug('Yay, work\'s over!')
        elif amIStillWorking=='Work!':
            logger.debug('Oh no, I\'m getting more work!')
            binWorker.work()
