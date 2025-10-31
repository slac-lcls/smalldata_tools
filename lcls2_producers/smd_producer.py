#!/usr/bin/env python

import numpy as np
import json
import psana
import time
from datetime import datetime

begin_job_time = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
start_job = time.time()
import argparse
import socket
import os
import logging
import requests
import sys
from glob import glob
from requests.auth import HTTPBasicAuth
from pathlib import Path
from importlib import import_module
from mpi4py import MPI

COMM = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

logger = logging.getLogger(__name__)
# logging.basicConfig(level=logging.DEBUG)
log_level = "INFO"
# log_level = DEBUG
log_format = "[ %(asctime)s | %(levelname)-3s | %(filename)s] %(message)s"
logging.basicConfig(format=log_format)
logger.setLevel(
    logging.INFO
)  # Set level here instead of in the basic config so other loggers are  not affected

if rank == 0:
    logger.info(f"MPI size: {size}")
    logger.info(
        "psana conda environment is {0}".format(os.environ["CONDA_DEFAULT_ENV"])
    )


# DEFINE DETECTOR AND ADD ANALYSIS FUNCTIONS
def define_dets(run, det_list):
    # Load DetObjectFunc parameters (if defined)
    # Assumes that the config file with function parameters definition
    # has been imported under "config"

    rois_args = []  # special implementation as a list to support multiple ROIs.
    dimgs_args = {}
    wfs_int_args = {}
    wfs_hitfinder_args = {}
    d2p_args = {}
    wfs_svd_args = {}
    droplet_args = {}
    azav_args = {}
    polynomial_args = {}
    sum_algo_args = {}

    # Get the functions arguments from the production config
    if "getROIs" in dir(config):
        rois_args = config.getROIs(run)
    if "getDetImages" in dir(config):
        dimgs_args = config.getDetImages(run)
    if "get_wf_integrate" in dir(config):
        wfs_int_args = config.get_wf_integrate(run)
    if "get_wf_hitfinder" in dir(config):
        wfs_hitfinder_args = config.get_wf_hitfinder(run)
    if "get_droplet2photon" in dir(config):
        d2p_args = config.get_droplet2photon(run)
    if "get_wf_svd" in dir(config):
        wfs_svd_args = config.get_wf_svd(run)
    if "get_droplet" in dir(config):
        droplet_args = config.get_droplet(run)
    if "get_azav" in dir(config):
        azav_args = config.get_azav(run)
    if "get_polynomial_correction" in dir(config):
        polynomial_args = config.get_polynomial_correction(run)
    if "get_sum_algos" in dir(config):
        sum_algo_args = config.get_sum_algos(run)

    dets = []

    for detname in det_list:
        havedet = detname in thisrun.detnames

        # Common mode (default: None)
        common_mode = None

        if not havedet:
            continue

        if detname.find("fim") >= 0 or detname.find("w8") >= 0:
            # why different name for w8? To differentiate full det from BLD
            det = DetObject(
                detname, thisrun, common_mode=common_mode, name=f"det_{detname}"
            )
        else:
            det = DetObject(detname, thisrun, common_mode=common_mode)
        logger.debug(f"Instantiated det {detname}: {det}")

        # HSD need special treatment due to their data structure.
        # TO REVIEW / REVISE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if detname.find("hsd") >= 0:  # and not args.nohsd:
            hsdsplit = hsdsplitFunc(writeHsd=False)
            if detname in rois_args:
                for sdetname in rois_args[detname]:
                    funcname = "%s__%s" % (sdetname, "ROI")
                    if rank == 0:
                        print(
                            f"sdetname {sdetname} args: {rois_args[detname][sdetname]}, func name: {funcname}"
                        )
                    RF = hsdROIFunc(
                        name=funcname,
                        writeArea=True,
                        ROI=rois_args[detname][sdetname],
                    )
                    hsdsplit.addFunc(RF)
            det.addFunc(hsdsplit)
            dets.append(det)
            continue

        ####################################
        ######## Standard detectors ########
        ####################################
        if detname in azav_args:
            det.addFunc(azimuthalBinning(**azav_args[detname]))

        if detname in rois_args:
            # ROI extraction
            for iROI, ROI in enumerate(rois_args[detname]):
                proj_ax = ROI.pop("proj_ax", None)

                thisROIFunc = ROIFunc(**ROI)
                if proj_ax is not None:
                    thisROIFunc.addFunc(projectionFunc(axis=proj_ax))
                det.addFunc(thisROIFunc)

        if detname in d2p_args:
            if rank == 0:
                print(f"\n\nSETTING UP DROPLET TO PHOTON FOR {det._name}")
                print(json.dumps(d2p_args[detname], indent=4))
                print("\n")
            if "nData" in d2p_args[detname]:  # defines how sparsifcation is saved
                nData = d2p_args[detname].pop("nData")
            else:
                nData = None
            if "get_photon_img" in d2p_args[detname]:  # save unsparsified image or not
                unsparsify = d2p_args[detname].pop("get_photon_img")
            else:
                unsparsify = False

            # Make individual analysis functions
            # (i) droplet
            droplet_dict = d2p_args[detname]["droplet"]
            dropfunc = dropletFunc(**droplet_dict)
            # (ii) droplet2Photon
            d2p_dict = d2p_args[detname]["d2p"]
            drop2phot_func = droplet2Photons(**d2p_dict)
            # (iii) sparsify
            sparsify = sparsifyFunc(nData=nData)

            # Assemble pipeline last to first
            if unsparsify:
                unsparsify_func = unsparsifyFunc()
                drop2phot_func.addFunc(unsparsify_func)
            drop2phot_func.addFunc(sparsify)
            dropfunc.addFunc(drop2phot_func)
            det.addFunc(dropfunc)

        # if detname in dimgs_args:
        #     # Detector image (det.raw.image())
        #     dimg_func = detImageFunc(**dimgs_args[detname])
        #     det.addFunc(dimg_func)

        if detname in wfs_int_args:
            # Waveform integration
            wfs_int_func = WfIntegration(**wfs_int_args[detname])
            det.addFunc(wfs_int_func)

        if detname in wfs_hitfinder_args:
            # Simple hit finder on waveform
            det.addFunc(SimpleHitFinder(**wfs_hitfinder_args[detname]))

        if detname in wfs_svd_args:
            if isinstance(
                wfs_svd_args[detname], list
            ):  # handle mutliple channel on the same det
                for i_ch, ch in enumerate(wfs_svd_args[detname]):
                    # Check that the basis file exists, skip instantiation if not
                    if not os.path.isfile(ch["basis_file"]):
                        logger.info(
                            f"SVD basis file for {detname} ch{i_ch} cannot be found. Will not instantiate SvdFit for it."
                        )
                        continue
                    this_func = SvdFit(**ch)
                    det.addFunc(this_func)
            else:
                det.addFunc(SvdFit(**wfs_svd_args[detname]))

        if detname in droplet_args:
            if "nData" in droplet_args[detname].keys():
                nData = droplet_args[detname].pop("nData")
            else:
                nData = None
            dFunc = dropletFunc(**droplet_args[detname])
            dFunc.addFunc(sparsifyFunc(nData=nData))
            det.addFunc(dFunc)

        if detname in polynomial_args:
            obj = det
            if "droplet" in polynomial_args[detname]:
                # If polynomial correction is applied to droplet data, we need to
                # add the droplet function first.
                droplet_args = polynomial_args[detname].pop("droplet")
                droplet_func = dropletFunc(**droplet_args)
                unsparsify = unsparsifyFunc()
                droplet_func.addFunc(unsparsify)
                det.addFunc(droplet_func)
                obj = unsparsify

            if "projection" in polynomial_args[detname]:
                # If projection is requested, add it to the polynomial correction.
                proj_args = polynomial_args[detname].pop("projection")
                proj_func = projectionFunc(**proj_args)

            poly_corr_func = PolynomialCurveCorrection(**polynomial_args[detname])
            if "proj_func" in locals():
                print("Add projection")
                poly_corr_func.addFunc(proj_func)
            obj.addFunc(poly_corr_func)

        if sum_algo_args == {}:
            # Store calib by default even if algos. dictionary not defined
            det.storeSum(sumAlgo="calib")
        else:
            # Add `all` algorithms
            if "all" in sum_algo_args:
                for algo in sum_algo_args["all"]:
                    det.storeSum(sumAlgo=algo)
            if detname in sum_algo_args:
                for algo in sum_algo_args[detname]:
                    det.storeSum(sumAlgo=algo)

        logger.debug(f"Rank {rank} Add det {detname}: {det}")
        dets.append(det)

    return dets


# General Workflow
# This is meant for arp which means we will always have an exp and run

fpath = os.path.dirname(os.path.abspath(__file__))
fpathup = "/".join(fpath.split("/")[:-1])
sys.path.append(fpathup)
if rank == 0:
    logger.info(fpathup)

from smalldata_tools.utilities import printMsg, checkDet
from smalldata_tools.common.detector_base import detData, getUserData, getUserEnvData
from smalldata_tools.lcls2.default_detectors import (
    detOnceData,
    epicsDetector,
    genericDetector,
)
from smalldata_tools.lcls2.hutch_default import defaultDetectors
from smalldata_tools.lcls2.DetObject import DetObject

from smalldata_tools.ana_funcs.roi_rebin import (
    ROIFunc,
    spectrumFunc,
    projectionFunc,
    imageFunc,
)
from smalldata_tools.ana_funcs.sparsifyFunc import sparsifyFunc, unsparsifyFunc
from smalldata_tools.ana_funcs.waveformFunc import WfIntegration, SimpleHitFinder
from smalldata_tools.ana_funcs.waveformFunc import getCMPeakFunc, templateFitFunc
from smalldata_tools.ana_funcs.waveformFunc import (
    hsdsplitFunc,
    hsdBaselineCorrectFunc,
    hitFinderCFDFunc,
    hsdROIFunc,
)
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.droplet2Photons import droplet2Photons
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning
from smalldata_tools.ana_funcs.smd_svd import SvdFit
from smalldata_tools.ana_funcs.detector_corrections import PolynomialCurveCorrection

import psplot
from psmon import publish


# Constants
HUTCHES = ["TMO", "RIX", "UED", "MFX"]

S3DF_BASE = Path("/sdf/data/lcls/ds/")
FFB_BASE = Path("/cds/data/drpsrcf/")
PSDM_BASE = Path(os.environ.get("SIT_PSDM_DATA", S3DF_BASE))
SD_EXT = Path("./hdf5/smalldata/")
logger.debug(f"PSDM_BASE={PSDM_BASE}")

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument(
    "--run", help="run", type=str, default=os.environ.get("RUN_NUM", "")
)
parser.add_argument(
    "--experiment",
    help="experiment name",
    type=str,
    default=os.environ.get("EXPERIMENT", ""),
)
parser.add_argument("--stn", help="hutch station", type=int, default=0)
parser.add_argument("--nevents", help="number of events", type=int, default=0)
parser.add_argument(
    "--directory", help="directory for output files (def <exp>/hdf5/smalldata)"
)
parser.add_argument("--gather_interval", help="gather interval", type=int, default=100)
parser.add_argument(
    "--norecorder", help="ignore recorder streams", action="store_true", default=False
)
parser.add_argument("--url", default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
parser.add_argument(
    "--epicsAll", help="store all epics PVs", action="store_true", default=False
)
parser.add_argument(
    "--full",
    help="store all data (please think before usig this)",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--default", help="store only minimal data", action="store_true", default=False
)
parser.add_argument(
    "--image",
    help="save everything as image (use with care)",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--tiff",
    help="save all images also as single tiff (use with even more care)",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--postRuntable",
    help="postTrigger for seconday jobs",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--wait", help="wait for a file to appear", action="store_true", default=False
)
parser.add_argument(
    "--rawFim", help="save raw Fim data", action="store_true", default=False
)
parser.add_argument(
    "--nohsd", help="dont save HSD data", action="store_true", default=False
)
parser.add_argument(
    "--nosum", help="dont save sums", action="store_true", default=False
)
parser.add_argument(
    "--noarch", help="dont use archiver data", action="store_true", default=False
)
parser.add_argument(
    "--psplot_live_mode",
    help="Run as a server for psplot live mode, i.e. no h5 file being written.",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--intg_delta_t",
    help="Offset for the integrating detector batch.",
    type=int,
    default=0,
)
parser.add_argument(
    "--config", help="Producer config file to use", default=None, type=str
)
parser.add_argument(
    "--all_events",
    help="Will write all events. If False, only writes h5 at slow detectors rate.",
    action="store_true",
    default=False,
)

args = parser.parse_args()

logger.debug("Args to be used for small data run: {0}".format(args))


###### Helper Functions ##########
def get_xtc_files(base, exp, run):
    """File all xtc files for given experiment and run"""
    run_format = "".join(["r", run.zfill(4)])
    data_dir = Path(base) / exp[:3] / exp / "xtc"
    xtc_files = list(data_dir.glob(f"*{run_format}*"))
    if rank == 0:
        logger.debug(f"xtc file list: {xtc_files}")
    return xtc_files


def get_sd_file(write_dir, exp, hutch):
    """Generate directory to write to, create file name"""
    if write_dir is None:
        if onS3DF:  # S3DF should now be the default
            write_dir = S3DF_BASE / hutch.lower() / exp / SD_EXT
        else:
            logger.error("On an unknown system, cannot figure where to save data.")
            logger.error("Please fix or pass a write_dir argument.")
            sys.exit()
    logger.debug(f"hdf5 directory: {write_dir}")

    write_dir = Path(write_dir)
    h5_f_name = write_dir / f"{exp}_Run{run.zfill(4)}.h5"
    if not write_dir.exists():
        if rank == 0:
            logger.info(f"{write_dir} does not exist, creating directory now.")
            try:
                write_dir.mkdir(parents=True)
            except (PermissionError, FileNotFoundError) as e:
                logger.error(
                    f"Unable to make directory {write_dir} for output"
                    f"exiting on error: {e}"
                )
                sys.exit()
    if rank == 0 and not args.psplot_live_mode:
        logger.info("Will write small data file to {0}".format(h5_f_name))
    elif rank == 0 and args.psplot_live_mode:
        logger.warning("Running in psplot_live mode, will not write any h5 file.")
    return h5_f_name


##### START SCRIPT ########
hostname = socket.gethostname()

# Parse hutch name from experiment and check it's a valid hutch
exp = args.experiment
run = args.run
station = args.stn
logger.debug("Analyzing data for EXP:{0} - RUN:{1}".format(args.experiment, args.run))

begin_prod_time = datetime.now().strftime("%m/%d/%Y %H:%M:%S")

hutch = exp[:3].upper()
if hutch not in HUTCHES:
    logger.error("Could not find {0} in list of available hutches".format(hutch))
    sys.exit()


# Get config file
def get_config_file(name, folder_path=Path(fpathup) / "lcls2_producers"):
    """Find config file using pathlib"""
    folder = Path(folder_path)
    target_file = folder / f"prod_config_{name}.py"

    if target_file.exists():
        return target_file.stem  # return the file name without extension
    else:
        if rank == 0:
            logger.error(f"Config file {target_file} not found.")
        sys.exit(1)


if args.config is None:
    prod_cfg = f"prod_config_{hutch.lower()}"
else:
    prod_cfg = get_config_file(args.config)
if rank == 0:
    logger.info(f"Producer cfg file: <{prod_cfg}>.")
config = import_module(prod_cfg)


# Figure out where we are running from and check that the
# xtc files are where we expect them.
onS3DF = False
useFFB = False
xtc_files = []

if hostname.find("sdf") >= 0:
    logger.debug("On S3DF")
    onS3DF = True
    if "ffb" in PSDM_BASE.as_posix():
        useFFB = True
        # wait for files to appear
        nFiles = 0
        n_wait = 0
        max_wait = 20  # 10s wait per cycle.
        waitFilesStart = datetime.now()
        while nFiles == 0:
            if n_wait > max_wait:
                raise RuntimeError(
                    "Waited {str(n_wait*10)}s, still no files available. Giving up."
                )
            xtc_files = get_xtc_files(PSDM_BASE, exp, run)
            nFiles = len(xtc_files)
            if nFiles == 0:
                if rank == 0:
                    print(
                        f"We have no xtc files for run {run} in {exp} in the FFB system, "
                        "we will wait for 10 second and check again."
                    )
                n_wait += 1
                time.sleep(10)
        waitFilesEnd = datetime.now()
        if rank == 0:
            print(f"Files appeared after {str(waitFilesEnd-waitFilesStart)} seconds")

    xtc_files = get_xtc_files(PSDM_BASE, exp, run)
    if len(xtc_files) == 0:
        raise RuntimeError(
            f"We have no xtc files for run {run} in {exp} in the offline system."
        )
else:
    logger.warning("On an unknow system, things may get weird.")

# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.directory, exp, hutch)

# Create data source.
datasource_args = {"exp": exp, "run": int(run)}
if useFFB:
    datasource_args["live"] = True

if rank == 0:
    logger.info("Opening the data source:")
if args.nevents != 0:
    datasource_args["max_events"] = args.nevents

# Setup if integrating detectors are requested.
if hasattr(config, "get_intg"):
    intg_main, intg_addl = config.get_intg(run)
    integrating_detectors = []
    skip_intg = False
else:
    intg_main, intg_addl = (None, None)
    skip_intg = True

if intg_main is not None and intg_main != "":
    ds = psana.DataSource(**datasource_args)
    thisrun = next(ds.runs())
    if not isinstance(thisrun, psana.psexp.null_ds.NullRun):
        detnames = thisrun.detnames
        if intg_main not in detnames:
            skip_intg = True  # skip integrating detector setup
            if rank == 0:
                logger.error(
                    f"Main integrating detector {intg_main} not found in the data."
                )
                logger.error("Skipping integrating detectors.")
                logger.error(
                    "Please check the integrating detector list in the config."
                )
        else:
            datasource_args["intg_det"] = intg_main
            datasource_args["intg_delta_t"] = args.intg_delta_t
            datasource_args["batch_size"] = 1
            os.environ["PS_SMD_N_EVENTS"] = (
                "1"  # must be 1 for any non-zero value of delta_t
            )
            integrating_detectors = [
                intg_main
            ] + intg_addl  # for the detector instantiation
    ds = None
    thisrun = None

if args.psplot_live_mode:
    datasource_args["psmon_publish"] = publish

ds = psana.DataSource(**datasource_args)

if ds.unique_user_rank():
    print("#### DATASOURCE AND PSANA ENV VAR INFO ####")
    print(f"Instantiated data source with arguments: {datasource_args}")
    print(f"MPI size: {size}")
    print(f"PS_EB_NODES={os.environ.get('PS_EB_NODES')}")
    print(f"PS_SRV_NODES={os.environ.get('PS_SRV_NODES')}")
    print(f"PS_SMD_N_EVENTS={os.environ.get('PS_SMD_N_EVENTS')}")  # defaults to 1000
    print(f"DS batchsize: {ds.batch_size}")
    print("#### END DATASOURCE AND PSANA ENV VAR INFO ####\n")
    logger.info(f"Unique user rank: {rank}")

thisrun = next(ds.runs())

# Generate smalldata object
if ds.unique_user_rank():
    logger.info(
        "Opening the h5file %s, gathering at %d" % (h5_f_name, args.gather_interval)
    )
if args.psplot_live_mode:
    if ds.unique_user_rank():
        logger.info("Setting up psplot_live plots.")
    psplot_configs = config.get_psplot_configs(int(run))
    psplot_callbacks = psplot.PsplotCallbacks()

    for key, item in psplot_configs.items():
        callback_func = item.pop("callback")

        psplot_callbacks.add_callback(callback_func(**item), name=key)
    small_data = ds.smalldata(
        filename=None, batch_size=args.gather_interval, callbacks=[psplot_callbacks.run]
    )
else:
    small_data = ds.smalldata(filename=h5_f_name, batch_size=args.gather_interval)
if ds.unique_user_rank():
    logger.info("smalldata file has been successfully created.")


##########################################################
##
## Setting up the default detectors
##
##########################################################
if not ds.is_srv():  # srv nodes do not have access to detectors.
    default_dets = defaultDetectors(hutch.lower(), thisrun)
    if ds.unique_user_rank():
        logger.info("Default detectors loaded:")
        for ddet in default_dets:
            logger.info(f"{ddet.name}: {ddet.detname}")

    #
    # add stuff here to save all EPICS PVs.
    #
    if args.full or args.epicsAll:
        epicsPV = [k[0] for k in thisrun.epicsinfo]
        if len(epicsPV) > 0:
            try:
                logger.info("epicsStore names for epicsAll", epicsPV)
            except:
                pass
            logger.info("adding all epicsPVs....")
            default_dets.append(
                epicsDetector(PVlist=epicsPV, name="epicsAll", run=thisrun)
            )
    elif hasattr(config, "epicsArchFilePV") and len(config.epicsArchFilePV) > 0:
        default_dets.append(
            epicsDetector(PVlist=config.epicsArchFilePV, name="epicsUser", run=thisrun)
        )

    if len(config.epicsOncePV) > 0:
        EODet = epicsDetector(PVlist=config.epicsOncePV, name="epicsOnce", run=thisrun)
    else:
        EODet = None
    EODetData = {"epicsOnce": {}}
    EODetTS = None

    default_det_aliases = [det.name for det in default_dets]

    if args.rawFim:
        for fim in default_det_aliases:
            if fim.find("fim") >= 0:  # are you really a FIM?
                default_dets.append(
                    genericDetector(fim, run=thisrun, h5name="%s_raw" % fim)
                )

    dets = []
    int_dets = []
    if not args.default:
        dets = define_dets(int(args.run), config.detectors)
        if not skip_intg:
            int_dets = define_dets(int(args.run), integrating_detectors)
    if ds.unique_user_rank():
        logger.info(f"Detectors: {[det._name for det in dets]}")
        logger.info(f"Integrating detectors: {[det._name for det in int_dets]}")
    logger.debug(f"Rank {rank} detectors: {[det._name for det in dets]}")
    logger.debug(
        f"Rank {rank} integrating detectors: {[det._name for det in int_dets]}"
    )

    det_presence = {}
    if args.full:
        try:
            aliases = [dn for dn in thisrun.detnames]
            vetoDets = ["epicsinfo"]  # at least for run 339 of rixx43518

            for alias in aliases:
                det_presence[alias] = 1
                if alias in default_det_aliases:
                    continue
                if alias in vetoDets:
                    continue
                try:
                    thisDet = DetObject(alias, thisrun)
                    if alias.find("hsd") >= 0 and not args.nohsd:
                        hsdsplit = hsdsplitFunc()
                        thisDet.addFunc(hsdsplit)
                    else:
                        fullROI = ROIFunc(writeArea=True)
                        thisDet.addFunc(fullROI)
                    if ds.unique_user_rank():
                        print("adding detector for %s" % alias)
                    dets.append(thisDet)
                except:
                    pass

        except:
            pass

    evt_num = (
        -1
    )  # set this to default until I have a useable rank for printing updates...
    if ds.unique_user_rank():
        logger.info("And now the event loop user....")

    normdict = {}
    for det in int_dets:
        normdict[det._name] = {"count": 0, "timestamp_min": 0, "timestamp_max": 0}

event_iter = thisrun.events()

for evt_num, evt in enumerate(event_iter):
    det_data = detData(default_dets, evt)

    # If we don't have the epics once data, try to get it!
    if EODet is not None and EODetData["epicsOnce"] == {}:
        EODetData = detData([EODet], evt)
        EODetTS = evt._seconds + 631152000  # Convert to linux time.

    # detector data using DetObject
    userDict = {}
    for det in dets:
        try:
            det.getData(evt)
            det.processFuncs()
            userDict[det._name] = getUserData(det)
            try:
                envData = getUserEnvData(det)
                if len(envData.keys()) > 0:
                    userDict[det._name + "_env"] = envData
            except:
                pass
            det.processSums()
            # print(userDict[det._name])
        except Exception as e:
            logger.warning(f"Failed analyzing det {det} on evt {evt_num}")
            print(e)
            pass

    # Combine default data & user data into single dict.
    det_data.update(userDict)

    # Integrating detectors

    if len(int_dets) > 0:
        userDictInt = {}
        # Get summed fast detectors' data for integrating detector event
        for det in int_dets:
            normdict[det._name]["count"] += 1

            # for now, sum up all default data....
            for k, v in det_data.items():
                if isinstance(v, dict):
                    for kk, vv in v.items():
                        sumkey = k + "_sum_" + kk
                        if k not in normdict[det._name]:
                            normdict[det._name][k] = {}
                            normdict[det._name][k][sumkey] = np.array(vv)
                        if sumkey in normdict[det._name][k].keys():
                            normdict[det._name][k][sumkey] += np.array(vv)
                        else:
                            normdict[det._name][k][sumkey] = np.array(vv)
                else:
                    sumkey = k + "_sum"
                    if sumkey in normdict[det._name]:
                        normdict[det._name][sumkey] += v
                    else:
                        normdict[det._name][sumkey] = v

            normdict[det._name]["timestamp_max"] = max(
                normdict[det._name]["timestamp_max"], evt.timestamp
            )
            if normdict[det._name]["timestamp_min"] == 0:
                normdict[det._name]["timestamp_min"] = evt.timestamp
            else:
                normdict[det._name]["timestamp_min"] = min(
                    normdict[det._name]["timestamp_min"], evt.timestamp
                )

            if evt.EndOfBatch():
                det.getData(evt)

                if det.evt.dat is None:
                    logger.info(
                        f"Rank {rank}: Integrating detector {det._name} has no data on evt {evt_num}"
                    )
                    userDictInt[det._name] = (
                        {}
                    )  # so we can still get the summed fast data

                else:
                    det.processFuncs()
                    userDictInt[det._name] = {}
                    tmpdict = getUserData(det)
                    for k, v in tmpdict.items():
                        userDictInt[det._name][k] = v

                    try:
                        envData = getUserEnvData(det)
                        if len(envData.keys()) > 0:
                            userDictInt[det._name + "_env"] = envData
                    except:
                        pass

                # save data in integrating det dictionary & reset norm dictionary
                for k, v in normdict[det._name].items():
                    if isinstance(v, dict):
                        for kk, vv in v.items():
                            userDictInt[det._name][kk] = vv
                            normdict[det._name][k][kk] = (
                                vv * 0
                            )  # may not work for arrays....
                    else:
                        userDictInt[det._name][k] = v
                        normdict[det._name][k] = v * 0  # may not work for arrays....
                # print(userDictInt)
                small_data.event(evt, userDictInt, align_group="intg")

    # store event-based data
    # if det_data is not None:
    # DO WE STILL WANT THAT???
    # #remove data fields from the save_def_dets list
    # if 'veto' in save_def_dets:
    #     for k in save_def_dets['veto']:
    #         v = det_data.pop(k, None)
    #     #for k,v in det_data.items():
    #     #    if k not in save_def_dets['veto']:
    #     #        save_det_data[k]=v
    # if 'save' in save_def_dets:
    #     save_det_data={}
    #     for k,v in det_data.items():
    #         if k in save_def_dets['save']:
    #             save_det_data[k]=v
    #     det_data = save_det_data
    # #save what was selected to be saved.
    # #print('SAVE ',det_data)

    if len(int_dets) == 0 or args.all_events:
        small_data.event(evt, det_data)
    else:
        scan_data = det_data.get("scan", {})
        timing_data = det_data.get("timing", {})
        data_for_smd = {"scan": scan_data, "timing": timing_data}
        small_data.event(evt, data_for_smd)

    # the ARP will pass run & exp via the environment, if I see that info, the post updates
    if (
        (evt_num < 10)
        or (evt_num < 100 and (evt_num % 10) == 0)
        or (evt_num < 1000 and evt_num % 100 == 0)
        or (evt_num % 1000 == 0)
    ):
        if (os.environ.get("ARP_JOB_ID", None)) is not None and ds.unique_user_rank():
            try:
                requests.post(
                    os.environ["JID_UPDATE_COUNTERS"],
                    json=[
                        {
                            "key": "<b>Current Event / rank </b>",
                            "value": evt_num + 1,
                        }
                    ],
                )
            except:
                print("ARP update post failed")
                pass
        elif ds.unique_user_rank():
            print("Processed evt %d" % evt_num)

if not ds.is_srv():
    sumDict = {"Sums": {}}
    for det in dets:
        for key in det.storeSum().keys():
            try:
                sumData = small_data.sum(det.storeSum()[key])
                sumDict["Sums"]["%s_%s" % (det._name, key)] = sumData
            except:
                print("Problem with data sum for %s and key %s" % (det._name, key))
    if len(sumDict["Sums"].keys()) > 0 and small_data.summary:
        small_data.save_summary(sumDict)

    if ds.unique_user_rank():
        logger.info("Saving detector configuration to UserDataCfg")
        userDataCfg = {}
        for det in default_dets:
            # Make a list of configs not to be saved as lists of strings don't work in ps-4.2.5
            noConfigSave = ["scan", "damage"]
            if det.name not in noConfigSave:
                userDataCfg[det.name] = det.params_as_dict()
        for det in dets:
            try:
                userDataCfg[det._name] = det.params_as_dict()
            except:
                userDataCfg[det.name] = det.params_as_dict()
        for det in int_dets:
            try:
                userDataCfg[det._name] = det.params_as_dict()
            except:
                userDataCfg[det.name] = det.params_as_dict()
        if EODet is not None:
            EODetData = detOnceData(EODet, EODetData, EODetTS, args.noarch)
            if EODetData["epicsOnce"] != {}:
                userDataCfg[EODet.name] = EODet.params_as_dict()
                Config = {"UserDataCfg": userDataCfg}
                Config.update(EODetData)
            else:
                Config = {"UserDataCfg": userDataCfg}
        else:
            Config = {"UserDataCfg": userDataCfg}
        if small_data.summary:
            small_data.save_summary(Config)  # this only works w/ 1 rank!

# Finishing up:
# The filesystem seems to make smalldata.done fail. Some dirty tricks
# are needed here.
# Hopefully this can be removed soon.
logger.debug(f"Smalldata type for rank {rank}: {small_data._type}")

logger.debug(f"smalldata.done() on rank {rank}")
small_data.done()


# Epics data from the archiver
h5_rank = None
if small_data._type == "client" and small_data._full_filename is not None:
    if small_data._client_comm.Get_rank() == 0:
        h5_rank = rank

if rank == h5_rank:
    logger.info(f"Getting epics data from Archiver (rank: {rank})")
    import asyncio
    import h5py
    from smalldata_tools.common.epicsarchive import EpicsArchive, ts_to_datetime

    h5_for_arch = h5py.File(h5_f_name, "a")

    ts = h5_for_arch["timestamp"]
    start = ts_to_datetime(
        min(ts[: int(1e5)])
    )  # assumes the early timestamps are in the first 10k events
    end = ts_to_datetime(
        max(ts[int(-1e5) :])
    )  # assumes the early timestamps are in the last 10k events

    epics_archive = EpicsArchive()
    loop = asyncio.get_event_loop()
    pvs = [pv[0] if isinstance(pv, tuple) else pv for pv in config.epicsPV]
    coroutines = [
        epics_archive.get_points(PV=pv, start=start, end=end, raw=True, useMS=True)
        for pv in pvs
    ]

    logger.debug(f"Run PV retrieval (async)")
    data = loop.run_until_complete(asyncio.gather(*coroutines))

    # Save to files
    h5_for_arch.create_group("epics_archiver")
    for pv_, data in zip(config.epicsPV, data):
        if isinstance(pv_, tuple):
            # In format ("PV", "alias")
            pv = pv_[0]
            alias = pv_[1]
        else:
            pv = pv_
            alias = None
        pv = pv.replace(":", "_")
        if data == []:
            continue
        data = np.asarray(data)
        dset_name = pv if alias is None else alias
        dset = h5_for_arch.create_dataset(f"epics_archiver/{dset_name}", data=data)
        logger.debug(f"Saved {pv} from archiver data.")


end_prod_time = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
end_job = time.time()
prod_time = (end_job - start_job) / 60
if ds.unique_user_rank():
    print("########## JOB TIME: {:03f} minutes ###########".format(prod_time))
logger.debug("rank {0} on {1} is finished".format(rank, hostname))


# if args.sort:
#    if ds.unique_user_rank():
#        import subprocess
#        cmd = ['timestamp_sort_h5', h5_f_name, h5_f_name]

if ds.unique_user_rank():
    if os.environ.get("ARP_JOB_ID", None) is not None:
        requests.post(
            os.environ["JID_UPDATE_COUNTERS"],
            json=[
                {
                    "key": "<b>Last Event</b>",
                    "value": "~ %d cores * %d evts" % (size, evt_num),
                }
            ],
        )
    else:
        print(f"Last Event: {evt_num}")


if args.postRuntable and ds.unique_user_rank():
    print("Posting to the run tables.")
    locStr = ""
    try:
        runtable_data = {
            "Prod%s_end" % locStr: end_prod_time,
            "Prod%s_start" % locStr: begin_prod_time,
            "Prod%s_jobstart" % locStr: begin_job_time,
            "Prod%s_ncores" % locStr: size,
        }
    except:
        runtable_data = {
            "Prod%s_end" % locStr: end_prod_time,
            "Prod%s_start" % locStr: begin_prod_time,
            "Prod%s_jobstart" % locStr: begin_job_time,
        }
    if args.default:
        runtable_data["SmallData%s" % locStr] = "default"
    else:
        runtable_data["SmallData%s" % locStr] = "done"

    time.sleep(5)

    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(args.experiment)
    logger.debug("URL: ", ws_url)
    user = (args.experiment[:3] + "opr").replace("dia", "mcc")
    if os.environ.get("ARP_LOCATION", None) == "S3DF":
        with open("/sdf/group/lcls/ds/tools/forElogPost.txt") as reader:
            answer = reader.readline()

        r = requests.post(
            ws_url,
            params={"run_num": args.run},
            json=runtable_data,
            auth=HTTPBasicAuth(args.experiment[:3] + "opr", answer[:-1]),
        )
        logger.debug(r)
    if det_presence != {}:
        rp = requests.post(
            ws_url,
            params={"run_num": args.run},
            json=det_presence,
            auth=HTTPBasicAuth(args.experiment[:3] + "opr", answer[:-1]),
        )
        logger.debug(rp)
