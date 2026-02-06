#!/usr/bin/env python

import json
import time
from datetime import datetime

import numpy as np
import psana

begin_job_time = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
start_job = time.time()

import argparse
import logging
import os
import socket

import requests

try:
    import psutil
except Exception:
    psutil = None
import sys
from importlib import import_module
from pathlib import Path

from mpi4py import MPI
from requests.auth import HTTPBasicAuth

COMM = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

COMM.Barrier()
job_perf_start = MPI.Wtime()
DEBUG_GLOBAL_TIMING = os.environ.get("SMD_DEBUG_TIMING", "0") != "0"
DEBUG_DET = os.environ.get("SMD_DEBUG_DET", "0") != "0"
DEBUG_EVT_TIMING = os.environ.get("SMD_DEBUG_EVT_TIMING", "0") != "0"
DEBUG_EVT_TIMING_EVERY = int(os.environ.get("SMD_DEBUG_EVT_TIMING_EVERY", "0"))
_last_mark = job_perf_start
_evt_last_mark = None
args_parsed_mark = None
config_import_mark = None
output_ready_mark = None


def _pss_cur_mb():
    if psutil is not None:
        try:
            return psutil.Process(os.getpid()).memory_full_info().pss / (1024**2)
        except Exception:
            pass
    try:
        with open("/proc/self/smaps_rollup", "r") as fh:
            for line in fh:
                if line.startswith("Pss:"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return float(parts[1]) / 1024.0
    except Exception:
        pass
    return -1.0


def _debug_timing_evt(evt_num):
    if not DEBUG_EVT_TIMING:
        return False
    if DEBUG_EVT_TIMING_EVERY <= 0:
        return evt_num == 0
    return evt_num % DEBUG_EVT_TIMING_EVERY == 0


def _debug_global_mark(label, now):
    global _last_mark
    if not DEBUG_GLOBAL_TIMING or rank != 0:
        return
    delta = now - _last_mark
    since = now - job_perf_start
    print(
        f"[DEBUG] {label} since_start={since:.6f}s delta={delta:.6f}s pss_mb={_pss_cur_mb():.2f}"
    )
    _last_mark = now


def _debug_evt_reset(now):
    global _evt_last_mark
    _evt_last_mark = now


def _debug_evt_mark(label, now):
    global _evt_last_mark
    if not DEBUG_EVT_TIMING:
        return
    if _evt_last_mark is None:
        _evt_last_mark = now
    delta = now - _evt_last_mark
    since = now - job_perf_start
    print(
        f"[DEBUG] rank {rank} {label} since_start={since:.6f}s delta={delta:.6f}s pss_mb={_pss_cur_mb():.2f}"
    )
    _evt_last_mark = now


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
    if os.getenv("CONDA_DEFAULT_ENV") is not None:
        logger.info(
            "psana conda environment is {0}".format(os.environ["CONDA_DEFAULT_ENV"])
        )
    else:
        spack_env = os.getenv("SPACK_ENV")
        if spack_env is not None:
            logger.info(f"psana spack environment is {spack_env}")
        else:
            logger.info("Could not determine what psana environment is in use.")


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
    azav_pyfai_args = {}
    polynomial_args = {}
    sum_algo_args = {}
    pressio_compression_args = {}

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
    if "get_azav_pyfai" in dir(config):
        azav_pyfai_args = config.get_azav_pyfai(run)
    if "get_polynomial_correction" in dir(config):
        polynomial_args = config.get_polynomial_correction(run)
    if "get_sum_algos" in dir(config):
        sum_algo_args = config.get_sum_algos(run)
    if "get_pressio_compression" in dir(config):
        pressio_compression_args = config.get_pressio_compression(run)

    dets = []

    for detname in det_list:
        det_start = time.perf_counter()
        det_last = det_start

        def _dd_mark(label, t0=det_start):
            nonlocal det_last
            if not DEBUG_DET:
                return
            now = time.perf_counter()
            print(
                "[DEBUG Rank{}] define_dets det={} {} delta={:.6f}s total={:.6f}s".format(
                    rank, detname, label, now - det_last, now - t0
                )
            )
            det_last = now

        havedet = detname in thisrun.detnames

        # Common mode (default: None)
        common_mode = None

        if not havedet:
            _dd_mark("skip_missing")
            continue

        if detname.find("fim") >= 0 or detname.find("w8") >= 0:
            # why different name for w8? To differentiate full det from BLD
            det = DetObject(
                detname, thisrun, common_mode=common_mode, name=f"det_{detname}"
            )
        else:
            det = DetObject(detname, thisrun, common_mode=common_mode)
        _dd_mark("DetObject")
        logger.debug(f"Instantiated det {detname}: {det}")

        #             **** Compression MUST be the first operation ****             #
        # It is "transparent" in that it just compresses then decompresses the data #
        if detname in pressio_compression_args:
            det.addFunc(pressioCompressDecompress(**pressio_compression_args[detname]))
            _dd_mark("pressio")

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
            _dd_mark("hsd_setup")
            dets.append(det)
            _dd_mark("append")
            continue

        ####################################
        ######## Standard detectors ########
        ####################################
        if detname in azav_args:
            det.addFunc(azimuthalBinning(**azav_args[detname]))
            _dd_mark("azav")

        if detname in azav_pyfai_args:
            det.addFunc(azav_pyfai(**azav_pyfai_args[detname]))
            _dd_mark("azav_pyfai")

        if detname in rois_args:
            # ROI extraction
            for iROI, ROI in enumerate(rois_args[detname]):
                proj_ax = ROI.pop("proj_ax", None)

                thisROIFunc = ROIFunc(**ROI)
                if proj_ax is not None:
                    thisROIFunc.addFunc(projectionFunc(axis=proj_ax))
            det.addFunc(thisROIFunc)
        _dd_mark("rois")

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
            _dd_mark("droplet2photon")

        # if detname in dimgs_args:
        #     # Detector image (det.raw.image())
        #     dimg_func = detImageFunc(**dimgs_args[detname])
        #     det.addFunc(dimg_func)

        if detname in wfs_int_args:
            # Waveform integration
            wfs_int_func = WfIntegration(**wfs_int_args[detname])
            det.addFunc(wfs_int_func)
            _dd_mark("wfs_int")

        if detname in wfs_hitfinder_args:
            # Simple hit finder on waveform
            det.addFunc(SimpleHitFinder(**wfs_hitfinder_args[detname]))
            _dd_mark("wfs_hitfinder")

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
            _dd_mark("wfs_svd")

        if detname in droplet_args:
            if "nData" in droplet_args[detname].keys():
                nData = droplet_args[detname].pop("nData")
            else:
                nData = None
            dFunc = dropletFunc(**droplet_args[detname])
            dFunc.addFunc(sparsifyFunc(nData=nData))
            det.addFunc(dFunc)
            _dd_mark("droplet")

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
            _dd_mark("poly_corr")

        if sum_algo_args == {}:
            # Store calib by default even if algos. dictionary not defined
            # NOTE: This can cause issues if you are relying on it for sums, and have
            #       added a waveform PV detector (like a QADC) to the detnames list.
            #       If this is causing issues, the simple fix is just to define the
            #       sum algorithms in the config file explicitly
            det.storeSum(sumAlgo="calib")
        else:
            # Add `all` algorithms
            if "all" in sum_algo_args:
                for algo in sum_algo_args["all"]:
                    det.storeSum(sumAlgo=algo)
            if detname in sum_algo_args:
                for algo in sum_algo_args[detname]:
                    det.storeSum(sumAlgo=algo)
        _dd_mark("sum_algos")

        logger.debug(f"Rank {rank} Add det {detname}: {det}")
        dets.append(det)
        _dd_mark("append")

    return dets


COMM.Barrier()
if DEBUG_GLOBAL_TIMING and rank == 0:
    _debug_global_mark("generic module imports", now=MPI.Wtime())

# General Workflow
# This is meant for arp which means we will always have an exp and run

fpath = os.path.dirname(os.path.abspath(__file__))
fpathup = "/".join(fpath.split("/")[:-1])
sys.path.append(fpathup)
if rank == 0:
    logger.info(fpathup)

import psplot
from psmon import publish

from smalldata_tools.ana_funcs.azav_pyfai import azav_pyfai
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning
from smalldata_tools.ana_funcs.compression import pressioCompressDecompress
from smalldata_tools.ana_funcs.detector_corrections import PolynomialCurveCorrection
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.droplet2Photons import droplet2Photons
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc, projectionFunc
from smalldata_tools.ana_funcs.smd_svd import SvdFit
from smalldata_tools.ana_funcs.sparsifyFunc import sparsifyFunc, unsparsifyFunc
from smalldata_tools.ana_funcs.waveformFunc import (
    SimpleHitFinder,
    WfIntegration,
    hsdROIFunc,
    hsdsplitFunc,
)
from smalldata_tools.common.detector_base import detData, getUserData, getUserEnvData
from smalldata_tools.lcls2.default_detectors import (
    detOnceData,
    epicsDetector,
    genericDetector,
)
from smalldata_tools.lcls2.DetObject import DetObject
from smalldata_tools.lcls2.hutch_default import defaultDetectors

# Moved the unused imports to comments to reduce confusion
# from smalldata_tools.ana_funcs.photons import photonFunc
# from smalldata_tools.ana_funcs.roi_rebin import (imageFunc, spectrumFunc)
# from smalldata_tools.ana_funcs.waveformFunc import (getCMPeakFunc,
#                                                    hitFinderCFDFunc,
#                                                    hsdBaselineCorrectFunc,
#                                                    templateFitFunc)
# from smalldata_tools.utilities import checkDet, printMsg

COMM.Barrier()
if DEBUG_GLOBAL_TIMING and rank == 0:
    _debug_global_mark("custom module imports", now=MPI.Wtime())

# Constants
HUTCHES = ["TMO", "RIX", "UED", "MFX", "XPP"]

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
parser.add_argument(
    "--psdm_dir",
    help="Override SIT_PSDM_DATA regardless of environment variables set.",
    type=str,
    default="",
)

args = parser.parse_args()


COMM.Barrier()
if DEBUG_GLOBAL_TIMING and rank == 0:
    _debug_global_mark("args parsed", now=MPI.Wtime())

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

det_presence = {}


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

if not args.psdm_dir:
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
                        f"Waited {str(n_wait*10)}s, still no files available. Giving up."
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
                print(
                    f"Files appeared after {str(waitFilesEnd-waitFilesStart)} seconds"
                )

        xtc_files = get_xtc_files(PSDM_BASE, exp, run)
        if len(xtc_files) == 0:
            raise RuntimeError(
                f"We have no xtc files for run {run} in {exp} in the offline system."
            )
    else:
        logger.warning("On an unknow system, things may get weird.")
else:
    logger.info(f"Requested data be found at: {args.psdm_dir}")

# Get output file, check if we can write to it
h5_f_name = get_sd_file(args.directory, exp, hutch)

# Create data source.
datasource_args = {
    "exp": exp,
    "run": int(run),
    "monitor": False,
    "log_level": "INFO",
}
if not args.psdm_dir:
    if useFFB:
        datasource_args["live"] = True

    if rank == 0:
        logger.info("Opening the data source:")
else:
    datasource_args["dir"] = args.psdm_dir
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
    intg_ds_start = time.perf_counter()
    ds = psana.DataSource(**datasource_args)
    intg_ds_end = time.perf_counter()
    intg_run_start = time.perf_counter()
    thisrun = next(ds.runs())
    intg_run_end = time.perf_counter()
    det_count = len(getattr(thisrun, "detnames", []))
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

if hasattr(config, "xdetectors") and config.xdetectors:
    datasource_args["xdetectors"] = config.xdetectors

COMM.Barrier()
if DEBUG_GLOBAL_TIMING and rank == 0:
    _debug_global_mark("ds args setup", now=MPI.Wtime())

ds = psana.DataSource(**datasource_args)
COMM.Barrier()
if DEBUG_GLOBAL_TIMING and rank == 0:
    _debug_global_mark("ds init", now=MPI.Wtime())

if ds.unique_user_rank():
    print("#### DATASOURCE AND PSANA ENV VAR INFO ####")
    print(f"Instantiated data source with arguments: {datasource_args}")
    print(f"MPI size: {size}")
    print(f"PS_EB_NODES={os.environ.get('PS_EB_NODES')}")
    print(f"PS_SRV_NODES={os.environ.get('PS_SRV_NODES')}")
    print(f"PS_SMD_N_EVENTS={os.environ.get('PS_SMD_N_EVENTS')}")  # defaults to 1000
    print(f"DS batchsize: {ds.batch_size}")
    print(f"SMD gather_interval: {args.gather_interval}")
    print("#### END DATASOURCE AND PSANA ENV VAR INFO ####\n")
    logger.info(f"Unique user rank: {rank}")


try:
    ps_srv_nodes = int(os.environ.get("PS_SRV_NODES", "0") or 0)
except ValueError:
    ps_srv_nodes = 0
if ps_srv_nodes > 0 and rank >= size - ps_srv_nodes:
    print(f"[INFO] srv node rank {rank} hostname {hostname}")
try:
    ps_eb_nodes = int(os.environ.get("PS_EB_NODES", "0") or 0)
except ValueError:
    ps_eb_nodes = 0
if ps_eb_nodes > 0 and 1 <= rank <= ps_eb_nodes:
    print(f"[INFO] eb node rank {rank} hostname {hostname}")

small_data = None
small_data_enabled = True
if ps_srv_nodes == 0:
    small_data_enabled = False
    if ds.unique_user_rank():
        logger.warning(
            "PS_SRV_NODES=0: skipping ds.smalldata and all small_data operations."
        )
thisrun = next(ds.runs())
COMM.Barrier()
run_init_mark = MPI.Wtime()
if DEBUG_GLOBAL_TIMING and rank == 0:
    _debug_global_mark("run init", now=run_init_mark)

# Generate smalldata object
if small_data_enabled:
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
    if args.psplot_live_mode:
        small_data = ds.smalldata(
            filename=None,
            batch_size=args.gather_interval,
            callbacks=[psplot_callbacks.run],
        )
    else:
        small_data = ds.smalldata(filename=h5_f_name, batch_size=args.gather_interval)
        # small_data = ds.smalldata(batch_size=args.gather_interval)
        # if ds.unique_user_rank():
        #    logger.info("Smalldata object created without filename.")
    if ds.unique_user_rank():
        logger.info("smalldata file has been successfully created.")
else:
    if ds.unique_user_rank():
        logger.warning("PS_SRV_NODES=0: small_data disabled; skipping ds.smalldata.")

# For debugging, bd timing
if DEBUG_GLOBAL_TIMING and ds.is_bd():
    bd_node_comm = thisrun.comms.get_bd_node_comm()
    bd_node_comm.Barrier()
    smd_init_mark = MPI.Wtime()
    print(
        f"[DEBUG] rank {rank} smd init since_start={smd_init_mark - job_perf_start:.6f}s delta={smd_init_mark - run_init_mark:.6f}s pss_mb={_pss_cur_mb():.2f}"
    )

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
    # add stuff here to save all EPICS PVss.
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
            intg_dets_start = time.perf_counter()
            int_dets = define_dets(int(args.run), integrating_detectors)
            intg_dets_end = time.perf_counter()
    if ds.unique_user_rank():
        logger.info(f"Detectors: {[det._name for det in dets]}")
        logger.info(f"Integrating detectors: {[det._name for det in int_dets]}")
    logger.debug(f"Rank {rank} detectors: {[det._name for det in dets]}")
    logger.debug(
        f"Rank {rank} integrating detectors: {[det._name for det in int_dets]}"
    )

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

if DEBUG_GLOBAL_TIMING and ds.is_bd():
    bd_node_comm = thisrun.comms.get_bd_node_comm()
    bd_node_comm.Barrier()
    default_det_setup_mark = MPI.Wtime()
    print(
        f"[DEBUG] rank {rank} default det setup since_start={default_det_setup_mark - job_perf_start:.6f}s delta={default_det_setup_mark - smd_init_mark:.6f}s pss_mb={_pss_cur_mb():.2f}"
    )

event_iter = thisrun.events()

cn_events = 0
sum_det_data_time = 0.0
first_evt_mark = None
smd_save_sum_mark = None
smd_save_cfg_mark = None

for evt_num, evt in enumerate(event_iter):
    if DEBUG_EVT_TIMING and evt_num == 0:
        _debug_evt_reset(MPI.Wtime())
    det_data = detData(default_dets, evt)
    if _debug_timing_evt(evt_num):
        _debug_evt_mark(f"evt {evt_num} default detData", now=MPI.Wtime())

    # If we don't have the epics once data, try to get it!
    if EODet is not None and EODetData["epicsOnce"] == {}:
        EODetData = detData([EODet], evt)
        EODetTS = evt._seconds + 631152000  # Convert to linux time.

    # detector data using DetObject
    userDict = {}
    for det in dets:
        try:
            det_start_mark = MPI.Wtime()
            det.getData(evt)
            det_getdata_mark = MPI.Wtime()
            if DEBUG_DET:
                print(
                    f"[DEBUG] rank {rank} evt {evt_num} det {det._name} getData time: {det_getdata_mark - det_start_mark:.6f}s"
                )
            det.processFuncs()
            det_process_mark = MPI.Wtime()
            if DEBUG_DET:
                print(
                    f"[DEBUG] rank {rank} evt {evt_num} det {det._name} processFuncs time: {det_process_mark - det_getdata_mark:.6f}s"
                )
            userDict[det._name] = getUserData(det)
            det_userdata_mark = MPI.Wtime()
            if DEBUG_DET:
                print(
                    f"[DEBUG] rank {rank} evt {evt_num} det {det._name} getUserData time: {det_userdata_mark - det_process_mark:.6f}s"
                )
            try:
                envData = getUserEnvData(det)
                if len(envData.keys()) > 0:
                    userDict[det._name + "_env"] = envData
            except:
                pass
            det_envdata_mark = MPI.Wtime()
            if DEBUG_DET:
                print(
                    f"[DEBUG] rank {rank} evt {evt_num} det {det._name} getUserEnvData time: {det_envdata_mark - det_userdata_mark:.6f}s"
                )
            det.processSums()
            det_processsums_mark = MPI.Wtime()
            if DEBUG_DET:
                print(
                    f"[DEBUG] rank {rank} evt {evt_num} det {det._name} processSums time: {det_processsums_mark - det_envdata_mark:.6f}s"
                )
            # print(userDict[det._name])
        except Exception as e:
            logger.warning(
                f"Failed analyzing det {det} on evt {evt_num}: {e}", exc_info=True
            )
            raise
            # logger.warning(f"Failed analyzing det {det} on evt {evt_num}")
            # print(e)
            # pass

    if _debug_timing_evt(evt_num):
        _debug_evt_mark(f"evt {evt_num} user dets", now=MPI.Wtime())
    # Combine default data & user data into single dict.
    det_data.update(userDict)

    # Integrating detectors

    if len(int_dets) > 0:
        intg_start = time.perf_counter()
        userDictInt = {}
        # Get summed fast detectors' data for integrating detector
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
                intg_evt_start = time.perf_counter()
                det.getData(evt)
                intg_get = time.perf_counter()
                intg_funcs = intg_get

                if det.evt.dat is None:
                    logger.info(
                        f"Rank {rank}: Integrating detector {det._name} has no data on evt {evt_num}"
                    )
                    userDictInt[det._name] = (
                        {}
                    )  # so we can still get the summed fast data

                else:
                    det.processFuncs()
                    intg_funcs = time.perf_counter()
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
                if small_data_enabled:
                    small_data.event(evt, userDictInt, align_group="intg")
                intg_write = time.perf_counter()
        if _debug_timing_evt(evt_num):
            _debug_evt_mark(f"evt {evt_num} intg dets", now=MPI.Wtime())

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

    if small_data_enabled:
        if len(int_dets) == 0 or args.all_events:
            small_data.event(evt, det_data)
        else:
            scan_data = det_data.get("scan", {})
            timing_data = det_data.get("timing", {})
            data_for_smd = {"scan": scan_data, "timing": timing_data}
            small_data.event(evt, data_for_smd)
    if _debug_timing_evt(evt_num):
        _debug_evt_mark(f"evt {evt_num} event store", now=MPI.Wtime())

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
                    timeout=10,
                )
            except:
                print("ARP update post failed")
                pass
        elif ds.unique_user_rank():
            print("Processed evt %d" % evt_num)
    cn_events += 1

    # Add first_evt_mark at the end of the loop to make sure fist-event
    # heavy processing time is captured.
    if evt_num == 0:
        first_evt_mark = MPI.Wtime()
    if DEBUG_GLOBAL_TIMING and evt_num == 0 and ds.is_bd():
        print(
            f"[DEBUG] rank {rank} first evt start since_start={first_evt_mark - job_perf_start:.6f}s delta={first_evt_mark - default_det_setup_mark:.6f}s pss_mb={_pss_cur_mb():.2f}"
        )

last_evt_mark = None
if (
    DEBUG_GLOBAL_TIMING
    and ds.is_bd()
    and not ds.is_srv()
    and first_evt_mark is not None
):
    last_evt_mark = MPI.Wtime()
    print(
        f"[DEBUG] rank {rank} last evt processing since_start={last_evt_mark - job_perf_start:.6f}s delta={last_evt_mark - first_evt_mark:.6f}s pss_mb={_pss_cur_mb():.2f}"
    )

if not ds.is_srv():
    sumDict = {"Sums": {}}
    if small_data_enabled:
        for det in dets:
            for key in det.storeSum().keys():
                try:
                    sumData = small_data.sum(det.storeSum()[key])
                    sumDict["Sums"]["%s_%s" % (det._name, key)] = sumData
                except:
                    print("Problem with data sum for %s and key %s" % (det._name, key))
        if len(sumDict["Sums"].keys()) > 0 and small_data.summary:
            small_data.save_summary(sumDict)

    if DEBUG_EVT_TIMING and ds.is_bd() and last_evt_mark is not None:
        smd_save_sum_mark = MPI.Wtime()
        print(
            f"[DEBUG] rank {rank} smd done save sum since_start={smd_save_sum_mark - job_perf_start:.6f}s delta={smd_save_sum_mark - last_evt_mark:.6f}s pss_mb={_pss_cur_mb():.2f}"
        )

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
        if small_data_enabled and small_data.summary:
            small_data.save_summary(Config)  # this only works w/ 1 rank!

    if DEBUG_EVT_TIMING and ds.is_bd() and smd_save_sum_mark is not None:
        smd_save_cfg_mark = MPI.Wtime()
        print(
            f"[DEBUG] rank {rank} smd done save cfg since_start={smd_save_cfg_mark - job_perf_start:.6f}s delta={smd_save_cfg_mark - smd_save_sum_mark:.6f}s pss_mb={_pss_cur_mb():.2f}"
        )

# Finishing up:
# The filesystem seems to make smalldata.done fail. Some dirty tricks
# are needed here.
# Hopefully this can be removed soon.
if small_data_enabled:
    logger.debug(f"Smalldata type for rank {rank}: {small_data._type}")

    logger.debug(f"smalldata.done() on rank {rank}")
    small_data.done()

if (
    DEBUG_EVT_TIMING
    and ds.is_bd()
    and not ds.is_srv()
    and smd_save_cfg_mark is not None
):
    smd_done_non_srv_mark = MPI.Wtime()
    print(
        f"[DEBUG] rank {rank} smd done all since_start={smd_done_non_srv_mark - job_perf_start:.6f}s delta={smd_done_non_srv_mark - smd_save_cfg_mark:.6f}s pss_mb={_pss_cur_mb():.2f}"
    )
# Special print only for srv nodes
if DEBUG_GLOBAL_TIMING and ds.is_srv():
    smd_done_mark = MPI.Wtime()
    print(
        f"[DEBUG] rank {rank} srv node done at since_start={smd_done_mark - job_perf_start:.6f}s delta={smd_done_mark - run_init_mark:.6f}s pss_mb={_pss_cur_mb():.2f}"
    )

# Epics data from the archiver
h5_rank = None
if small_data_enabled:
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

    if hasattr(config, "epicsPV") and config.epicsPV:
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
    h5_for_arch.close()


MPI.COMM_WORLD.Barrier()
if DEBUG_GLOBAL_TIMING and rank == 0:
    _debug_global_mark("psana all done", now=MPI.Wtime())

if ds.unique_user_rank():
    import h5py

    h5_for_arch = h5py.File(h5_f_name, "a")

    def write_config(file_handle, base_name, cfg_dict):
        for key, val in cfg_dict.items():
            if isinstance(val, dict):
                file_handle.create_group(f"{base_name}/{key}")
                write_config(file_handle, f"{base_name}/{key}", val)
            else:
                file_handle.create_dataset(f"{base_name}/{key}", data=val)

    if "UserDataCfg" not in h5_for_arch.keys():
        # This guard is necessary for serial datasource when running interactively
        h5_for_arch.create_group("UserDataCfg")
        write_config(h5_for_arch, "UserDataCfg", Config["UserDataCfg"])

    h5_for_arch.close()


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
            timeout=10,
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
            timeout=10,
        )
        logger.debug(r)
    if det_presence != {} and os.environ.get("ARP_LOCATION", None) == "S3DF":
        rp = requests.post(
            ws_url,
            params={"run_num": args.run},
            json=det_presence,
            auth=HTTPBasicAuth(args.experiment[:3] + "opr", answer[:-1]),
            timeout=10,
        )
        logger.debug(rp)
