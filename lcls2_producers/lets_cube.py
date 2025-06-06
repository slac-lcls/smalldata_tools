import os
import sys
import time
import numpy as np
import psana
from mpi4py import MPI

COMM = MPI.COMM_WORLD
rank = COMM.Get_rank()
size = COMM.Get_size()

import logging

from pathlib import Path
from importlib import import_module

smdata_tools_path = Path(__file__).parent.parent
sys.path.insert(0, str(smdata_tools_path))

import smalldata_tools
from smalldata_tools.lcls2.DetObject import NullDetObject

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
log_format = smalldata_tools.LogConfig.FullFormat

handler = logging.StreamHandler()
formatter = logging.Formatter(log_format)
handler.setFormatter(formatter)
logger.addHandler(handler)

# from smalldata_tools.lcls2.cube.utils import PsanaNode


def parse_args():  # move to another file?
    import argparse

    parser = argparse.ArgumentParser(description="Arguments for LCLS2 cube processing.")
    parser.add_argument(
        "-r",
        "--run",
        type=int,
        default=int(os.environ.get("RUN_NUM", "")),
        help="Run number to process",
    )
    parser.add_argument(
        "-e",
        "--experiment",
        type=str,
        default=os.environ.get("EXPERIMENT", ""),
        help="Experiment name",
    )
    parser.add_argument(
        "-C",
        "--config",
        type=str,
        default=None,
        required=True,
        help="Configuration file name (without the cube_config and .py suffix)",
    )
    parser.add_argument(
        "--batchsize", type=int, default=1024, help="Batch size for processing"
    )
    parser.add_argument("--nevents", help="number of events", type=int, default=-1)
    parser.add_argument(
        "--directory", help="directory to save the output", type=str, default=None
    )
    return parser.parse_args()


start_time = time.time()

args = parse_args()
exp = args.experiment
hutch = exp[:3]
run = args.run

# Setup datasource
os.environ["PS_SMD_N_EVENTS"] = f"{args.batchsize}"  # SMD0 to EB batches

datasource_args = {"exp": exp, "run": run}
datasource_args["batch_size"] = args.batchsize  # EB to DB batches
if args.nevents != -1:
    datasource_args["max_events"] = args.nevents
# datasource_args['smd_callback'] = binning_obj.smd_callback
# datasource_args['smd_callback'] = callbacks.simple_callback
# datasource_args['small_xtc'] = ['xgmd']

ds = psana.DataSource(**datasource_args)
myrun = next(ds.runs())

if ds.unique_user_rank():
    print("\n#### DATASOURCE AND PSANA ENV VAR INFO ####")
    print(f"Instantiated data source with arguments: {datasource_args}")
    print(f"MPI size: {size}")
    print(f"PS_EB_NODES={os.environ.get('PS_EB_NODES')}")
    print(f"PS_SRV_NODES={os.environ.get('PS_SRV_NODES')}")
    print(f"PS_SMD_N_EVENTS={os.environ.get('PS_SMD_N_EVENTS')}")  # defaults to 1000
    print(f"DS batchsize: {ds.batch_size}")
    print("#### END DATASOURCE AND PSANA ENV VAR INFO ####\n")

""" Main processing """
if ds.is_srv():
    from smalldata_tools.lcls2.cube.srv import CubeSrv

    cube_srv = CubeSrv()
    cube_srv.set_file_handle(
        exp=exp,
        run_num=run,
        filepath=args.directory,
    )
    cube_srv.run()

else:
    import smalldata_tools.lcls2.cube as cube
    import smalldata_tools.lcls2.cube.event_engine as event_engine
    from smalldata_tools.lcls2.cube.utils import get_config_file

    config_module = get_config_file(args.config, smdata_tools_path / "lcls2_producers")
    config = import_module(config_module)

    cube_obj = cube.cube.get_cube(myrun, engine=event_engine.smalldata_tools_engine)

    processors = config.detectors(myrun)
    # weed out NullDetObjects:
    temp = []
    for detector in processors:
        if isinstance(detector, NullDetObject):
            if ds.unique_user_rank():

                logger.warning(f"Detector {detector._name} is not defined.")
        else:
            temp.append(detector)
    processors = temp
    cube_obj.add_processors(processors)

    screener = config.screener(myrun)
    cube_obj.set_event_screener(screener)

    cube_obj.run()


COMM.Barrier()
if ds.unique_user_rank():
    print(f"### Job done; total time = {(time.time() - start_time)/60:.2f} minutes.")
