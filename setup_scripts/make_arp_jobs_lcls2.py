import sys
import requests
import json
import shutil
from krtc import KerberosTicket
import argparse
from pathlib import Path
import os

parser = argparse.ArgumentParser()
parser.add_argument(
    "--experiment",
    help="experiment name",
    type=str,
    default=os.environ.get("EXPERIMENT", ""),
)
parser.add_argument(
    "--queue", help="Partition on which to run the jobs", type=str, default="milano"
)
parser.add_argument(
    "--psplot_live", help="Make psplot-live job as well", action="store_true"
)
parser.add_argument(
    "--config", help="Config file for the producer", type=str, default=None
)
parser.add_argument(
    "--cube", help="Cube config to make cube job for.", type=str, default=None
)
parser.add_argument(
    "--all_events",
    help="Save all events, not just integrating detectors",
    action="store_true",
)
args = parser.parse_args()

exp = args.experiment
hutch = exp[:3].lower()
SDF_BASE = Path(f"/sdf/data/lcls/ds/{hutch}/{exp}")

# job arguments
if "milano" in args.queue:
    location = "S3DF"
    nodes = 2
    account = f"lcls:{exp}"

    # smd
    executable_smd = str(
        SDF_BASE / "results/smalldata_tools/arp_scripts/submit_smd2.sh"
    )
    trigger = "START_OF_RUN"

    # cube
    executable_cube = str(
        SDF_BASE / "results/smalldata_tools/arp_scripts/lets_cube_lclc2.sh"
    )
    trigger_cube = "MANUAL"

    # summaries
    executable_summaries = str(
        SDF_BASE / "results/smalldata_tools/arp_scripts/submit_plots.sh"
    )
    queue_summaries = args.queue
    run_param_name = "SmallData"
    args_summaries = ""
else:
    print("No known system related to this queue. Exit now.")
    sys.exit(0)


job_defs = []

# smallData job
parameters = f"--partition {args.queue} --postRuntable --nodes {nodes} --wait"
if args.config is not None:
    parameters += f" --config {args.config}"

smd_job_def = {
    "name": "smd",
    "executable": executable_smd,
    "trigger": trigger,
    "location": location,
    "parameters": parameters,
}

job_defs.append(smd_job_def)

if args.all_events:
    smd_job_def["parameters"] += " --all_events"
    smd_job_def["trigger"] = "MANUAL"
    job_defs.append(smd_job_def)

# psplot_live job
if args.psplot_live:
    parameters = f"--partition {args.queue} --nodes {nodes} --wait --psplot_live"
    if args.config is not None:
        parameters += f" --config {args.config}"

    job_defs.append(
        {
            "name": "smd_psplot_live",
            "executable": executable,
            "trigger": trigger,
            "location": location,
            "parameters": parameters,
        }
    )

# cube job
if args.cube is not None:
    job_defs.append(
        {
            "name": "cube",
            "executable": executable_cube,
            "trigger": trigger_cube,
            "location": location,
            "parameters": f"--config {args.cube}",
        }
    )

# summary plots job
if hutch in ["rix", "xcs", "xpp", "mfx"]:
    job_defs.append(
        {
            "name": "run_summary",
            "executable": executable_summaries,
            "trigger": "RUN_PARAM_IS_VALUE",
            "run_param_name": run_param_name,
            "run_param_value": "done",
            "location": location,
            "parameters": f"--queue {queue_summaries} {args_summaries}",
        }
    )


for job_def in job_defs:
    krbheaders = KerberosTicket("HTTP@pswww.slac.stanford.edu").getAuthHeaders()
    ws_url = "https://pswww.slac.stanford.edu/ws-kerb/lgbk/lgbk/{}/ws/create_update_workflow_def".format(
        exp
    )
    r = requests.post(ws_url, headers=krbheaders, json=job_def)
    r.raise_for_status()
    print("\nJOB CREATION LOG: {}".format(r.json()))
