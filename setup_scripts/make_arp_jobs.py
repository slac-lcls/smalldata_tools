import sys
import requests
import json
import shutil
from krtc import KerberosTicket
import argparse
from pathlib import Path
import os

parser = argparse.ArgumentParser()
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--queue', help='Queue on which to run the jobs', type=str, default='psanaq')
parser.add_argument('--cube', help='Make cube job as well', type=int, default=0)
args = parser.parse_args()

exp = args.experiment
hutch = exp[:3].lower()
FFB_BASE = Path("/cds/data/drpsrcf/{hutch}/{exp}/scratch")
PSANA_BASE = Path("/cds/data/psdm/{hutch}/{exp}")
SDF_BASE = Path(f"/sdf/data/lcls/ds/{hutch}/{exp}")

# job arguments
if 'milano' in args.queue:
    location = 'S3DF'
    cores = 120
    account = f"lcls:{exp}"

    # smd
    executable = str(SDF_BASE / 'results/smalldata_tools/arp_scripts/submit_smd.sh')
    trigger = 'START_OF_RUN'
    
    # cube
    executable_cube = str(SDF_BASE / 'results/smalldata_tools/arp_scripts/cubeRun.sh')
    run_param_name = 'SmallData'
    args_cube = f''
    
    # summaries
    executable_summaries = str(SDF_BASE / 'results/smalldata_tools/arp_scripts/submit_plots.sh')
    queue_summaries = args.queue
    args_summaries = ''
else:
    print('No known system related to this queue. Exit now.')
    sys.exit(0)
    


    
job_defs = []

# smallData job
job_defs.append( {
    'name': 'smd',
    'executable': executable,
    'trigger': trigger,
    'location': location,
    'parameters': f'--queue {args.queue} --norecorder --postRuntable --cores {cores} --wait'
    } )

# cube job
if args.cube>0:
    job_defs.append( {
        'name': 'cube',
        'executable': executable_cube,
        'trigger': 'RUN_PARAM_IS_VALUE',
        'run_param_name' : run_param_name,
        'run_param_value' : 'done',
        'location': location,
        'parameters': f'--queue {args.queue} --cores {cores} {args_cube}'
        } )

# summary plots job
if hutch in ['rix','xcs','xpp']:
    job_defs.append( {
        'name': 'run_summary',
        'executable': executable_summaries,
        'trigger': 'RUN_PARAM_IS_VALUE',
        'run_param_name' : run_param_name,
        'run_param_value' : 'done',
        'location': location,
        'parameters': f'--queue {queue_summaries} {args_summaries}'
    } )



for job_def in job_defs:
    krbheaders = KerberosTicket('HTTP@pswww.slac.stanford.edu').getAuthHeaders()
    ws_url = 'https://pswww.slac.stanford.edu/ws-kerb/lgbk/lgbk/{}/ws/create_update_workflow_def'.format(exp)
    r = requests.post(ws_url, headers=krbheaders, json=job_def)
    r.raise_for_status()
    print('\nJOB CREATION LOG: {}'.format(r.json()))
