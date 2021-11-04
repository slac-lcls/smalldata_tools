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
FFB_BASE = Path("/cds/data/drpsrcf/{}/{}/scratch".format(hutch, exp))
PSANA_BASE = Path("/cds/data/psdm/{}/{}".format(hutch, exp))

if 'ffb' in args.queue:
    location = 'SRCF_FFB'
    cores = 60
    executable = str(FFB_BASE / 'smalldata_tools/arp_scripts/submit_smd.sh')
    executable_cube = str(FFB_BASE / 'smalldata_tools/arp_scripts/cubeRun.sh')
    trigger = 'START_OF_RUN'
else:
    location = 'SLAC'
    cores = 12
    executable = str(PSANA_BASE / 'results/smalldata_tools/arp_scripts/submit_smd.sh')
    executable_cube = str(PSANA_BASE / 'results/smalldata_tools/arp_scripts/cubeRun.sh')
    trigger = 'ALL_FILES_TRANSFERRED'

job_defs = []
job_defs.append( {
    'name': 'smd',
    'executable': executable,
    'trigger': trigger,
    'location': location,
    'parameters': '--queue {} --norecorder --postRuntable --cores {} --wait'.format(args.queue, cores)
    } )

if args.cube>0:
    job_defs.append( {
        'name': 'cube',
        'executable': executable_cube,
        'trigger': 'MANUAL',
        'location': location,
        'parameters': '--queue {} --norecorder --postRuntable --cores {} --wait'.format(args.queue, cores)
        } )

for job_def in job_defs:
    krbheaders = KerberosTicket('HTTP@pswww.slac.stanford.edu').getAuthHeaders()
    ws_url = 'https://pswww.slac.stanford.edu/ws-kerb/lgbk/lgbk/{}/ws/create_update_workflow_def'.format(exp)
    r = requests.post(ws_url, headers=krbheaders, json=job_def)
    r.raise_for_status()
    print('\nJOB CREATION LOG: {}'.format(r.json()))
