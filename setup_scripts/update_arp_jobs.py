import requests
import json
import shutil
from krtc import KerberosTicket
import argparse
from pathlib import Path
import os

parser = argparse.ArgumentParser()
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
args = parser.parse_args()

exp = args.experiment
hutch = exp[:3].lower()
FFB_BASE = Path("/cds/data/drpsrcf/{}/{}/scratch".format(hutch, exp))
PSANA_BASE = Path("/cds/data/psdm/{}/{}".format(hutch, exp))

krbheaders = KerberosTicket('HTTP@pswww.slac.stanford.edu').getAuthHeaders()
ws_url = 'https://pswww.slac.stanford.edu/ws-kerb/lgbk/lgbk/{}/ws/workflow_definitions'.format(exp)
r = requests.get(ws_url, headers=krbheaders)
job_defs = r.json()['value']

for job in job_defs:
    if job['name']=='smd' and job['location']=='SRCF_FFB':
        print('MODIFYING JOB {}'.format(job))
        id = job['_id']
        job_def = {
            '_id': id,
            'name': 'smd',
            'executable': str(PSANA_BASE / 'results/smalldata_tools/arp_scripts/submit_smd.sh'),
            'trigger': 'MANUAL',
            'location': 'SLAC',
            'parameters': '--queue psanaq --norecorder --postRuntable --cores 12 --wait' 
        }
        ws_url = 'https://pswww.slac.stanford.edu/ws-kerb/lgbk/lgbk/{}/ws/create_update_workflow_def'.format(exp)
        r = requests.post(ws_url, headers=krbheaders, json=job_def)
        r.raise_for_status()
        print('\nJOB MODIFICATION LOG: {}'.format(r.json()))

    elif job['name']=='cube' and job['location']=='SRCF_FFB':
        print('MODIFYING JOB {}'.format(job))
        id = job['_id']
        job_def = {
            '_id': id,
            'name': 'cube',
            'executable': str(PSANA_BASE / 'results/smalldata_tools/arp_scripts/cubeRun.sh'),
            'trigger': 'MANUAL',
            'location': 'SLAC',
            'parameters': '--queue psanaq --norecorder --postRuntable --cores 12 --wait' 
        }
        ws_url = 'https://pswww.slac.stanford.edu/ws-kerb/lgbk/lgbk/{}/ws/create_update_workflow_def'.format(exp)
        r = requests.post(ws_url, headers=krbheaders, json=job_def)
        r.raise_for_status()
        print('\nJOB MODIFICATION LOG: {}'.format(r.json()))
