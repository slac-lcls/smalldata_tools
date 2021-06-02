import requests
import json
import shutil
from krtc import KerberosTicket
import argparse
from pathlib import Path
import os

parser = argparse.ArgumentParser()
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--psana', help='Setup psana job', type=bool, default=0)
parser.add_argument('--ffb', help='Setup ffb job', type=bool, default=0)
parser.add_argument('--queue', help='Queue on the ffb', type=str, default='ffbh3q')
args = parser.parse_args()

exp = args.experiment
hutch = exp[:3].lower()
FFB_BASE = Path("/cds/data/drpsrcf/{}/{}/scratch".format(hutch, exp))
PSANA_BASE = Path("/cds/data/psdm/{}/{}".format(hutch, exp))

if args.psana and not args.ffb:
    job_def = {
        'name': 'smd',
        'executable': str(PSANA_BASE / 'results/smalldata_tools/examples/smalldata_producer.py'),
        'trigger': 'MANUAL',
        'location': 'SLAC',
        'parameters': '--queue psanaq --norecorder --postRuntable --cores 12 --wait' 
        }
elif args.ffb:
    job_def = {
        'name': 'smd',
        'executable': str(FFB_BASE / '/smalldata_tools/examples/smalldata_producer.py'),
        'trigger': 'START_OF_RUN',
        'location': 'SRCF_FFB',
        'parameters': '--queue {} --norecorder --postRuntable --cores 60 --wait'.format(args.queue)
        }
else:
    job_def = None

if job_def is not None:
    krbheaders = KerberosTicket('HTTP@pswww.slac.stanford.edu').getAuthHeaders()
    ws_url = 'https://pswww.slac.stanford.edu/ws-kerb/lgbk/lgbk/{}/ws/create_update_workflow_def'.format(exp)
    r = requests.post(ws_url, headers=krbheaders, json=job_def)
    r.raise_for_status()
    print('\nJOB CREATION LOG: {}'.format(r.json()))