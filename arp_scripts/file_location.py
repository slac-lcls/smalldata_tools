#!/bin/python

import requests
import json
import argparse
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

exp = 'xppc00121'
run = 120

parser = argparse.ArgumentParser()
parser.add_argument('--run', '-r', type=int, action='store', help="Run number")
parser.add_argument('--experiment', '-e', type=str, action='store', help="Experiment name")

def get_file_location(exp, run):
    ws_url = f"https://pswww.slac.stanford.edu/ws/lgbk/lgbk/{exp}/ws/{run}/files_for_live_mode_at_location"

    # check sdf location
    r = requests.get(ws_url, params={"location": "S3DF"})
    r.raise_for_status()
    on_sdf = r.json()['value']['all_present']
    
    if on_sdf:
        logger.info("Data files found on s3df.")
        print("/sdf/data/lcls/ds/")
    else:
        # test ffb location
        r = requests.get(ws_url, params={"location": "SRCF_FFB"})
        r.raise_for_status()
        on_ffb = r.json()['value']['all_present']
        if on_ffb:
            logger.info('Data files found on ffb.')
            print("/sdf/data/lcls/drpsrcf/ffb")
        else:
            logger.info("Elog is not aware of data neither on S3DF or FFB.")
    return 0


if __name__ == "__main__":
    args = parser.parse_args()
    get_file_location(args.experiment, args.run)

