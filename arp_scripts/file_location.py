#!/usr/bin/env python

import requests
import json
import argparse
import logging
import sys


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def get_file_location(exp, run, location):
    ws_url = f"https://pswww.slac.stanford.edu/ws/lgbk/lgbk/{exp}/ws/{run}/files_for_live_mode_at_location"

    # check sdf location
    r = requests.get(ws_url, params={"location": location})
    r.raise_for_status()
    on_sdf = r.json()['value']['all_present']
    
    if on_sdf:
        logger.info("Data files found on s3df.")
        print("/sdf/data/lcls/ds/")
    else:
        print('/sdf/data/lcls/drpsrcf/ffb')
        # test ffb location (comment for now, we might not need that)
        #r = requests.get(ws_url, params={"location": "SRCF_FFB"})
        #r.raise_for_status()
        #on_ffb = r.json()['value']['all_present']
        #if on_ffb:
        #    logger.info('Data files found on ffb.')
        #    print("/sdf/data/lcls/drpsrcf/ffb")
        #else:
        #    logger.info("Elog is not aware of data neither on S3DF or FFB.")
        #    logger.info("Try by ignoring stream 80.")
        #    ignore_stream(r.json())
    return 0


def ignore_stream(r, stream=80):
    s80_present = True
    count = 0 # count files that are present on the drp
    for file in r['value']['files']:
        p = file['path']
        is_present = file['is_present']
        if not is_present:
            if not f"s{stream}" in p:
                count+=1
            else:
                s80_present = False
    # if only the s80 is missing, go to ffb location
    if count==0 and s80_present is False:
        print("/sdf/data/lcls/drpsrcf/ffb")
                
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', '-r', type=int, action='store', help="Run number")
    parser.add_argument('--experiment', '-e', type=str, action='store', help="Experiment name")
    parser.add_argument('--location', '-l', help="The DM location name", default="S3DF")
    parser.add_argument('--verbose', '-v', action='store_true', help="Turn on verbose logging")
    args = parser.parse_args()

    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    get_file_location(args.experiment, args.run, args.location)

