#!/usr/bin/env python

"""
Post run params to the elog
"""
import os
import sys
import json
import argparse
import requests
from krtc import KerberosTicket
from urllib.parse import urlparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Example for posting run time parameters to the elog')
    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('--url', default="https://pswww.slac.stanford.edu/ws-kerb/lgbk/")
    parser.add_argument('experiment_name', help="The name of the experiment")
    parser.add_argument('run_number', help="run number of the experiment")
    args = parser.parse_args()
    ws_url = args.url + "/run_control/{0}/ws/add_run_params".format(args.experiment_name)
    krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
    runtable_data = {
      "usrData_1": "1.0",
      "usrData_2": "2.0",
      "usrData_3": "3.0",
      "usrData_4": "4.0"
    }
    r = requests.post(ws_url, headers=krbheaders, params={"run_num": args.run_number}, json=runtable_data)
    print(r)


