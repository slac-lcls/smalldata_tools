#!/usr/bin/env python

import argparse 
import os
from datetime import datetime
from smalldata_tools.epicsarchive import EpicsArchive
import logging
import requests
import json
from krtc import KerberosTicket

def PV2JSON(n):
    return n.replace(".", u"\uFF0E").replace("$", u"\uFF04")

def add_run_to_list(v, rlist, allruns, exp):
    v = int(v)
    if v in allruns:
        rlist.append(v)
    else:
        logger.warn("%d is not a run in experiment %s!" % (v, exp))

parser = argparse.ArgumentParser()
parser.add_argument('--runs', 
                    help='comma-separated run number or ranges, or ALL', 
                    type=str, 
                    default='ALL')
parser.add_argument('--experiment', 
                    help='experiment name', 
                    type=str, 
                    default=os.environ.get('EXPERIMENT', ''))
parser.add_argument('--pv', 
                    help='PV list', 
                    type=str, 
                    default="")
parser.add_argument('--file', 
                    help='PV file (one PV per line)', 
                    type=str, 
                    default="")
parser.add_argument('--url',
                    type=str,
                    default="https://pswww.slac.stanford.edu/ws-kerb/lgbk")

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    args = parser.parse_args()
    logger.debug('Args to be used for update run tables: {0}'.format(args))

    # Build the PV list
    if args.pv != "":
        pvlist = [PV2JSON(l.strip()) for l in args.pv.split(",")]
    else:
        pvlist = []
    if args.file != "":
        with open(args.file) as f:
            pvlist = pvlist + [PV2JSON(l.strip()) for l in f.readlines()]

    logger.debug('pvlist: {0}'.format(pvlist))

    archive = EpicsArchive()
    krbheaders = KerberosTicket("HTTP@pswww.slac.stanford.edu").getAuthHeaders()
    rundata = {}
    allruns = []
    for r in requests.get("%s/lgbk/%s/ws/runs" % (args.url, args.experiment), 
                          headers=krbheaders).json()['value']:
        n = int(r['num'])
        rundata[n] = r
        allruns.append(n)

    if args.runs == "ALL":
        rlist = allruns
    else:
        rlist = []
        r = args.runs.split(',')
        for v in r:
            if not "-" in v:
                add_run_to_list(v, rlist, allruns, args.experiment)
            else:
                vl = v.split("-")
                if len(vl) != 2:
                    raise Exception("Bad syntax for runs: %s" % args.runs)
                for v in range(int(vl[0]), int(vl[1])+1):
                    add_run_to_list(v, rlist, allruns, args.experiment)
    logger.debug('rlist: {0}'.format(rlist))

    for run in rlist:
        old_params = requests.get("%s/lgbk/%s/ws/runs/%d" % (args.url, args.experiment, run), 
                                  headers=krbheaders).json()['value']['params']

        start  = datetime.fromisoformat(rundata[run]['begin_time'])
        logger.debug("Run %d: start %s (%.9f)" % (run, rundata[run]['begin_time'], start.timestamp()))
        #print("old_params=%s" % old_params)

        new_params = {}
        for p in pvlist:
            if p in old_params.keys():
                logger.debug("%s is already a run param in run %d, skipping..." % (p, run))
                continue
            v = archive.get_point(PV=p, when=start, value_only=True)
            if v is not None:
                new_params[PV2JSON(p)] = v

        #print("new_params=%s" % new_params)
        if len(new_params) > 0:
            r = requests.post("%s/run_control/%s/ws/add_run_params" % (args.url, args.experiment),
                              params={"run_num": str(run)}, json=new_params, headers=krbheaders)
            print(r)
