"""
Getting data from the elog run tables
This example returns all the timetool calibration runs, checks if the files are on tape 
    and submits smalldata jobs if the files are not present
"""
import requests
import json
from krtc import KerberosTicket
from urllib.parse import urlparse
import glob
import os

ws_url = "https://pswww.slac.stanford.edu/ws-kerb/lgbk"
krbheaders = KerberosTicket("HTTP@" + urlparse(ws_url).hostname).getAuthHeaders()
r = requests.get("%s/lgbk/ws/experiments"%ws_url, headers=krbheaders)
all_experiments = r.json()["value"]

instrument = "XPP"
#instrument = "XCS"
paramList=["%s:SCAN:SCANVAR00"%instrument]

instrument_experiments = [ x for x in all_experiments  if x["instrument"] == instrument]
instrument_experiment_names = [ x['name'] for x in instrument_experiments ]

for exp_name in instrument_experiment_names:
    rr = requests.get("%s/lgbk/%s/ws/get_run_params_for_all_runs"%(ws_url, exp_name), params={"param_names":",".join(paramList) }, headers=krbheaders)
    run_param_values = rr.json()['value']
    for rv in run_param_values:
        runNum = rv["num"]
        scanVar = rv["params"].get(paramList[0],"N/A")
        if scanVar == "txt":
            presentXtc=glob.glob('/reg/d/psdm/%s/%s/xtc/*-r%04d*xtc'%(exp_name[:3], exp_name, int(runNum)))
            if len(presentXtc)==0:
                #print("{0}".format(exp_name))
                print("{0}, {1} - need to restore".format(exp_name, runNum))
                break
            else:
                #print("{0}, {1} - on disk".format(exp_name, runNum))
                h5Name='/reg/d/psdm/xpp/xpptut15/scratch/timetool_calib/%s_Run%03d.h5'%(exp_name, int(runNum))
                haveh5=glob.glob(h5Name)
                if len(haveh5)==0:
                    print("{0}, {1} - on disk, no h5".format(exp_name, runNum))
                    #os.system('./examples/smallDataRun -d timetool_calib/ -e %s -r %d -j 12'%(exp_name, int(runNum)))
                else:
                    print("{0}, {1} - on disk, have h5".format(exp_name, runNum))

            #break
        #print("{0}, {1}, {2}".format(exp_name, rv["num"], rv["params"].get(paramList[0],"N/A")))
