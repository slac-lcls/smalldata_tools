import sys
import numpy as np
import argparse
from matplotlib import pyplot as plt
import socket
import os
import psana 

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run",type=int)
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--dir", help="output directory")
parser.add_argument("--file", help="file name")
args = parser.parse_args()

if not args.run:
    run=int(raw_input('provide a run number of experiment %s:'%expname))
else:
    run=args.run

hutches=['amo','sxr','xpp','xcs','mfx','cxi','mec']
if not args.exp:
    hostname=socket.gethostname()
    foundHutch=False
    for ihutch in hutches:
        if hostname.find(ihutch)>=0:
            hutch=ihutch.upper()
            foundHutch=True
            break
    if not foundHutch and hostname.find('psusr')>=0:
        if hostname.find('psusr11')>=0:
            hutch='AMO'
        elif hostname.find('psusr12')>=0:
            hutch='SXR'
        elif hostname.find('psusr13')>=0:
            hutch='XPP'
        elif hostname.find('psusr21')>=0:
            hutch='XCS'
        elif hostname.find('psusr22')>=0:
            hutch='CXI'
        elif hostname.find('psusr23')>=0:
            hutch='MEC'
        elif hostname.find('psusr24')>=0:
            hutch='MFX'
        if hutch!='':
            foundHutch=True
    else:
        #then check current path
        path=os.getcwd()
        for ihutch in hutches:
            if path.find(ihutch)>=0:
                hutch=ihutch.upper()
                foundHutch=True
                break
        if not foundHutch:
            print 'I cannot figure out which hutch we are in, so cannot determine the current experiment'
            sys.exit()

    try:
        import RegDB.experiment_info
        expname=RegDB.experiment_info.active_experiment(hutch)[1]
    except:
        print 'could not determine experiment name, will quit'
        sys.exit()
else:
    expname=args.exp

dirname=''
if args.dir:
    dirname = args.dir
    if dirname.find('/')<=0:
        dirname+='/'
fname=''
if args.file:
    fname = args.file

#from smalldata_tools import SmallDataAna
sys.path.append('./smalldata_tools')

ds = psana.DataSource('exp=%s:run=%s'%(expname, run))
det = psana.Detector('epix10ka2m')
ecfg = ds.env().configStore().get(....)

