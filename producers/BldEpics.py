#!/usr/bin/env python
import tables
import numpy as np
import holoviews as hv
hv.extension('bokeh')
import panel as pn
pn.extension()
import os
import argparse
import sys
import logging
try:
    basestring
except NameError:
    basestring = str
fpath=os.path.dirname(os.path.abspath(__file__))
fpathup = '/'.join(fpath.split('/')[:-1])
try:
    fpath = os.environ.get('MYDIR', fpathup).replace('/arp_scripts','')
except:
    fpath = fpathup
sys.path.append(fpath)
from smalldata_tools.utilities import rebin
from smalldata_tools.utilities import image_from_dxy

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument('--run', help='run', type=str, default=os.environ.get('RUN_NUM', ''))
parser.add_argument('--experiment', help='experiment name', type=str, default=os.environ.get('EXPERIMENT', 'meclx3519'))
parser.add_argument('--epicslist', default=None, help="path&name of text file with list of PVs (one/line)")
args = parser.parse_args()
logger.debug('Args to be used for Bld&EPICS extraction: {0}'.format(args))

expname = args.experiment
run = int(args.run)

dat = tables.open_file('/cds/data/drpsrcf/%s/%s/scratch/hdf5/smalldata/%s_Run%04d.h5'%(expname[:3],expname,expname,run)).root

tiffdir = '/cds/data/psdm/%s/%s/scratch/run%d'%(expname[:3],expname,run)

def SafeFloat(self, string, precision=None):
    try:
        floatstr = float(string)
        if precision:
            if precision == 20:
                res = '%.20f' % floatstr
            elif precision == 15:
                res = '%.15f' % floatstr
            else:
                res = '%.10f' % floatstr
        else:
            res = '%.5f' % floatstr
    except (ValueError, TypeError):
        if string == '':
            res = '_'
        else:
            res = string.replace(' ','_')
        
    return res

def get_bld_data(key):
    if key == 'FEEGasDetEnergy':
        return 0.25*(dat.gas_detector.f_11_ENRC.read()+dat.gas_detector.f_12_ENRC.read()+dat.gas_detector.f_21_ENRC.read()+dat.gas_detector.f_22_ENRC.read())
    elif key == 'MEC-XT2-BMMON-02':    
        return dat.xt2_ipm2.sum.read()
    elif key == 'MEC-XT2-BMMON-03':    
        return dat.xt2_ipm3.sum.read()
    elif key == 'event':
        return np.arange(dat.fiducials.shape[0])
    else: 
        return None



###
# code for bld
###
mecblds = ['event', 'MEC-XT2-BMMON-02', 'MEC-XT2-BMMON-03', 'FEEGasDetEnergy']
bldfilename = tiffdir + '/run%d_Bld.txt' % run
if not os.path.isdir(tiffdir):
    os.mkdir(tiffdir)
bldfile = open(bldfilename, 'w')
        
blddict = {k: get_bld_data(k) for k in mecblds if get_bld_data(k) is not None }
bldfile.write(' '.join([k for k in blddict])+'\n')

for ievt,evt in enumerate(blddict['event']):
    line=''
    for k in blddict:
        line+='%g '%blddict[k][ievt]
    print(line)
    bldfile.write(line+'\n')
bldfile.close()


###
#epics section
###
epicsfilename = tiffdir + '/run%s_epicsarch.txt' % run

if args.epicslist is None:
    pvs = ['SIOC:SYS0:ML00:AO627',    # XFEL Energy
           'MEC:NOTE:LAS:NSDELAY',    # NS Laser Delay
           'MEC:NS1:MMS:02.RBV',      # AB Waveplate position
           'MEC:NS1:MMS:01.RBV',      # EF Waveplate position
           'MEC:LAS:MMN:30.RBV',      # GH Waveplate position
           'MEC:LAS:MMN:29.RBV',      # IJ Waveplate position
           'MEC:PPL:MMN:09.RBV',      # Target Y position
           'MEC:USR:MMS:17.RBV',      # Target X position
           'MEC:LAS:DDG:05:aDelaySI', # Streak 1 delay
           'MEC:LAS:DDG:05:eDelaySI'] # Streak 2 delay
else:
    with open(args.epicslist, 'r') as file:
        pvs = [ line[:-1] for line in file]

epData = getattr(dat,'epicsAll',None)
epDataCfg = getattr(dat.UserDataCfg,'epicsAll',None)
if epData is None:
    epData = getattr(dat,'epics',None)
    epDataCfg = getattr(dat.UserDataCfg,'epics',None)
if epData is None:
    print('no epics in h5')
    sys.exit()

aliases = [alias.decode('utf-8') for alias in getattr(epDataCfg,'PVlist',None).read()]
PVnames = getattr(epDataCfg,'PVlist_PV',None)
if PVnames is not None:
    PVnames = [pv.decode('utf-8') for pv in PVnames.read()]
header = [pv for pv in pvs if (pv in aliases or pv in PVnames) ]

data={tpv: getattr(epData, tpv).read() if tpv in aliases else getattr(epData, aliases[PVnames.index(tpv)]).read() for tpv in header}

epicsfile = open(epicsfilename, 'w')
epicsfile.write(' '.join(header)+'\n')

for ievt in np.arange(dat.fiducials.shape[0]):
    epicsfile.write(' '.join(['%g'%val[ievt] for key,val in data.items()])+'\n')
epicsfile.close()
