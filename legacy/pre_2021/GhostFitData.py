import sys
import numpy as np
import argparse
from matplotlib import pyplot as plt
import socket
import os
import time
import itertools
import h5py
from mpi4py import MPI
import psana
import tables

#sys.path.insert(1,'./smalldata_tools')
from smalldata_tools.SmallDataUtils import getUserData
from smalldata_tools.roi_rebin import ROIFunc
from smalldata_tools.DetObject import DetObject

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

if not args.exp:
    expname=(raw_input('provide the experiment name:'))
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


smdFile=tables.open_file('/reg/d/psdm/xcs/%s/hdf5/smalldata/%s_Run%03d.h5'%(expname, expname, int(run))).root
xon=smdFile.lightStatus.xray.read()
allfids = smdFile.fiducials.read()
alltimes = smdFile.event_time.read()
try:
    allepix10k = np.nansum(np.nansum(smdFile.epix10ka2m.azav.read(),axis=1),axis=1)
except:
    allepix10k = np.zeros_like(allfids)
allipm5=smdFile.ipm5.sum.read()
allevr162=smdFile.evr.code_162.read()

drop=(xon==0)
drop[0]=0 #e cannot get the pre-shot for this one.
predrop=np.append(drop[1:],[0]).astype(bool)

offfids = allfids[drop] #this look ok
offtimes = alltimes[drop] 
preofffids = allfids[predrop] #this really does not.
preofftimes = alltimes[predrop]

#if expname=='xcsx35617':
#    Ped_epix10k = getPed_epix10k_xcsx35617(int(run))
#elif expname=='xcslt8717':
#    Ped_epix10k = getPed_epix10k_xcslt8717(int(run))

dsenv = (psana.DataSource('exp=%s:run=%d:smd'%(expname, int(run)))).env()
runIdx = psana.DataSource('exp=%s:run=%d:idx'%(expname, int(run))).runs().next()
#the ghost fitting should actually make a great pedestal!
if expname=='xcslt8717':
    epixDet = DetObject('epix10ka2m', dsenv, int(run), common_mode=80)
    epixDetPre = DetObject('epix10ka2m', dsenv, int(run), common_mode=80)
    #epixDet = DetObject('epix10ka2m', dsenv, int(run), common_mode=81)
    #epixDetPre = DetObject('epix10ka2m', dsenv, int(run), common_mode=81)
    #if Ped_epix10k is not None:
    #    epixDet.setPed(Ped_epix10k)
    #    epixDetPre.setPed(Ped_epix10k)
elif expname=='xcslt4017':
    #I'm not so sure this is what we want here. Maybe we did used method 81.
    epixDet = DetObject('epix10ka2m', dsenv, int(run), common_mode=80)
    epixDetPre = DetObject('epix10ka2m', dsenv, int(run), common_mode=80)
else:
    epixDet = DetObject('epix10ka2m', dsenv, int(run), common_mode=80)
    epixDetPre = DetObject('epix10ka2m', dsenv, int(run), common_mode=80)
    #epixDet = DetObject('epix10ka2m', dsenv, int(run), common_mode=81)
    #epixDetPre = DetObject('epix10ka2m', dsenv, int(run), common_mode=81)
    #if Ped_epix10k is not None:
    #    epixDet.setPed(Ped_epix10k)
    #    epixDetPre.setPed(Ped_epix10k)
epixDet.addFunc(ROIFunc(ROI=[0,1e6], writeArea=True, name='full')) #works.
epixDetPre.addFunc(ROIFunc(ROI=[0,1e6], writeArea=True, name='full')) #works.

#ped=np.ones([16,300,300])
ped = epixDet.ped
#if we are using official corrections, then we need 1st of pedestals to get correct shape
if ped.ndim == 4: 
    ped = ped[0]
moduleShape=[offfids.shape[0]]
for dimShp in ped[0].shape: #get module shape as we will create one dataset per module
    moduleShape.append(dimShp)
moduleShape = tuple(moduleShape)

#calculate information for whole events to be added to hdf5
prepredrop=np.append(predrop[1:],[0]).astype(bool)
epix10k_drop = allepix10k[drop]
epix10k_pre = allepix10k[predrop]
epix10k_prepre = allepix10k[prepredrop]
xon_drop = xon[drop]
xon_pre = xon[predrop]
xon_prepre = xon[prepredrop]
ipm5_drop = allipm5[drop]
ipm5_pre = allipm5[predrop]
ipm5_prepre = allipm5[prepredrop]
evr162_drop = allevr162[drop]
evr162_pre = allevr162[predrop]
evr162_prepre = allevr162[prepredrop]

if dirname!='':
    outFileName='%s/GhostEvents_%s_Run%03d.h5'%(dirname,expname,run)
else:
    outFileName='/reg/d/psdm/xcs/%s/scratch/GhostEvents_%s_Run%03d.h5'%(expname,expname,run)

fout = h5py.File(outFileName, "w",driver='mpio',comm=MPI.COMM_WORLD)
fout.create_dataset('epix10k_off',data=epix10k_drop)
fout.create_dataset('epix10k_pre',data=epix10k_pre)
fout.create_dataset('epix10k_prepre',data=epix10k_prepre)
fout.create_dataset('xon_off',data=xon_drop)
fout.create_dataset('xon_pre',data=xon_pre)
fout.create_dataset('xon_prepre',data=xon_prepre)
fout.create_dataset('ipm5_off',data=ipm5_drop)
fout.create_dataset('ipm5_pre',data=ipm5_pre)
fout.create_dataset('ipm5_prepre',data=ipm5_prepre)
fout.create_dataset('evr162_off',data=evr162_drop)
fout.create_dataset('evr162_pre',data=evr162_pre)
fout.create_dataset('evr162_prepre',data=evr162_prepre)

offModuleDatasets=[]
preModuleDatasets=[]
for imodule,module in enumerate(ped):
    offModuleDatasets.append(fout.create_dataset('epix10ka2m_off_%02d'%imodule,moduleShape))
    preModuleDatasets.append(fout.create_dataset('epix10ka2m_pre_%02d'%imodule,moduleShape))

print('we have %d events for ghost fitting'%offfids.shape[0])
for ievt, evtfid, evttime, preevtfid, preevttime in itertools.izip(itertools.count(), offfids, offtimes, preofffids, preofftimes):
    print('get events...', evttime, evtfid, '  ---  ', preevttime, preevtfid)
    evtt = psana.EventTime(evttime,evtfid)
    evt = runIdx.event(evtt)
    preevtt = psana.EventTime(preevttime,preevtfid)
    preevt = runIdx.event(preevtt)
    try:
        epixDet.getData(evt)
        epixDet.processFuncs()
        epixDict=getUserData(epixDet)

        epixDetPre.getData(preevt)
        epixDetPre.processFuncs()
        epixDictPre=getUserData(epixDetPre)
        
        for offDs, preDs, module, modulePre in zip( offModuleDatasets, preModuleDatasets, epixDict['full_area'], epixDictPre['full_area']):
            offDs[ievt] = module
            preDs[ievt] = modulePre
        #print 'data off: ',offModuleDatasets[0][ievt][5][5],' pre ',preModuleDatasets[0][ievt][5][5]
    except:
        pass

fout.close()
#t1 = time.time()
#print 'this took %g seconds '%(t1-t0)
