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
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import TheilSenRegressor
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import HuberRegressor
from future.utils import iteritems
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run",type=int)
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--dir", help="output directory")
parser.add_argument("--stat", help="keep stat", action='store_true')
parser.add_argument("--Huber135", help="also use Huber w/default settings", action='store_true')
parser.add_argument("--Huber", help="also use Huber epsilon 2.5", action='store_true')
parser.add_argument("--RANSAC", help="also use RANSAC",type=int)
parser.add_argument("--TheilSen", help="also use TheilSen",type=int)
parser.add_argument("--TheilSen150", help="also use TheilSen150",type=int)
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

#ghostDir='/reg/neh/home/snelson/gitMaster_smalldata_tools/smalldata_tools/output'
ghostDir='/reg/d/psdm/xcs/%s/scratch/GhostImageFiles/'%(expname)
ghostFile=tables.open_file('%s/GhostEvents_%s_Run%03d.h5'%(ghostDir,expname, int(run)))#.root
smdFile=tables.open_file('/reg/d/psdm/xcs/%s/hdf5/smalldata/%s_Run%03d.h5'%(expname, expname, int(run))).root

if dirname!='':
    outFileName='%s/GhostFits_%s_Run%03d.h5'%(dirname,expname,run)
else:
    outFileName='/reg/d/psdm/xcs/%s/scratch/GhostFits_%s_Run%03d.h5'%(expname,expname,run)
#fout = h5py.File(outFileName, "w")
fout = h5py.File(outFileName, "w",driver='mpio',comm=MPI.COMM_WORLD)
print('opened output file: ',outFileName)

resultsDict={'times':{}, 'score':{},'p0':{},'p1':{}}
p0_dset = {}
p1_dset = {}
time_dset = {}
score_dset = {}
p0_array = {}
p1_array = {}
time_array = {}
score_array = {}

oData = ghostFile.get_node('/epix10ka2m_off_00')
nRows=oData.shape[1]
nCols=oData.shape[2]

estimators={'OLS': LinearRegression(fit_intercept=True)}
if args.Huber:
    estimators['Huber135']=HuberRegressor(epsilon=2.5) #default setting
if args.Huber135:
    estimators['Huber135']=HuberRegressor(epsilon=1.35) #default setting
#force RANSAC to keep high fraction of events
if args.RANSAC:
    estimators['RANSAC']=RANSACRegressor(random_state=args.RANSAC, min_samples=0.7) 
#calculate TheilSen on a larger number of events, limit the # of iterations to make it run ok.
if args.TheilSen:
    estimators['TheilSen']=TheilSenRegressor(random_state=args.TheilSen, n_subsamples=int(oData.shape[0]*0.9), max_subpopulation=200)
    estimators['TheilSen100']=TheilSenRegressor(random_state=args.TheilSen, n_subsamples=int(oData.shape[0]*0.9), max_subpopulation=100)
if args.TheilSen150:
    estimators['TheilSen150']=TheilSenRegressor(random_state=args.TheilSen, n_subsamples=int(oData.shape[0]*0.9), max_subpopulation=150)

#maxModules=1
#maxRows=10
#maxCols=10
maxModules=16
maxRows=nRows
maxCols=nCols
keepStat=False
if args.stat:
    keepStat=args.stat

nanArray=np.empty((16,nRows,nCols))
nanArray.fill(np.nan)
p0ModuleArray={}
for name,estimator in iteritems(estimators):
    p0_dset[name] = fout.create_dataset('epix10ka2m_%s_p0'%name,data=nanArray.copy())
    p1_dset[name] = fout.create_dataset('epix10ka2m_%s_p1'%name,data=nanArray.copy())
    p0_array[name] = nanArray[0].copy()
    p1_array[name] = nanArray[0].copy()
    if keepStat:
        time_dset[name] = fout.create_dataset('epix10ka2m_%s_time'%name,data=nanArray.copy())
        score_dset[name] = fout.create_dataset('epix10ka2m_%s_score'%name,data=nanArray.copy())
        time_array[name] = nanArray[0].copy()
        score_array[name] = nanArray[0].copy()

try:
    pixelGain = smdFile.UserDataCfg.epix10ka2m.pixelGain.read()
except:
    pixelGain = np.ones_like(nanArray)

for iModule in range(max(maxModules, size)):
    #use MPI to treat modules in parallel.
    if maxModules<=size and not (iModule==rank):
        continue

    offData = ghostFile.get_node('/epix10ka2m_off_%02d'%iModule).read()
    preData = ghostFile.get_node('/epix10ka2m_pre_%02d'%iModule).read()
    moduleGain = pixelGain[iModule]

    for iRow in range(0,maxRows):
        if iRow%5==0: print('fitting pixels in row %d of module %d '%(iRow,iModule))
        for iCol in range(0,maxCols):
            gain = moduleGain[iRow,iCol]
            poffData = offData[:,iRow,iCol]/gain
            ppreData = preData[:,iRow,iCol]/gain

            if (~(np.isnan(poffData))).astype(int).sum()==0:
                continue

            PppreData=ppreData[:,np.newaxis]
            line_x = np.array([0, np.nanmax(ppreData)])

            for name,estimator in iteritems(estimators):
                try:
                    t0=time.time()
                    estimator.fit(PppreData, poffData)
                    if keepStat:
                        time_array[name][iRow][iCol]=time.time()-t0
                    try:
                        p0_array[name][iRow][iCol]=estimator.intercept_*gain
                        p1_array[name][iRow][iCol]=estimator.coef_[0]
                    except:
                        try:
                            p0_array[name][iRow][iCol]=estimator.estimator_.intercept_*gain
                            p1_array[name][iRow][iCol]=estimator.estimator_.coef_[0]
                        except:
                            pass
                        pass
                    if keepStat:
                        try:
                            score_array[name][iRow][iCol]=estimator.score(PppreData, poffData)
                        except:
                            pass
                except:
                    pass

    for name,estimator in iteritems(estimators):
        p0_dset[name][iModule] = p0_array[name]
        p1_dset[name][iModule] = p1_array[name]
        if keepStat:
            time_dset[name][iModule] = time_array[name]
            score_dset[name][iModule] = score_array[name]
    
fout.close()
