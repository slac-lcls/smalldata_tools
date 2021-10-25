import numpy as np
import logging
import time
import itertools

import psana
from smalldata_tools.DetObject import DetObject
from smalldata_tools.SmallDataUtils import getUserData
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

logger = logging.getLogger()


DUMMY = False
def lengthy_computation(*args):
    time.sleep(10+np.random.randint(10))
    return 1


def bin_distribution(bins_info, func=None):
    """ Distribute bin analysis to worker on a per bin basis. Generally should be run 
    from rank 0 in MPI.
    Workflow: if a worker is not busy, it gets the first bin with a 'not done'
    or 'not in progress' status. Workers' status is stored in array 'working' and the bin 
    processing status is in 'bin_done'.
    
    Args:
        bins_info: bin-dependent parameter to pass to the worker. Typically a list of idx
            corresponding to the bin. For cube: list of fiducials and event time array.
        func: what to do with returned results. Must accept bin_idx and data as kw argument. Generally
            a function to write data to file.
    """
    nBins = len(bins_info)
    logger.info('Total number of bins: {}'.format(nBins))
    working = np.zeros((size-1, 2)) # rank currently working and bin_idx on which it is working
    bin_status = np.zeros(nBins) # whether bin has been processed or not. 0.5=in progress, 1=done
    while(not np.all(bin_status==1)):
        for worker_id, at_work in enumerate(working):
            if at_work[0]==0:
                try:
                    bin_idx = np.argwhere(bin_status==0)[0][0]
                except Exception as e:
                    continue
                job_info = {'bin_idx': bin_idx, 'info': []}
                for info in bins_info[bin_idx]:
                    job_info['info'].append(info)
                logger.debug('Send job to rank {}: {}'.format(worker_id+1, len(job_info)))
                bin_status[bin_idx] = 0.5
                comm.send(job_info, dest=worker_id+1) # rank 0 does not do jobs
                working[worker_id,:] = [1, bin_idx]
#         logger.debug('Working: {}'.format(working))
        
        t1 = time.time()
        job_done = comm.recv(source=MPI.ANY_SOURCE)
        t2 = time.time()
        bin_idx = job_done['idx']
        logger.info('Bin {} received on rank 0 after {}s.'.format(bin_idx, t2-t1))
        if func is not None:
            output = func(data=job_done['data'], bin_idx=bin_idx)
        
        bin_status[bin_idx] = 1
        worker_id = np.argwhere(working[:,1]==bin_idx)[0][0]
        working[worker_id,0] = 0
        logger.debug('New bin status: {}'.format(bin_status))
    logger.info('**** DONE ****')
    
    """ Closing worker's while loop """
    for worker_id, at_work in enumerate(working):
        if at_work[0]==0:
            comm.send('done', dest=worker_id+1)
    return


class BinWorker(object):
    def __init__(self, run, expname):
        """ Make index-based datasource and get area detector info from root rank.
        """
        logger.debug('Initializing worker {}.'.format(rank))
        self.run = int(run)
        self.expname = expname
        bcast_var = None
        dsname = comm.bcast(bcast_var, root=0)
        print(dsname)
        print('********** Start setup.')
        t0 = time.time()
        self.dsIdx = psana.DataSource(str(dsname))
        print('********** Datasource on rank {}: {}'.format(rank, time.time()-t0))
        self.dsIdxRun = next(self.dsIdx.runs())
        print('********** Setup on rank {}: {}'.format(rank, time.time()-t0))
        bcast_var = None
        self.targetVarsXtc = comm.bcast(bcast_var, root=0) # get det info from rank 0
        logger.debug('Detectors info received on rank {}. {}'.format(rank, self.targetVarsXtc))
        self.dets = []
        for det_info in self.targetVarsXtc:
            try:
                detname = det_info['source']
                cm = det_info.get('common_monde',None)
                det = DetObject(detname, self.dsIdx.env(), self.run, common_mode=cm, name=detname)
                self.dets.append(det)
                # add detector analysis function (only support full image for now)
                det.addFunc(ROIFunc(writeArea=True))
            except Exception as e:
                print('Could not make detector {}. Abort'.format(detname))
                print(e)
                comm.Abort()
        logger.info('Rank {} has datasource and detectors.'.format(rank))
        return
    
    def work(self):
        done = 0
        while(not done):
            job_info = comm.recv(source=0)
            if job_info=='done':
                logger.debug('Rank {} done.'.format(rank))
                done = 1
                continue
            logger.info('Rank {} got some job. Shots in bin #{}: {}'\
                        .format(rank, job_info['bin_idx'], len(job_info['info'][0])))
            if DUMMY: # just to test things
                out = lengthy_computation()
            else:
                out = self.process_bin(job_info['info'], job_info['bin_idx']) # out[0]: summed_data, out[1]: n_in_bin
            job_done = {'idx': job_info['bin_idx'], 'data': out}
            comm.send(job_done, dest=0)
        logger.debug('Rank {} out of while loop.'.format(rank))
        return
    
    def process_bin(self, bin_info, bin_idx):
        """ bin_info[0]: fiducials, bin_info[1]: evttime
        """
        summed_data = {}
        n_in_bin = {}
        # prepare data for each det
        for det in self.dets:
            summed_data[det._name] = 0
            n_in_bin[det._name] = 0
        # get detectors data for the bin
        for ievt, evtfid, evttime in zip(itertools.count(), bin_info[0], bin_info[1]):
            evtt = psana.EventTime(int(evttime),int(evtfid))
            evt = self.dsIdxRun.event(evtt)
            data = self.process_area_dets(evt)
            for detname in data.keys():
                n_in_bin[detname]+=1
                summed_data[detname]+=data[detname]
        # make full image for each det
        for det, thisDetDict in zip(self.dets, self.targetVarsXtc):
            if isinstance(summed_data[det._name],int):
                logger.info('No data in bin {}.'.format(bin_idx))
                summed_data[det._name] = np.nan
                continue
            if hasattr(det, 'x'):
                summed_data[det._name] = det.det.image(self.run, summed_data[det._name])
#             if 'image' in thisDetDict:
#                 if thisDetDict['image']==True:
#                     summed_data[thisDetDict['source']] = \
#                         det.det.image(self.run, summed_data[thisDetDict['source']])
        return summed_data, n_in_bin
        
    def process_area_dets(self, evt):
        """ Process area detectors defined in self.targetVarXtc for a given events. """
        det_data = {}
        for det, thisDetDict in zip(self.dets, self.targetVarsXtc):
            try:
                det.getData(evt)
                det.processFuncs()
                thisDetDataDict = getUserData(det)
                img = None
                for key in thisDetDataDict.keys():
                    if key=='full_area':
                        img = thisDetDataDict[key]
                    elif key.find('ROI')>=0:
                        img = thisDetDataDict[key]
                if img is None:
                    print('Problem with getting detector area data.')
                    continue
                if 'thresADU' in thisDetDict:
                    img[img<thisDetDict['thresADU']] = 0
                elif 'thresRms' in thisDetDict:
                    img[img<thisDetDict['thresRms']*det.rms] = 0

                det_data[det._name] = img # can onky handle full area ROIFunc for now
                
                # calculate variance (see ... for ref)
                
#                     if not (key=='full_area' or key.find('ROI')>=0 or key.find('photon_img')>=0):
#                         continue
#                     if (key=='full_area' or key.find('ROI')>=0):
#                         if 'thresADU' in thisDetDict:
#                             thisDetDataDict[key][thisDetDataDict[key]<thisDetDict['thresADU']]=0
#                         elif 'thresRms' in thisDetDict:
#                             thisDetDataDict[key][thisDetDataDict[key]<thisDetDict['thresRms']*det.rms]=0
#                         dArray[ib%bins_per_job]=dArray[ib%bins_per_job]+thisDetDataDict[key]
#                     else: #if key.find('photon_img')
#                         dIArray[ib%bins_per_job]=dIArray[ib%bins_per_job]+thisDetDataDict[key]

#                     x = thisDetDataDict[key]
#                     oldM = dMArray
#                     dMArray = dMArray + (x-dMArray)/(ievt+1)
#                     dSArray = dSArray + (x-dMArray)*(x-oldM)
            except Exception as e:
                print('Failed to get data for this event for det {}.\n{}'.format(det._name, e))
        return det_data
    
    
class mpi_bin_data(object):
    def __init__(self, bin_idx):
        self.bin_idx = bin_idx