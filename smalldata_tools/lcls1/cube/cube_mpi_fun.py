import numpy as np
import logging
import time
import itertools

import psana
from smalldata_tools.DetObject import DetObject, DetObjectFunc
from smalldata_tools.SmallDataUtils import getUserData
from smalldata_tools.ana_funcs.roi_rebin import ROIFunc
from smalldata_tools.ana_funcs.photons import photonFunc
from smalldata_tools.ana_funcs.droplet import dropletFunc
from smalldata_tools.ana_funcs.droplet2Photons import droplet2Photons
from smalldata_tools.ana_funcs.azimuthalBinning import azimuthalBinning
from smalldata_tools.ana_funcs.azav_pyfai import azav_pyfai

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

logger = logging.getLogger()


DUMMY = False
def lengthy_computation(*args):
    time.sleep(10+np.random.randint(10))
    return 1


class Bin_distribution(object):
    """ Handles the distribution of the bin analysis to worker on a per bin basis. Generally should be run 
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
    
    def __init__(self, bins_info, file):
        self.bins_info = bins_info
        self.nBins = len(self.bins_info)
        logger.info('Total number of bins: {}'.format(self.nBins))
        self.working = np.zeros((size-1, 2)) # rank currently working and bin_idx on which it is working
        self.bin_status = np.zeros(self.nBins) # whether bin has been processed or not. 0.5=in progress, 1=done
        self.file = file
        return
    
    
    def distribute(self):
        while(not np.all(self.bin_status==1)):
            for worker_id, at_work in enumerate(self.working):
                if at_work[0]==0:
                    try:
                        bin_idx = np.argwhere(self.bin_status==0)[0][0]
                    except Exception as e:
                        continue
                    job_info = {'bin_idx': bin_idx, 'info': []}
                    for info in self.bins_info[bin_idx]:
                        job_info['info'].append(info)
                    logger.debug('Send job to rank {}: {}'.format(worker_id+1, len(job_info)))
                    self.bin_status[bin_idx] = 0.5
                    comm.send(job_info, dest=worker_id+1) # rank 0 does not do jobs
                    self.working[worker_id,:] = [1, bin_idx]
            # logger.debug('Working: {}'.format(self.working))

            t1 = time.time()
            job_done = comm.recv(source=MPI.ANY_SOURCE)
            t2 = time.time()
            bin_idx = job_done['idx']
            logger.info('Bin {} received on rank 0 after {}s.'.format(bin_idx, t2-t1))
            
            self.save_bin_to_h5(bin_idx, job_done['data'])
            
            self.bin_status[bin_idx] = 1
            worker_id = np.argwhere(self.working[:,1]==bin_idx)[0][0]
            self.working[worker_id,0] = 0
            logger.debug('New bin status: {}'.format(self.bin_status))
        logger.info('**** DONE ****')

        """ Send stop signal to workers """
        for worker_id, at_work in enumerate(self.working):
            if at_work[0]==0:
                comm.send('done', dest=worker_id+1)
        return
    
    
    def save_bin_to_h5(self, bin_idx, data):
        """ data[0]: summed_data, data[1]: n_in_bin, data[2]: proc_data
        """
        for detname in data[0].keys():
            # data
            dset_name = f'{detname}_data'
            dset = self.file[dset_name]
            dset[bin_idx] = data[0][detname]
            # n_in_bin
            dset_name = f'{detname}_nEntries'
            dset = self.file[dset_name]
            dset[bin_idx] = data[1][detname]
            # proc data
            if data[2]:
                for key in data[2][detname].keys():
                    dset_name = f'{detname}_{key}'
                    dat = data[2][detname][key]
                    shape = (self.nBins,)+dat.shape
                    dset = self.file.require_dataset(dset_name, shape=shape, dtype=float)
                    dset[bin_idx] = dat
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
        logger.info('********** Datasource on rank {}: {}s'.format(rank, time.time()-t0))
        self.dsIdxRun = next(self.dsIdx.runs())
        self.parse_detectors()
        logger.info('Rank {} has datasource and detectors.'.format(rank))
        print('********** Setup on rank {}: {}s'.format(rank, time.time()-t0))
        return
    
    
    def parse_detectors(self):
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
                # add full ROI analysis
                det.addFunc(ROIFunc(writeArea=True)) # add full image (always)
                
                # add other to process on the image of each bin
                if 'det_proc' in det_info.keys():
                    for func_args in det_info['det_proc']:
                        fname = func_args.pop('name')
                        print(f'Add DetObjectFunc {fname} as post-process.')
                        func = globals()[fname]
                        func = func(**func_args)
                        func._proc = False # not process every event
                        det.addFunc(func)
            except Exception as e:
                print('Could not make detector {}. Abort'.format(detname))
                print(e)
                comm.Abort()
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
                out = self.process_bin(job_info['info'], job_info['bin_idx']) 
                # INFO: out[0]: summed_data, out[1]: n_in_bin, out[3]: proc_data
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
        
        # get summed detectors data for the bin
        for ievt, evtfid, evttime in zip(itertools.count(), bin_info[0], bin_info[1]):
            evtt = psana.EventTime(int(evttime),int(evtfid))
            evt = self.dsIdxRun.event(evtt)
            data = self.process_event(evt)
            for detname, dat in data.items():
                if dat is None:
                    continue
                else:
                    n_in_bin[detname]+=1
                    summed_data[detname]+=dat
        
        # post-process on the summed data in the bin
        summed_data, proc_data = self.process_summed_bin(summed_data, bin_idx)
        return summed_data, n_in_bin, proc_data
    
    
    def process_event(self, evt):
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
                        img = thisDetDataDict[key].astype(float) # needed for detectors whose data are uint16 (Rayonix)
                    elif key.find('ROI')>=0:
                        img = thisDetDataDict[key].astype(float)
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
                det_data[det._name] = None
        return det_data
    
    
    def process_summed_bin(self, summed_data, bin_idx):
        proc_data = {}
        for det, thisDetDict in zip(self.dets, self.targetVarsXtc):
            if isinstance(summed_data[det._name],int):
                logger.info('No data in bin {}.'.format(bin_idx))
                summed_data[det._name] = np.nan
                continue
            
            # make full image for each det if requested
            if hasattr(det, 'x') and thisDetDict['image']==1:
                summed_data[det._name] = det.det.image(self.run, summed_data[det._name])
            
            # process the additional functions
            for func in [det.__dict__[k] for k in det.__dict__ if isinstance(det.__dict__[k], DetObjectFunc)]:
                if isinstance(func, (ROIFunc, photonFunc, dropletFunc, droplet2Photons)):
                    continue
                proc_data[det._name] = func.process(summed_data[det._name])
                logger.debug(f'Processed data keys: {proc_data[det._name].keys()}')
        return summed_data, proc_data

    
    
class mpi_bin_data(object):
    def __init__(self, bin_idx):
        self.bin_idx = bin_idx
