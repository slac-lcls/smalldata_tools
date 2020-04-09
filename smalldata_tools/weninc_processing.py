from __future__ import absolute_import, division, print_function

import os
import time
import psana
import tables
import numpy as np
from mpi4py import MPI

class Batch():
    def __init__(self, batch_size):
        self.data = {}
        self.counter = 0
        self.batch_size = batch_size
        self.comm = MPI.COMM_WORLD
        
    def add_data(self, key, value, kind='var'):
        if key in self.data:
            l = self.data[key][0]
            assert type(value) == type(l[-1])
            l.append(value)
        else:
            self.data[key] = ([value], kind)   
            
    def event_complete(self):
        self.counter += 1
        if self.counter == self.batch_size:
            self.send()
            self.counter = 0
                
    def send(self):
        d = {}
        for key, value in self.data.items():
            data, kind = value
            t = type(data[0])
            if t is np.ndarray:
                if kind == 'var':
                    d[key] = np.concatenate(data)
                else:
                    d[key] = np.array(data)
            elif t is int:
                d[key] = np.array(data, dtype=np.int64)
            elif t is float:
                d[key] = np.array(data, dtype=np.float32)
            else:
                raise TypeError('value must be np.ndarray or int or float')
        self.comm.send(d, dest=0, tag=1)
        self.data.clear()

def events(ds, batch):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank() - 1
    size = comm.Get_size()
    nworkers = size - 1
    fee_det = psana.Detector('FEEGasDetEnergy')
    for i, evt in enumerate(ds.events()):
        if i%nworkers != rank:
            continue 
        evt_id = evt.get(psana.EventId)
        if evt_id is None:
            print('Missing event id in event %d' %i)
            continue
        batch.add_data('time', evt_id.time()[0] << 32 | evt_id.time()[1])
        fee_gas = fee_det.get(evt)
        if fee_gas is None:
            pulse_energy = -1.0
        else:
            pulse_energy = 0.5*(fee_gas.f_21_ENRC() + fee_gas.f_22_ENRC())
        batch.add_data('pulse_energy', pulse_energy)
        yield evt
        batch.event_complete()
    batch.send()
    comm.send('worker finished', dest=0, tag=0)

def master(output_file):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    nworkers = size - 1
    fh = tables.open_file(output_file, 'w')
    dsets = {}
    active_workers = nworkers
    status = MPI.Status()
    nevents = 0
    print_events = 0
    start = time.time()
    while active_workers > 0:
        data = comm.recv(source=MPI.ANY_SOURCE, status=status)
        tag = status.Get_tag()
        if tag == 0:
            active_workers -= 1
            print('Worker %d finished' %status.Get_source())
        else:
            if 'time' in data:
                nevents += data['time'].size
            else:
                print('no data in batch')
            if nevents >= print_events:
                print_events += 1000
                print('%s %d events processed' %(time.ctime(), nevents))
            for key, value in data.items():
                if key not in dsets:
                    a = tables.Atom.from_dtype(value.dtype)
                    print(key, value.shape)
                    if len(value.shape) > 1:
                        shape = (0,) + value.shape[1:]
                    else:
                        shape = (0,)
                    dsets[key] = fh.create_earray(fh.root, key, a, shape)
                dsets[key].append(value)
    fh.close()
    end = time.time()
    rate = nevents / (end - start)
    print('Processed %d events with %.1f Hz' %(nevents, rate))