import os
import copy
import numpy as np

from future.utils import iteritems
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
try:
    basestring
except NameError:
    basestring = str

    
#epix10k: thermisotr value to temp in C
def getThermistorTemp(x):
    if x==0: return 0.
    u = x/16383.0 * 2.5
    i = u / 100000
    r = (2.5 - u)/i
    try:
        l = np.log(r/10000)
        t = 1.0 / (3.3538646E-03 + 2.5654090E-04 * l + 1.9243889E-06 * (l*l) + 1.0969244E-07 * (l*l*l))
        return t - 273.15
    except:
        return np.nan
