#!/usr/bin/env python

from datetime import datetime
import numpy as np
import psana
import time
import argparse
import socket
import os
import logging 
import requests
import sys
from glob import glob
from PIL import Image
from requests.auth import HTTPBasicAuth

fpath=os.path.dirname(os.path.abspath(__file__))
fpathup = '/'.join(fpath.split('/')[:-1])
sys.path.append(fpathup)
print(fpathup)

from smalldata_tools.utilities import printMsg, checkDet
from smalldata_tools.SmallDataUtils import setParameter, getUserData, getUserEnvData
from smalldata_tools.SmallDataUtils import defaultDetectors, detData
from smalldata_tools.SmallDataDefaultDetector import epicsDetector, eorbitsDetector
from smalldata_tools.SmallDataDefaultDetector import bmmonDetector, ipmDetector
from smalldata_tools.SmallDataDefaultDetector import encoderDetector, adcDetector
from smalldata_tools.DetObject import DetObject
from smalldata_tools.roi_rebin import ROIFunc, spectrumFunc, projectionFunc, sparsifyFunc, imageFunc
from smalldata_tools.waveformFunc import getCMPeakFunc, templateFitFunc
from smalldata_tools.droplet import dropletFunc
from smalldata_tools.photons import photonFunc
from smalldata_tools.azimuthalBinning import azimuthalBinning
