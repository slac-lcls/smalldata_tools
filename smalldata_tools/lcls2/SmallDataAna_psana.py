import numpy as np
import logging
import socket
import psana

from smalldata_tools.BaseSmallDataAna_psana import BaseSmallDataAna_psana
from smalldata_tools.DetObject_lcls2 import DetObject_lcls2

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SmallDataAna_psana_lcls2(BaseSmallDataAna_psana):
    def __init__(self, expname='', run=-1, dirname='', filename='', plotWith='matplotlib'):
        super().__init__(expname=expname, run=run, dirname=dirname, filename=filename, plotWith=plotWith)
        
        xtcdirname = None
        hostname = socket.gethostname()
        if hostname.find('drp-srcf')>=0:
            xtcdirname='/cds/data/drpsrcf/{self.hutch.lower()/{expname}/xtc'
        self.ds = psana.DataSource(exp=expname, run=run, dir=xtcdirname)
        return