from abc import ABCMeta, abstractmethod
import numpy as np
import logging

from psmon import publish
from psmon.plots import Image, XYPlot


logger = logging.getLogger(__name__)

class PsplotCallbacks:
    def __init__(self):
        self.callbacks = {}

    def run(self, data_dict):
        for name, item in self.callbacks.items():
            if item.has_needed(data_dict):
                item.process(data_dict)
    
    def add_callback(self, item, name):
        self.callbacks[name] = item
        return


class PsplotCallback(metaclass=ABCMeta):
    def __init__(self, det_dict):
        self.det_dict = det_dict
        logger.info(f"Instantiating {type(self).__name__} with det_dict {det_dict}")
        self._needed = self.need()
    
    @abstractmethod
    def process(self, data_dict):
        pass

    def has_needed(self, data_dict):
        return all( map(lambda x: x in data_dict.keys(), self._needed) )

    def need(self):
        """
        Return a list of field needed in the data dict to run the function
        """
        # Can probably be done better
        needed = []
        for key, val in self.det_dict.items():
            if isinstance(val, str):
                needed.append(val)
        return needed



class PsplotScan(PsplotCallback, metaclass=ABCMeta):
    """ Psplot callback base class for plotting scans """
    def __init__(self, det_dict):
        """
        """
        print(det_dict)
        super().__init__(det_dict)
        self.det_binned = None
        self.i0_binned = None
        self.counts = None
        self.bins = det_dict.get('bins') #, default=None)
        self.nbins = self.bins.size + 1
        return

    @abstractmethod
    def process(self, data_dict):
        pass


class SpectrumScan(PsplotScan):
    def __init__(self, det_dict):
        super().__init__(det_dict)
        return

    def process(self, data_dict):
        det_dat = data_dict[self.det_dict['data']]
        fim1 = data_dict[self.det_dict['norm']]
        count = data_dict[self.det_dict['count']]
        scan_var = data_dict[self.det_dict['scan']] / count

        det_dat = self.hitfinder(det_dat)
        det_dat = det_dat[900:1050]
        i0_fim1 = np.nansum(fim1)

        if self.det_binned is None:
            self.det_binned = np.zeros((self.nbins, *det_dat.shape))
            self.i0_binned = np.zeros((self.nbins))[:, np.newaxis]
            self.counts = np.zeros((self.nbins))[:, np.newaxis]
            print(f"SHAPE: {self.det_binned.shape}")

        bin_idx = np.digitize(scan_var, self.bins)
        self.det_binned[bin_idx] += det_dat
        self.i0_binned[bin_idx, :] += i0_fim1
        self.counts[bin_idx, :] += count

        im = self.det_binned / self.i0_binned
        im /= np.nanmean(im)  # make it a nicer range to work with in the plot
        im[np.where(im == 0)] = np.nan

        scan_plot = Image(0, 'Andor VLS scan', im, aspect_lock=False)
        publish.send('andor_vls_scan', scan_plot)

        y = det_dat
        x = np.arange(y.size)
        spectrum = XYPlot(0, 'Andor VLS', x, [y, y+0.5])
        publish.send('andor_vls', spectrum)

    @staticmethod
    def hitfinder(data, bkg_roi=(500,800), threshold=3500, std_threshold=4):
        amean = data[bkg_roi[0]:bkg_roi[1]].mean()
        astd = data[bkg_roi[0]:bkg_roi[1]].std()

        data[data>threshold] = 0
        data[data < (amean + std_threshold*astd)] = 0
        data[data > (amean + std_threshold*astd)] = 1
        return data

# DOCUMENTATION
# XYPlot(
#     ts,
#     title,
#     xdata,
#     ydata,
#     xlabel=None,
#     ylabel=None,
#     leg_label=None,
#     leg_offset=None,
#     formats='-',
#     xdate=False,
#     ydate=False,
# )
# 
# Image(
#     ts,
#     title,
#     image,
#     xlabel=None,
#     ylabel=None,
#     aspect_ratio=None,
#     aspect_lock=True,
#     pos=None,
#     scale=None,
# )
# 
# publish.send
# Signature: publish.send(topic, data)
# Docstring:
# Publishes a data object to all clients suscribed to the topic.

# Arguments
#  - topic: The name of the topic to which the data is being published.
#  - data: The data object to be published to suscribers.
