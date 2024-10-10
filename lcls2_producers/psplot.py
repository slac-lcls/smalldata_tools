from abc import ABCMeta, abstractmethod
import numpy as np
import logging

from psmon import publish
from psmon.plots import Image, XYPlot, MultiPlot


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
    def __init__(self, **kwargs):
        logger.info(f"Instantiating {type(self).__name__} with arguments {kwargs}")
        self.data_fields = kwargs['data_fields']
        self._needed = self.need()

    @abstractmethod
    def process(self, data_dict):
        pass

    def has_needed(self, data_dict):
        return all(map(lambda x: x in data_dict.keys(), self._needed))

    def need(self):
        """
        Return a list of field needed in the data dict to run the function
        """
        # Can probably be done better
        needed = []
        for key, val in self.data_fields.items():
            if isinstance(val, str):
                needed.append(val)
        return needed


class PsplotScan(PsplotCallback, metaclass=ABCMeta):
    """Psplot callback base class for plotting scans"""

    def __init__(self, bins=None, **kwargs):
        """ """
        super().__init__(**kwargs)
        self.det_binned = None
        self.i0_binned = None
        self.counts = None
        self.bins = bins
        self.nbins = self.bins.size + 1
        return

    @abstractmethod
    def process(self, data_dict):
        pass


class SpectrumScan(PsplotScan):
    def __init__(self, bins=None, **kwargs):
        print(f"kwargs: {kwargs}")
        self._range = kwargs.pop("spectrum_range", None)
        self._lineouts_idx = kwargs.pop("lineouts_idx", None)
        super().__init__(bins=bins, **kwargs)
        return

    def process(self, data_dict):
        det_dat = data_dict[self.data_fields["data"]]
        i0_dat = data_dict[self.data_fields["norm"]]
        count = data_dict[self.data_fields["count"]]
        scan_var = data_dict[self.data_fields["scan"]] / count

        if self._range is not None:
            det_dat = det_dat[self._range[0] : self._range[1]]
        i0_dat = np.nansum(i0_dat)

        if self.det_binned is None:
            self.det_binned = np.zeros((self.nbins, *det_dat.shape))
            self.i0_binned = np.zeros((self.nbins))[:, np.newaxis]
            self.counts = np.zeros((self.nbins))[:, np.newaxis]
            print(f"SHAPE: {self.det_binned.shape}")

        bin_idx = np.digitize(scan_var, self.bins)
        self.det_binned[bin_idx] += det_dat
        self.i0_binned[bin_idx, :] += i0_dat
        self.counts[bin_idx, :] += count

        im = self.det_binned / self.i0_binned
        im /= np.nanmean(im)  # make it a nicer range to work with in the plot
        im[np.where(im == 0)] = np.nan

        scan_plot = Image(0, "Spectrum Scan", im, aspect_lock=False)
        publish.send("spectrum_scan", scan_plot)

        ys = []
        for ii, idx in enumerate(self._lineouts_idx):
            x = self.bins
            y = im[:-1, idx]  # cut last bin to match shape
            ys.append(y)  # lineouts at idx
        lineouts = XYPlot(0, f"Spectrum Lineout", x, ys, leg_label=self._lineouts_idx)
        publish.send(f"spectrum_lineout", lineouts)


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
