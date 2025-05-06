from abc import ABCMeta, abstractmethod
import numpy as np
import logging

from psmon import publish
from psmon.plots import Image, XYPlot, MultiPlot
from scipy.ndimage import rotate 


logger = logging.getLogger(__name__)

def process_VLS_2D(vls, bg_roi = [[120,-1],[20,100]],threshold_min=20,threshold_max=4000,photon_count = False,rotate = False):
    #bgs = vls[bg_roi[0][0]:bg_roi[0][1],bg_roi[1][0]:bg_roi[1][1]].mean(axis = (1,2))
    #im_nbg = vls - bgs[:,np.newaxis,np.newaxis]
    #astd,amean = im_nbg[bg_roi[0][0]:bg_roi[0][1],bg_roi[1][0]:bg_roi[1][1]].std(),im_nbg[bg_roi[0][0]:bg_roi[0][1],bg_roi[1][0]:bg_roi[1][1]].mean()
    im_nbg = vls-717
    im_proc = im_nbg.copy()
    im_proc[im_proc>threshold_max] = 0
    im_proc[im_proc<threshold_min] = 0
    if photon_count:
        im_proc[im_proc>threshold_min] = 1
    if rotate:
        # Rotate the frames by X degrees. The second input is the angle in degrees. 
        im_proc = rot(im_proc,4.5,axes=(2,1))
        im_proc = im_proc[:,12:127,:]
    return im_proc

    
def process_ANDOR_Basic(andor, bg_roi_x = [0,100],bg_roi_y=[0,100],\
                        absolute_theshold = False,threshold = False,photon_count = False,\
                        threshold_min=4, threshold_max=4000):
    # 1D processing
    if len(andor.shape) ==2:
        andor_nbg = subtract_bg(andor,bg_roi_x)
        andor_proc = andor_nbg.copy()
        if threshold:
            astd,amean = andor_nbg[:,bg_roi_x[0]:bg_roi_x[1]].std(),andor_nbg[:,bg_roi_x[0]:bg_roi_x[1]].mean()
            andor_proc[andor_proc>threshold_max] = 0
            if absolute_threshold:
                andor_proc[andor_proc<(threshold_min)] = 0
            else:
                andor_proc[andor_proc<(threshold_min*astd)] = 0
            if photon_count:
                andor_proc[andor_proc>threshold_min] = 1
    elif len(andor.shape)==3:
        # Shape of full image is (2048 x 512)
        bgs = andor[:,bg_roi_y[0]:bg_roi_y[1],bg_roi_x[0]:bg_roi_x[1]].mean(axis = (1,2))
        andor_nbg = andor - bgs[:,np.newaxis,np.newaxis]
        andor_proc = andor_nbg.copy()
        if threshold:
            astd = andor_nbg[:,bg_roi_y[0]:bg_roi_y[1],bg_roi_x[0]:bg_roi_x[1]].std()
            amean = andor_nbg[:,bg_roi_y[0]:bg_roi_y[1],bg_roi_x[0]:bg_roi_x[1]].mean()
            andor_proc[andor_proc>threshold_max] = 0
            if absolute_threshold:
                andor_proc[andor_proc<(threshold_min)] = 0
            else:
                andor_proc[andor_proc<(threshold_min*astd)] = 0
            if photon_count:
                andor_proc[andor_proc>threshold_min] = 1
    return andor_proc
def sum_ANDOR(andor_ims):
    if len(andor_ims.shape) ==2:
        return np.sum(andor_ims,axis = 1)
    elif len(andor_ims.shape) ==3:
        return np.sum(andor_ims,axis = (1,2))

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
            elif isinstance(val, list):
                [needed.append(v) for v in val]
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
        if isinstance(self.bins, list):
            self.nbins = [b.size+1 for b in self.bins]
        else:
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
        self._labels = kwargs.pop('labels', None)
        super().__init__(bins=bins, **kwargs)
        return

    def process(self, data_dict):
        evc_laser_off = 273
        evc_laser_on = 273

        # VLS processing 
        im = data_dict[self.data_fields["data"]]
        #vls_bg = np.nanmean(im[480:510,20:80],axis = (0,1))
        #vls_background_subtracted = im - vls_bg[np.newaxis,np.newaxis]
        #vls_rotated = rotate(vls_background_subtracted,4.5,order=3,axes = (1,0),cval=0)
        #vls_cropped = vls_rotated[11:510,:]
        det_dat =process_ANDOR_Basic(im)


        
        i0_dat = data_dict[self.data_fields["norm"]]
        count = data_dict[self.data_fields["count"]]
        scan_vars = [data_dict[ss]/count for ss in self.data_fields["scan_var"]]
        evc = data_dict[self.data_fields["evc"]]

        # Processing of the detectors
        if self._range is not None:
            det_dat = det_dat[self._range[0] : self._range[1]]
        i0_dat = np.nansum(i0_dat[4:])

        if self.det_binned is None:
            self.det_binned = [np.zeros((nb, *det_dat.shape)) for nb in self.nbins]
            self.i0_binned = [np.zeros((nb))[:, np.newaxis] for nb in self.nbins]
            self.counts = [np.zeros((nb))[:, np.newaxis] for nb in self.nbins]

            self.det_binned_on = [np.zeros((nb, *det_dat.shape)) for nb in self.nbins]
            self.i0_binned_on = [np.zeros((nb))[:, np.newaxis] for nb in self.nbins]
            self.counts_on = [np.zeros((nb))[:, np.newaxis] for nb in self.nbins]

            self.det_binned_off = [np.zeros((nb, *det_dat.shape)) for nb in self.nbins]
            self.i0_binned_off = [np.zeros((nb))[:, np.newaxis] for nb in self.nbins]
            self.counts_off = [np.zeros((nb))[:, np.newaxis] for nb in self.nbins]
            print([f"SHAPE: {db.shape}" for db in self.det_binned])

            self.det_binned_diff = [np.zeros((nb, *det_dat.shape)) for nb in self.nbins]

        for ii, (scan_var, bins, nbs) in enumerate(zip(scan_vars, self.bins, self.nbins)):
            print(scan_var,bins,np.digitize(scan_var, bins))
            bin_idx = np.digitize(scan_var, bins)

            self.det_binned[ii][bin_idx] += det_dat
            self.i0_binned[ii][bin_idx, :] += i0_dat
            self.counts[ii][bin_idx, :] += count
            im = self.det_binned[ii] / self.i0_binned[ii]
            im /= np.nanmean(im)  # make it a nicer range to work with in the plot
            im[np.where(im == 0)] = np.nan
            scan_plot = Image(0, f"{self._labels[ii]}", im, aspect_lock=False)
            publish.send(f"{self._labels[ii]}_scan_all", scan_plot)
            y = np.nanmean(im[:-1,self._lineouts_idx[0]:self._lineouts_idx[1]],axis = 1)
            x = bins
            pfy_all = XYPlot(0, f"{self._labels[ii]} Lineout", x, y)
            publish.send(f"{self._labels[ii]}_scan_pfy_all", pfy_all)


            if (evc[evc_laser_off]/count)<0.5:
                # Pumped analysis
                self.det_binned_on[ii][bin_idx] += det_dat
                self.i0_binned_on[ii][bin_idx, :] += i0_dat
                self.counts_on[ii][bin_idx, :] += count

                #scan_plot_on = Image(0, f"{self._labels[ii]}", im_on, aspect_lock=False)

            else:
                # Unumped analysis
                self.det_binned_off[ii][bin_idx] += det_dat
                self.i0_binned_off[ii][bin_idx, :] += i0_dat
                self.counts_off[ii][bin_idx, :] += count
                
                #scan_plot_off = Image(0, f"{self._labels[ii]}", im_off, aspect_lock=False)

            # Calculate things to plot
            im_on = self.det_binned_on[ii] / self.i0_binned_on[ii]
            im_on /= np.nanmean(im_on)  # make it a nicer range to work with in the plot
            im_on[np.where(im_on == 0)] = np.nan
            pfy_on = np.nanmean(im_on[:-1,self._lineouts_idx[0]:self._lineouts_idx[1]],axis = 1)
            em_on = np.nanmean(im_on,axis = 0)

            im_off = self.det_binned_off[ii] / self.i0_binned_off[ii]
            im_off /= np.nanmean(im_off)  # make it a nicer range to work with in the plot
            im_off[np.where(im_off == 0)] = np.nan
            pfy_off = np.nanmean(im_off[:-1,self._lineouts_idx[0]:self._lineouts_idx[1]],axis = 1)
            em_off = np.nanmean(im_off,axis = 0)

            self.det_binned_diff[ii] = (self.det_binned_on[ii] / self.i0_binned_on[ii]) - (self.det_binned_off[ii] / self.i0_binned_off[ii])
            im_diff = self.det_binned_diff[ii]
            im_diff /= np.nanmean(im_diff)
            pfy_diff = np.nanmean(im_diff[:-1,self._lineouts_idx[0]:self._lineouts_idx[1]],axis = 1)
            em_diff = np.nanmean(im_diff,axis = 0)

            scan_plot_on = Image(0, f"{self._labels[ii]} scan pumped", im_on, aspect_lock=False)
            scan_plot_off = Image(0, f"{self._labels[ii]} scan unpumped", im_off, aspect_lock=False)
            scan_plot_diff = Image(0, f"{self._labels[ii]} scan difference", im_diff, aspect_lock=False)
            publish.send(f"{self._labels[ii]}_scan_diff", scan_plot_diff)
            
            multiplot_RIXS = MultiPlot(0,'RIXS',ncols = 2)
            multiplot_RIXS.add(scan_plot_on)
            multiplot_RIXS.add(scan_plot_off) 
            publish.send(f"{self._labels[ii]}_scan_onoff",multiplot_RIXS)
            
            pfy_diff_plot = XYPlot(0, f"{self._labels[ii]} scan PFY difference", x, pfy_diff)
            pfy_on_off_plot = XYPlot(0, f"{self._labels[ii]}  scan PFY", [x,x], [pfy_on,pfy_off],formats = ['.-','x-'],leg_label = ['Pumped','Unpumped'])
            multiplot_PFY = MultiPlot(0,'PFY difference',ncols = 2)
            multiplot_PFY.add(pfy_on_off_plot)
            multiplot_PFY.add(pfy_diff_plot)
            publish.send(f"{self._labels[ii]}_scan_pfy_diff",multiplot_PFY)

            em_diff_plot = XYPlot(0, f"{self._labels[ii]} scan emission difference", x, em_diff)
            em_on_off_plot = XYPlot(0, f"{self._labels[ii]}  scan emission", [x,x], [em_on,em_off],formats = ['.-','x-'],leg_label = ['Pumped','Unpumped'])
            multiplot_em = MultiPlot(0,'Emission difference',ncols = 2)
            multiplot_em.add(em_on_off_plot)
            multiplot_em.add(em_diff_plot)
            publish.send(f"{self._labels[ii]}_scan_emission",multiplot_em)

            

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
