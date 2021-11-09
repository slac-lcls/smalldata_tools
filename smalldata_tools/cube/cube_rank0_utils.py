import logging
import numpy as np
from pathlib import Path
from bokeh.plotting import figure
from bokeh.models import Span
import panel as pn
import holoviews as hv
hv.extension('bokeh')

logger = logging.getLogger()

TOTAL_WIDTH = 850


def make_report(anaps, hist_list, filters, varList, exp, run):
    # HISTOGRAMS
    hists_tab = pn.GridSpec(height=len(hist_list)*220, width=TOTAL_WIDTH, name='Histograms')
    
    filts = []
    for ii,(varName, hist_param) in enumerate(hist_list.items()):
        var = anaps.sda.getVar(varName)
        low = [filt[1] for filt in filters if filt[0]==varName]
        high = [filt[2] for filt in filters if filt[0]==varName]
        
        if low and high:
            filts.append( np.logical_and.reduce([var>low[0], var<high[0]]) )
        else:
            filts.append(None)
            
    for ii,(varName, hist_param) in enumerate(hist_list.items()):
        var = anaps.sda.getVar(varName)
        filt_var = np.logical_and.reduce([f for j,f in enumerate(filts) if (f is not None) and (j!=ii)])
        var_filtered = var[filt_var]
        low = [filt[1] for filt in filters if filt[0]==varName]
        high = [filt[2] for filt in filters if filt[0]==varName]
        
        if low and high:
            ratio = np.logical_or.reduce([var<low[0], var>high[0]]).sum()/var.size
            ratio_filtered = np.logical_or.reduce([var_filtered<low[0], var_filtered>high[0]]).sum()/var_filtered.size
            title = '{}. Filter: {} - {}. Rejected: {:.1%} / {:.1%}'.format(varName, low, high, ratio, ratio_filtered)
        else:
            title = '{}'.format(varName)
        if hist_param is None:
            bins = np.linspace(np.percentile(var,0.01), np.percentile(var,99.9), 60)
        else:
            bins = np.linspace(hist_param[0], hist_param[1], hist_param[2])
        counts, edges = np.histogram(var, bins=bins)
        counts_filtered, edges = np.histogram(var_filtered, bins=bins)
        fig = figure(title=title, x_axis_label=varName, y_axis_label='Count', background_fill_color="#fafafa")
        fig.quad(top=counts, bottom=0, left=edges[:-1], right=edges[1:], line_color='black', 
                 fill_color='#000080')
        fig.quad(top=counts_filtered, bottom=0, left=edges[:-1], right=edges[1:], line_color='black', 
                 fill_color='#59CBE8', fill_alpha=0.7)
        if low and high:
            l1 = Span(location=low[0], dimension='height', line_color='orange', line_width=2)
            l2 = Span(location=high[0], dimension='height', line_color='orange', line_width=2)
            fig.renderers.extend([l1,l2])
        l = ii%2 # is left
        hists_tab[((ii+1)//2-l):((ii+1)//2-l+1),l:l+1] = fig
    tabs = pn.Tabs(hists_tab)
    
    # DETECTORS DIAGNOTICS
    dets = []
    for var in varList:
        if isinstance(var, dict):
            if var['source']:
                dets.append(var)
        elif var not in anaps.sda._fields:
            dets.append(var)
    dets_tab = pn.GridSpec(height=500*len(dets), width=TOTAL_WIDTH, name='Detectors')
    for ii,det in enumerate(dets):
        threshold = None
        if isinstance(det, dict):
            detname = det['source']
            if 'thresADU' in det:
                threshold = det['thresADU']
        else:
            detname = det
        counts, bin_centers, bin_edges = anaps.pixelHistogram(detname=detname, numEvts=100, nSkip=0, common_mode=None, nBins=150)
    #     fig = hv.Histogram((bin_edges, counts))
    #     fig.opts(title='jungfrau1M', fill_color='#000080', logy=True)
        fig = figure(title=detname, x_axis_label='Intensity (ADU)', y_axis_label='Count', 
                     y_axis_type="log", background_fill_color="#fafafa")
        fig.quad(top=counts, bottom=0.1, left=bin_edges[:-1], right=bin_edges[1:], line_color='black', fill_color='#000080')
        if threshold is not None:
            l1 = Span(location=threshold, dimension='height', line_color='orange', line_width=2)
            fig.renderers.extend([l1])
        dets_tab[ii:ii+1,0] = fig
    
#         anaps.AvImage(detname=det)
#         detname, img, av_img = anaps.getAvImage(detname=detname)
#         img = anaps.__dict__[det].det.image(run, img)
#         xrange, yrange = [0,img.shape[0]], [0,img.shape[1]]
#         fig = hv.Image(img, bounds=(xrange[0],yrange[0], xrange[1], yrange[1])).opts(cmap='viridis',clim=(0,10))
#         dets[1:2,0] = fig
    
    tabs.append(dets_tab)
    
    # save to stats dir
    reports_dir = Path(''.join(['/cds/data/psdm/', exp[:3], '/', exp, '/stats/summary/Cube_', str(run)]))
    if not reports_dir.exists():
        reports_dir.mkdir()
    report_file = reports_dir / 'report.html'
    logger.info('Saving report to {}'.format(report_file))
    print('Saving report to {}'.format(report_file))
    tabs.save(report_file)
    
    return tabs