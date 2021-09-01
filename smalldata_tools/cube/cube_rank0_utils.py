import logging
import numpy as np
from bokeh.plotting import figure
from bokeh.models import Span
import panel as pn

logger = logging.getLogger()


def make_report(ana, hist_list, filters, exp, run):
#     gspec = pn.GridSpec(sizing_mode='stretch_both')

    # HISTOGRAMS
    gspec = pn.GridSpec(height=len(hist_list)*600, width=1000, name='Histograms')
    size = 3
    for ii,varName in enumerate(hist_list):
        var = ana.getVar(varName)
        low = [filt[1] for filt in filters if filt[0]==varName]
        high = [filt[2] for filt in filters if filt[0]==varName]

        title = '{}. Limits: {} - {}'.format(varName, low, high)
        counts, edges = np.histogram(var, np.linspace(np.percentile(var,0.1), np.percentile(var,99.9), 60))
        fig = figure(title=title, x_axis_label=varName, y_axis_label='Count', background_fill_color="#fafafa")
        fig.quad(top=counts, bottom=0, left=edges[:-1], right=edges[1:], line_color='black', fill_color='#000080')
        if low and high:
            l1 = Span(location=low[0], dimension='height', line_color='orange', line_width=2)
            l2 = Span(location=high[0], dimension='height', line_color='orange', line_width=2)
            fig.renderers.extend([l1,l2])
        l = ii%2 # is left
        gspec[size*(ii-l):size*(ii-l+1)-1,l:l+1] = fig
    tabs = pn.Tabs(gspec)
    
    # save to stats dir
    reports_dir = results_dir = ''.join(['/cds/data/psdm/', exp[:3], '/', exp, '/stats/summary/Cube_', run])
    report_file = ''.join([reports_dir, '/report.html'])
    logger.info('Saving report to {}'.format(report_file))
    tabs.save(report_file)
    