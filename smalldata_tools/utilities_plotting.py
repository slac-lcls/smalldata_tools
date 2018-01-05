import os
import json
import copy
import numpy as np
import time
import h5py
import tables
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
import resource

import bokeh
import bokeh.plotting as bp
from bokeh.models import PanTool, SaveTool, HoverTool, ResetTool, ResizeTool
from bokeh.models import WheelZoomTool, BoxZoomTool

from collections import deque
from itertools import islice
from bisect import insort, bisect_left
 
import sys

try:
    sys.path.append('/reg/neh/home/snelson/gitMaster_smalldata_tools/')
    import bokeh_utils
    cmaps = bokeh_utils.get_all_mpl_palettes(allmaps=['cool','gray','jet','spectral'])
except:
  print 'no bokeh utils'
  pass

def create_range_slider(vabsmin,vabsmax,vmin,vmax,im,step=0.1):
    JS_code_slider = """                                                                            
        var vmin = rslider.range[0];
        var vmax = rslider.range[1];
        im.glyph.color_mapper.high = vmax;
        im.glyph.color_mapper.low = vmin; 
        im.data_source.trigger('change');
    """
    callback_slider = bokeh.models.CustomJS(args=dict( im=im),
                                        code=JS_code_slider)

    try:
        rslider = bokeh.models.RangeSlider(title="low/high limit",start=vabsmin,end=vabsmax,step=step, range=[vmin,vmax],callback=callback_slider,orientation="horizontal")
    except:
        rslider = bokeh.models.RangeSlider(title="low/high limit",start=vabsmin,end=vabsmax,step=step, value=[vmin,vmax],callback=callback_slider,orientation="horizontal")

    callback_slider.args['rslider'] = rslider

    return rslider

def plotImageBokeh(data, plotWidth=600, plotHeight=400, xRange=None, yRange=None, plot_title='', dateTime=False, plotMinP=None, plotMaxP=None, plotMin="auto", plotMax="auto", initial_cmap='jet', output_quad=False, tools=None):
    if plotMinP is None: 
        plotMinP = np.nanpercentile(data, 5)
    if plotMaxP is None: 
        plotMaxP = np.nanpercentile(data, 95)
    if xRange is None:
        xRange=(0,data.shape[0])
    if yRange is None:
        yRange=(0,data.shape[1])
    if tools is None:
        if bokeh.__version__=='0.12.6':
            tools = 'pan, wheel_zoom, box_zoom, resize, reset'
        else:
            tools = 'pan, wheel_zoom, box_zoom, reset'
    
    if dateTime:
        created_map = bokeh_utils.create_map(data2plot=data,palette_name=initial_cmap,
                                            fig_width_pxls=plotWidth, fig_height_pxls=plotHeight,x_range=xRange,
                                            y_range=yRange,title=plot_title,x_axis_type="datetime", tools=tools,
                                            cmaps=cmaps,vmin=plotMinP,vmax=plotMaxP, output_quad=output_quad)
    else: 
        created_map = bokeh_utils.create_map(data2plot=data,palette_name=initial_cmap,
                                            fig_width_pxls=plotWidth, fig_height_pxls=plotHeight,x_range=xRange,
                                            y_range=yRange,title=plot_title, tools=tools,
                                            cmaps=cmaps,vmin=plotMinP,vmax=plotMaxP, output_quad=output_quad)

    if output_quad:
        p, im, q = created_map
    else:
        p, im = created_map
        
    # Controls
    vabsmin,vabsmax = plotMin, plotMax
    step=(plotMaxP-plotMinP)/50.
    range_slider = create_range_slider(np.nanmin(data),np.nanmax(data),plotMinP,plotMaxP,im=im,step=step)
    select_cm = bokeh_utils.create_cmap_selection(im,cmaps=cmaps, value=initial_cmap)

    # Layout
    layout = bp.gridplot([[p],[range_slider,select_cm]])

    if output_quad:
        return layout,p,im,q
    else:
        return layout,p,im


def getColors(nColors=None, maxCol=2.6, minCol=0.25):
    """
    class to return a list of colors names to be used in plots.
    Parameters:
    nColors: number of colors to return, def: return all passing
             brightness selection (not too black or white)
    maxCol: maximum brightness (sum red+blue+yellow), def 2.6
    minCol: minimum brightness (sum red+blue+yellow), def 0.25
    """
    colors = mcolors.cnames
    #sort them
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.hex2color(color)[:3])),name) for name, color in colors.items() if (np.array(mcolors.hex2color(color)[:3]).sum()<maxCol and np.array(mcolors.hex2color(color)[:3]).sum()>minCol))
    sorted_names = [name for hsv, name in by_hsv]
    if nColors is None:
        return sorted_names
    else:
        dcol = int(len(sorted_names)/nColors)
        color_names=[]
        for i in range(nColors+1):
            color_names.append(sorted_names[i*dcol])
        return color_names

def plot1d(data, **kwargs):
    """
    class to plot 1-d data (histograms, scans).
    Input are xData (1-d array, optionally list of arrays)

    Parameters
    ----------
    xData: vector of x-values for the data, needs to have the same shape
    xLabel: string to be used for labelling a X-axis
    yLabel: string to be used for labelling the Y-axis
    plotWith: method for plotting, options are matplotlib, 
              bokeh_notebook or bokeh w/ output as a html file
    fig: if using matplotlib, a figure can be passed in which the new data 
         will be plotted
    runLabel: label to indicate the run ('Run001'), def 'Run', can be list
              is passing array
    plotTitle: title to show on plot
    plotDirname: directory to which bokeh html files will be written (def ./)
    marker: pass marker style (def o, options: *,^,+,x,s,D)
    markersize: marker size, default is 5
    line_dash: connect markers with line style as defined (default: none)
    tools: for bokeh plotting, tools can be passed. Default tools will be used
           otherwise
    """
    if isinstance(data, np.ndarray):
        data = data.tolist()
    if isinstance(data[0], int) or isinstance(data[0], float):
        data = [data]
    xData = kwargs.pop("xData", [])
    if xData == []:
        for yData in data:
            xData.append(np.arange(np.array(yData).shape[0]).tolist())
    if isinstance(xData[0], int) or isinstance(xData[0], float):
        xData = [xData]

    xLabel = kwargs.pop("xLabel", "")
    yLabel = kwargs.pop("yLabel", "")
    plotWith = kwargs.pop("plotWith", "matplotlib")
    RunLabel = kwargs.pop("runLabel", "Run")
    if isinstance(RunLabel, list):
        runLabel = RunLabel[0]
    else:
        runLabel = RunLabel
    plotTitle = kwargs.pop("plotTitle","%s histogram for %s"%(xLabel, runLabel))
    marker = kwargs.pop("marker", "o")
    line_dash = kwargs.pop("line_dash", None)
    markersize = kwargs.pop("markersize", 5)
    tools = kwargs.pop("tools", None)
    plotDirname = kwargs.pop("plotDirname", "./") 
    plotFilename =  kwargs.pop("plotFilename", '%s/%s_%s_histogram'%(plotDirname,runLabel, xLabel.replace('/','_')))    
    if len(data)>1:
        sorted_colors = getColors(nColors=len(data))
    if plotWith.find('matplotlib')>=0:
        if 'fig' in kwargs.keys():
            fig = kwargs.pop("fig")
        else:
            plt.figure(figsize=(8,5))
        if len(kwargs)>0:
            print 'found unexpected parameters to plot1d, will ignore', kwargs
        if line_dash is None:
            line_dash = 'none'
        for ic,(idata,ixdata) in enumerate(zip(data, xData)):
            if ic > 0:
                plt.plot(ixdata,idata,marker=marker, color=sorted_colors[ic],linestyle=line_dash)
            else:
                plt.plot(ixdata,idata,marker=marker,linestyle=line_dash,markersize=markersize)
        plt.xlabel(xLabel)
        plt.ylabel(yLabel)
        plt.title(plotTitle)                  
        if plotWith.find('file')>=0:
            plt.savefig('%s.jpeg'%(plotFilename))

    elif plotWith.find('bokeh')>=0:
        #IMPLEMENT ME
        #bokeh, use BoxSelectTool to get selected region
        #https://stackoverflow.com/questions/34164587/get-selected-data-contained-within-box-select-tool-in-bokeh
        #set cuts does not need to be specified explicly, can simple be gotten from callback of tool.....
        if tools is None:
            pan=PanTool()
            wheel_zoom=WheelZoomTool()
            box_zoom=BoxZoomTool()
            save=SaveTool()
            reset=ResetTool()
            hover=HoverTool(tooltips=[
                ("(x,y)","($x, $y)")
            ])
            tools = [pan, wheel_zoom,box_zoom,save,hover,reset]
            if bokeh.__version__=='0.12.6':
                resize=ResizeTool()
                tools = [pan, wheel_zoom,box_zoom,resize,save,hover,reset]
        p = bp.figure(title=plotTitle, x_axis_label=xLabel, y_axis_label=yLabel,tools=tools)
        for ic,(idata,ixdata) in enumerate(zip(data, xData)):
            if len(data)==1:
                color='black'
            else:
                color=sorted_colors[ic]
            if isinstance(RunLabel, list):
                iRunLabel=RunLabel[ic]
            else:
                iRunLabel=runLabel
            if marker=='o':
                p.circle(ixdata, idata, legend=iRunLabel, size=markersize, color=color)
            elif marker=='*':
                p.asterisk(ixdata, idata, legend=iRunLabel, size=markersize, color=color)
            elif marker=='+':
                p.cross(ixdata, idata, legend=iRunLabel, size=markersize, color=color)
            elif marker=='x':
                p.x(ixdata, idata, legend=iRunLabel, size=markersize, color=color)
            elif marker=='s':
                p.square(ixdata, idata, legend=iRunLabel, size=markersize, color=color)
            elif marker=='^':
                p.triangle(ixdata, idata, legend=iRunLabel, size=markersize, color=color)
            elif marker=='D':
                p.diamond(ixdata, idata, legend=iRunLabel, size=markersize, color=color)
            else:
                print 'unsupported marker option'
            if line_dash is not None:
                p.line(ixdata, idata, line_dash=line_dash, line_color=color)
        if plotWith=='bokeh_notebook':
            if len(kwargs)>0:
                print 'found unexpected parameters to plot1d, will ignore', kwargs
            bp.output_notebook()
            bp.show(p)
        else:
            if len(kwargs)>0:
                print 'found unexpected parameters to plot1d, will ignore', kwargs
            bp.output_file('%s.html'%(plotFilename))
            bp.save(p)
                
    elif plotWith != 'no_plot':
        print 'plotting using %s is not implemented yet, options are matplotlib, bokeh_notebook, bokeh_html or no_plot'%plotWith
            
    return


