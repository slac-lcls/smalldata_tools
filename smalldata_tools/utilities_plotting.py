import os
import json
import copy
import numpy as np
import time
import h5py
import tables
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import gridspec
import matplotlib.cm
import resource

import bokeh
import bokeh.plotting as bp
from bokeh.models import PanTool, SaveTool, HoverTool, ResetTool
from bokeh.models import WheelZoomTool, BoxZoomTool, BoxSelectTool
from bokeh.models import ColumnDataSource
try:
    from bokeh.models import ResizeTool
except:
    pass

import sys

def valid_mpl_palettes():
    return sorted(matplotlib.cm.cmaps_listed.keys() + matplotlib.cm.datad.keys())

def get_all_mpl_palettes(allmaps=None):
    """
    map_names, palettes_dict, palette_nums_dict = get_all_palettes()
    """
    if type(allmaps)==type(None):
        # Getting all the map names available in matplotlib
        allmaps = valid_mpl_palettes()
    palette_nums = []
    palettes = []
    map_names = []
    for cmap_name in allmaps:
        palette_nums += [[]]
        palettes += [[]]
        cmap = matplotlib.cm.get_cmap(cmap_name)
        if cmap.N < 256:
            pts = np.linspace(0,255,cmap.N).round()
            rr = np.interp(np.arange(256),pts,np.array(cmap.colors)[:,0])
            gg = np.interp(np.arange(256),pts,np.array(cmap.colors)[:,1])
            bb = np.interp(np.arange(256),pts,np.array(cmap.colors)[:,2])
            palette = []
            for i in range(256):
                palette += [matplotlib.colors.rgb2hex((rr[i],gg[i],bb[i]))]
            palettes[-1] += palette
            palette_nums[-1] += [int(m[1:],16) for m in palette]
            map_names += [cmap_name+"_smooth"]
            palette_nums += [[]]
            palettes += [[]]
        cmap = matplotlib.cm.get_cmap(cmap_name,256)
        palette = [matplotlib.colors.rgb2hex(m) for m in cmap(np.arange(cmap.N))]  # always 256 pts!!!
        palettes[-1] += palette
        palette_nums[-1] += [int(m[1:],16) for m in palette]
        map_names += [cmap_name]
    palettes_dict = dict(zip(map_names, palettes))
    palette_nums_dict = dict(zip(map_names, palette_nums))
    return {"map_names":map_names,"palettes_dict":palettes_dict,"palette_nums_dict":palette_nums_dict}

#
# define cmaps.
#
cmaps = get_all_mpl_palettes(allmaps=['jet','gray','cool','hot','seismic','CMRmap','nipy_spectral'])


def create_cmap_selection(im,cmaps=None,value="jet"):
    if type(cmaps)==type(None):
        cmaps = get_all_palettes()
    bmem_cmaps = bokeh.models.ColumnDataSource(data=cmaps["palettes_dict"])
    call_cm = bokeh.models.CustomJS(args=dict(im=im,bmem_cmaps=bmem_cmaps),code="""
        sel_map = select_cm.value; // select_cm will be added to the args dictionary after creating the widget!!
        numpalette = bmem_cmaps.data[sel_map];  // needs to be added if a list of list

        var numglyph_palette1 = im.glyph.color_mapper["palette"];
        for (var i=0; i<numpalette.length; i++){
            numglyph_palette1[i] = numpalette[i];   // changing values of the palette with those of selected one
        }
        // this is needed to refresh the colorbar:
        var old = im.glyph.color_mapper.high;
        im.glyph.color_mapper.high = old+1;
        im.data_source.trigger.change.emit();
        im.glyph.color_mapper.high = old;
        im.data_source.trigger.change.emit();
        """)
    select_cm = bokeh.models.Select(title="Color Map:", value =value, options=cmaps["map_names"],callback=call_cm)
    call_cm.args.update(dict(select_cm=select_cm))

    return select_cm

def create_map(data2plot=None,palette_name="jet",fig_width_pxls=800,
               fig_height_pxls=500,x_range=(0,1),y_range=(0,1),
               title="MAP",x_axis_type="linear", cmaps=None,cb_title="",
               create_colorbar=True, min_border_left=20,min_border_right=10,
               min_border_top=30, min_border_bottom=10,title_font_size="12pt",
               title_align="center",plotMin="auto",plotMax="auto",
               output_quad=False,
               tools= ["box_zoom,wheel_zoom,pan,reset,previewsave,resize"]):
    """
    x_axis_type: "linear", "log", "datetime", "auto"
    """
    if type(cmaps)==type(None):
        cmaps = get_all_palettes()
    if plotMin=="auto":
        plotMin = np.nanmin(data2plot)
    if plotMax=="auto":
        plotMax = np.nanmax(data2plot)
    try:
        x0=x_range[0]
        x1=x_range[1]
    except:
        x0 = x_range.start
        x1 = x_range.end
    try:
        y0=y_range[0]
        y1=y_range[1]
    except:
        y0 = y_range.start
        y1 = y_range.end
                    
    p = bp.figure(x_range=x_range, y_range=y_range,x_axis_type=x_axis_type,
                  plot_width=fig_width_pxls, plot_height=fig_height_pxls, 
                  min_border_left=min_border_left,
                  min_border_right=min_border_right,
                  min_border_top=min_border_top,
                  min_border_bottom=min_border_bottom,
                  tools= tools,title=title)
    p.title.text_font_size = title_font_size
    p.title.align = title_align
    im = p.image(image=[data2plot],dw=[x1-x0],dh=[y1-y0],x=[x0],y=[y0],palette=cmaps["palettes_dict"][palette_name])
    im.glyph.color_mapper.high = plotMax
    im.glyph.color_mapper.low = plotMin
    imquad = p.quad(top=[y1], bottom=[y0], left=[x0], right=[x1],alpha=0) # This is used for hover and taptool
    if create_colorbar:
        color_bar = bokeh.models.ColorBar(color_mapper=im.glyph.color_mapper, label_standoff=12, location=(0,0))
        p.add_layout(color_bar, 'right')
    if output_quad:
        return p,im,imquad
    else:
        return p,im


def create_range_input_button(plotMin,plotMax,im):
    JS_code_plotMin_plotMax = """
        var plotMin = parseFloat(input_plotMin.value);
        var plotMax = parseFloat(input_plotMax.value);
        im.glyph.color_mapper.high = plotMax;
        im.glyph.color_mapper.low = plotMin;
        //im.data_source.trigger('change');
        im.data_source.change.emit();
    """
    callback_update = bokeh.models.CustomJS(code=JS_code_plotMin_plotMax)
    input_plotMin = bokeh.models.widgets.TextInput(value="%g"%plotMin,title="plotMin")
    input_plotMax = bokeh.models.widgets.TextInput(value="%g"%plotMax,title="plotMax")
    button_update =  bokeh.models.Button(label="Update plotMin/plotMax", callback=callback_update)

    callback_update.args['input_plotMin'] = input_plotMin
    callback_update.args['input_plotMax'] = input_plotMax
    callback_update.args['im'] = im

    return input_plotMin,input_plotMax,button_update

def create_range_slider(vabsmin,vabsmax,plotMin,plotMax,im,step=0.1):
    JS_code_slider = """                                                                            
        var plotMin = rslider.range[0];
        var plotMax = rslider.range[1];
        im.glyph.color_mapper.high = plotMax;
        im.glyph.color_mapper.low = plotMin; 
        im.data_source.trigger('change');
    """
    JS_code_slider_b0215 = """                                                                            
        var plotMin = rslider.value[0];
        var plotMax = rslider.value[1];
        im.glyph.color_mapper.high = plotMax;
        im.glyph.color_mapper.low = plotMin; 
        im.data_source.change.emit();
    """

    try:
        callback_slider = bokeh.models.CustomJS(args=dict( im=im),
                                                code=JS_code_slider)
        rslider = bokeh.models.RangeSlider(title="low/high limit",start=vabsmin,end=vabsmax,step=step, range=[plotMin,plotMax],callback=callback_slider,orientation="horizontal")
    except:
        callback_slider = bokeh.models.CustomJS(args=dict( im=im),
                                                code=JS_code_slider_b0215)
        rslider = bokeh.models.RangeSlider(title="low/high limit",start=vabsmin,end=vabsmax,step=step, value=[plotMin,plotMax],callback=callback_slider,orientation="horizontal")

    callback_slider.args['rslider'] = rslider

    return rslider

def create_img_slider_scale(im, imOrg,valStart=0,coords=None,p = None):
    JS_code_slider = """                                                                            
        var img_idx_val = rslider.value; 
        var cmin=1000000000000
        var cminIdx=-1        

        for (var i=0; i < sourceCoords.data.coords.length; i++){
            if (Math.abs(sourceCoords.data.coords[i] - img_idx_val) < cmin){
                cmin = Math.abs(sourceCoords.data.coords[i] - img_idx_val);
                cminIdx=i;
            }
        }
        var arShp = sourceShp.data.imShp[1]*sourceShp.data.imShp[2]
        var d3data_flat = source.data.im3d
        var newArr=[]
        for(var i = 0; i < d3data_flat.length; i += arShp) {
            newArr.push(d3data_flat.slice(i, i + arShp));
        }
        var selArr = newArr[cminIdx]
        im.data_source.data['image'][0] = selArr;  
        //im.data_source.trigger('change')
        im.data_source.change.emit()
        p.title.text = "bin value: "+img_idx_val+", idx "+cminIdx
        p.update()
    """
    if coords is not None and isinstance(coords, np.ndarray):
        coords = coords.tolist()
    if coords is not None and not isinstance(coords, list):
        print 'if coords is given it needs to be a list.'
        coords = None
    if coords is None:
        coords=np.arange(imOrg.shape[0]).tolist()
        coords.append(imOrg.shape[0])
        step=1
    else:
        step = (max(coords)-min(coords))/(imOrg.shape[0]-1)

    #is valStart is not index, calculate.
    if isinstance(valStart, float) or valStart>=imOrg.shape[0]:
        valStart = np.argmin(abs(coords-np.ones_like(coords)*valStart))

    source = ColumnDataSource(data=dict(im3d=imOrg))
    sourceShp = ColumnDataSource(data=dict(imShp = imOrg.shape))
    sourceCoords = ColumnDataSource(data=dict(coords = coords))

    callback_slider = bokeh.models.CustomJS(args=dict(im=im,p=p,source=source, sourceShp=sourceShp, sourceCoords=sourceCoords), 
                                            code=JS_code_slider)

    rslider = bokeh.models.Slider(title="image",start=coords[0],end=coords[-1],step=step,
                                  value=valStart, callback=callback_slider, orientation="horizontal")
    callback_slider.args['rslider'] = rslider

    return rslider
 
def plotImageBokeh(data, **kwargs):
    """
    function to plot images with bokeh
    plot will have colormap selection and scale selection
    default bokeh tools include zoom.
    Input are data (2-d array)

    Parameters
    ----------
    plotTitle: title to show on plot, default ''
    plotLog: color scale relating to log of Z-values, default False
    plotWidth: width of figure in pxls, default 600
    plotHeight: height of figure in pxls, default 400    
    plotMinP: minimum for color scale to start, default "auto" (uses np.nanpercentile,1)
    plotMaxP: maximum for color scale to start, default "auto" (uses np.nanpercentile,99)
    plotMin: minimum for color scale, default "auto" (uses np.nanmin)
    plotMax: maximum for color scale, default "auto" (uses np.nanmax)
    xRange: values for x-axis, defaults to range(image size)
    yRange: values for y-axis, defaults to range(image size)
    rangeInput: add tools to enter exact ranges of color scale in addition to slider, default True
    tools: for bokeh plotting, tools can be passed. Default tools will be used otherwise
    palette name: name of default color palette, default is jet
    datetime: x-axis data in form of date time, default False
    """
    plot_title = kwargs.pop("plot_title","")
    plotLog = kwargs.pop("plotLog",False)
    plotWidth = kwargs.pop("plotWidth",600)
    plotHeight = kwargs.pop("plotHeight",400)
    rangeInput = kwargs.pop("rangeInput",True)
    palette_name = kwargs.pop("palette_name","jet")
    plotMin = kwargs.pop("plotMin","auto")
    plotMax = kwargs.pop("plotMax","auto")
    plotMinP = kwargs.pop("plotMinP","auto")
    plotMaxP = kwargs.pop("plotMaxP","auto")
    xRange = kwargs.pop("xRange",None)
    yRange = kwargs.pop("yRange",None)
    tools = kwargs.pop("tools",None)
    dateTime = kwargs.pop("dateTime",None)
    output_quad = kwargs.pop("output_quad",False)
    if len(kwargs)>0:
        print 'found unexpected parameters to plotMarker, will ignore', kwargs

    if plotLog: data = np.log(data)
    try:
        data[np.isinf(data)]=np.nan 
        data[np.isneginf(data)]=np.nan 
    except:
        print 'could not set inf data to nan. '
    if plotMin=="auto":
        plotMin = np.nanmin(data)
    if plotMax=="auto":
        plotMax = np.nanmax(data)
    if plotMinP=="auto": 
        plotMinP = np.nanpercentile(data, 5)
    if plotMaxP=="auto": 
        plotMaxP = np.nanpercentile(data, 95)
    if xRange is None:
        xRange=(0,data.shape[0])
    if yRange is None:
        yRange=(0,data.shape[1])
    if tools is None:
        pan=PanTool()
        wheel_zoom=WheelZoomTool()
        box_zoom=BoxZoomTool()
        save=SaveTool()
        reset=ResetTool()
        hover=HoverTool(tooltips=[
            ("(x,y)","($x, $y)")
        ])
        box_select=BoxSelectTool()
        tools = [pan, wheel_zoom,box_zoom,save,hover,box_select,reset]
        if bokeh.__version__=='0.12.6':
            resize=ResizeTool()
            tools = [pan, wheel_zoom,box_zoom,save,hover,box_select,reset,resize]
    
    if dateTime:
        created_map = create_map(data2plot=data,palette_name=palette_name,
                                 fig_width_pxls=plotWidth, 
                                 fig_height_pxls=plotHeight,x_range=xRange,
                                 y_range=yRange,title=plot_title,
                                 x_axis_type="datetime", tools=tools,
                                 cmaps=cmaps,plotMin=plotMinP,plotMax=plotMaxP,
                                 output_quad=output_quad)
    else: 
        created_map = create_map(data2plot=data,palette_name=palette_name,
                                 fig_width_pxls=plotWidth, 
                                 fig_height_pxls=plotHeight,x_range=xRange,
                                 y_range=yRange,title=plot_title, tools=tools,
                                 cmaps=cmaps,plotMin=plotMinP,plotMax=plotMaxP,
                                 output_quad=output_quad)

    if output_quad:
        p, im, q = created_map
    else:
        p, im = created_map
        
    # Controls
    #vabsmin,vabsmax = plotMin, plotMax
    step=(plotMaxP-plotMinP)/100.
    range_slider = create_range_slider(plotMin,plotMax,plotMinP,plotMaxP,im=im,step=step)
    select_cm = create_cmap_selection(im,cmaps=cmaps, value=palette_name)
    if rangeInput:
        range_input = create_range_input_button(plotMin,plotMax,im=im)
        # Layout
        layout = bp.gridplot([[p],[range_slider,select_cm], [range_input[2], range_input[0], range_input[1]]])
    else:
        layout = bp.gridplot([[p],[range_slider,select_cm]])

    if output_quad:
        return layout,p,im,q
    else:
        return layout,p,im

def plot2d_from3d(data2plot, **kwargs):
    """
    class to plot images with bokeh w.
    plot will have colormap selection and scale selection
    another slider allows to select the image out of a 3-d cube
    default bokeh tools include zoom.
    Input are data (3-d array)

    Parameters
    ----------
    title: title to show on plot, default 'Binned Image Data'
    coords: coordinates on Z-axis, default None (range of cube depth)
    init_plot: value of coordinate for initial image, default None (summed image)
    plotLog: color scale relating to log of Z-values, default False
    rangeInput: add tools to enter exact ranges of color scale in addition to slider, default True
    cmaps: list of color maps, default None, meaning all defined.
    plotWidth: width of figure in pxls, default 600
    plotHeight: height of figure in pxls, default 400    
    plotMinP: minimum for color scale to start, default "auto" (uses np.nanpercentile,1)
    plotMaxP: maximum for color scale to start, default "auto" (uses np.nanpercentile,99)
    plotMin: minimum for color scale, default "auto" (uses np.nanmin)
    plotMax: maximum for color scale, default "auto" (uses np.nanmax)
    x_range: values for x-axis, defaults to range(image size)
    y_range: values for y-axis, defaults to range(image size)
    tools: for bokeh plotting, tools can be passed. Default tools will be used otherwise
    palette name: name of default color palette, default is jet
    x_axis_type: x-axis type, "linear", "log", "datetime", "auto", default linear
    min_border_left: figure borders: left, default 20
    min_border_right: figure borders: right, default 10
    min_border_top: figure borders: top, default 30
    min_border_bottom: figure borders: bottom, default 10
    title_font_size: font size of title, default 12pt
    title_align: alignment of plot title, default center
    """
    title = kwargs.pop("title","")
    plotLog = kwargs.pop("plotLog",False)
    plotWidth = kwargs.pop("plotWidth",600)
    plotHeight = kwargs.pop("plotHeight",400)
    rangeInput = kwargs.pop("rangeInput",True)
    palette_name = kwargs.pop("palette_name","jet")
    plotMin = kwargs.pop("plotMin","auto")
    plotMax = kwargs.pop("plotMax","auto")
    plotMinP = kwargs.pop("plotMinP","auto")
    plotMaxP = kwargs.pop("plotMaxP","auto")
    x_range = kwargs.pop("x_range",None)
    y_range = kwargs.pop("y_range",None)
    tools = kwargs.pop("tools",None)
    dateTime = kwargs.pop("dateTime",None)
    cmaps = kwargs.pop("cmaps",None)
    coord = kwargs.pop("coord",None)
    init_plot = kwargs.pop("init_plot",None)
    output_quad = kwargs.pop("output_quad",False)
    min_border_right = kwargs.pop("min_border_right",10)
    min_border_left = kwargs.pop("min_border_left",20)
    min_border_top = kwargs.pop("min_border_top",30)
    min_border_bottom = kwargs.pop("min_border_bottom",10)
    title_font_size = kwargs.pop("title_font_size","12pt")
    title_align = kwargs.pop("title_align","center")
    x_axis_type = kwargs.pop("x_axis_type","linear")
    if len(kwargs)>0:
        print 'found unexpected parameters to plotMarker, will ignore', kwargs

    if type(cmaps)==type(None):
        cmaps = get_all_mpl_palettes(allmaps=['jet','gray','cool','hot','seismic','CMRmap','nipy_spectral'])

    if plotLog:
        data2plot = np.log(data2plot)
        data2plot[np.isinf(data2plot)]=np.nan 
        data2plot[np.isneginf(data2plot)]=np.nan 

    #get auto scale. Use full set of images.
    if plotMin=="auto":
        plotMin = np.nanmin(data2plot)
    if plotMax=="auto":
        plotMax = np.nanmax(data2plot)
    if plotMinP=="auto": 
        plotMinP = np.nanpercentile(data2plot, 5)
    if plotMaxP=="auto": 
        plotMaxP = np.nanpercentile(data2plot, 95)

    if tools is None:
        pan=PanTool()
        wheel_zoom=WheelZoomTool()
        box_zoom=BoxZoomTool()
        save=SaveTool()
        reset=ResetTool()
        hover=HoverTool(tooltips=[
            ("(x,y)","($x, $y)")
        ])
        box_select=BoxSelectTool()
        tools = [pan, wheel_zoom,box_zoom,save,hover,box_select,reset]
        if bokeh.__version__=='0.12.6':
            resize=ResizeTool()
            tools = [pan, wheel_zoom,box_zoom,save,hover,box_select,reset,resize]

    #deal with getting (initial) 2-d image to plot
    if len(data2plot.shape)<2:
        print 'data to plot has less than 2 dims'
        return 
    elif len(data2plot.shape)>3:
        print 'data to plot has more than 3 dims'
        return
    elif len(data2plot.shape)==2:
        init_plot=-1
        init_dat = data2plot
    else:
        if init_plot is not None:
            initIdx=None
            if coord is not None:
                if isinstance(init_plot, int):
                    if (init_plot>=data2plot.shape[0] or init_plot<0):
                        initIdx=0
                    else:
                        initIdx = init_plot
                else:
                    initIdx = np.argmin(abs(np.array(coord)-init_plot))
            if not isinstance(initIdx, int):
                if not isinstance(init_plot, int):
                    print 'init_plot needs to be integer, using z-axis to be implemented later, will start at first image'
                    initIdx = 0
                else:
                    initIdx = init_plot
        else:
            initIdx = 0
        init_dat = data2plot[initIdx]
        
    if coord is None or (len(data2plot.shape)==3 and data2plot.shape[0]!=len(coord)):
        coord = np.arange(0,data2plot.shape[0])              
            
    #plot range X&Y
    if x_range is None:
        x0=0; x1=init_dat.shape[0]
    else:
        x_range = np.array(x_range)
        if x_range.shape[0]==1:
            x0=min(x_range,0); x1 = max(x_range,0)
        else:
            x0=min(x_range); x1 = max(x_range)

    if y_range is None:
        y0=0
        y1=init_dat.shape[1]
    else:
        y_range = np.array(y_range)
        if y_range.shape[0]==1:
            y0=min(y_range,0); y1 = max(y_range,0)
        else:
            y0=min(y_range); y1 = max(y_range)

    imgSource = ColumnDataSource(
            {'value': init_dat})
    #create figure.
    p = bp.figure(x_range=(x0, x1), y_range=(y0, y1),x_axis_type=x_axis_type,
                  plot_width=plotWidth,plot_height=plotHeight, 
                  min_border_left=min_border_left,min_border_right=min_border_right,
                  title=title,min_border_top=min_border_top,min_border_bottom=min_border_bottom,
                  tools= tools)
    p.title.text_font_size = title_font_size
    p.title.align = title_align
    im = p.image(image=[init_dat],dw=[x1-x0],dh=[y1-y0],x=[x0],y=[y0],palette=cmaps["palettes_dict"][palette_name])

    im.glyph.color_mapper.high = plotMax
    im.glyph.color_mapper.low = plotMin
    imquad = p.quad(top=[y1], bottom=[y0], left=[x0], right=[x1],alpha=0) # This is used for hover and taptool
    
    color_bar = bokeh.models.ColorBar(color_mapper=im.glyph.color_mapper, label_standoff=12, location=(0,0))
    p.add_layout(color_bar, 'right')

    #create more tools.
    #colormap selection
    select_cm = create_cmap_selection(im,cmaps=cmaps, value=palette_name)

    #range slider.
    step=(plotMaxP-plotMinP)/50.
    #vabsmin,vabsmax = plotMin, plotMax
    range_slider = create_range_slider(np.nanmin(data2plot),np.nanmax(data2plot),plotMinP,plotMaxP,im=im,step=step)

    if rangeInput:
        range_input = create_range_input_button(plotMin,plotMax,im=im)

    if len(data2plot.shape)==3:
        img_slider = create_img_slider_scale(im, data2plot,valStart=init_plot,coords=coord,p=p)
        #put them togeter.
        if rangeInput:
            layout = bp.gridplot([[p],[range_slider,select_cm,img_slider], [range_input[2], range_input[0], range_input[1]]])
        else:
            layout = bp.gridplot([[p],[range_slider,select_cm,img_slider]])
    else:
        if rangeInput:
            layout = bp.gridplot([[p],[range_slider,select_cm], [range_input[2], range_input[0], range_input[1]]])
        else:
            layout = bp.gridplot([[p],[range_slider,select_cm]])

    if output_quad:
        return layout, p,im,imquad
    else:
        return layout, p,im


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
            color_names.append(sorted_names[min(i*dcol, len(sorted_names)-1)])
        return color_names

def plotMarker(data, **kwargs):
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
    fig: if using matplotlib, you can pass a gridspec for the figure
         in which the data will be plotted
         if passing anything when using bokeh, the figure will be returned
    runLabel: label to indicate the run ('Run001'), def 'Run', can be list
              is passing array
    plotTitle: title to show on plot
    plotDirname: directory to which bokeh html files will be written (def ./)
    marker: pass marker style (def o, options: *,^,+,x,s,D)
    markersize: marker size, default is 5
    line_dash: connect markers with line style as defined (default: none)
    color: pass (list of) colors, otherwise a set of colors will be automatically picked
    tools: for bokeh plotting, tools can be passed. Default tools will be used
           otherwise
    ylim: limits for y-axis, automatically chosen by matplotlib/bokeh otherwise
    width_height: measurements for to be created figures
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
    color = kwargs.pop("color", None)
    if color is not None and len(color) != len(data):
        color = None
    if color is None and len(data)>1:
        print 'ndataset ',len(data)
        sorted_colors = getColors(nColors=len(data))
    ylim = kwargs.pop("ylim", None)
    width_height = kwargs.pop("width_height", None)
    if len(kwargs)>1 or (len(kwargs)==1 and kwargs.keys()[0]!='fig'):
        print 'found unexpected parameters to plotMarker, will ignore', kwargs

    if plotWith.find('matplotlib')>=0:
        if 'fig' in kwargs.keys():
            fig = kwargs.pop("fig")
        else:
            if width_height is None: width_height = (8,5)
            plt.figure(figsize=width_height)
            fig = (gridspec.GridSpec(1,1)[0])
        if line_dash is None:
            line_dash = 'none'

        for ic,(idata,ixdata) in enumerate(zip(data, xData)):
            if isinstance(marker, list) and len(marker)==len(data):
                plotMarker = marker[ic]
            else:
                plotMarker = marker
            if color is None and ic == 0:
                plt.subplot(fig).plot(ixdata,idata,marker=plotMarker,linestyle=line_dash,markersize=markersize)
            elif color is None:
                plt.subplot(fig).plot(ixdata,idata,marker=plotMarker, color=sorted_colors[ic],linestyle=line_dash)
            else:
                plt.subplot(fig).plot(ixdata,idata,marker=plotMarker, color=color[ic],linestyle=line_dash)

        plt.subplot(fig).set_xlabel(xLabel)
        plt.subplot(fig).set_ylabel(yLabel)
        plt.subplot(fig).set_title(plotTitle)
        if ylim is not None:
            plt.subplot(fig).set_ylim(np.nanmin(np.array(ylim)),np.nanmax(np.array(ylim)))
        if plotWith.find('file')>=0:
            plt.savefig('%s.jpeg'%(plotFilename))

    elif plotWith.find('bokeh')>=0:
        #IMPLEMENT ME
        #bokeh, use BoxSelectTool to get selected region
        #https://stackoverflow.com/questions/34164587/get-selected-data-contained-within-box-select-tool-in-bokeh
        #set cuts does not need to be specified explicly, can simple be gotten from callback of tool.....
        if 'fig' in kwargs.keys():
            fig = kwargs.pop("fig")
        else:
            fig = None
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
        if width_height is None:
            p = bp.figure(title=plotTitle, x_axis_label=xLabel, y_axis_label=yLabel,tools=tools)
        else:
            p = bp.figure(title=plotTitle, x_axis_label=xLabel, y_axis_label=yLabel,tools=tools, width=width_height[0], height=width_height[1])
        for ic,(idata,ixdata) in enumerate(zip(data, xData)):
            if color is not None:
                thiscolor=color[ic]
            elif len(data)==1:
                thiscolor='black'
            else:
                thiscolor=sorted_colors[ic]
            if isinstance(RunLabel, list):
                iRunLabel=RunLabel[ic]
            else:
                iRunLabel=runLabel
                if len(data)==1:
                    iRunLabel=None
            if isinstance(marker, list) and len(marker)==len(data):
                plotMarker = marker[ic]
            else:
                plotMarker = marker
            if plotMarker=='o':
                p.circle(ixdata, idata, legend=iRunLabel, size=markersize, color=thiscolor, name='p%d'%ic)
            elif plotMarker=='*':
                p.asterisk(ixdata, idata, legend=iRunLabel, size=markersize, color=thiscolor, name='p%d'%ic)
            elif plotMarker=='+':
                p.cross(ixdata, idata, legend=iRunLabel, size=markersize, color=thiscolor, name='p%d'%ic)
            elif plotMarker=='x':
                p.x(ixdata, idata, legend=iRunLabel, size=markersize, color=thiscolor, name='p%d'%ic)
            elif plotMarker=='s':
                p.square(ixdata, idata, legend=iRunLabel, size=markersize, color=thiscolor, name='p%d'%ic)
            elif plotMarker=='^':
                p.triangle(ixdata, idata, legend=iRunLabel, size=markersize, color=thiscolor, name='p%d'%ic)
            elif plotMarker=='D':
                p.diamond(ixdata, idata, legend=iRunLabel, size=markersize, color=thiscolor, name='p%d'%ic)
            else:
                print 'unsupported marker option'
            if line_dash is not None:
                p.line(ixdata, idata, line_dash=line_dash, line_color=thiscolor, name='p%d'%ic)
        if plotWith=='bokeh_notebook':
            bp.output_notebook()
            if fig is None:
                bp.show(p)
        else:
            bp.output_file('%s.html'%(plotFilename))
            bp.save(p)
        if fig is not None:
            return p
                
    elif plotWith != 'no_plot':
        print 'plotting using %s is not implemented yet, options are matplotlib, bokeh_notebook, bokeh_html or no_plot'%plotWith
            
    return



def plotImage(image, **kwargs):
    """
    class to plot 1-d data (histograms, scans).
    Input are xData (1-d array, optionally list of arrays)

    Parameters
    ----------
    xLabel: string to be used for labelling a X-axis
    yLabel: string to be used for labelling the Y-axis
    extend: x,y ranges for plot, form: x0, y0, x1, y1. Use bin counter as default
    plotWith: method for plotting, options are matplotlib, 
              bokeh_notebook or bokeh w/ output as a html file
    fig: if using matplotlib, you can pass a gridspec for the figure
         in which the data will be plotted
         if passing anything when using bokeh, the figure will be returned
    runLabel: label to indicate the run ('Run001'), def 'Run', can be list
              is passing array
    plotTitle: title to show on plot
    plotDirname: directory to which bokeh html files will be written (def ./)
    tools: for bokeh plotting, tools can be passed. Default tools will be used
           otherwise
    ylim: limits for z-axis, automatically chosen by matplotlib/bokeh otherwise
    width_height: measurements for to be created figures
    dateTime: default False, True will make the plot use time coords on x-axis, bokeh only (?)
    output_quad: default False, bokeh only, determines what of layout is return
    plotLog: default False, bokeh only
    """
    
    if not isinstance(image, np.ndarray):
        image = np.array(image)
    if len(image.shape)!=2:
        print 'plotImage need to be given a 2-d array, got this shape',image.shape,' instead, will quit'
        return

    xLabel = kwargs.pop("xLabel", "")
    yLabel = kwargs.pop("yLabel", "")
    extent = kwargs.pop("extent", [0,image.shape[0],0,image.shape[1]])
    plotWith = kwargs.pop("plotWith", "matplotlib")
    RunLabel = kwargs.pop("runLabel", "Run")
    if isinstance(RunLabel, list):
        runLabel = RunLabel[0]
    else:
        runLabel = RunLabel
    plotTitle = kwargs.pop("plotTitle","%s vs %s for %s"%(xLabel, yLabel, runLabel))
    tools = kwargs.pop("tools", None)
    plotDirname = kwargs.pop("plotDirname", "./") 
    plotFilename =  kwargs.pop("plotFilename", '%s/%s_%s_vs_%s'%(plotDirname,runLabel, xLabel.replace('/','_'), yLabel.replace('/','_')))
    ylim = kwargs.pop("ylim", None)
    width_height = kwargs.pop("width_height", None)
    dateTime = kwargs.pop("dateTime", False)
    output_quad = kwargs.pop("output_quad", False)
    plotLog = kwargs.pop("plotLog", False)
    rangeInput = kwargs.pop("rangeInput", True)
    if len(kwargs)>1 or (len(kwargs)==1 and kwargs.keys()[0]!='fig'):
        print 'found unexpected parameters to plotMarker, will ignore', kwargs

    if plotWith.find('matplotlib')>=0:
        fig = None
        if 'fig' in kwargs.keys():
            fig = kwargs.pop("fig")
        if fig is None:
            if width_height is None: width_height = (8,5)
            plt.figure(figsize=width_height)
            fig = (gridspec.GridSpec(1,1)[0])
        if ylim is None:
            ylim = [np.nanpercentile(image,1), np.nanpercentile(image,99.5)]
        plt.subplot(fig).imshow(image, clim=ylim, extent=extent, aspect='auto', interpolation='none',origin='lower')
        plt.subplot(fig).set_xlabel(xLabel)
        plt.subplot(fig).set_ylabel(yLabel)
        plt.subplot(fig).set_title(plotTitle)

    elif plotWith.find('bokeh')>=0:
        if width_height is None: width_height = (600,400)
        if ylim is None:
            ylim=['auto', 'auto']
        fig = kwargs.pop("fig", None)
        #retValues are layout,p,im[,q] (latter if output_quad is True)
        retValues = plotImageBokeh(image, plotWidth=width_height[0], plotHeight=width_height[1], palette_name='jet', output_quad=output_quad, tools=tools,  plot_title=plotTitle, dateTime=dateTime, xRange=(extent[0], extent[1]), yRange=(extent[2], extent[3]), plotMinP="auto", plotMaxP="auto", plotMin=ylim[0], plotMax=ylim[1],plotLog=plotLog,rangeInput=rangeInput)

        if fig is not None:
            return retValues
        if output_quad:
            layout,p,im,q = retValues
        else:
            layout,p,im = retValues
        if plotWith=='bokeh_notebook':
            bp.output_notebook()
            #bp.show(p)
            bp.show(layout)
        else:
            bp.output_file('%s/%s_%s_vs_%s.html'%(self.plot_dirname,self.runLabel, plotvars[0].replace('/','_'), plotvars[1].replace('/','_')))
            #bp.save(p)
            bp.save(layout)

    elif plotWith != 'no_plot':
        print 'plotting using %s is not implemented yet, options are matplotlib, bokeh_notebook, bokeh_html or no_plot'%plotWith
            
    return


def plot3d_img_time(data2plot, **kwargs):
    """
    class to plot images with bokeh w. slider on 1st axis (for cube display)
    plot ROI against coords.
    Input are data (3-d array)

    Parameters
    ----------
    title: title to show on plot, default 'Binned Image Data'
    coords: coordinates on Z-axis, default None (range of cube depth)
    init_plot: value of coordinate for initial image, default None (summed image)
    plotLog: color scale relating to log of Z-values, default False
    rangeInput: add tools to enter exact ranges of color scale in addition to slider, default True
    cmaps: list of color maps, default None, meaning all defined.
    plotWidth: width of figure in pxls, default 600
    plotHeight: height of figure in pxls, default 400    
    plotMinP: minimum for color scale to start, default "auto" (uses np.nanpercentile,1)
    plotMaxP: maximum for color scale to start, default "auto" (uses np.nanpercentile,99)
    plotMin: minimum for color scale, default "auto" (uses np.nanmin)
    plotMax: maximum for color scale, default "auto" (uses np.nanmax)
    x_range: values for x-axis, defaults to range(image size)
    y_range: values for y-axis, defaults to range(image size)
    tools: for bokeh plotting, tools can be passed. Default tools will be used otherwise
    palette name: name of default color palette, default is jet
    x_axis_type: x-axis type, "linear", "log", "datetime", "auto", default linear
    min_border_left: figure borders: left, default 20
    min_border_right: figure borders: right, default 10
    min_border_top: figure borders: top, default 30
    min_border_bottom: figure borders: bottom, default 10
    title_font_size: font size of title, default 12pt
    title_align: alignment of plot title, default center
    """
    title = kwargs.pop("title","")
    plotLog = kwargs.pop("plotLog",False)
    plotWidth = kwargs.pop("plotWidth",600)
    plotHeight = kwargs.pop("plotHeight",400)
    rangeInput = kwargs.pop("rangeInput",True)
    palette_name = kwargs.pop("palette_name","jet")
    plotMin = kwargs.pop("plotMin","auto")
    plotMax = kwargs.pop("plotMax","auto")
    plotMinP = kwargs.pop("plotMinP","auto")
    plotMaxP = kwargs.pop("plotMaxP","auto")
    x_range = kwargs.pop("x_range",None)
    y_range = kwargs.pop("y_range",None)
    tools = kwargs.pop("tools",None)
    dateTime = kwargs.pop("dateTime",None)
    cmaps = kwargs.pop("cmaps",None)
    coord = kwargs.pop("coord",None)
    init_plot = kwargs.pop("init_plot",None)
    output_quad = kwargs.pop("output_quad",False)
    min_border_right = kwargs.pop("min_border_right",10)
    min_border_left = kwargs.pop("min_border_left",20)
    min_border_top = kwargs.pop("min_border_top",30)
    min_border_bottom = kwargs.pop("min_border_bottom",10)
    title_font_size = kwargs.pop("title_font_size","12pt")
    title_align = kwargs.pop("title_align","center")
    x_axis_type = kwargs.pop("x_axis_type","linear")
    if len(kwargs)>0:
        print 'found unexpected parameters to plotMarker, will ignore', kwargs

    if type(cmaps)==type(None):
        cmaps = get_all_mpl_palettes(allmaps=['jet','gray','cool','hot','seismic','CMRmap','nipy_spectral'])

    if tools is None:
        pan=PanTool()
        wheel_zoom=WheelZoomTool()
        box_zoom=BoxZoomTool()
        save=SaveTool()
        reset=ResetTool()
        hover=HoverTool(tooltips=[
            ("(x,y)","($x, $y)")
        ])
        box_select=BoxSelectTool()
        tools = [pan, wheel_zoom,box_zoom,save,hover,box_select,reset]
        if bokeh.__version__=='0.12.6':
            resize=ResizeTool()
            tools = [pan, wheel_zoom,box_zoom,save,hover,box_select,reset,resize]

    #get auto scale. Use full set of images.
    if plotMin=="auto":
        plotMin = np.nanmin(data2plot)
    if plotMax=="auto":
        plotMax = np.nanmax(data2plot)

    #deal with getting (initial) 2-d image to plot
    if len(data2plot.shape)<2:
        print 'data to plot has less than 2 dims'
        return 
    elif len(data2plot.shape)>3:
        print 'data to plot has more than 3 dims'
        return
    elif len(data2plot.shape)==2:
        init_plot=-1
        init_dat = data2plot
    else:
        if init_plot is not None:
            if not isinstance(init_plot, int):
                print 'init_plot needs to be integer, using z-axis to be implemented later, will start at first image'
                init_plot = 0
        else:
            init_plot = 0
        print 'DEBUG: start w/ plot ',init_plot
        init_dat = data2plot[init_plot]
        
    if coord is None or (len(data2plot.shape)==3 and data2plot.shape[0]!=len(coord)):
        coord = np.arange(0,data2plot.shape[0])              
            
    #plot range X&Y
    if x_range is None:
        x0=0; x1=init_dat.shape[0]
    else:
        x_range = np.array(x_range)
        if x_range.shape[0]==1:
            x0=min(x_range,0); x1 = max(x_range,0)
        else:
            x0=min(x_range); x1 = max(x_range)

    if y_range is None:
        y0=0
        y1=init_dat.shape[1]
    else:
        y_range = np.array(y_range)
        if y_range.shape[0]==1:
            y0=min(y_range,0); y1 = max(y_range,0)
        else:
            y0=min(y_range); y1 = max(y_range)

    imgSource = ColumnDataSource(
            {'value': init_dat})

    data1d = data2plot.copy()
    while len(data1d.shape)>1:
        data1d=np.nanmean(data1d,axis=1)
    bin1d = coord #np.arange(data1d.shape[0])

    #create figure.
    p = bp.figure(x_range=(x0, x1), y_range=(y0, y1),x_axis_type=x_axis_type,
                              plot_width=plotWidth,plot_height=plotHeight, 
                              min_border_left=min_border_left,min_border_right=min_border_right,
                              title=title,min_border_top=min_border_top,min_border_bottom=min_border_bottom,
                              tools= tools)
    p.title.text_font_size = title_font_size
    p.title.align = title_align
    im = p.image(image=[init_dat],dw=[x1-x0],dh=[y1-y0],x=[x0],y=[y0],palette=cmaps["palettes_dict"][palette_name])
    p1d = plotMarker(data1d, xData=bin1d, fig='return_me',plotWith='bokeh_notebook',
                     width_height=(plotWidth,plotHeight),plotTitle='ROI vs scan')

    im.glyph.color_mapper.high = plotMax
    im.glyph.color_mapper.low = plotMin
    imquad = p.quad(top=[y1], bottom=[y0], left=[x0], right=[x1],alpha=0) # This is used for hover and taptool
    
    color_bar = bokeh.models.ColorBar(color_mapper=im.glyph.color_mapper, label_standoff=12, location=(0,0))
    p.add_layout(color_bar, 'right')

    #p.add_tools(box_select)
    source = ColumnDataSource(data=dict(im3d=data2plot))
    sourceShp = ColumnDataSource(data=dict(imShp = data2plot.shape))

    callback = bokeh.models.CustomJS(args=dict(p1d=p1d, sourceShp=sourceShp, source=source), code="""
        /// get BoxSelectTool dimensions from cb_data parameter of Callback
        var geometry = cb_data['geometry'];
        var arShp = sourceShp.data.imShp[1]*sourceShp.data.imShp[2]
        var d3data_flat = source.data.im3d
        nPixROI = (Math.round(geometry['y1']) - Math.round(geometry['y0']))*(Math.round(geometry['x1']) - Math.round(geometry['x0']))
        ROIsums=[]

        var ir = 0;
        for(var i = 0; i < d3data_flat.length; i += arShp) {
            loc_image = d3data_flat.slice(i, i + arShp)
            var sumTest = loc_image.reduce(function(a, b) { return a + b; }, 0);

            ir = 0; //reset on each (flattened) image
            var debug_ir_lastRow = -1
            var local_ROIsum=0

            for(var irc = 0; irc < loc_image.length; irc += sourceShp.data.imShp[1]) { //split flat ar into rows
                if ( ( ir >= Math.round(geometry['y0'])) && (ir < Math.round(geometry['y1'])) ) { //select rows in ROI
                    debug_ir_lastRow = ir
                    thisSlice = loc_image.slice(irc+Math.round(geometry['x0']), irc+Math.round(geometry['x1']))
                    var sum = thisSlice.reduce(function(a, b) { return a + b; }, 0);

                    local_ROIsum += sum
                }
                ir = ir + 1
 
            }
            ROIsums.push(local_ROIsum/nPixROI)
        }
        p1d.select(name='p0')[0].data_source.data['y'] = ROIsums
        p1d.select(name='p0')[0].data_source.change.emit()
        p1d.title.text = "ROI: X \"+Math.round(geometry['x0'])+":"+Math.round(geometry['x1'])+" Y "+Math.round(geometry['y0'])+":"+Math.round(geometry['y1'])
        p1d.update()
        """)

    box_select = BoxSelectTool(callback=callback)
    p.add_tools(box_select)

    #colormap selection
    select_cm = create_cmap_selection(im,cmaps=cmaps, value=palette_name)

    #range slider.
    if plotMinP=="auto": 
        plotMinP = np.nanpercentile(data2plot, 5)
    if plotMaxP=="auto": 
        plotMaxP = np.nanpercentile(data2plot, 95)
    plotMin = np.nanmin(data2plot)
    plotMax = np.nanmax(data2plot)
    step=(plotMaxP-plotMinP)/50.
    range_slider = create_range_slider(plotMin, plotMax,plotMinP,plotMaxP,im=im,step=step)
    if rangeInput:
        range_input = create_range_input_button(plotMin,plotMax,im=im)

    if len(data2plot.shape)==3:
        img_slider = create_img_slider_scale(im, data2plot,valStart=init_plot,coords=coord,p=p)
        if rangeInput:
            layout = bp.gridplot([[p],[range_slider,select_cm,img_slider],[range_input[2], range_input[0], range_input[1]],[p1d]])
        else:
            layout = bp.gridplot([[p],[range_slider,select_cm,img_slider],[p1d]])
    else:
        if rangeInput:
            layout = bp.gridplot([[p],[range_slider,select_cm],[range_input[2], range_input[0], range_input[1]],[p1d]])
        else:
            layout = bp.gridplot([[p],[range_slider,select_cm],[p1d]])

    if output_quad:
        return layout, p,im,p1d,imquad
    else:
        return layout, p,im,p1d
