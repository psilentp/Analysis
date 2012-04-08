__author__ = 'psilentp'

from numpy import array
from heka_io import HekaIO,gbi
from chaco.example_support import COLOR_PALETTE
from enable.example_support import DemoFrame, demo_main
# Enthought library imports
from enable.api import Window, Component, ComponentEditor
from traits.api import HasTraits, Instance
from enthought.traits.ui.api import Item, Group, View
# Chaco imports
from chaco.api import create_line_plot, add_default_axes,\
    add_default_grids, OverlayPlotContainer, PlotLabel,\
    create_scatter_plot, Legend
from chaco.tools.api import PanTool, ZoomTool, LegendTool,\
    TraitsTool, DragZoom
from read_heka import *

#===============================================================================
# # Create the Chaco plot.
#===============================================================================
def _create_plot_component():
    container = OverlayPlotContainer(padding = 50, fill_padding = True,
        bgcolor = "lightgray", use_backbuffer=True)

    fprefix = './test_data/CEN184/'
    #filename1 = fprefix + 'THL_2012-03-21_18-40-42_000.dat'
    #filename2  = fprefix + 'THL_2012-03-21_18-44-42_000.dat'
    filename = './test_data/CEN111/THL_2011-07-09_15-02-54_000.dat'
    ioreader = HekaIO(filename)
    #read a block
    blo = ioreader.read_block(group = 4)

    f = open(filename)
    head = BundleHeader(f)
    head.load(f)
    bi = head.oBundleItems[2]
    pgf = PGFFile(f,bi)

    value_mapper = None
    index_mapper = None
    plots = {}

    firstplot = True
    for seg in blo.segments:
        #prococol building
        #print seg.annotations
        for a_sig in seg.analogsignals:
            x = array(a_sig.times)
            y = array(a_sig)
            ch = int(a_sig.annotations['trSourceChannel'])
            plot = create_line_plot((x,y), width=0.5,color=tuple(COLOR_PALETTE[ch]))
            plot.index.sort_order = "ascending"
            plot.bgcolor = "white"
            plot.border_visible = True

            #code for protocols
            print "###########################"
            print pgf.tree['children'][a_sig.annotations['pgf_index']]['children'][0]['children'][1]['contents'].seVoltage
            for key in ['pgf_index','trTraceCount','trAdcChannel','trSourceChannel','swStimCount']:
                print "%s:%s"%(key,a_sig.annotations[key])

            if not firstplot:
                plot.value_mapper = value_mapper
                value_mapper.range.add(plot.value)
                plot.index_mapper = index_mapper
                index_mapper.range.add(plot.index)
            else:
                value_mapper = plot.value_mapper
                index_mapper = plot.index_mapper
                add_default_grids(plot)
                add_default_axes(plot)
                plot.index_range.tight_bounds = False
                plot.index_range.refresh()
                plot.value_range.tight_bounds = False
                plot.value_range.refresh()
                plot.tools.append(PanTool(plot))
                # The ZoomTool tool is stateful and allows drawing a zoom
                # box to select a zoom region.
                zoom = ZoomTool(plot, tool_mode="box", always_on=False)
                plot.overlays.append(zoom)
                # The DragZoom tool just zooms in and out as the user drags
                # the mouse vertically.
                dragzoom = DragZoom(plot, drag_button="right")
                plot.tools.append(dragzoom)
                # Add a legend in the upper right corner, and make it relocatable
                legend = Legend(component=plot, padding=10, align="ur")
                legend.tools.append(LegendTool(legend, drag_button="right"))
                #print a_sig.annotations
                plot.overlays.append(legend)
                firstplot = False
            container.add(plot)
            plots["sweep %s"%a_sig.annotations['trLabel'][:4]] = plot
            # Set the list of plots on the legend
    legend.plots = plots
    # Add the title at the top
    container.overlays.append(PlotLabel(blo.annotations['grLabel'],
        component=container,
        font = "swiss 16",
        overlay_position="top"))
    # Add the traits inspector tool to the container
    container.tools.append(TraitsTool(container))

    return container
#===============================================================================
# Attributes to use for the plot view.
size=(800,700)
title="Simple Line Plot"
#===============================================================================
# # Demo class that is used by the demo.py application.
#===============================================================================
class Demo(HasTraits):
    plot = Instance(Component)

    traits_view = View(
        Group(
            Item('plot', editor=ComponentEditor(size=size),
                show_label=False),
            orientation = "vertical"),
        resizable=True, title=title,
        width=size[0], height=size[1]
    )

    def _plot_default(self):
        return _create_plot_component()

demo = Demo()
#===============================================================================
# Stand-alone frame to display the plot.
#===============================================================================
class PlotFrame(DemoFrame):
    def _create_window(self):
        # Return a window containing our plots
        return Window(self, -1, component=_create_plot_component())

if __name__ == "__main__":
    demo_main(PlotFrame, size=size, title=title)
    #print gbi()
    # EOF