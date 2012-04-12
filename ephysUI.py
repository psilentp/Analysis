from traits.api import HasTraits, Instance
from traitsui.api import View, Item
from chaco.tools.api import PanTool, ZoomTool
from chaco.api import Plot, ArrayPlotData, VPlotContainer,GridPlotContainer
from enable.component_editor import ComponentEditor
from numpy import cos, linspace, sin
from heka_io import *
import numpy as np


###In order to meld the axonio with the hekaio I will need to make get_block return
###data structured in the same way for both file-formats. Maybe it would be better to create an
###Irregularly Sampled signal to represent the protocol.
def get_block(block_num = 1):
    def get_stimtrace(epochs,channel):
        times = []
        vms = []
        for ep in epochs:
            if ep.annotations['chDacChannel'] == channel:
                times.append(float(ep.time))
                vms.append(float(ep.annotations['value']))
                times.append(float(ep.time)+float(ep.duration))
                vms.append(ep.annotations['value'])
        return {'x':np.array(times),'y':np.array(vms)}
    #filename = './test_data/CEN184/THL_2012-03-21_18-40-42_000.dat'
    filename = './test_data/CEN184/THL_2012-03-21_18-44-42_000.dat'
    #filename = './test_data/CEN111/THL_2011-07-09_15-02-54_000.dat'
    ioreader = HekaIO(filename)
    blo = ioreader.read_block(group = block_num)

    protocol_list = list()
    sweep_list = list()
    max_channels = 0
    for seg in blo.segments:
        chnl_map = dict()
        for ep_index, ep in enumerate(seg.epochs):
            if ep.label == 'protocol_epoch':
                chnl_map.update({ep.annotations['chDacChannel']:
                                get_stimtrace(seg.epochs,ep.annotations['chDacChannel'])})
            if len(chnl_map.keys()) > max_channels:
                max_channels = len(chnl_map.keys())
        protocol_list.append(chnl_map)
        for a_sig in seg.analogsignals:
            sweep_list.append({'x':np.array(a_sig.times),'y':np.array(a_sig)})
    return(sweep_list,protocol_list,max_channels)

class BlockPlot(HasTraits):
    protocols_plots = Instance(VPlotContainer)
    signals_plots = Instance(VPlotContainer)

    traits_view = View(
        Item('protocols_plots',editor=ComponentEditor(size=(300,20))),
        Item('signals_plots',editor=ComponentEditor(size=(300,400))),
            width=700, height=700, resizable=True, title="Chaco Plot")


    def __init__(self):
        super(BlockPlot, self).__init__()
        sweeps, protocols,num_channels = get_block(block_num=3)
        datasourses = list()
        ### create the sweeps plot
        x = sweeps[0]['x']
        sweep_data = {'x':x}
        sweep_keys = list()
        for index,item in enumerate(sweeps):
            sweep_data.update({'x'+str(index):item['x']})
            sweep_data.update({'y'+str(index):item['y']})
            sweep_keys.append(('x'+str(index),'y'+str(index)))
        sweep_plot_data = ArrayPlotData(**sweep_data)
        datasourses.append(sweep_plot_data)
        sweeps = Plot(datasourses[0],padding = 30)
        sweeps.padding_bottom=60
        for key in sweep_keys:
            if not key == 'x':
                sweeps.plot(key, type="line", color="blue")
        sweeps.border_visible = False
        zoom = ZoomTool(sweeps,tool_mode="box", always_on=False)
        sweeps.tools.append(zoom)
        sweeps.overlays.append(zoom)
        sweeps.tools.append(PanTool(sweeps))
        #plot_list = [sweeps]
        plot_list = list()
        ### create the protocol plot... list of plots for multiple stimulus channels?

        for chnum in range(num_channels):
            chkey = protocols[0].keys()[chnum]
            x = protocols[0][chkey]['x']
            prot_data = {'x':x}
            for index, item in enumerate(protocols):
                prot_data.update({str(index):item[chkey]['y']})
            prot_plot_data = ArrayPlotData(**prot_data)
            datasourses.append(prot_plot_data)
            prots = Plot(datasourses[chnum+1],padding = 30 )
            prots.padding_bottom=5
            for key in prot_data.keys():
                if not key == 'x':
                    prots.plot(('x',key),type = 'line',color = 'blue')
            prots.border_visible = False
            prots.x_axis.visible = False
            prots.index_range = sweeps.index_range

            plot_list.append(prots)
        self.signals_plots = VPlotContainer(sweeps)
        self.protocols_plots = VPlotContainer(*plot_list)


if __name__ == "__main__":
    BlockPlot().configure_traits()
