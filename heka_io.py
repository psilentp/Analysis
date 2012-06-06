from __future__ import absolute_import
import neo
from neo.io.baseio import BaseIO
from neo.core import Block, Segment, AnalogSignal, EventArray
from neo.io.tools import create_many_to_one_relationship
import numpy as np
import quantities as pq
from read_heka import *

def gbi():
    def get_stimtrace(epochs):
        times = []
        vms = []
        for ep in epochs:
            if ep.annotations['channel'] == 3:
                times.append(float(ep.time))
                vms.append(float(ep.annotations['value']))
                times.append(float(ep.time)+float(ep.duration))
                vms.append(ep.annotations['value'])
        return (np.array(times),np.array(vms))
    import pylab as plb
    #filename = './test_data/CEN184/THL_2012-03-21_18-40-42_000.dat'
    #filename = './test_data/CEN184/THL_2012-03-21_18-44-42_000.dat'
    #filename = './test_data/CEN111/THL_2011-07-09_15-02-54_000.dat'
    filename = './test_data/CEN189/THL_2012-05-05_03-47-10_000.dat'
    ioreader = HekaIO(filename)
    blo = ioreader.read_block(group = 2)
    for seg in blo.segments:
        ax1 = plb.subplot(2,1,1)
        x,y = get_stimtrace(seg.epochs)
        plb.plot(x,y,'o-')
        plb.subplot(2,1,2,sharex = ax1)
        for a_sig in seg.analogsignals:
            x = np.array(a_sig.times)
            y = np.array(a_sig)
            plb.plot(x,y,)
    plb.show()

class HekaIO(BaseIO):
    is_readable = True
    is_writable = False
    supported_objects = [Segment,AnalogSignal,Block,EventArray]
    readable_objects = [Segment,AnalogSignal,Block]
    writeable_objects = []
    has_header = True
    is_streameable = False
    writeable_objects = []

    read_params = {
        Segment : [
            ('segment_duration',
                {'value': 15.0, 'label' : 'Segment size (s.)'}),
            ('num_analogsignal',
                {'value': 8, 'label': 'Number of recording points'})
        ]
    }
    write_params = None
    name = 'example'
    extentions = ['nof']
    mode = 'file'

    def __init__(self,filename = './test_data/CEN184/THL_2012-03-21_18-40-42_000.dat'):
        BaseIO.__init__(self)
        self.filename = filename
        f = open(filename)
        #create a bundle header object with the file
        head = BundleHeader(f)
        #load the file
        head.load(f)
        #get the .pgf and .pul items in the file
        for bi in head.oBundleItems:
            if str(bi.oExtension)[0:4] == '.pgf':
                self.pgf = PGFFile(f,bi)
            if str(bi.oExtension)[0:4] == '.pul':
                self.pul = PULFile(f,bi)
        f.close()

    def read_block(self,
                   lazy = False,
                   cascade = True,
                   group = 0):
        blo = Block(name = 'test')
        if cascade:
            tree = getbyroute(self.pul.tree,[0,group])
            for i,child in enumerate(tree['children']):
                blo.segments.append(self.read_segment(group=group,series = i))
            annotations = tree['contents'].__dict__.keys()
            annotations.remove('readlist')
            for a in annotations:
                d = {a:str(tree['contents'].__dict__[a])}
                blo.annotate(**d)
        create_many_to_one_relationship(blo)
        return blo

    def read_segment(self,
                    lazy = False,
                    cascade = True,
                    group = 0,
                    series = 0):
        seg = Segment( name = 'test')
        if cascade:
            tree = getbyroute(self.pul.tree,[0,group,series])
            for sw,sweep in enumerate(tree['children']):
                if sw == 0:
                    starttime = pq.Quantity(float(sweep['contents'].swTimer),'s')
                for ch,channel in enumerate(sweep['children']):
                    sig = self.read_analogsignal(group=group,
                                            series=series,
                                            sweep=sw,
                                            channel = ch)
                    annotations = sweep['contents'].__dict__.keys()
                    annotations.remove('readlist')
                    for a in annotations:
                        d = {a:str(sweep['contents'].__dict__[a])}
                        sig.annotate(**d)
                    sig.t_start = pq.Quantity(float(sig.annotations['swTimer']),'s') - starttime
                    seg.analogsignals.append(sig)
            annotations = tree['contents'].__dict__.keys()
            annotations.remove('readlist')
            for a in annotations:
                d = {a:str(tree['contents'].__dict__[a])}
                seg.annotate(**d)
        create_many_to_one_relationship(seg)
        ### add protocols to signals
        for sig_index,sig in enumerate(seg.analogsignals):
            pgf_index = sig.annotations['pgf_index']
            st_rec = self.pgf.tree['children'][pgf_index]['contents']
            chnls = [ch for ch in self.pgf.tree['children'][pgf_index]['children']]
            for ch_index, chnl in enumerate(chnls):
                ep_start = sig.t_start
                for se_epoch_index, se_epoch in enumerate(chnl['children']):
                    se_rec = se_epoch['contents']
                    se_duration = pq.Quantity(float(se_rec.seDuration),'s')
                    if not(int(se_rec.seVoltageSource)):
                        se_voltage = pq.Quantity(float(se_rec.seVoltage),'V')
                    else:
                        se_voltage = pq.Quantity(float(chnl['contents'].chHolding),'V')
                    epoch = neo.Epoch(ep_start,se_duration,'protocol_epoch',value=se_voltage,channel_index=ch_index)
                    fully_annototate(chnl,epoch)
                    epoch.annotations['sig_index'] = sig_index
                    ep_start = ep_start + se_duration
                    seg.epochs.append(epoch)
        return seg

    def read_analogsignal(self,
                        lazy = False,
                        cascade = True,
                        group = 0,
                        series = 0,
                        sweep = 0,
                        channel = 0):
        tree = getbyroute(self.pul.tree,[0,group,series,sweep,channel])
        f = open(self.filename)
        sig = gettrace(tree['contents'],f)
        f.close()
        annotations = tree['contents'].__dict__.keys()
        annotations.remove('readlist')
        for a in annotations:
            d = {a:str(tree['contents'].__dict__[a])}
            sig.annotate(**d)
        pgf_index = series_count(self.pul,[0,group,series])
        sig.annotate(pgf_index = pgf_index)
        return sig

def fully_annototate(heka_tree,neo_object):
    annotations = heka_tree['contents'].__dict__.keys()
    annotations.remove('readlist')
    for a in annotations:
        d = {a:str(heka_tree['contents'].__dict__[a])}
        neo_object.annotate(**d)

def getleafs(tree_obj,f):
    if isinstance(tree_obj['contents'], TraceRecord):
        tr = [gettrace(tree_obj['contents'],f)]
        return copy.copy(tr)
    else:
        leaflist = list()
        for child in tree_obj['children']:
            [leaflist.append(leaf) for leaf in getleafs(child,f)]
        return leaflist

def protocol_signal(seg,sig_index,channel_index):
    """return a protocol signal for the apropriate trace"""
    criteria = lambda ep: ep.annotations['sig_index'] == sig_index and ep.annotations['channel_index'] == channel_index
    epochs = filter(criteria,seg.epochs)
    #print [ep.annotations['value'] for ep in epochs]
    total_duration = seg.analogsignals[sig_index].times[-1]
    #reduce(lambda x,y: x + y.time,epochs,pq.Quantity(0,'s'))
    sampling_period = seg.analogsignals[sig_index].sampling_period
    units = seg.epochs[0].annotations['value'].units
    prot_array = np.array(seg.analogsignals[0])#np.zeros(int(total_duration/sampling_period))
    prot_array[:] = epochs[-1].annotations['value']
    left_index = 0
    for epoch in epochs:
        right_index = int(epoch.duration/sampling_period) + left_index
        #print right_index,left_index
        prot_array[left_index:right_index] = epoch.annotations['value']
        left_index = right_index
    return AnalogSignal(prot_array,sampling_period = sampling_period,units = units)
                 
def gettrace(trec,f):
    import numpy as np
    format_type_lenghts = [2,4,4,8]
    format_type = [np.int16,np.int32,np.float32,np.float64]
    pointsize = format_type_lenghts[int(trec.trDataFormat)]
    dtype = format_type[int(trec.trDataFormat)]
    f.seek(int(trec.trData))
    byte_string = f.read(int(trec.trDataPoints)*pointsize)
    import numpy as np
    ydata = np.fromstring(byte_string,dtype = dtype)
    tunit = pq.Quantity(1,str(trec.trXUnit))
    yunit = pq.Quantity(1,str(trec.trYUnit))
    sig = AnalogSignal(ydata*float(trec.trDataScaler)*yunit,
            sampling_period=float(trec.trXInterval)*tunit,
            units = trec.trYUnit[0])
    annotations = trec.__dict__.keys()
    annotations.remove('readlist')
    for a in annotations:
        d = {a:str(trec.__dict__[a])}
        sig.annotate(**d)
    return sig
