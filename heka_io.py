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
            times.append(float(ep.time))
            vms.append(float(ep.label))
            times.append(float(ep.time)+float(ep.duration))
            vms.append(ep.label)
        return (np.array(times),np.array(vms))
    import pylab as plb
    #filename = './test_data/CEN184/THL_2012-03-21_18-40-42_000.dat'
    #filename = './test_data/CEN184/THL_2012-03-21_18-44-42_000.dat'
    filename = './test_data/CEN111/THL_2011-07-09_15-02-54_000.dat'
    ioreader = HekaIO(filename)
    blo = ioreader.read_block(group = 0)
    for seg in blo.segments:
        ax1 = plb.subplot(2,1,1)
        x,y = get_stimtrace(seg.epochs)
        plb.plot(x,y,'o-')
        plb.subplot(2,1,2,sharex = ax1)
        for a_sig in seg.analogsignals:
            for key in ioreader.pgf.tree['children'][a_sig.annotations['pgf_index']]['contents'].__dict__:
                print key + ':' + str(ioreader.pgf.tree['children'][a_sig.annotations['pgf_index']]['contents'].__dict__[key])
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
        #print 'here'
        BaseIO.__init__(self)
        self.filename = filename
        f = open(filename)
        #create a bundle header object with the file
        head = BundleHeader(f)
        #load the file
        head.load(f)
        #print head
        #get the .pgf and .pul items in the file
        for bi in head.oBundleItems:
            if str(bi.oExtension)[0:4] == '.pgf':
                self.pgf = PGFFile(f,bi)
            if str(bi.oExtension)[0:4] == '.pul':
                #print 'here'
                self.pul = PULFile(f,bi)
        f.close()

    def read_block(self,
                   lazy = False,
                   cascade = True,
                   group = 0):
        blo = Block(name = 'test')
        if cascade:
            tree = getbyroute(self.pul.tree,[0,group])
            #print tree['contents'].__dict__
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
        #ep_start = pq.Quantity(0,'s')
        for sig in seg.analogsignals:
            pgf_index = sig.annotations['pgf_index']
            ep_start = sig.t_start
            st_rec = self.pgf.tree['children'][pgf_index]['contents']
            for se_epoch in self.pgf.tree['children'][pgf_index]['children'][0]['children']:
                se_rec = se_epoch['contents']
                se_duration = pq.Quantity(float(se_rec.seDuration),'s')
                se_voltage = se_rec.seVoltage
                if not se_rec.seVoltage.data:
                    se_voltage = self.pgf.tree['children'][pgf_index]['children'][0]['contents'].chHolding
                epoch = neo.Epoch(ep_start,se_duration,se_voltage)
                ep_start = ep_start + se_duration
                #print se_duration
                #print ep_start
                seg.epochs.append(epoch)
            #final_duration = pq.Quantity(float(st_rec.stSweepInterval),'s') - ep_start + sig.t_start
            #print ep_start
            #print final_duration
            #epoch = neo.Epoch(ep_start,final_duration,-0.04)
            #seg.epochs.append(epoch)
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

def getleafs(tree_obj,f):
    if isinstance(tree_obj['contents'], TraceRecord):
        tr = [gettrace(tree_obj['contents'],f)]
        #print type(tr)
        return copy.copy(tr)
    else:
        leaflist = list()
        for child in tree_obj['children']:
            [leaflist.append(leaf) for leaf in getleafs(child,f)]
        return leaflist

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
