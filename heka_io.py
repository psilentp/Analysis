import neo
from __future__ import absolute_import
from neo.io.baseio import BaseIO
from neo.core import Block, Segment, AnalogSignal, EventArray
from neo.tools import create_many_to_one_relationship
import numpy as np
import quantites as pq

class ExampleIO(BaseIO):
    is_readable = True
    is_writable = False
    supported_objects = [Segment,AnalogSignal,EventArray]
    readable_objects = [Segment,AnalogSignal]
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

    def __init_(self,filename = None):
        BaseIO.__init__(self)
        self.filename = filename
        f = open('./test_data/CEN184/THL_2012-03-21_18-40-42_000.dat')
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

    def read_segment(self,
                    lazy = False,
                    cascade = True,
                    path):
        seg = Segment( name = 'test')
        if cascade:
            seg.analogsignals = getleafs(self.pul.tree)


    def read_analogsignal(self,
                        lazy = False,
                        cascade = True):
        pass

def pu(str):
    from quantities import ms,kHz,pA,mV
    if str == 'pA':
        return pA
    if str == 'mV':
        return mV
    if str == 'ms':
        return ms

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
    from quantities import ms,kHz,nA,mV
    import ephys as ep
    import numpy as np
    format_type_lenghts = [2,4,4,8]
    format_type = [np.int16,np.int32,np.float32,np.float64]
    pointsize = format_type_lenghts[int(trec.trDataFormat)]
    dtype = format_type[int(trec.trDataFormat)]
    f.seek(int(trec.trData))
    byte_string = f.read(int(trec.trDataPoints)*pointsize)
    import numpy as np
    ydata = np.fromstring(byte_string,dtype = dtype)
    return AnalogSignal(ydata*float(trec.trDataScaler)*pu(trec.trXunit[0]),
                        sampling_period=float(trec.trXInterval)*pu(trec.XUnit[0]),
                        units = trec.trXUnit[0])
