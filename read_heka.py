__author__ = 'psilentp'

""" script to parse a packed HEKA data file"""

from struct import *
from collections import namedtuple
from heka_structs import *
import ephys as ep
import json

class TreeFile(object):
    """Abstract class to read and parse a Tree-format type subfile.
    It is constructed by passing in the packed HEKA file and the bundle
    item object of interest---usually a member of the BundleHeader class"""
    def __init__(self,b_file,b_item): #bundle file, bundle item objects
        self.b_item = b_item
        b_file.seek(int(b_item.oStart))
        mnumber = unpack('4s',b_file.read(4))
        assert(mnumber[0] in "eerT") #  assert the magic number
        self.num_levels = unpack('L',b_file.read(4))[0] #get the number of lev
        self.l_sizes = list() #list to hold the number of levels
        #(self.num_levels):
        for l in range(self.num_levels):
            self.l_sizes.append(unpack('L',b_file.read(4))[0])


    def load_level(self,level,b_file):
        #pack_str = str(self.pack_strs[level]) + 'L'
        #unpack the data into contents - change this so that contents is a class
        #that knows how to read in it's own data something like
        #contents = self.record_classes[level](self.l_sizes[level])
        #pass the level size to the constructor to ensure that the level size is
        #as expected, then call self.contents.read(b_file)
        #print level
                #print len(self.record_classes)
        #print len(self.l_sizes)

        contents = self.record_classes[level](self.l_sizes[level])
        contents.load(b_file)
        #print self.l_sizes[level]
        #print contents.read_size()-4
        #print contents

        #print contents.__dict__
        #print 'reading %s bytes'%self.l_sizes[level]
        #bytestr = repr(b_file.read(self.l_sizes[level]))
        #print bytestr
        children = list()
        node = {'contents':contents,'children':children}
        #num_children = unpack('L',b_file.read(4))[0]
        #print 'num_chil: ' + str(num_children)
        #for child in range(num_children):
        #    self.load_level(level+1,b_file)
        for child in range(contents.num_children):
            node['children'].append(self.load_level(level+1,b_file))
        return node

    def disp_lev(self,node,depth,rstr):
        rstr += depth*'\t'
        for line in str(node['contents']).split('\n'):
            rstr += line + ' '
        for child in node['children']:
            rstr += self.disp_lev(child,depth+1,'\n')
        return rstr

    def __str__(self):
        return self.disp_lev(self.tree,0,' ')


class PGFFile(TreeFile):
    """subclass of TreeFile, used to parse the PGF file. Just gets a file
    object and the BundleHeader for that file"""
    def __init__(self,b_file,b_item):
        TreeFile.__init__(self,b_file,b_item)
        RootRecSize             = 584
        StimulationRecSize      = 248
        ChannelRecSize          = 400
        StimSegmentRecSize      = 80
        #print self.l_sizes
        assert([RootRecSize,StimulationRecSize,
                ChannelRecSize,StimSegmentRecSize] == self.l_sizes)
        self.record_classes = (RootRecord_PGF,
                               StimulationRecord,
                               ChannelRecord,
                               StimSegmentRecord)
        self.tree = self.load_level(0,b_file)

class PULFile(TreeFile):
    """subclass of TreeFile, used to parse the PGF file. Just gets a file
    object and the BundleHeader for that file"""
    def __init__(self,b_file,b_item):
        TreeFile.__init__(self,b_file,b_item)
        RootRecSize          = 544
        GroupRecSize         = 128
        SeriesRecSize        = 1408
        SweepRecSize         = 160
        TraceRecSize         = 408
        assert([RootRecSize,GroupRecSize,SeriesRecSize,
                SweepRecSize,TraceRecSize] == self.l_sizes)
        self.record_classes = (RootRecord_PUL,
                               GroupRecord,
                               SeriesRecord,
                               SweepRecord,
                               TraceRecord)
        self.tree = self.load_level(0,b_file)

def make_dict(tree,pth):
    #print(pth)
    nodedict = dict()
    #add the children of a record
    for i,child in enumerate(tree['children']):
        p = copy.copy(pth)
        p.append(i)
        nodedict.update({"%03d"%(i) + '_' + child['contents'].dictLabel().rstrip() : make_dict(child,p)})
        #add the contents of the record
    """
    if isinstance(tree['contents'], heka_structs.TraceRecord):
        nodedict.update({'path':str(pth)})
        for key in tree['contents'].repdict().keys():
            if not key in ['readlist']:
                nodedict.update({key:tree['contents'].__dict__[key]})
    """
    return nodedict

class HEKAEncoder(json.JSONEncoder):
    def default(self,obj):
        if isinstance(obj,HEKA_type):
            return str(obj.data)
        return json.JSONEncoder.default(self,obj)

class TraceEncoder(json.JSONEncoder):
    import ephys
    def default(self,obj):
        if isinstance(obj,ep.Trace):
            #print 'here'
            maxsize = 2000
            downrate = 1
            if len(obj.y) > maxsize:
                downrate = len(obj.y)/(maxsize/2)
                #print downrate
            dec = obj[::downrate]
            dec.set_yunits('pA')
            dec.dt.units = 'ms'
            y = list(dec.y.__array__())
            x = list(dec.get_t_array(timebase = 'global').__array__())
            rlist = list()
            for xi,yi in zip(x,y):
                rlist.append({'x':xi,'y':yi})
            return rlist
        return json.JSONEncoder.default(self,obj)

def getbyroute(tree_obj,route):
    if len(route) == 1:
        return tree_obj['contents']
    for turn in route[1:]:
        tree_obj = tree_obj['children']
        tree_obj = tree_obj[turn]
    return tree_obj

def series_count(pul,route):
    count = 0
    for i in range(route[1]):
        tree = pul.tree['children'][i]
        count += int(tree['contents'].num_children)
    count += route[2]
    return count

def gettrace(trec,f):
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
    return ep.Trace(y = ydata*float(trec.trDataScaler),
        dt=float(trec.trXInterval),
        yunits = trec.trYUnit[0],
        xunits = trec.trXUnit[0])

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

def plotleafs(tree_obj,f):
    if isinstance(tree_obj['contents'], TraceRecord):
        gettrace(tree_obj['contents'],f).plot(color = 'k')
    else:
        for child in tree_obj['children']:
            plotleafs(child,f)

