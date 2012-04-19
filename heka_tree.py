__author__ = 'psilentp'

from heka_structs import *
from heka_io import *
from read_heka import *
from traits.api import HasTraits,Instance
from traitsui.api import View, Item,ValueEditor

class ViewablePGF(HasTraits):
    def __init__(self,pgf):
        self.pgf = pgf

default_view = View(Item(name='pgf',editor = ValueEditor()),height=100,width=100,resizable = True)

filename = './test_data/CEN111/THL_2011-07-09_15-02-54_000.dat'
ioreader = HekaIO(filename)
viewer = ViewablePGF(ioreader.pgf.tree)
viewer.configure_traits(view=default_view)
