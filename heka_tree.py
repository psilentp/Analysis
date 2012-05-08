__author__ = 'psilentp'

from heka_structs import *
from heka_io import *
from read_heka import *
from traits.api import HasTraits,Instance,Int,Str,List,File

from traitsui.api import View, Item,ValueEditor,TreeEditor,TreeNode,Group,Handler
no_view = View()


class ViewablePGF(HasTraits):
    def __init__(self,pgf):
        self.pgf = pgf
"""
default_view = View(Item(name='pgf',editor = ValueEditor()),height=200,width=200,resizable = True)
start_globs = set(globals().keys())
print start_globs
filename = './test_data/CEN111/THL_2011-07-09_15-02-54_000.dat'
mylist = [1,2,3,4,5]
ioreader = HekaIO(filename)
#viewer = ViewablePGF(ioreader.pgf.tree)
new_globs = set(globals().keys())
root = globals()
[root.pop(key) for key in start_globs.difference(new_globs)]
viewer = ViewablePGF(root)
viewer.configure_traits(view=default_view)
change
"""

import os
cell_directory = './test_data/'
dir_contents = os.listdir(cell_directory)
cen_strs = [x for x in dir_contents if x[:3] == 'CEN']

class Trial(HasTraits):
    name = Str('<unknown>')

class CeNeuron(HasTraits):
    number = Int()
    cellid = Str('<unknown>')
    files = List(File)
    trials = List(Trial)
    snum = Str()
    def __init__(self,number,cellid,trials):
        self.number = number
        self.cellid = cellid
        self.trials = trials
        self.snum = str(number)
        self.files = [x for x in os.listdir(cell_directory + 'CEN' + str(number) + '/')]

class CellBrowser(HasTraits):
    name = Str('C elegans neurons')
    cens = List(CeNeuron)

class Root(HasTraits):
    browser = Instance(CellBrowser)

trial1 = Trial(name = 'trial1')
trial2 = Trial(name = 'trial2')
trial3 = Trial(name = 'trial3')
trial4 = Trial(name = 'trial4')
trial5 = Trial(name = 'trial5')

cells = [CeNeuron(number = int(x[3:]), cellid = 'AVB', trials = []) for x in cen_strs]
"""
cell1 = CeNeuron(number = 1,
                cellid ='AVA',
                trials = [trial1,trial2]
                )

cell2 = CeNeuron(number = 2,
                cellid = 'AVB',
                trials = [trial3,trial4,trial5]
                )
"""
tree_root = Root(browser = CellBrowser(cens=cells))

tree_editor = TreeEditor(
    nodes = [
        TreeNode( node_for = [CellBrowser],
            auto_open = True,
            children = 'cens',
            label = 'name',
            view = no_view ),
        TreeNode( node_for = [CeNeuron],
            auto_open = True,
            children = 'trials',
            label = 'snum',
            view = View(Group(
                        Item('cellid'),
                        Item('files'), orientation = 'vertical',show_left =True)) ),
        TreeNode( node_for = [Trial],
            auto_open = False,
            children = '',
            label = 'name',
            view = no_view)
    ]
)

view = View(
    Group(
        Item(
            name = 'browser',
            id = 'cen',
            editor = tree_editor,
            resizable = True ),
        orientation = 'vertical',
        show_labels = True,
        show_left = True, ),
    title = 'Cells',
    id ='traitsui.tests.tree_editor_test',
    dock = 'horizontal',
    buttons = [ 'Undo', 'OK', 'Cancel' ],
    resizable = True,
    width = .3,
    height = .3 )

if __name__ == '__main__':
    tree_root.configure_traits( view = view )