#
#  process_optical.py
#  
#
#  Created by Theodore Lindsay on 12/9/09.
#  Copyright (c) 2009 University of Oregon. All rights reserved.
#

import matplotlib
#matplotlib.use('pdf')
import cenfunctions as cnf
import abfloader as abl
import numpy as np

import cendb as cdb

import psilentplib as psl
import pylab as plb
import os
#import gc

#import the datamap dicts to map data.
#from datafiles import CEN_data_map
#from datafiles import c_list

def write_sorted_plist(d,f):
    from collections import OrderedDict
    keylist = d.keys()
    def custom_comp(k):
        try:
            return int(k)
        except ValueError:
            return -1
    keylist.sort(key = custom_comp)
    d2 = OrderedDict()
    for key in keylist:
        d2.update({key:d[key]})
    print d2.keys()
    writer = SortedPlistWriter(f)
    writer.writeDict(d2)
    
import plistlib
class SortedPlistWriter(plistlib.PlistWriter):
    def writeDict(self,d):
        self.beginElement("dict")
        items = d.items()
        for key, value in items:
            assert isinstance(key, (str,unicode)),"keys must be strings"
            self.simpleElement("key",key)
            self.writeValue(value)
        self.endElement("dict")

def make_cell(cennumber,celldb,cdm):
    #construct a cell and add it to the celldb
    celldb.update({cennumber:cdb.Cell(cennum = cennumber,cdm = cdm)})
    #transaction.commit()
    
def make_plots(CEN,cell):
    dirstr = 'cellfigs/CEN' + str(CEN)
    try:
        os.mkdir(dirstr)
    except OSError:
        print "removing directory " + dirstr
        import shutil
        shutil.rmtree(dirstr)
        os.mkdir(dirstr)
        
    exp_dict = {
    'fam':" subtracted i-fam",
    'captrans':" captrans",
    'iv':" IV_data",
    'presynaptic_ic':" Presynaptic Optical DRC (IC)",
    'presynaptic_vc':" Presynaptic Optical DRC (VC)",
    'postsynaptic_ic':" Presynaptic Synaptic Potential (IC)",
    'synaptic_stim':" Postsynaptic Synaptic Currents (VC)",
    'long_ifam':" Long i-fam",
    'long_iv':" IV curve taken from long-ifam",
    's_slope':" tau of early currents",
    'best_n':" IV curve with steepest negative slope",
    'rev_fam':" Reversal Potential",
    'amp_curve':" Isyn vs Vm",
    'blue_pulse':" long pulse",
    'ic_bluepulse_long': " ic long blue pulse",
    'vc_bluepulse_long': " vc long blue pulse",
    'synaptic_stim_low': " Postsynaptic Currents (VC) low power",
    'amp_curve': " evoked I vs V curve",
    }
    
    for attr in exp_dict.keys():
        try:
            p = cell.__getattribute__(attr)
        except AttributeError:
            print ("CEN%s is missing %s " % (CEN,attr))
            continue
        p.plot()
        plb.suptitle("CEN:" + str(CEN) + exp_dict[attr],fontsize = 20)
        plb.savefig("%s%s%s%s%s" % (dirstr,'/',str(CEN),'_',attr))
        plb.close()
        
def just_plot():
    celldb = dba.get_db()
    explist = CEN_data_map.keys()
    for CEN in explist:
        make_plots(CEN,celldb[CEN])
        

def multiload(CEN):
    from ZODB import FileStorage, DB
    
    storage = FileStorage.FileStorage('optical_db.fs')
    db = DB(storage)
    conn = db.open()
    
    dbroot = conn.root()
    celldb = dbroot['celldb']
    
    if not(CEN in celldb):
            print "loading CEN" + str(CEN)
            try:
                make_cell(CEN,celldb,CEN_data_map[CEN])
                make_plots(CEN,celldb[CEN])
                #gc.collect()
            except MemoryError:
                dba.close_db()
                print("Memory Error!")
    import transaction
    transaction.commit()

def add_new_cells(explist):
    """update the database with some new cells,
    passed in as a list of cennums."""
    
    #Load the experimental conditions property list
    import plistlib
    plistfile = open('tfile.plist','r')
    d = plistlib.readPlist(plistfile)
    #create a list of cells
    c_list = list()
    for i in d:
        if 'cennum' in d[i]:
            c_list.append(int(d[i]['cennum']))
    #create the cell data map dictionary (keyed by integers)
    CEN_data_map = dict()
    [CEN_data_map.update({i:d[str(i)]}) for i in c_list]
    plistfile.close()
    
    for CEN in explist:
        import persistent
        from ZODB import FileStorage, DB
        from BTrees.IOBTree import IOBTree
        import transaction
        storage = FileStorage.FileStorage('optical_db.fs')
        #create a db object using the file
        db = DB(storage,large_record_size=150000000)
        #open a connection to the db object
        connection = db.open()
        #get the root of the connection
        dbroot = connection.root()
        #get a reference to celldb
        try:
            celldb = dbroot['celldb']
        except KeyError:
            dbroot['celldb'] = IOBTree()
            celldb = dbroot['celldb']
        #celldb = dba.get_db()
        #load and store the data for a cell
        if not(CEN in celldb):
            print "loading CEN" + str(CEN)
            try:
                make_cell(CEN,celldb,CEN_data_map[CEN])
                make_plots(CEN,celldb[CEN])
                #gc.collect()
            except MemoryError:
                dba.close_db()
                print("Memory Error!")
                break
        import transaction
        transaction.commit()
        storage.close()
    import transaction
    transaction.commit()
    storage.close()

def main():
    #Load the experimental conditions property list
    import plistlib
    plistfile = open('tfile.plist','r')
    d = plistlib.readPlist(plistfile)
    #create a list of cells
    c_list = list()
    for i in d:
        if 'cennum' in d[i]:
            c_list.append(int(d[i]['cennum']))
    #create the cell data map dictionary (keyed by integers)
    CEN_data_map = dict()
    [CEN_data_map.update({i:d[str(i)]}) for i in c_list]
    
    #use explist to subindex the cells that will be processed in 
    #current analysis
    #explist = [1212,1213,1214,1215,1216,1217,1218,1219,1220,1221,1222,1223,1224] #explist for ga_blue_pulse
    """
    explist =   [1175,1176,1177,1178,1179,1180,1181,1182,1183,1184,
                1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,
                1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,
                1205,1206,1207,1208,1209,1210,1211]
                
    explist += [759,761,764,765,767,769,770,790,791,792]
    """
    #explist = [759]
    #explist = [1188]
    #explist = c_list
    #c_list.remove(1172)
    explist = [1224]
    r = 0
    for CEN in explist:
        import persistent
        from ZODB import FileStorage, DB
        from BTrees.IOBTree import IOBTree
        import transaction
        storage = FileStorage.FileStorage('optical_db.fs')
        #create a db object using the file
        db = DB(storage)
        #open a connection to the db object
        connection = db.open()
        #get the root of the connection
        dbroot = connection.root()
        #get a reference to celldb
        celldb = dbroot['celldb']
        #celldb = dba.get_db()
        #load and store the data for a cell
        if not(CEN in celldb):
            print "loading CEN" + str(CEN)
            try:
                make_cell(CEN,celldb,CEN_data_map[CEN])
                make_plots(CEN,celldb[CEN])
                #gc.collect()
            except MemoryError:
                dba.close_db()
                print("Memory Error!")
                break
            if r > 5:
                break
            r += 1
        import transaction
        transaction.commit()
        storage.close()
    import transaction
    transaction.commit()
    storage.close()




if __name__ == '__main__':

    main()
