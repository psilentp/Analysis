#
#  db_access.py
#  
#
#  Created by Theodore Lindsay on 12/21/09.
#  Copyright (c) 2009 University of Oregon. All rights reserved.
#

import persistent
from ZODB import FileStorage, DB
from BTrees.IOBTree import IOBTree
import transaction
    
#storage = FileStorage.FileStorage('/Users/psilentp/Desktop/Analysis/optical_db.fs')

storage = FileStorage.FileStorage('/Volumes/Storage/psilentp/Desktop/Ephys/Analysis/optical_db.fs')

def get_db():
    #open up the file storage
    #create a db object using the file
    db = DB(storage)
    #open a connection to the db object
    connection = db.open()
    #get the root of the connection
    dbroot = connection.root()
    #get a reference to celldb
    celldb = dbroot['celldb']
    return celldb
    
def make_db():
    #create a db object using the file
    db = DB(storage)
    #open a connection to the db object
    connection = db.open()
    #get the root of the connection
    dbroot = connection.root()
    #create a 'celldb' in root 
    dbroot['celldb'] = IOBTree()
    #get a reference to celldb
    celldb = dbroot['celldb']
    return celldb
    
def close_db():
    transaction.commit()
    storage.close()
