import csv
f = open('db_entries_Sheet.csv')
reader = csv.reader(f)
data = [row for row in reader]
newdict = dict()
for cellrow in data[1:74]:
    celldict = dict()
    filename = cellrow[6]
    celldict.update({'cennum':int(cellrow[5][3:])+1224})
    celldict.update({'hekaID':int(cellrow[5][3:])})
    if not(cellrow[7] == ''):
        celldict.update({'on_cell_captrans':{'files':filename,'path':[int(cellrow[7])-1,0]}})
    if not(cellrow[8] == ''):
        celldict.update({'on_cell_ifam':{'files':filename,'path':[int(cellrow[8])-1]}})
    if not(cellrow[9] == ''):
        celldict.update({'in_cell_captrans':{'files':filename,'path':[int(cellrow[9])-1,0]}})
    if not(cellrow[10] == ''):
        celldict.update({'in_cell_ifam':{'files':filename,'path':[int(cellrow[10])-1]}})
    if not(cellrow[11] == ''):
        celldict.update({'long_ifam':{'files':filename,'path':[int(cellrow[11])-1]}})
    if not(cellrow[12] == ''):
        rpot_trials = list()
        for trial in cellrow[12].split(','):
            try:
                rpot_trials.append({'files':filename,'path':[int(trial)-1]})
            except ValueError:
                    pass
        celldict.update({'rev_pot':{'trials':rpot_trials}})
    cennum =int(cellrow[5][3:])+1224
    newdict.update({str(cennum):celldict})

f.close()

import plistlib
f = open('_tfile.plist')
olddict = plistlib.readPlist(f)
for key in newdict.keys():
    olddict.update({key:newdict[key]})
f.close()

from process_optical import write_sorted_plist

f = open('tfile.plist','w')
write_sorted_plist(olddict,f)
f.close()