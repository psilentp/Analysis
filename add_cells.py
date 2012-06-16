#cells = [1322,10000,1323,1324]

import csv
f = open('db_entriesSheet.csv')
reader = csv.reader(f)
cells = [int(hid)+1224 for hid in [line[5][3:] for line in reader][1:]]
cells = cells[:cells.index(1333)+1]


import subprocess
for cell in cells:
    subprocess.call(['./add_cell.py',str(cell)])