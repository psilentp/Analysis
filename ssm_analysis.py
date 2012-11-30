"""analysis of ssm data"""

#CEN184/THL_2012-03-21_18-44-42_000.dat group 1 usable.
#CEN183/THL_2012-03-19_19-43-43_000.dat has small ~8mV jumps.
#CEN172/THL_2012-03-09_15-31-14_000.dat groups 7&8 have good jumpers.
#CEN17
"""
import heka_io as hio
import pylab as plb
file_prfx = "/Volumes/BigGuy/CELLS/"
file_pofix = "CEN172/THL_2012-03-09_15-31-14_000.dat"
file_name = file_prfx + file_pofix
reader = hio.HekaIO(file_name)
for i in range(8)[6:]:
	blo = reader.read_block(group = i)
	seg = blo.segments[0]
	[plb.plot(sig.times,sig) for sig in seg.analogsignals[::2]]
plb.show()
"""


plistfile = open('tfile.plist','r')
import plistlib
d = plistlib.readPlist(plistfile)
for cen in d.keys():
	try:
		if d[cen]['sID'].upper() in ['SRA-6+','ASH(L)','ASH(R)','ASH']:
			print cen + " " + d[cen]['sID']
	except (KeyError,TypeError):
		pass
		#print "CEN %s has no sID"%(cen)

#import db_access as dba
#db = dba.get_db()

file_prfx = "/Volumes/WD_PASSPORT/CENs/CEN"
file_pofix = "1162/2010_01_12_0166.abf"
file_name = file_prfx + file_pofix

import abfloader as abl
data = abl.load_data(file_name)
print data['signals'].keys()
