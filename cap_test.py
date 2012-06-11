import heka_io
from heka_io import HekaIO
import cenfunctions as cnf
def get_captrace(group = 0):
    #h_reader = HekaIO('/Volumes/Storage/psilentp/Desktop/Ephys/Analysis/test_data/CEN111/THL_2011-07-09_15-02-54_000.dat')
    #h_reader = HekaIO('/Volumes/Storage/psilentp/Desktop/Ephys/Analysis/test_data/CEN189/THL_2012-05-05_03-47-10_000.dat')
    h_reader = HekaIO('/Users/psilentp/Documents/Projects/Analysis/test_data/CEN189/THL_2012-05-05_03-47-10_000.dat')
    blo = h_reader.read_block(group = group)
    seg = blo.segments[0]
    tmsirs = [cnf.TimeSeries(sig) for sig in seg.analogsignals]
    sum = reduce(lambda x,y:x+y,tmsirs)
    average = sum/len(tmsirs)
    protocol = cnf.TimeSeries(heka_io.protocol_signal(seg,0,0))
    return (average,protocol)


on_cell,stim = get_captrace(group = 0)
in_cell,dummy = get_captrace(group = 2)
subt = in_cell - on_cell

on_cell.set_tunits('s')
in_cell.set_tunits('s')
stim.set_tunits('s')
subt.set_tunits('s')
crecord = cnf.CapRecord(on_cell,in_cell,subt,stim)
crecord.plot()
import pylab as plb
plb.show()

#ax = plb.subplot(2,1,1)
#protocol.plot(marker = 'o')
#plb.subplot(2,1,2,sharex = ax)
#average.plot(marker = 'o')