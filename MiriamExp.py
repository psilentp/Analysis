__author__ = 'psilentp'

from heka_io import *
import quantities as pq
import pylab as plb
import numpy as np

def ts(sig,*key,**kwargs):
#time slice - cut up a wave given start and end times tunits corresponds
#to the units for the passed start and stop params
    tunits = kwargs.pop('tunits',sig.sampling_period.units)
    timebase = kwargs.pop('timebase','local')
    dtunits = sig.sampling_period.units
    dt = sig.sampling_period
    ind = [None,None,None]
    if np.size(key) == 3:
        ind[0] = int(pq.Quantity(key[0],tunits).rescale(dtunits)/dt)
        ind[1] = int(pq.Quantity(key[1],tunits).rescale(dtunits)/dt)
        ind[2] = int(pq.Quantity(key[2],tunits).rescale(dtunits)/dt)
    elif np.size(key) == 2:
        ind[0] = int(pq.Quantity(key[0],tunits).rescale(dtunits)/dt)
        ind[1] = int(pq.Quantity(key[1],tunits).rescale(dtunits)/dt)
        ind[2] = 1
    elif np.size(key) == 1:
        ind[0] = int(pq.Quantity(key[0],tunits).rescale(dtunits)/dt)
        ind[1] = int(pq.Quantity(key[0],tunits).rescale(dtunits)/dt) + 1
        ind[2] = 1
    selection = slice(ind[0],ind[1],ind[2])
    ob = sig.__getitem__(selection)
    if(timebase == 'global'):
        ob.tstart +=sig.tstart.rescale(tunits) + pq.Quantity(key[0],tunits).rescale(dtunits)
        #return np.array(ar)
    return ob


filename = './test_data/CEN184/THL_2012-03-21_18-44-42_000.dat'
ioreader = HekaIO(filename)
block_num = 3
blo = ioreader.read_block(group = block_num)
#stimulus power is stored in the second epoch in which annotations['chDacChannel'] == '2'
input_output = list()
for seg in blo.segments:
    sig = seg.analogsignals[0]
    eps = [x for x in seg.epochs if x.annotations['chDacChannel'] == '2']
    chunk = ts(sig,eps[1].time,eps[1].time+eps[1].duration)
    input_output.append({'in':eps[1].annotations['value'],'out':np.mean(chunk)})

x = [io['in'] for io in input_output]
y = [io['out'] for io in input_output]
plb.plot(x,y,'-o')
print input_output
plb.show()