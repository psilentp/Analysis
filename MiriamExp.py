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


fprefix = '/Volumes/Data/CENs/CEN'

jefferson = [[''
            ]
washington = [['184','THL_2012-03-21_18-40-42_000.dat','THL_2012-03-21_18-44-42_000.dat'],
              ['181','THL_2012-03-15_17-49-05_000.dat','THL_2012-03-15_17-54-17_000.dat'],
              ['180','THL_2012-03-15_15-05-04_000.dat','THL_2012-03-15_15-08-40_000.dat'],
              ['174','THL_2012-03-12_18-34-23_000.dat','THL_2012-03-12_18-38-41_000.dat'],
              ['173','THL_2012-03-12_17-28-39_000.dat','THL_2012-03-12_17-32-19_000.dat'],
              ['172','THL_2012-03-09_15-26-56_000.dat','THL_2012-03-09_15-31-14_000.dat']]


washington = [[fprefix+c+'/'+x,fprefix+c+'/'+y] for c,x,y in washington]

jefferson = []
for i,block_num in enumerate([1,2,5,1,2,4]):
    filename = washington[i][1]#'./test_data/CEN184/THL_2012-03-21_18-44-42_000.dat'
    ioreader = HekaIO(filename)
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
    x = [0.01, 0.03, 0.1, 0.3, 1, 3]
    x = x[:len(y)]
    print len(y)
    #print y
    plb.plot(x,y,'-o')
    plb.gca().set_xscale('log')
    print input_output
plb.show()
