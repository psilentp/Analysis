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

errors = [['182','THL_2012-03-19_17-40-54_000.dat','THL_2012-03-19_17-49-46_000.dat'],
    ['179','THL_2012-03-14_20-02-26_000.dat','THL_2012-03-14_20-06-51_000.dat'],
    ['177','THL_2012-03-14_16-35-53_000.dat','THL_2012-03-14_16-39-22_000.dat'],
    ['176','THL_2012-03-13_19-37-19_000.dat','THL_2012-03-13_19-40-47_000.dat'],
    ['175','THL_2012-03-13_15-34-06_000.dat','THL_2012-03-13_15-38-02_000.dat']]

jefferson = [['182','THL_2012-03-19_17-40-54_000.dat','THL_2012-03-19_17-49-46_000.dat'],
             ['179','THL_2012-03-14_20-02-26_000.dat','THL_2012-03-14_20-06-51_000.dat'],
             ['176','THL_2012-03-13_19-37-19_000.dat','THL_2012-03-13_19-40-47_000.dat'],
             ['183','THL_2012-03-19_19-39-39_000.dat','THL_2012-03-19_19-43-43_000.dat'],
             ['178','THL_2012-03-14_17-31-43_000.dat','THL_2012-03-14_17-36-20_000.dat'],
             ]

washington = [['184','THL_2012-03-21_18-40-42_000.dat','THL_2012-03-21_18-44-42_000.dat'],
              ['181','THL_2012-03-15_17-49-05_000.dat','THL_2012-03-15_17-54-17_000.dat'],
              ['180','THL_2012-03-15_15-05-04_000.dat','THL_2012-03-15_15-08-40_000.dat'],
              ['174','THL_2012-03-12_18-34-23_000.dat','THL_2012-03-12_18-38-41_000.dat'],
              ['173','THL_2012-03-12_17-28-39_000.dat','THL_2012-03-12_17-32-19_000.dat'],
              ['172','THL_2012-03-09_15-26-56_000.dat','THL_2012-03-09_15-31-14_000.dat']]
[3,
 1,
 1]

washington = [[fprefix+c+'/'+x,fprefix+c+'/'+y] for c,x,y in washington]
jefferson = [[fprefix+c+'/'+x,fprefix+c+'/'+y] for c,x,y in jefferson]
errors = [[fprefix+c+'/'+x,fprefix+c+'/'+y] for c,x,y in errors]

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
    #print len(y)
    #print y
    plb.plot(x,y,'-',color = 'g',lw=3,alpha=0.5)
    #plb.gca().set_xscale('log')
    #print input_output

for i,block_num in enumerate([1,2,1,3,1]):
    filename = jefferson[i][1]#'./test_data/CEN184/THL_2012-03-21_18-44-42_000.dat'
    ioreader = HekaIO(filename)
    blo = ioreader.read_block(group = block_num)
    #stimulus power is stored in the second epoch in which annotations['chDacChannel'] == '2'
    input_output = list()
    for seg in blo.segments[:6]:
        sig = seg.analogsignals[0]
        eps = [x for x in seg.epochs if x.annotations['chDacChannel'] == '2']
        chunk = ts(sig,eps[1].time,eps[1].time+eps[1].duration)
        input_output.append({'in':eps[1].annotations['value'],'out':np.mean(chunk)})

    x = [io['in'] for io in input_output]
    y = [io['out'] for io in input_output]
    x = [0.01, 0.03, 0.1, 0.3, 1, 3]
    x = x[:len(y)]
    #print len(y)
    #print y
    plb.plot(x,y,'-',color = 'b',lw=3,alpha=0.5)
    plb.gca().set_xscale('log')

    #print input_output
#plb.legend()


def get_block(filename =  filename,block_num = 1):
    def get_stimtrace(epochs,channel):
        times = []
        vms = []
        for ep in epochs:
            if ep.annotations['chDacChannel'] == channel:
                times.append(float(ep.time))
                vms.append(float(ep.annotations['value']))
                times.append(float(ep.time)+float(ep.duration))
                vms.append(ep.annotations['value'])
        return {'x':np.array(times),'y':np.array(vms)}
        #filename = './test_data/CEN184/THL_2012-03-21_18-40-42_000.dat'
    #filename = './test_data/CEN184/THL_2012-03-21_18-44-42_000.dat'
    #filename = './test_data/CEN111/THL_2011-07-09_15-02-54_000.dat'
    ioreader = HekaIO(filename)
    blo = ioreader.read_block(group = block_num)

    protocol_list = list()
    sweep_list = list()
    max_channels = 0
    for seg in blo.segments:
        chnl_map = dict()
        for ep_index, ep in enumerate(seg.epochs):
            if ep.label == 'protocol_epoch':
                chnl_map.update({ep.annotations['chDacChannel']:
                                     get_stimtrace(seg.epochs,ep.annotations['chDacChannel'])})
            if len(chnl_map.keys()) > max_channels:
                max_channels = len(chnl_map.keys())
        protocol_list.append(chnl_map)
        for a_sig in seg.analogsignals:
            sweep_list.append({'x':np.array(a_sig.times),'y':np.array(a_sig)})
    return(sweep_list,protocol_list,max_channels)


# plot example trace
plb.figure()
filename = '/Volumes/Data/CENs/CEN176/THL_2012-03-13_19-40-47_000.dat'
data = get_block(filename=filename,block_num=1)
ax = plb.subplot(2,1,1)
for prot in data[1][:6]:
    print prot.keys()
    plb.plot(prot['2']['x'],prot['2']['y'],color ='k')
plb.subplot(2,1,2,sharex = ax)
for sig in data[0]:
    plb.plot(sig['x'],sig['y'], color = 'k')
plb.show()