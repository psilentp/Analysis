import quantities as pq
import pylab as plb
import numpy as np

def get_low_filter(signal,Hz):
    #return a low-pass-filterd version of the time series
    from scipy.signal import filtfilt
    filter_order = 3
    sampfreq = signal.sampling_rate
    filtfreq = pq.Quantity(Hz,'Hz')
    passband_peram = float((filtfreq/sampfreq).simplified)
    from scipy.signal import butter
    [b,a]=butter(filter_order,passband_peram)
    filterd = pq.Quantity(filtfilt(b,a,signal),signal.units)
    return filterd

def sub_baseline(signal,start,stop,**kwargs):
    #baseline subtract the Trace, start and stop define the interval,
    #in object time units
    tunits = kwargs.pop('tunits',signal.sampling_period.units)
    timebase = kwargs.pop('timebase','local')
    #print 'here'
    #print timebase
    #print timebase == 'global'
    if timebase == 'global':
        #print 'global'
        start -= signal.t_start.magnitude
        stop -= signal.t_stop.magnitude
    baseline = pq.Quantity(np.mean(ts(signal,start,stop,tunits = tunits)),signal.units)
    return (signal - baseline)

def ts(signal,*key,**kwargs):
    #time slice - cut up a wave given start and end times tunits corresponds
    #to the units for the passed start and stop params
    tunits = kwargs.pop('tunits',signal.sampling_period.units)
    timebase = kwargs.pop('timebase','local')
    dtunits = signal.sampling_period.units
    dt = signal.sampling_period
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
    ob = signal.__getitem__(selection)
    ob.t_start = pq.Quantity(0,tunits)
    if(timebase == 'global'):
        ob.t_start+=signal.t_start.rescale(tunits) + pq.Quantity(key[0],tunits).rescale(dtunits)
    #return np.array(ar)
    return ob



groups = [4,8,15]
#stimtimes = [0,30,30,60,60,90]
stimtimes = [0,1,1,2,2,2]
stimvals =  [50,50,0,0,50,50]
stimax = plb.subplot(2,1,1)
plb.plot(stimtimes,stimvals)
plb.subplot(2,1,2,sharex = stimax)
import heka_io
filename = './test_data/CEN189/THL_2012-05-05_03-47-10_000.dat'
ioreader = heka_io.HekaIO(filename)
for group in groups:
    block = ioreader.read_block(group=group)
    seg = block.segments[0]
    signal = seg.analogsignals[0]
    #filtered = signal
    #choped = signal
    choped = ts(signal,29,31,timebase = 'local')
    #choped = sub_baseline(choped,0,4)
    filtered = get_low_filter(choped,100)
    #plb.plot(signal.times[::10],signal[::10],color = 'k' )
    plb.plot(choped.times[::10],filtered[::10],color = 'k' )

plb.show()