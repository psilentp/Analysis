import pylab as plb
import quantities as pq
import numpy as np
from neo import io


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

dfiles = {897:'SSM_data/CEN897_p_potential.abf',
          828:'SSM_data/CEN828_p_potential.abf',
          799:'SSM_data/CEN799_p_potential.abf'}

tslices = {897:(572,740),828:(70,180),799:(550,900)}
signals = dict()
ax = plb.subplot(3,1,1)
for i,CEN in enumerate(dfiles.keys()):
    print CEN
    reader = io.AxonIO(filename = dfiles[CEN])
    cell = reader.read()
    seg = cell.segments[0]
    signal1 = ts(seg.analogsignals[0],tslices[CEN][0],tslices[CEN][1])
    if CEN==897:
        signal1 -= pq.Quantity(20,'mV')
    signals.update({CEN:signal1})
    plb.subplot(3,1,i+1,sharex = ax)
    plb.plot(signal1.times[::20],signal1[::20])
    plb.gca().set_ybound([0,-45])
plb.show()

    #signal2 = seg.analogsignals[1]
    #plb.subplot(2,1,1)
    #plb.plot(signal1.times[::20],signal1[::20])
    #plb.subplot(2,1,2)
    #plb.plot(signal2.times[::20],signal2[::20])
    #plb.show()