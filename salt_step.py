import quantities as pq
import pylab as plb

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

groups = [4,8,15]
stimtimes = [0,30,30,60,60,90]
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
    filtered = get_low_filter(signal,100)
    plb.plot(signal.times[::10],filtered[::10],color = 'k' )

plb.show()