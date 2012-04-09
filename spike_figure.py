from neo.io import axonio
import numpy as np
import pylab as plb
import neo
import quantities as pq
""" code for Shawn's review paper"""

def fig_1():
    filename = './test_data/CEN914/2008_12_04_0011.abf'
    ioreader = axonio.AxonIO(filename)
    blo = ioreader.read()
    seg = blo.segments[0]
    times = seg.analogsignals[0].times
    sind = np.nonzero(times > 314.5)[0][0]
    eind = np.nonzero(times > 317.5)[0][0]
    signal_chunk = seg.analogsignals[0][sind:eind]
    stim_chunk = seg.analogsignals[2][sind:eind]

    pulse_start = 10400
    pulse_end = 21600

    stim_times =  [times[sind],
                   signal_chunk.times[pulse_start],
                   signal_chunk.times[pulse_start],
                   signal_chunk.times[pulse_end],
                   signal_chunk.times[pulse_end],
                   times[eind]]

    stim = [0,0,3,3,0,0]
    stim_times = [x-stim_times[0] for x in stim_times]

    stim_ax = plb.axes((0.1,0.8,0.8,0.1))
    plb.plot(stim_times,stim)

    sig_ax = plb.axes((0.1,0.1,0.8,0.6),sharex = stim_ax)
    plb.plot(signal_chunk.times-signal_chunk.times[0],signal_chunk)
    #plb.plot(signal_chunk)
    #plb.plot(stim_chunk.times-stim_chunk.times[0],stim_chunk)



######################
def fig_2():
    plb.figure()


    filenames = ['./test_data/CEN914/2008_12_04_0013.abf','./test_data/CEN914/2008_12_04_0020.abf']

    filename = filenames[0]
    ioreader = axonio.AxonIO(filename)

    stim_ax = plb.axes((0.1,0.8,0.8,0.1))
    prots = ioreader.read_protocol()
    print prots
    plb.plot(prots[0].analogsignals[0].times,prots[0].analogsignals[0])

    sig_ax = plb.axes((0.1,0.1,0.8,0.6),sharex = stim_ax)
    blo = ioreader.read()
    indxs = range(len(blo.segments))
    #indxs = [0,1,2]
    for seg in map(blo.segments.__getitem__,indxs):
        sig = seg.analogsignals[0]
        sig.t_start
        times = sig.times - sig.times[0]
        plb.plot(times,sig,'k')

    filename = filenames[1]
    ioreader = axonio.AxonIO(filename)
    blo = ioreader.read()
    indxs = range(len(blo.segments))
    #indxs = [1,2]
    for seg in map(blo.segments.__getitem__,indxs):
        sig = seg.analogsignals[0]
        sig.t_start
        times = sig.times - sig.times[0]
        plb.plot(times,sig,'k')


fig_1()
#plb.gca().set_xbound((0.52,0.65))
plb.show()