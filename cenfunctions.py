#
#  cenfunctions.py
#  
#
#  Created by Theodore Lindsay on 12/15/09.
#  Copyright (c) 2009 University of Oregon. All rights reserved.
#

import weakref

from pylab import *
import lsqfit as lsq
import quantities as quan
import linecache
import numpy as np
import psilentplib as psl
import tables as pyt
import copy as cp
import persistent
import transaction

from scipy import integrate
#from scipy import mean
    
class SubPanel(object):
    """Utility class for plotting, contains the bounds of the sub-pannel in 
    the range of 0-1 and a function to remap itself within some other 
    sub-pannel."""
    def __init__(self,bounds):
        self.bounds = bounds

    def rect_remap(self,re,bounds = None):
        """return the coordinates of the sub-pannel remamped within the bounds 
        of another global pannel."""
        if not bounds is None:
            self.bounds = bounds
            
        bound_l = self.bounds[0]
        bound_b = self.bounds[1]
        bound_w = self.bounds[2]
        bound_h = self.bounds[3]

        re_l = re[0]
        re_b = re[1]
        re_w = re[2]
        re_h = re[3]

        ret_l = bound_w*re_l + bound_l
        ret_b = bound_h*re_b + bound_b
        ret_w = bound_w*re_w
        ret_h = bound_h*re_h

        return [ret_l,ret_b,ret_w,ret_h]

class Event(object):
    """Abstract class to hold the data of an arbitrary time domain event"""
    def __init__(self, **kwargs):
        event_start = kwargs.pop('event_start',
                        quan.Quantity(0,'s',dtype = 'Float32')) 
        event_duration = kwargs.pop('event_duration', 
                        quan.Quantity(0,'s',dtype = 'float32'))
        event_type = kwargs.pop('event_type', '')
        event_name = kwargs.pop('event_name','')
        tunits = kwargs.pop('tunits','s')
        yunits = kwargs.pop('yunits','V')
        
        self.event_start = self.toquant(event_start,tunits)
        self.event_duration = self.toquant(event_duration,tunits)
        self.tunits = tunits
        self.yunits = yunits
        self.event_name = event_name
        self.event_type = event_type
        self.tunits_plot = 's'
        
    def set_type(self,typestr):
        self.event_type = typestr

    def set_name(self,namestr):
        self.event_name = namestr

    def get_type(self):
        return self.event_type

    def toquant(self,item,typestr=''):
        if type(item) == quan.Quantity:
            return item
        else:
            return quan.Quantity(item,typestr,dtype = 'float32')

    def ts(self,*keys):
        pass
        
#class CompoundEvent(Event):
#    """multiple simultaneously occuring events"""
#    def __init__(self,events,**kwargs):
#        Event.__init__(self,**kwargs)
#        self.events = cp.deepcopy(events)
        
class TimeData(Event):
    """class for generalized set of time-ordered data, accepts uneven sampling"""
    
    def __init__(self, event_start = quan.Quantity(0,'s'), 
                    event_name = '',
                    event_type = '',
                    event_duration = quan.Quantity(0,'s', dtype = 'float32'),
                    t = quan.Quantity([],'s',dtype = 'float32'), 
                    y = quan.Quantity([],'V',dtype = 'float32'),
                    tunits = 's',
                    yunits = 'V'):

        if (size(t) == 0):
            t = arange(0,size(self.y)*dt, self.dt)
            event_duration = t[-1] - t[0]
        elif not(size(t) == size(y)):
            raise EventException("X and Y must be of identical size")

        Event.__init__(self,
                        event_start = event_start,
                        event_duration = event_duration, 
                        event_name = event_name, 
                        event_type = event_type,
                        tunits = tunits,
                        yunits = yunits)

        self.t = self.toquant(t,tunits)
        self.y = self.toquant(y,yunits)
        
        self.interp = np.array(zeros(size(self.t)), dtype = bool)
        self.event_type = 'TimeData'
        
        self.tunits_plot = 's'
        
    def __add__(self,c):
        pass_t = c.t
        pass_y = c.y
        
        pass_y.units = self.y.units
        return TimeData(event_start = self.event_start,
                        event_duration = self.event_duration,
                        t = self.t, 
                        y = (self.y + pass_y))

    def __sub__(self,c):
        pass_t = c.t
        pass_y = c.y

        pass_y.units = self.y.units
        return TimeData(event_start = self.event_start,
                        event_duration = self.event_duration,
                        t = self.t, 
                        y = (self.y - pass_y))

    def __mul__(self,c):
        pass_t = c.t
        pass_y = c.y
        
        pass_y.units = self.y.units
        
        return TimeData(event_start = self.event_start,
                        event_duration = self.event_duration,
                        t = self.t, 
                        y = ((self.y * pass_y)))

    def __getitem__(self,key):
        if type(key) == slice:

            if (key.step == None):
                sample_stride = 1
            else:
                sample_stride = key.step
            
            t_index = arange(int(key.start),int(key.stop),int(sample_stride))
            return TimeData(event_start = self.event_start + key.start,
                            event_duration = self.t[t_index[-1]] 
                            - self.t[t_index[0]], 
                            y = self.y[ind])
        #if key is simply an index
        return TimeData(event_start = key, 
                        event_duration = self.dt, 
                        y = self.y[key])

    def ts(self,*key):
        if size(key) > 1:
            return "slice"
        else:
            return "key"

class TimeSeries(Event):
    """a class for generalized evenly sampled time series data"""
    def __init__(self,*args,**kwargs):
    #parse the **kwargs
        import neo
        tunits = kwargs.pop('tunits', 's')
        yunits = kwargs.pop('yunits', 'V')
        hdf5_key = kwargs.pop('hdf5_key', '')
        event_start = cp.copy(kwargs.pop('event_start', quan.Quantity(0,tunits))) 
        event_name = kwargs.pop('event_name','unnamed')
        event_type = kwargs.pop('event_type', 'time_series')
        event_duration = cp.copy(kwargs.pop('event_duration', quan.Quantity([0],tunits,dtype = 'float32')))
        t = cp.copy(kwargs.pop('t', quan.Quantity([],tunits,dtype = 'float32')))
        y = cp.copy(kwargs.pop('y', quan.Quantity([],yunits,dtype = 'float32'))) 
        dt = cp.copy(kwargs.pop('dt', quan.Quantity(1,tunits)))
        
        if len(args) > 0 and (type(args[0]) is neo.core.analogsignal.AnalogSignal):
            sig = args[0]
            tunits = sig.times.units
            yunits = sig.units
            event_start = sig.t_start
            event_name = sig.name
            t = cp.copy(sig.times)
            y = array(sig)
            dt = sig.sampling_period
                 
    #ensure that the three type parameters (t,y,and dt) are of Quantity type,
    #otherwise, convert them using tunits and yunits
    
        if not(type(t) == quan.Quantity):
            t = quan.Quantity(t,tunits)
        if not(type(y) == quan.Quantity):
            y = quan.Quantity(y,yunits)
        if not(type(dt) == quan.Quantity):
            dt = quan.Quantity(dt,tunits)
        if not(type(event_duration) == quan.Quantity):
            event_duration = quan.Quantity(event_duration,tunits)
        if not(type(event_start) == quan.Quantity):
            event_start = quan.Quantity(event_start,tunits)
        #print 'tunits ' + str(tunits)
        #print 'yunits ' + str(yunits)
        #print 'y.units ' + str(y.units)
    #call the super-class constructor
        Event.__init__(self,
                        event_start = event_start,
                        event_duration = event_duration,
                        event_name = event_name, 
                        event_type = event_type, 
                        tunits = tunits, 
                        yunits = yunits)
                        
    #assign the data to class members
        self.dt = dt
        self.y = y
        self.t = t
        
        self.hdf5_key = hdf5_key
        
    
        
    def __add__(self,c):
    #overload + operator so that only y is added.
        try: pass_y = cp.copy(c.y)
        except AttributeError:
            pass_y = cp.copy(c)
        ob = cp.copy(self)
        ob.y = (self.y + pass_y)
        return ob
        
    def __sub__(self,c):
    #overload - operator so that only y is subtracted.
        try: pass_y = c.y#cp.copy(c.y)
        except AttributeError:
            pass_y = c#cp.copy(c)
        ob = cp.copy(self)
        ob.y = (self.y - pass_y)
        return ob
        
    def __mul__(self,c):
    #overload * operator so that the y's are multiplied.
        try: pass_y = cp.copy(c.y)
        except AttributeError:
            pass_y = cp.copy(c)
        ob = cp.copy(self)
        ob.y = (self.y*pass_y)
        return ob
        
    def __div__(self,c):
    #overload * operator so that the y's are divided.
        try: pass_y = cp.copy(c.y)
        except AttributeError:
            pass_y = cp.copy(c)
        ob = cp.copy(self)
        ob.y = (self.y/pass_y)
        return ob
        
    def __pow__(self,c):
    #overload ** operator so that only y's are powered.
        try: pass_y = cp.copy(c.y)
        except AttributeError:
            pass_y = cp.copy(c)
        ob = cp.copy(self)
        ob.y = (self.y**pass_y)
        return ob
        
    def __getitem__(self,key):
    #overload slice operator cut the ts up by point index. Adjust event_start
    #as an offset from self.event_start.
        if type(key) == slice:
            if (key.step == None):
                sample_stride = 1
            else:
                sample_stride = key.step
            ob = cp.copy(self)
            #ob.event_start = quan.Quantity(self.event_start + (self.dt * key.start),self.tunits)
            #ob.event_duration = self.event_duration - ob.event_start
            #ob.y = self.y.__getitem__(key)
            ob.dt = cp.copy(self.dt*sample_stride)
            
            ob.y = cp.copy(type(self.y).__getitem__(self.y,key))
            return ob
            
    def __str__(self):
    #convert to string
        t = self.get_t_array()
        return str([self.y,t])

    def __array__(self):
    #cast as array, just return y data
        return np.array(self.y)
        
    def plot(self,**kwargs):
    #plot- requires matplotlib
        #,color = 'k'
        timebase = kwargs.pop('timebase','local')
        tunits_plot = kwargs.pop('tunits',self.tunits_plot)
        plot_type = kwargs.pop('plot_type','standard')
        tpoints = quan.Quantity(self.get_t_array(timebase = timebase),self.dt.units)
        if plot_type == 'standard':
            plot(tpoints.rescale(tunits_plot),self.y,**kwargs)
            return
        if plot_type == 'errorbar':
            errorbar(tpoints.rescale(tunits_plot),self.y,**kwargs)
            return
        
    def sub_baseline(self,start,stop,**kwargs):
    #baseline subtract the TimeSeries, start and stop define the interval,
    #in object time units
        tunits = kwargs.pop('tunits',self.dt.units)
        timebase = kwargs.pop('timebase','local')
        #print 'here'
        #print timebase
        #print timebase == 'global'
        if timebase == 'global':
            print 'global'
            start -= self.event_start.magnitude
            stop -= self.event_start.magnitude
        baseline = quan.Quantity(mean(self.ts(start,stop,tunits = tunits)),self.yunits)
        return (self - baseline)
        
    def integrate(self,*args,**kwargs):
    #Integrate the TimeSeries.
        tunits = kwargs.pop('tunits',self.dt.units)
        start, stop = 0,inf
        try:
            start = args[0]
            stop = args[1]
        except IndexError:
            pass
        units = self.y.units * self.dt.units
        return(quan.Quantity(integrate.trapz(self.ts(start,stop,tunits = tunits).y,dx = self.dt),units))
    
    #def mean(self,*args):
    #Get the mean of the TimeSeries.
#        start, stop = 0,inf
#        try:
#            start = args[0]
#            stop = args[1]
#        except IndexError:
#            pass
#        units = self.y.units
#        return(quan.Quantity(mean(self.ts(start,stop).y),units))

    def get_ylbl(self):
    #return a y label for plotting
        return "Current (" + str(self.yunits).split()[1] + ")"
    
    def get_tlbl(self):
    #return a t label for plotting
        return "Time (" + str(self.tunits).split()[1] + ")"
        
    def get_duration(self):
    #return the duration
        return (size(self.y)-1)*self.dt

    def get_t_array(self,timebase = 'local',**kwargs):
        tunits = kwargs.pop('tunits',self.dt.units)
    #return an array of t points corresponding to each y point.
        start = 0
        stop = float(self.get_duration().rescale(tunits))
        ar = linspace(start, stop, num = size(self.y))
        if(timebase == 'local'):
            return array(ar)
            #ar = arange(quan.Quantity(0,
            #            self.dt.units),
            #            self.get_duration(),
            #            self.dt,
            #            dtype = float64)
            #start = float64(self.event_start)
        if(timebase == 'global'):
            ar+=self.event_start.rescale(tunits)
            return array(ar)
        return ar
    
    def get_low_filter(self,Hz):
    #return a low-pass-filterd version of the time series
        import filtfilt
        filter_order = 3
        sampfreq = 1/self.dt.rescale('sec')
        filtfreq = quan.Quantity(Hz,'Hz')
        passband_peram = float((filtfreq/sampfreq).simplified)
        from scipy.signal import butter
        [b,a]=butter(filter_order,passband_peram)
        ob = cp.copy(self)
        ob.y = quan.Quantity(filtfilt.filtfilt(b,a,self.y),self.yunits)
        return ob
    
    def get_normalized(self):
    #normalize so that y goes from 0 to 1
        n = self.y-np.amin(self.y)
        n /= np.amax(n)
        return array(n)
        
    def get_thresh_crossings(self,thresh,tunits = 's'):
    #find the timepoints that cross a threshold
        mask = (self.y > thresh).astype(int8)
        d = np.diff(mask)
        idx, = d.nonzero()
        if tunits == 'points':
            return idx
        else:
            return self.get_t_array(tunits = tunits)[idx]
    
    def low_filter(self,Hz):
    #filter signal in place
        self.y = self.get_low_filter(Hz).y
    
    def get_trend(self,ksize = 151):
        ob = cp.copy(self)
        from scipy.signal import medfilt
        #print ksize
        ob.y = medfilt(self.y,kernel_size = ksize)
        ob.y = quan.Quantity(ob.y,self.yunits)
        return ob
        
    def get_detrend(self,ksize = 299):
        return self - self.get_trend(ksize = ksize)

    def ts(self,*key,**kwargs):
    #time slice - cut up a wave given start and end times tunits corresponds
    #to the units for the passed start and stop params
        tunits = kwargs.pop('tunits',self.dt.units)
        dtunits = self.dt.units
        ind = [None,None,None]
        if size(key) == 3:
            ind[0] = int(quan.Quantity(key[0],tunits).rescale(dtunits)/self.dt)
            ind[1] = int(quan.Quantity(key[1],tunits).rescale(dtunits)/self.dt)
            ind[2] = int(quan.Quantity(key[2],tunits).rescale(dtunits)/self.dt)
        elif size(key) == 2:
            ind[0] = int(quan.Quantity(key[0],tunits).rescale(dtunits)/self.dt)
            ind[1] = int(quan.Quantity(key[1],tunits).rescale(dtunits)/self.dt)
            ind[2] = 1
        elif size(key) == 1:
            ind[0] = int(quan.Quantity(key[0],tunits).rescale(dtunits)/self.dt)
            ind[1] = int(quan.Quantity(key[0],tunits).rescale(dtunits)/self.dt) + 1
            ind[2] = 1
        selection = slice(ind[0],ind[1],ind[2])
        ob = self.__getitem__(selection)
        return ob
        
    def tr(self,*key,**kwargs):
    #get an array of timepoints given start and stop times, again tuints corresponds
    #to the units of the passed params
        tunits = kwargs.pop('tunits',self.dt.units)
        dtunits = dtunits = self.dt.units
        k0 = float(quan.Quantity(key[0],tunits).rescale(dtunits))
        k1 = float(quan.Quantity(key[1],tunits).rescale(dtunits))
        retarray = quan.Quantity(np.arange(k0,k1,float(self.dt)),tunits)
        return np.array(retarray)
        
    def set_start(self,starttime,**kwargs):
        #set the global starttime
        tunits = kwargs.pop('tunits',self.dt.units)
        self.event_start = quan.Quantity(starttime,tunits) 
    
    def set_yunits(self,units):
        self.yunits = units
        self.y.units = units

    def set_tunits(self,units):
        self.tuints = units
        self.dt.units = units
        self.t.units = units
        
    def animate(self,n,speed,box):
        import time
        ion()
        #speed = 100
        tstart = time.time()
        y = self.ts(0,box,tunits = 's')
        x = y.get_t_array()
        pnts = len(x)
        line, = plot(x,np.ma.masked_where(x>0,np.array(y.y)))
        ax = gca()
        ax.set_ybound([self.y.min(),self.y.max()])
        ax.set_xbound([0,x[-1]*1.1])
        for i in arange(1,len(x),speed):
            line.set_ydata(np.ma.masked_where(x/self.dt>i,np.array(y.y)))
            draw()
        for i in arange(1,n):
            line.set_ydata(np.array(self.y[i*speed:i*speed+pnts]))
            draw()
        print 'FPS:' , (n+(len(x)/speed))/(time.time()-tstart)
            
        
                      
class StimSig(TimeSeries):
#    all_ts = {}
#    created = 0
    
    def __init__(self,events,**kwargs):
#        print StimSig.created
#        StimSig.all_ts[id(self)] = StimSig.created
#        StimSig.created += 1
        Event.__init__(self,**kwargs)
        self.events = events.copy() #cp.copy(events)
        self.events['signals'] = cp.copy(events['signals'])
        self.events['stims'] = cp.copy(events['stims'])
        #self.events['signals'] = events['signals']
        #self.events['stims'] = events['stims']
        #TimeSeries.__init__(self,**kwargs)
        kwargs.pop('event_duration','')
        kwargs.pop('event_start','')
    
#    def __del__(self):
#        StimSig.all_ts[id(self)] = 0
        
    def __copy__(self):
        #ob = StimSig.__new__(StimSig)
        #ob.__init__(self.events,event_start = self.event_start,event_duration = self.event_duration)
        ob = StimSig(self.events,event_start = self.event_start,event_duration = self.event_duration)
        #ob.__init__(cp.copy(self.events),event_start = cp.copy(self.event_start),event_duration = cp.copy(self.event_duration))
        #ob.events['signals'] = cp.copy(self.events['signals'])
        #ob.events['stims'] = cp.copy(self.events['stims'])
        return ob
        
    def __setattr__(self,name,value):
        if name == 'y':
            self.events['signals'].y = cp.copy(value)#cp.copy(value)#quan.Quantity(value)
        else:
            object.__setattr__(self, name, value)
            
    def __getattribute__(self,name):
        if name in ['y','dt','event_duration','event_start','tunits','yunits']:
            return self.events['signals'].__getattribute__(name)
        else:
            return object.__getattribute__(self,name)
    
    def __getitem__(self,key):
        if type(key) == slice:
            if (key.step == None):
                sample_stride = 1
            else:
                sample_stride = key.step
        ob = cp.copy(self)
        #ob.events['signals'] = cp.copy(self.events['signals'].__getitem__(key))
        #ob.events['stims'] = cp.copy(self.events['stims'].__getitem__(key))
        ob.events['signals'] = self.events['signals'].__getitem__(key)
        ob.events['stims'] = self.events['stims'].__getitem__(key)
        return ob
        
    def stim_plot():
        subplot(2,1,1)
        self.events['stims'].plot()
        self.events['signals'].plot()
        
class EventRecord(list):
    """Holds multiple events aligned on a common time-base"""
    def __init__(self,data,**kwargs):
        self.max_dt = kwargs.pop('max_dt',1)
        self.record_dt = kwargs.pop('record_dt', quan.Quantity(1,'s'))
        self.plot_panel = kwargs.pop('panel',SubPanel([0,0,1,1]))
        self.record_name = kwargs.pop('record_name','')
        self.record_type = kwargs.pop('record_type','')
        self.event_types = set('empty')
        self.record_duration = quan.Quantity(0,'s')
        list.__init__([])
        if (data != None):
            [self.append(x) for x in data]
    
    def __copy__(self):
#        ob = list.__copy__(self) #type(self).__new__(type(self))
        return EventRecord(self,max_dt = self.max_dt,plot_panel = self.plot_panel, record_name = self.record_name,event_types =self.event_types,record_duration = self.record_duration)
        #ob.max_dt = self.max_dt
        #ob.record_dt = self.record_dt
        #ob.plot_panel = self.plot_panel
        #ob.record_name = self.record_name
        #ob.event_types = cp.copy(self.event_types)
        #ob.record_duration = self.record_duration
        #return ob
    
    #def __deepcopy__(self):
    #    ob = EventRecord(self)
    #    ob.max_dt = self.max_dt
    #    ob.record_dt = self.record_dt
    #    ob.plot_panel = self.plot_panel
    #    ob.record_name = self.record_name
    #    ob.event_types = cp.deepcopy(self.event_types)
    #    ob.record_duration = self.record_duration
    #    return ob
    
    def append(self,event):
        list.append(self,event)
        self.event_types.add(event.event_type)
        if((event.event_start + event.event_duration) > self.record_duration):
            self.duration = event.event_start + event.event_duration

    def num_items(self):
        return size(self)
        
    def hdf5_write(self,node):
        """Write event record data to a pytable HDF5 file node. If the node is a leaf: update. If the node is a group: create a new leaf and write. Return a reference to the existing/new leaf"""
        pass
    
    def hdf5_load(self,node):
        """construct the event record from a given HDF5 file leaf. Check that it conforms to the correct type."""
        pass
        
#class ChannelRecord(EventRecord):
#    """Holds data from a given DAC or ADC channel - forces the contained TimeSeries objects to maintain identical units.  Constructor can accept a simple matrix. Convenience method for generating a group of TimeSeries data"""
    
#    def __init__(self,data,**kwargs):
    
#        self.yunits = kwargs.pop('yunits')#,quan.Quantity(1,'V'))
        #if type(self.yunits) != quan.Quantity:
        #    self.yunits = quan.Quantity(self.yunits[0],self.yunits[1])
    
#        self.tunits = kwargs.pop('tunits')#,quan.Quantity(1,'s'))
        #if type(self.tunits) != quan.Quantity:
        #    self.tunits = quan.Quantity(self.tunits[0],self.tunits[1])
    
#        self.dt = kwargs.pop('dt')#,quan.Quantity(1,'s'))
        #if type(self.dt) != quan.Quantity:
        #    self.dt = quan.Quantity(self.dt[0],self.dt[1])
        
#        self.channel_name = kwargs.pop('channel_name', 'data_channel')
    
#        data = [TimeSeries(y = row, dt = cp.copy(self.dt), event_type = self.channel_name, tunits = self.tunits, yunits = self.yunits) for row in data] 
#        EventRecord.__init__(self,data = data,**kwargs) 
    
class EpisodicRecord(EventRecord):
    """Parent class for episodic stimulation protocols"""
    def __init__(self, signals, stims, **kwargs):
        self.signals = signals #list of list of signals
        self.stims = stims #list of list of stimulation protocols
        tracedicts = [StimSig({'signals':x,'stims':y}) 
                      for x,y in zip(signals,stims)]
        #tracedicts = [StimSig({'signals':cp.copy(x),'stims':cp.copy(y)}) 
        #              for x,y in zip(signals,stims)]
        EventRecord.__init__(self,tracedicts)
   
    def __copy__(self):
        #ob = EventRecord.__copy__(self)
        ob =  type(self).__new__(type(self))
        ob.__init__([x.events['signals'] for x in self],[x.events['stims'] for x in self])
        #[ob.append(cp.copy(x)) for x in self]
        return ob
        
    def __add__(self,c):
        #ob = EventRecord.__copy__(self)
        #ob = cp.copy(self)
        ob = type(self).__new__(type(self))
        if isinstance(c,EpisodicRecord):
             ob.__init__([x.events['signals']+y.events['signals'] for x,y in zip(self,c)],[x.events['stims'] for x in self])
            #[ob.append(x+y) for x,y in zip(self,c)]
        else:
            ob.__init__([x.events['signals']+c for x in self],[x.events['stims'] for x in self])
            #[ob.append(x+c) for x in self]
        return ob
    
    def __sub__(self,c):
        ob = self[0:0]#cp.copy(self)
        #ob = EventRecord.__copy__(self)
        #ob = EpisodicRecord.__new__(EpisodicRecord)
        ob = type(self).__new__(type(self))
        if isinstance(c,EpisodicRecord):
            ob.__init__([x.events['signals']-y.events['signals'] for x,y in zip(self,c)],[x.events['stims'] for x in self])
            #[ob.append(x-y) for x,y in zip(self,c)]
        else:
            ob.__init__([x.events['signals']-c for x in self],[x.events['stims'] for x in self])
            #[ob.append(x-c) for x in self]
        return ob
        
    def __mul__(self,c):
        #ob = EventRecord.__copy__(self)
        #ob = cp.copy(self)
        ob = type(self).__new__(type(self))
        if isinstance(c,EpisodicRecord):
            ob.__init__([x.events['signals']*y.events['signals'] for x,y in zip(self,c)],[x.events['stims'] for x in self])
            #[ob.append(x*y) for x,y in zip(self,c)]
        else:
            ob.__init__([x.events['signals']*c for x in self],[x.events['stims'] for x in self])
            #[ob.append(x*c) for x in self]
        return ob
        
    def __div__(self,c):
        #ob = EventRecord.__copy__(self)
        #ob = cp.copy(self)
        ob = type(self).__new__(type(self))
        if isinstance(c,EpisodicRecord):
            ob.__init__([x.events['signals']/y.events['signals'] for x,y in zip(self,c)],[x.events['stims'] for x in self])
            #[ob.append(x/y) for x,y in zip(self,c)]
        else:
            ob.__init__([x.events['signals']/c for x in self],[x.events['stims'] for x in self])
        return ob
        
class FamilyRecord(EpisodicRecord):
    """Container class with ploting methods for the data from a voltage or 
    current family"""
    def __init__(self,
                response_list = [],
                cmnd_list = [],
                record_name = 'unnamed', 
                record_type= 'Family'):
        if not len(response_list) == len(cmnd_list):
            raise EventException(
                    "response_list and cmnd_list: must be of identical size")
        
        EpisodicRecord.__init__(self,
                                response_list,
                                cmnd_list,
                                record_name = record_name, 
                                record_type = record_type)

    def plot(self):
        self.panel_plot([0.1,0.05,.8,.85])
    
    def get_panel(self,bounds,new_panel = None):
        if not(new_panel == None):
            panel = new_panel
        else:
            panel = self.plot_panel
            
        ax1_rect = [0.0,0.1,1,0.9]
        ax2_rect = [0.0,0.0,1,0.09]

        ax1_rect = panel.rect_remap(ax1_rect,bounds = bounds)
        ax2_rect = panel.rect_remap(ax2_rect,bounds = bounds)
        
        ax1 = axes(ax1_rect)
        ax2 = axes(ax2_rect,sharex = ax1)
        return([ax1,ax2])
        
    def hdf5_write(self,group):
        """Write family record data to a pytable HDF5 file node. If the node is 
        a leaf: update. If the node is a group: create a new leaf and write. 
        Return a reference to the existing/new leaf"""
        #get the dimentionality of the event matrix
        rows = len(self.event_list)
        columns = size(self.event_list[0].y)
        print("rows:" + str(rows) + "  columns:" + str(columns))
        #get the name of the event
        r_name = self.record_name
        r_title = self.record_type
        #create a subgroup in group to hold the data. 
        #Set the attributes of the group to describe the event.
        wfile = group._v_file
        try:
            wfile.removeNode(group,name = r_name)
        except pyt.NodeError:
           pass
           
        subgroup = wfile.createGroup(group,r_name)
            
    def hdf5_load(self,node):
        """construct the family record from a given HDF5 file leaf. Check that 
        it conforms to the correct type."""
        pass
        
    
class IVCurve(object):
    def __init__(self,i,v):
        self.i = i
        self.v = v
        self.fit,self.i_c = self.linear_leak_transform()
        
    def linear_leak_transform(self):
        a,b = np.polyfit(self.v[0:3],self.i[0:3],1)
        y = np.polyval([a,b],self.v)
        return y,np.array(self.i)-y
        
    def full_regress(self):
        a,b = np.polyfit(self.v,self.i,1)
        y = np.polyval([a,b],self.v)
        return y
        
    def get_ohmic_fit(self):
        return np.polyfit(self.v,self.i,1)
        
    def plot(self,**kwargs):
        mode = kwargs.pop('mode','normal')
        ax = axes()
        if mode == 'normal':
            plot(self.v,self.i,'-o',label = "ljp & SR corrected IV")
            plot(self.v,self.fit, label = "Passive-leak estimation")
            plot(self.v,self.i_c,'-o', label = "Passive-leak-subtracted IV")
        if mode == 'rp':
            plot([v-quan.Quantity(15,'mV') for v in self.v],self.i,'o',ms=8,label = "Peak amplitude of synaptic current")
            plot([v-quan.Quantity(15,'mV') for v in self.v],self.full_regress(),lw=2, label = "Linear Fit")
            ax.annotate("E_Cl:-47mV", xy=(-47, 2),  xycoords='data',
                xytext = (-47,4),textcoords = 'data',
                horizontalalignment='center', verticalalignment='center',
                arrowprops = dict(arrowstyle ='->'))
               
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        legend(loc = 0)
        override = {'fontsize': 'large','position': (0,0)}
        xlabel(unicode("Vm (" + self.v[0].dimensionality.unicode + ")"),override)
        override = {'fontsize': 'large','position': (0,0),'rotation':'horizontal'}
        if mode == 'normal':
            ylabel(unicode("Im (" + self.i[0].dimensionality.unicode + ")"),override)
        if mode == 'rp':
            ylabel(unicode("Isyn (" + self.i[0].dimensionality.unicode + ")"),override)
        
        
    def __add__(self,c):
        return IVCurve([s+o for s,o in zip(self.i,c.i)],self.v)
    
    def __sub__(self,c):
        return IVCurve([s-o for s,o in zip(self.i,c.i)],self.v)
                    
#        ax = axes()
#        plot(iv.v,iv.i,'-o',label = "ljp & SR corrected IV")
#        ax.spines['left'].set_position('zero')
#        ax.spines['bottom'].set_position('zero')
#        ax.spines['top'].set_visible(False)
#        ax.spines['right'].set_visible(False)
#       ax.yaxis.set_ticks_position('left')
#        ax.xaxis.set_ticks_position('bottom')
#        legend(loc = 0)
#        override = {'fontsize': 'large','position': (0,0)}
#        xlabel(unicode("Vm (" + iv.v[0].dimensionality.unicode + ")"),override)
#        override = {'fontsize': 'large','position': (0,0),'rotation':'horizontal'}
#        ylabel(unicode("Im (" + iv.i[0].dimensionality.unicode + ")"),override)
        
class CurrentFamily(FamilyRecord):
    """Holds the data from a Famly of Currents evoked from series of voltage 
    pulses"""
    
    def __init__(self,
                 response_list = [],
                 cmnd_list = [],
                 record_name = 'unnamed', 
                 record_type = 'current_family'):
        FamilyRecord.__init__(self,
                              response_list = response_list,
                              cmnd_list = cmnd_list,
                              record_name = record_name, 
                              record_type = record_type)
    
    def choose_iv_nslope(self, *args,**kwargs):
        """ finds the i_v with the greatest n-ness in the slope
            pass, start, stop and number of bins, optionally pass the time units"""
        tunits = kwargs.pop('tunits','s')
        #when to look
        try:
            start = args[0]
            stop = args[1]
            bins = args[2]
        except IndexError:
            start, stop, bins = 0.02,0.1,10
        i_list = list()
        p = np.linspace(start,stop,bins)
        startpoints = p[:-1]
        stoppoints = p[1:]
        for sta,sto in zip(startpoints,stoppoints):
            i = [ssig.events['signals'].ts(sta,sto,tunits = tunits).y.mean() for ssig in self]
            i_list.append((i,sta,sto))
        def ilist_comp(x,y):
            dx = np.diff(np.array(x[0]))
            dy = np.diff(np.array(y[0]))
            if dx[np.argmin(dx)] < dy[np.argmin(dy)]:
                return 1
            elif dx[np.argmin(dx)] == dy[np.argmin(dy)]:
                return 0
            else:
                return -1
        #print i_list
        i_list.sort(ilist_comp)
        i = i_list[-1][0]
        v = [ ssig.events['stims'].ts(start,stop,tunits = tunits).y.mean() for ssig in self]
        rcurve = IVCurve(i,v)
        rcurve.sampled_at = i_list[-1][1]
        return rcurve
                     
        
    def get_iv(self,*args,**kwargs):
        """Return an IV curve - can pass start and stop time, defaults to sec"""
        tunits = kwargs.pop('tunits','s')
        mode = kwargs.pop('mode','mean')
        try:
            start = args[0]
            stop = args[1]
        except IndexError:
            start, stop = 0.02,0.1
        if mode is 'mean':
            i = [ ssig.events['signals'].ts(start,stop,tunits = tunits).y.mean() for ssig in self]
        else:
            i = [ np.abs(ssig.events['signals'].ts(start,stop,tunits = tunits).y).argmax() for ssig in self]
            i = [ ssig.events['signals'].ts(start,stop,tunits = tunits).y[ind] for ind, ssig in zip(i,self)]
        v = [ ssig.events['stims'].ts(start,stop,tunits = tunits).y.mean() for ssig in self]
        return IVCurve(i,v)
        
    def sr_ljp_correct(self,ra,ljp):
        """series resistance compensation"""
        signals = [ssig.events['signals'] for ssig in self]
        errors = [sig*ra for sig in signals]
        corrected = [stms.events['stims']-e-ljp for e,stms in zip(errors,self)]
        for i in corrected:
            i.y.units = 'mV'
        return CurrentFamily(response_list = signals,cmnd_list = corrected)
                       
    def panel_plot(self,bounds,**kwargs):
        new_panel = kwargs.pop('new_panel',None)
        [ax1, ax2] = FamilyRecord.get_panel(self,bounds,new_panel = new_panel)
        
        axes(ax1)
        [x.events['signals'].set_yunits('pA') for x in self]
        [x.events['signals'].plot(**cp.copy(kwargs)) for x in self]
        #formatting
        xlbs = ax1.get_xticklabels()
        [l.set_visible(False) for l in xlbs]
        ylabel(unicode("Im (" + self[0].events['signals'].y.dimensionality.unicode + ")",'utf-8'))
        ####
        
        axes(ax2)
        [x.events['stims'].plot(**cp.copy(kwargs)) for x in self]
        
        #formatting
        ylbs = ax2.get_yticklabels()
        [l.set_visible(False) for l in ylbs]
        [l.set_visible(True) for l in ylbs[::2]]
        ylabel(unicode("Vm (" + self[0].events['stims'].y.dimensionality.unicode + ")",'utf-8'))
        xlabel(unicode("Time (" + self[0].events['stims'].tunits_plot + ")"))
        ####
        
        ax1.axis('tight')
        ax2.axis('tight')
        kill_spines(ax1)
        kill_spines(ax2)
    
    def plot(self,**kwargs):
        self.panel_plot([0.1,0.1,.8,.85],**kwargs)

class StimCurrentFam(CurrentFamily):
    
    def __init__(self,
                response_list = [],
                cmnd_list = [],
                stim = None,
                record_name = 'unnamed', 
                record_type = 'current_family'):
        CurrentFamily.__init__(self,
                            response_list = response_list,
                            cmnd_list = cmnd_list,
                            record_name = record_name,
                            record_type = record_type)
        self.stimtrace = stim
        
    def sr_ljp_correct(self,ra,ljp):
        """series resistance compensation"""
        signals = [ssig.events['signals'] for ssig in self]
        errors = [sig*ra for sig in signals]
        corrected = [stms.events['stims']-e-ljp for e,stms in zip(errors,self)]
        for i in corrected:
            i.y.units = 'mV'
        return StimCurrentFam(response_list = signals,cmnd_list = corrected,stim = cp.copy(self.stimtrace))
        
    
    def get_panel(self,bounds,new_panel = None):
        if not(new_panel == None):
            panel = new_panel
        else:
            panel = self.plot_panel
        ax0_rect = [0.0,0.85,1,0.1]
        ax1_rect = [0.0,0.1,1,0.7]
        ax2_rect = [0.0,0.0,1,0.09]
        
        ax0_rect = panel.rect_remap(ax0_rect,bounds = bounds)
        ax1_rect = panel.rect_remap(ax1_rect,bounds = bounds)
        ax2_rect = panel.rect_remap(ax2_rect,bounds = bounds)
        
        
        ax1 = axes(ax1_rect)
        ax2 = axes(ax2_rect,sharex = ax1)
        ax0 = axes(ax0_rect,sharex = ax2)
        return([ax0,ax1,ax2])

    def panel_plot(self,bounds,**kwargs):
        new_panel = kwargs.pop('new_panel',None)
        [ax0,ax1,ax2] = self.get_panel(bounds)
        axes(ax0)
        kwa = cp.copy(kwargs)
        kwa.update({'color':'b'})
        self.stimtrace.plot(**kwa)
        
        axes(ax1)
        #print 'here'
        [x.events['signals'].set_yunits('pA') for x in self]
        [x.events['signals'].plot(**cp.copy(kwargs)) for x in self]
        #formatting
        xlbs = ax1.get_xticklabels()
        [l.set_visible(False) for l in xlbs]
        ylabel(unicode("Im (" + self[0].events['signals'].y.dimensionality.unicode + ")",'utf-8'))
        ####
        
        axes(ax2)
        [(x.events['stims']-quan.Quantity(15,'mV')).plot(**cp.copy(kwargs)) for x in self]
        
        #formatting
        ylbs = ax2.get_yticklabels()
        [l.set_visible(False) for l in ylbs]
        [l.set_visible(True) for l in ylbs[::2]]
        ylabel(unicode("Vm (" + self[0].events['stims'].y.dimensionality.unicode + ")",'utf-8'))
        xlabel(unicode("Time (" + self[0].events['stims'].tunits_plot + ")"))
        ####
        ax0.axis('tight')
        ax1.axis('tight')
        ax2.axis('tight')
        ax1.set_ylim([-70,100])
        ax2.set_xlim([0,3])
        kill_spines(ax0)
        kill_spines(ax1)
        kill_spines(ax2)
        ax0.spines['bottom'].set_visible(False)
        ax0.spines['left'].set_visible(False)
        setp(ax0.get_xticklabels(), visible=False)
        setp(ax0.get_xticklines(), visible=False)
        #ax0.set_xticks([])
        ax0.set_yticks([])
        ax2.set_yticks(ax2.get_yticks()[::3])
    
    def plot(self,**kwargs):
        self.panel_plot([0.1,0.1,.8,.85],**kwargs)
        
    def get_amp_curve(self,**kwargs):
        """Return an evoked-current vs V curve - can pass start and stop time, defaults to sec"""
        tunits = kwargs.pop('tunits','s')
        base_epoch = kwargs.pop('base_epoch',[0.88,0.92])
        stim_epoch =  kwargs.pop('stim_epoch',[0.97,1.04])
        baseline = self.get_iv(base_epoch[0],base_epoch[1],tunits = 's')
        peak = self.get_iv(stim_epoch[0],stim_epoch[1],tunits = 's',mode = 'peak')
        #print "here"
        return peak-baseline
    
class VoltageFamily(FamilyRecord): 
    """Holds the data from a Famly of Voltages evoked from series of current 
    pulses"""
    def __init__(self,
                 response_list = [],
                 cmnd_list = [],
                 record_name = 'unnamed', 
                 record_type = 'Voltage_family'):
        FamilyRecord.__init__(self,
                              response_list = response_list,
                              cmnd_list = cmnd_list,
                              record_name = record_name, 
                              record_type = record_type)
        
    def panel_plot(self,bounds,**kwargs):
        new_panel = kwargs.pop('new_panel',None)
        [ax1, ax2] = FamilyRecord.get_panel(self,bounds,new_panel = new_panel)
        
        axes(ax1)
        [x.events['signals'].plot() for x in self]
        xlbs = ax1.get_xticklabels()
        [l.set_visible(False) for l in xlbs]
        
        axes(ax2)
        [x.events['stims'].plot() for x in self]
        
        for trace in ax1.get_lines():
            trace.set_color('black')

        for trace in ax2.get_lines():
            trace.set_color('black')

        ax1.axis('tight')
        ax2.axis('tight')
    
    def plot(self,**kwargs):
        self.panel_plot([0.1,0.05,.8,.85],**kwargs)

class OpticalFamily(FamilyRecord):
    def __init__(self,
                 response_list = [],
                 cmnd_list = [],
                 record_name = 'unnamed', 
                 record_type = 'Optical_family'):
        
        FamilyRecord.__init__(self,
                              response_list = response_list,
                              cmnd_list = cmnd_list,
                              record_name = record_name, 
                              record_type = record_type)
    
    def err_plot(self,bounds,err_fam,**kwargs):
        new_panel = kwargs.pop('new_panel',None)
        c_map = kwargs.pop('c_map',cm.hot)
        color_start,color_stop = kwargs.pop('c_range',(0,0.7))
        [ax1, ax2] = FamilyRecord.get_panel(self,bounds,new_panel = new_panel)
        
        axes(ax1)
        
        #Get a range of colormap indices to use
        l = len(self)
        
        for i,(e,s) in enumerate(zip(reversed(err_fam),reversed(self))):
            c_val = (color_stop/(l+1)) * (i+1) + color_start
            pse = s + e
            mse = s - e
            xe,ye = poly_between(s.get_t_array(tunits = 's'),pse.y,mse.y)
            ax1.fill(xe,ye,fc = c_map(c_val),ec = 'None', alpha = 0.2)
            s.plot(color = c_map(c_val),lw=1)
        axes(ax2)
        
        [x.events['stims'].plot() for x in reversed(self)]
        ax2.set_yscale('log')
        
        
        #for i,trace in enumerate(ax1.get_lines()):
        #store the colormap index for a line and put it in val
        #    val = (color_stop/(l+1)) * (i+1)
        #    trace.set_color(cm.hot(val))
        #   trace.set_linewidth(3)

        for i,trace in enumerate(ax2.get_lines()):
            c_val = (color_stop/(l+1)) * (i+1) + color_start
            trace.set_color(c_map(c_val))

        ax1.axis('tight')
        ax2.axis('tight')
                
    def panel_plot(self,bounds,**kwargs):
        new_panel = kwargs.pop('new_panel',None)
        
        c_map = kwargs.pop('c_map',cm.hot)
        color_start,color_stop = kwargs.pop('c_range',(0,0.7))

        [ax1, ax2] = FamilyRecord.get_panel(self,bounds,new_panel = new_panel)
        
        axes(ax1)
        [x.events['signals'].plot() for x in reversed(self)]
        #formatting
        xlbs = ax1.get_xticklabels()
        [l.set_visible(False) for l in xlbs]
        lb = self.get_label_str(self[0].events['signals'].y.dimensionality.string)
        ylabel(lb)
        ####
        
        axes(ax2)
        [x.events['stims'].plot() for x in reversed(self)]
        
        #formatting
        ylbs = ax2.get_yticklabels()
        [l.set_visible(False) for l in ylbs]
        [l.set_visible(True) for l in ylbs[::2]]
        
        lb = self.get_label_str(self[0].events['stims'].y.dimensionality.string)
        
        ylabel(unicode(lb))
        
        lb = self.get_label_str(self[0].events['stims'].tunits_plot)
        xlabel(unicode(lb))
        #xlabel(unicode("Time (" + self[0].events['stims'].tunits_plot + ")"))
        
        ####

        ax2.set_yscale('log')
        #Get a range of colormap indices to use
        l = len(ax1.get_lines())
        
        for i,trace in enumerate(ax1.get_lines()):
            #stor the colormap index for a line and put it in val
            c_val = (color_stop/(l+1)) * (i+1)
            trace.set_color(c_map(c_val))
            trace.set_linewidth(1)

        for i,trace in enumerate(ax2.get_lines()):
            c_val = (color_stop/(l+1)) * (i+1)
            trace.set_color(c_map(c_val))

        ax1.axis('tight')
        ax2.axis('tight')
        kill_spines(ax1)
        kill_spines(ax2)
        
    def plot_stims(self,**kwargs):
        """call the plot command on all the stims"""
        c_map = kwargs.pop('c_map',cm.hot)
        color_start,color_stop = kwargs.pop('c_range',(0,0.7))
        [x.events['stims'].plot() for x in reversed(self)]
        ax1 = gca()
        l = len(ax1.get_lines())
        for i,trace in enumerate(ax1.get_lines()):
            #store the colormap index for a line and put it in val
            c_val = (color_stop/(l+1)) * (i+1)
            trace.set_color(c_map(c_val))
            trace.set_linewidth(1)

    def plot_signals(self,**kwargs):
        """call the plot command on all the signals"""
        c_map = kwargs.pop('c_map',cm.hot)
        color_start,color_stop = kwargs.pop('c_range',(0,0.7))
        [x.events['signals'].plot() for x in reversed(self)]
        ax1 = gca()
        l = len(ax1.get_lines())
        for i,trace in enumerate(ax1.get_lines()):
            #store the colormap index for a line and put it in val
            c_val = (color_stop/(l+1)) * (i+1)
            trace.set_color(c_map(c_val))
            trace.set_linewidth(1)
    
    def plot(self,**kwargs):
        try:
            err_fam = kwargs.pop('err_fam')
            self.err_plot([0.1,0.05,.8,.85],err_fam,**kwargs)
            return
        except KeyError:
            pass
        try:
            fig_type = kwargs.pop('subtype')
            if fig_type == 'stims':
                self.plot_stims(**kwargs)
            elif fig_type == 'signals':
                self.plot_signals(**kwargs)
        except KeyError:
            self.panel_plot([0.1,0.1,.8,.85],**kwargs)
        
    def get_label_str(self,unit_str):
        a = quan.Quantity('1',unit_str)
        measuredic = {'Vm':'V','Im':'A','Time':'s','Irradiance\n':'mW/(mm*mm)'}
        [measuredic.update({s:quan.Quantity(1,measuredic[s])}) for s in measuredic.keys()]
        for q in measuredic.keys():
            if measuredic[q].dimensionality.simplified == a.dimensionality.simplified:
                return unicode(q + "(" + a.dimensionality.unicode + ")",'utf-8')
        
class StartSlope():
    def __init__(self,fam):
        start = 2.2
        stop = 5
        self.trace = cp.copy(fam[-1])
        self.tstart = quan.Quantity(start,'ms')
        self.tstart.units = self.trace.dt.units
        self.ydat = array(self.trace.ts(start,stop,tunits ='ms').y)
        self.xdat = array(self.trace.ts(start,stop,tunits ='ms').get_t_array())
        self.fit_slope()
        self.fit_tau()
        
    def fit_slope(self):
        self.slope,self.intercept = np.polyfit(self.xdat,self.ydat,1)
        self.fit = np.polyval([self.slope,self.intercept],self.xdat)
        
    def fit_tau(self):
         #get the initial guess by linear regression on log transformed data
        lg,offset = psl.lntransform(self.ydat) 
        m,b = polyfit(self.xdat,lg,1)

        ##----reduced model= single exponential----###
        #set perameters
        self.i_not = psl.Parameter(e**b + offset,"i_not")
        self.tau = psl.Parameter(abs(1/m),"tau")
        self.i_sted = psl.Parameter(offset,"i_sted")
        self.reduced = psl.FitModel([self.i_not,self.tau,self.i_sted],psl.single_exp,self.xdat,self.ydat,model_name = "single_exp")
    
    def plot(self):
        self.trace.ts(0,6,tunits = 'ms').plot(tunits = 'us',label = 'trace')
        xcords = self.xdat+float(self.tstart)
        #plot(xcords,self.fit,lw = 6, alpha = 0.6,label = 'fit')
        expfit = self.reduced.fx(self.xdat)
        plot(xcords,expfit, lw = 6, alpha = 0.6, label = 'expfit') 
        legend()
        ylabel(unicode("Im (" + self.trace.y.dimensionality.unicode + ")",'utf-8'))
        xlabel(unicode("Time (" + self.trace.dt.dimensionality.unicode + ")"))
        kill_spines(gca())
               
        
class CapRecord(dict):
    """Holds the data for a set of capacitive transients and implements some 
    plotting methods"""
    
    def __init__(self,
                 oncell_cap,
                 incell_cap,
                 sub_cap,
                 prot,
                 panel = SubPanel([0,0,1,1,])):
        """Initialize with 3 TimeSeries objects: the oncell capacitive transient
         object, the incell capacitive transient object and the subtracted 
         capacitive transient and the step protocol.""" 
        dict.__init__(self)
        oncell_cap.set_yunits('pA')
        incell_cap.set_yunits('pA')
        sub_cap.set_yunits('pA')
        
        self.event_keys= ['oncell_cap',
                          'incell_cap',
                          'sub_cap']
                          
        self['oncell_cap'] = oncell_cap
        self['incell_cap'] = incell_cap
        self['sub_cap'] = sub_cap
        self['prot'] = prot
        self.calc_capacitance()
        #dictionary to hold the string representation of 
        #each of the three traces.
        self.str_items = {'oncell':False,'incell':False,'wholecell':True} 
        self.plot_panel = panel
        
    def __str__(self):
        """Returns calculated parameters for the oncell, 
        incell and whole_cell transients. Future versions should utilize a 
        settable class atribute that indicates which sets of parameters are 
        returned -- the whole_cell data being the most relevent and probably 
        the default"""
        retstr = ''
        if self.str_items['oncell']:
            retstr += "------oncell parameters-----\n"
            retstr += "Rt:"
            retstr += str(self.oc_rt)[0:4] + ' ' + str(self.oc_rt)[-3:]
            retstr += "\nRa:"
            retstr += str(self.oc_ra)[0:4] + ' ' + str(self.oc_ra)[-3:]
            retstr += "\nRm:"
            retstr += str(self.oc_rm)[0:4] + ' ' + str(self.oc_rm)[-3:]
            retstr += "\nCm:"
            retstr += str(self.oc_cm)[0:4] + ' ' + str(self.oc_cm)[-3:]
            retstr += '\n'
        if self.str_items['incell']:
            retstr += "------incell parameters-----\n"
            retstr += "Rt:"
            retstr += str(self.ic_rt)[0:4] + ' ' + str(self.ic_rt)[-3:]
            retstr += "\nRa:"
            retstr += str(self.ic_ra)[0:4] + ' ' + str(self.ic_ra)[-3:]
            retstr += "\nRm:"
            retstr += str(self.ic_rm)[0:4] + ' ' + str(self.ic_rm)[-3:]
            retstr += "\nCm:"
            retstr += str(self.ic_cm)[0:4] + ' ' + str(self.ic_cm)[-3:]
            retstr += '\n'
        if self.str_items['wholecell']:
            retstr += "------wholecell parameters-----\n"
            retstr += "          Rt:"
            retstr += str(self.wc_rt)[0:4] + ' ' + str(self.wc_rt)[-3:]
            retstr += "\n          Ra:"
            retstr += str(self.wc_ra)[0:4] + ' ' + str(self.wc_ra)[-3:]
            retstr += "\n          Rm:"
            retstr += str(self.wc_rm)[0:4] + ' ' + str(self.wc_rm)[-3:]
            retstr += "\n          Cm:"
            retstr += str(self.wc_cm)[0:4] + ' ' + str(self.wc_cm)[-3:]
            retstr += '\n'
        return retstr
        
    def oncell_str(self,Val):
        self.str_items['oncell'] = Val
        
    def incell_str(self,Val):
        self.str_items['incell'] = Val
        
    def wholecell_str(self,Val):
        self.str_items['wholecell'] = Val
    
    def __repr__(self):
        retstr = ''
        retstr += "-------------------- \nTraces Loaded: \n-------------------- \n"
        retstr += self.keys()
        
    def plot(self):
        """Wrapper for panel_plot. Allow more options for diferent styles of 
        plots in future versions ... i.e. plot_fits."""
        self.panel_plot([0.1,0.1,0.8,0.85])
    
    def panel_plot(self,bounds,new_panel = None):
        """ Nice little summary of the cap data - includes the positive going 
        cap current trace (oncell,incell and whole cell) with the fits 
        overlayed. Future versions should include a key and print 
        the whole cell parameters."""
        
        if not(new_panel == None):
            panel = new_panel
        else:
            panel = self.plot_panel
            
        ax1_rect = [0.0,0.2,1,.8]
        ax2_rect = [0.0,0.0,1,0.15]

        ax1_rect = panel.rect_remap(ax1_rect,bounds = bounds)
        ax2_rect = panel.rect_remap(ax2_rect,bounds = bounds)
        
        ax1 = axes(ax1_rect)
        ax1.set_color_cycle(['m','m','c','c','g','g'])
        alpha_list = [0.2,0.2,0.7] #alpha values for the display fit curves
        plotrange = [0.00102, 0.00125] #time range to plot the fits
        
        events = [self[k] for k in self.event_keys]
        line_labels = self.event_keys
        for e,f,a,l in zip(events[0:3],self.fit_list[0:3],alpha_list,line_labels):
            dat = e.ts(*plotrange,tunits = 's') #a slice (plotrange) of cap traces data.
            xpnts = e.ts(*self.pfitrange,tunits = 's').get_t_array()
            #xpnts -= xpnts[2]
            fx = f.fx(xpnts)# The fit Points
            dat.plot(marker = 'o',label = l) #plot the sampled data points.
            start = float(self.pfitrange[0] - plotrange[0])#-xpnts[1]
            stop = float(self.pfitrange[1] - plotrange[0]) 
            xpnts = quan.Quantity(xpnts,e.dt.units)
            xpnts += quan.Quantity(start,'s')
            plot(xpnts.rescale('s'),fx, lw = 6, alpha = a) #plot the fit curves
            

        xlbs = ax1.get_xticklabels()
        [l.set_visible(False) for l in xlbs]
        ylabel(unicode("Im (" + e.y.dimensionality.unicode + ")",'utf-8'))
        legend()
        
        ax1.axis('tight')
        
        figtext(0.3,0.7,unicode(str(self),'utf-8'))
        #legend(("oncell_dat","_nolegend_","incell_data","_nolegend_","sub_data","_nolegend_"))
        ax2 = axes(ax2_rect, sharex = ax1)
        self['prot'].ts(*plotrange,tunits = 's').plot(color = 'k')
        
        ylbs = ax2.get_yticklabels()
        [l.set_visible(False) for l in ylbs]
        [l.set_visible(True) for l in ylbs[::3]]
        
        ylabel(unicode("Vm (" + self['prot'].y.dimensionality.unicode + ")"))
        xlabel(unicode("Time (" + dat.tunits_plot + ")"))
        ax2.axis('tight')
        kill_spines(ax1)
        kill_spines(ax2)
        

    def plot_fits(self):
        """a three pannel subplot showing the data with the 
        fitted data and the single exponential fits."""
        ax = subplot(3,1,1)
        max_bound = [0,0]
        axlist = list()
        for i,f in enumerate(self.fit_list[:-3]):
            a = subplot(3,1,i+1, sharey = ax)
            axlist.append(a)
            f.plot()
            psl.reduce_xticks(3)
            axis('tight')
            current_ybound = a.get_ylim()
            current_xbound = a.get_xlim()
            if (current_ybound[1] - current_ybound[0]) > (max_bound[1] - max_bound[0]):
                max_bound = list(current_ybound)
        ax.set_ybound(max_bound)
        xtext = (current_xbound[1]-current_xbound[0])*0.5+current_xbound[0]
        ytext = (max_bound[1]-max_bound[0])*0.5+max_bound[0]

        for i,f in enumerate(self.fit_list[:-3]):
            axes(axlist[i])
            text(xtext,ytext,str(f),verticalalignment = 'center', horizontalalignment = 'center')

        figure()
        ax = subplot(3,1,1)
        axlist = list()
        for i,f in enumerate(self.fit_list[3:]):
            a = subplot(3,1,i+1, sharey = ax)
            axlist.append(a)
            f.plot()
            psl.reduce_xticks(3)
            axis('tight')
            current_ybound = a.get_ylim()
            current_xbound = a.get_xlim()
            if (abs(current_ybound[1] - current_ybound[0])) > abs((max_bound[1] - max_bound[0])):
                max_bound = list(current_ybound)
        ax.set_ybound(max_bound)

        for i,f in enumerate(self.fit_list[3:]):
            axes(axlist[i])
            text(xtext,ytext*-1,str(f),verticalalignment = 'center', horizontalalignment = 'center')

    def plot_summary(self):
        range = [0.001,0.0015]
        self.event_list[0].ts(*range).plot()
        self.event_list[0].ts(*range).plot()
        self.event_list[0].ts(*range).plot()

    def calc_capacitance(self):
        """does the work of calculating the capacitance for all three traces"""
        #trange_pos = [0.00106, 0.0014] #time range to use to calculate capacitance 
        trange_pos = [0.00106, 0.00116]
        self.pfitrange = trange_pos #set the fit range
        self.fit_list = [] #initialize a empty list to hold the fits
        for k in self.event_keys: #itt through all the traces except the last one
            i = self[k] # get the trace by key
            y = i.ts(*trange_pos,tunits = 's') # get the y data in the time slice 
            x = i.ts(*trange_pos,tunits = 's').get_t_array() #get the x data in the time slice
            y = array(y.y) #convert the y data into an array for the fits
            self.fit_list.append(self.__make_fits__(x,y)) # make the fit using x and y array
        trange_neg = [0.00006, 0.00016]#same as above but for the negative going captran
        for k in self.event_keys:
            i = self[k]
            y = i.ts(*trange_neg,tunits = 's')
            x = i.ts(*trange_neg,tunits = 's').get_t_array()    
            y = array(y.y)
            self.fit_list.append(self.__make_fits__(x,y))
            
        #voltage cmnd during negative going step
        v_neg = mean(self['prot'].ts(*trange_neg,tunits = 's'))
        #voltage cmnd during postitive going step
        v_pos = mean(self['prot'].ts(*trange_pos,tunits = 's'))
        #calculate delta v
        delta_v = v_pos - v_neg
        
        #Make bool array with status of V step
        v_up = np.array(self['prot']) > -45
        
        [oc_rt,oc_ra,oc_rm,oc_qt,oc_cm] = self.__p__(0,3,'oncell_cap',delta_v,v_up)
        [ic_rt,ic_ra,ic_rm,ic_qt,ic_cm] = self.__p__(1,4,'incell_cap',delta_v,v_up)
        [wc_rt,wc_ra,wc_rm,wc_qt,wc_cm] = self.__p__(2,5,'sub_cap',delta_v,v_up)
        
#------------------set units--------------------------#
        v_units = self['prot'].y.units
        i_units = self['incell_cap'].y.units
        r_units = v_units / i_units
        r_units.units = quan.ohm
        q_units = self['incell_cap'].dt.units * self['incell_cap'].y.units
        q_units.units = quan.coulomb
        c_units = q_units / v_units
        c_units.units = quan.farad

        giga_ohm_s = u"G\u03A9"
        mega_ohm_s = u"M\u03A9"
        gos = giga_ohm_s.encode("utf-8")
        mos = mega_ohm_s.encode("utf-8")
        GOhm = quan.UnitQuantity('gigaohm', quan.ohm*1e9, symbol= gos)
        MOhm = quan.UnitQuantity('gigaohm', quan.ohm*1e6, symbol= mos)
        pF = quan.UnitQuantity('picofarad', quan.farad*1e-12, symbol= "pF")
        
        self.oc_ra = quan.Quantity(oc_ra,r_units.units)* r_units.magnitude
        self.oc_rt = quan.Quantity(oc_rt,r_units.units)* r_units.magnitude
        self.oc_rm = quan.Quantity(oc_rm,r_units.units)* r_units.magnitude
        self.oc_qt = quan.Quantity(oc_qt,q_units.units)* q_units.magnitude
        self.oc_cm = quan.Quantity(oc_cm,c_units.units)* c_units.magnitude

        self.oc_ra.units = MOhm
        self.oc_rt.units = GOhm
        self.oc_rm.units = GOhm
        self.oc_cm.units = pF
        
        self.ic_ra = quan.Quantity(ic_ra,r_units.units)* r_units.magnitude
        self.ic_rt = quan.Quantity(ic_rt,r_units.units)* r_units.magnitude
        self.ic_rm = quan.Quantity(ic_rm,r_units.units)* r_units.magnitude
        self.ic_qt = quan.Quantity(ic_qt,q_units.units)* q_units.magnitude
        self.ic_cm = quan.Quantity(ic_cm,c_units.units)* c_units.magnitude

        self.ic_ra.units = MOhm
        self.ic_rt.units = GOhm
        self.ic_rm.units = GOhm
        self.ic_cm.units = pF

        self.wc_ra = quan.Quantity(wc_ra,r_units.units)* r_units.magnitude
        self.wc_rt = quan.Quantity(wc_rt,r_units.units)* r_units.magnitude
        self.wc_rm = quan.Quantity(wc_rm,r_units.units)* r_units.magnitude
        self.wc_qt = quan.Quantity(wc_qt,q_units.units)* q_units.magnitude
        self.wc_cm = quan.Quantity(wc_cm,c_units.units)* c_units.magnitude

        self.wc_ra.units = MOhm
        self.wc_rt.units = GOhm
        self.wc_rm.units = GOhm
        self.wc_cm.units = pF

    def __p__(self,pos_fit_ind,neg_fit_ind,trace_ind,delta_v,v_up):
        """calculate whole cell parameters and return them in a list. Pass the 
        positive and negative fit objects from the capacitive transients"""
    
        #whole cell change in current
        i1 = self.fit_list[pos_fit_ind].prm_byname('i_sted')
        i2 = self.fit_list[neg_fit_ind].prm_byname('i_sted')
        delta_i = i1 - i2
        #print delta_i
        ##use current during the 
        trace = self[trace_ind]
        i1 = np.mean(trace.ts(0.00205,0.00215,tunits = 's'))
        i2 = np.mean(trace.ts(0.00090,0.00100,tunits = 's'))
        
        
        delta_i = i1 - i2
        #print delta_i
        #whole cell tau
        tau = self.fit_list[pos_fit_ind].prm_byname('tau')

        #calculate wc charge
        intrange = [0.00100, 0.0014]
        #array to integrate to calculate charge transfer during 
        
        integ_array = TimeSeries(y = np.where(v_up,
                                              np.array(self[trace_ind]) - i1,
                                              v_up),
                                 dt = self[trace_ind].dt)
        
        q1 = integrate.trapz(np.array(integ_array.ts(*intrange,tunits = 's')),dx=float(integ_array.dt))        
        q2 = delta_i * tau
        rt = delta_v / delta_i 
        qt = q1 + q2
        
        ra = tau * delta_v / qt
        rm = rt - ra 
        cm = qt * rt / (delta_v * rm)
        
        return [rt,ra,rm,qt,cm]

    def __make_fits__(self,x,y):
        """does the work of fitting the exponential model to the data. 
        Gets passed an array of time(x) and current(y) and appends the 
        FitModel to the CapRecord object. Future versions should search for the 
        greatest positive and negative going point and set the fit range 
        accordingly. Perhaps this should be the work of the calling function"""

        #get the initial guess by linear regression on log transformed data
        lg,offset = psl.lntransform(y) 
        m,b = polyfit(x,lg,1)

        ##----reduced model= single exponential----###
        #set perameters
        i_not = psl.Parameter(e**b + offset,"i_not")
        tau = psl.Parameter(abs(1/m),"tau")
        i_sted = psl.Parameter(offset,"i_sted")
        #model
        #ired = self.single_exp#lambda plist,x: plist[0]*e**(-1*x/plist[1])+ plist[2]
        #fit the model
        reduced = psl.FitModel([i_not,tau,i_sted],psl.single_exp,x,y,model_name = "single_exp")
        
        ##----full model= double exponential----###
        i_not_1 = psl.Parameter(e**b + offset,"i_not_1")
        tau_1 = psl.Parameter(abs(1/m),"tau_1")
        i_not_2 = psl.Parameter(0,"i_not_2")
        tau_2 = psl.Parameter(abs(tau_1()),"tau_2")
        i_sted = psl.Parameter(0,"i_sted")
        #model
        #ifull = self.double_exp#lambda plist,x: plist[0]*e**(-1*x/plist[1]) + plist[2]*e**(-1*x/plist[3]) + plist[4]
        #fit the model
        full = psl.FitModel([i_not_1,tau_1,i_not_2,tau_2,i_sted],psl.double_exp,x,y,model_name = "double_exp")
        
        #test models against each other for explained variance
        sig = psl.ftest_models(full,reduced)
        
        #set an atribute on the reduced model to indicate whether a double is
        #significant or not.
        if sig['pval'] > 0.05:
            reduced.double_sig = False
        else:
            reduced.double_sig = True
        #return just the reduced model
        return reduced
    
class EventException(Exception):
    def __init__(self, value):
        self.parameter = value
    def __str__(self):
        return repr(self.parameter)

def kill_spines(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

def load_signals_atf(hdr,signal_str = '', type_str = ''):
    
    f_data = loadtxt(hdr.file_path, skiprows=11)
    path = hdr.file_path

    indices = hdr.get_signal_cols(signal_str)
    
    return_list = list()
    
    for i in indices:
        tcol = hdr.col_num('Time')
        ycol = i
        ## get units
        tunits = hdr.sgnl_units('Time')
        yunits = hdr.sgnl_units(hdr.signal_list[i])
        
        return_list.append(TimeSeries(t = f_data[:,tcol],y = f_data[:,ycol], event_type = type_str, tunits = tunits,yunits = yunits, event_name = 'signal_' + str(i)))
    
    return return_list


###------------------------test code---------------------------###
def main():
    cell = CeNeuron(915)
    figure()
    cell.cap_records.plot()
    a = cell.testarr
    h5file = pyt.openFile("trial.h5", mode = 'w')
    group = h5file.root
    group2 = h5file.createGroup(h5file.root,"ivgroup")
    figure()
    a.plot()
    ar = a.hdf5_write(group)
    a.hdf5_load(ar)
    figure()
    a.plot()
    h5file.close()

#figure()
#cell.cap_records.plot_fits()
    print(cell.cap_records)
    Show()
    
if __name__ == '__main__':
    main()


