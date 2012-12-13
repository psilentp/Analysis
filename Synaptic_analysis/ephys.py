import pylab as plb
import numpy as np
import quantities as quan
import copy as cp


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

        if (np.size(t) == 0):
            t = arange(0,np.size(self.y)*dt, self.dt)
            event_duration = t[-1] - t[0]
        elif not(np.size(t) == np.size(y)):
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
        if np.size(key) > 1:
            return "slice"
        else:
            return "key"

class Trace(Event):
    """a class for generalized evenly sampled time series data"""
    def __init__(self,**kwargs):
    #parse the **kwargs
        tunits = kwargs.pop('tunits', 's')
        yunits = kwargs.pop('yunits', 'V')
        hdf5_key = kwargs.pop('hdf5_key', '')
        event_start = cp.copy(kwargs.pop('event_start', quan.Quantity(0.0,tunits))) 
        event_name = kwargs.pop('event_name','unnamed')
        event_type = kwargs.pop('event_type', 'time_series')
        event_duration = cp.copy(kwargs.pop('event_duration', quan.Quantity([0],tunits,dtype = 'float32')))
        t = cp.copy(kwargs.pop('t', quan.Quantity([],tunits,dtype = 'float32')))
        y = cp.copy(kwargs.pop('y', quan.Quantity([],yunits,dtype = 'float32'))) 
        dt = cp.copy(kwargs.pop('dt', quan.Quantity(1,tunits)))
                 
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
            ob = Trace(
                tunits =            cp.copy(self.tunits),
                yunits =            cp.copy(self.yunits),
                event_start =       cp.copy(self.event_start),
                event_name =        cp.copy(self.event_name),
                event_type =        cp.copy(self.event_type),
                event_duration =    cp.copy(self.event_duration),
                t =                 cp.copy(self.t),
                y =                 cp.copy(type(self.y).__getitem__(self.y,key)),
                dt =                cp.copy(self.dt*sample_stride))
            #ob.event_start = quan.Quantity(self.event_start + (self.dt * key.start),self.tunits)
            #ob.event_duration = self.event_duration - ob.event_start
            #ob.y = self.y.__getitem__(key)
            #ob.dt = cp.copy(self.dt*sample_stride)
            #print quan.version.__version__
            #ob.y = cp.copy(type(self.y).__getitem__(self.y,key))
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
            plb.plot(tpoints.rescale(tunits_plot),self.y,**kwargs)
            return
        if plot_type == 'errorbar':
            plb.errorbar(tpoints.rescale(tunits_plot),self.y,**kwargs)
            return
        
    def sub_baseline(self,start,stop,**kwargs):
    #baseline subtract the Trace, start and stop define the interval,
    #in object time units
        tunits = kwargs.pop('tunits',self.dt.units)
        timebase = kwargs.pop('timebase','local')
        #print 'here'
        #print timebase
        #print timebase == 'global'
        if timebase == 'global':
            #print 'global'
            start -= self.event_start.magnitude
            stop -= self.event_start.magnitude
        baseline = quan.Quantity(np.mean(self.ts(start,stop,tunits = tunits)),self.yunits)
        return (self - baseline)
        
    def integrate(self,*args,**kwargs):
    #Integrate the Trace.
        from scipy import integrate
        tunits = kwargs.pop('tunits',self.dt.units)
        start, stop = 0,float('inf')
        try:
            start = args[0]
            stop = args[1]
        except IndexError:
            pass
        units = self.y.units * self.dt.units
        return(quan.Quantity(integrate.trapz(self.ts(start,
                                                    stop,tunits = tunits).y,
                                                    dx = self.dt)
                            ,units))
    
    #def mean(self,*args):
    #Get the mean of the Trace.
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
        return (np.size(self.y)-1)*self.dt

    def get_t_array(self,timebase = 'local',**kwargs):
        tunits = kwargs.pop('tunits',self.dt.units)
    #return an array of t points corresponding to each y point.
        start = 0
        stop = float(self.get_duration().rescale(tunits))
        ar = np.linspace(start, stop, num = np.size(self.y))
        if(timebase == 'local'):
            return np.array(ar)
        if(timebase == 'global'):
            ar+=self.event_start.rescale(tunits)
            return np.array(ar)
        return ar
    
    def get_low_filter(self,Hz):
    #return a low-pass-filterd version of the time series
        import FiltFilt
        filter_order = 3
        sampfreq = 1/self.dt.rescale('sec')
        filtfreq = quan.Quantity(Hz,'Hz')
        passband_peram = float((filtfreq/sampfreq).simplified)
        from scipy.signal import butter
        [b,a]=butter(filter_order,passband_peram)
        ob = cp.copy(self)
        ob.y = quan.Quantity(FiltFilt.filtfilt(b,a,self.y),self.yunits)
        return ob
    
    def get_normalized(self):
    #normalize so that y goes from 0 to 1
        n = self.y-np.amin(self.y)
        n /= np.amax(n)
        return np.array(n)
        
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
        timebase = kwargs.pop('timebase','local')
        dtunits = self.dt.units
        ind = [None,None,None]
        if np.size(key) == 3:
            ind[0] = int(quan.Quantity(key[0],tunits).rescale(dtunits)/self.dt)
            ind[1] = int(quan.Quantity(key[1],tunits).rescale(dtunits)/self.dt)
            ind[2] = int(quan.Quantity(key[2],tunits).rescale(dtunits)/self.dt)
        elif np.size(key) == 2:
            ind[0] = int(quan.Quantity(key[0],tunits).rescale(dtunits)/self.dt)
            ind[1] = int(quan.Quantity(key[1],tunits).rescale(dtunits)/self.dt)
            ind[2] = 1
        elif np.size(key) == 1:
            ind[0] = int(quan.Quantity(key[0],tunits).rescale(dtunits)/self.dt)
            ind[1] = int(quan.Quantity(key[0],tunits).rescale(dtunits)/self.dt) + 1
            ind[2] = 1
        selection = slice(ind[0],ind[1],ind[2])
        ob = self.__getitem__(selection)
        if(timebase == 'global'):
            ob.event_start+=self.event_start.rescale(tunits) + quan.Quantity(key[0],tunits).rescale(dtunits)
            #return np.array(ar)
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
        
    def animate(self,n,speed,box):
        import time
        ion()
        #speed = 100
        tstart = time.time()
        y = self.ts(0,box,tunits = 's')
        x = y.get_t_array()
        pnts = len(x)
        line, = plb.plot(x,np.ma.masked_where(x>0,np.array(y.y)))
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