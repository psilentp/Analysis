from pylab import *
from scipy import optimize
from scipy import stats
from numpy import *
import copy as cp

#import copy_reg
#import pickle, marshal, types

def lntransform(x):
    offset = abs(min(x)) + 1
    return [log(x + offset), offset]

#
# register a pickle handler for code objects

#def fun_unpickler(data):
#    return marshal.loads(data)

#def code_pickler(code):
#    return code_unpickler, (marshal.dumps(code),)

#copy_reg.pickle(types.CodeType, code_pickler, code_unpickler)

def single_exp(plist,x):
    return plist[0]*e**(-1*x/plist[1])+ plist[2]
    
def double_exp(plist,x):
    return plist[0]*e**(-1*x/plist[1]) + plist[2]*e**(-1*x/plist[3]) + plist[4]

class FitModel:
    def __init__(self,params,function,x,y,model_name = ""):
        self.params = params
        self.function = function
        self.x = x
        self.y = y
        p0 = fit(self.params,self.function,self.x,self.y)
        self.paramVals = [i() for i in self.params]
        self.ss = sum((y-function(self.paramVals,x))**2)
        self.df = size(y) - size(params)
        self.model_name = model_name
        
    
    def __str__(self):
        buildstr = self.model_name + ": \n" 
        for i in self.params:
            buildstr = buildstr + str(i) + "\n"
        return buildstr

    def __repr__(self):
        return self.__str__()
    
    def prm_byname(self,param_name):
        l = [x.get_name() for x in self.params]
        i = l.index(param_name)
        return self.params[i]()

    def plot(self):
        plot(self.x,self.y,'bo')
        plot(self.x,self.function(self.paramVals,self.x),'red')

    def fx(self,tpoints):
        fx = self.function(self.paramVals,tpoints)
        return fx

def ftest_models(full,reduced):
    F_null = ((reduced.ss - full.ss)/(reduced.df - full.df))/(full.ss/full.df)
    pval = float(stats.f.sf(F_null, (reduced.df - full.df), full.df))
    return {'F_null':F_null,'pval':pval}

class Parameter:
    def __init__(self, value, name):
        self.value = value
        self.name = name

    def set(self, value):
        self.value = value

    def get_name(self):
        return self.name
    
    def __mul__(self,value):
        return self.value*value

    def __call__(self):
        return self.value

    def __str__(self):
        return self.name + "=" + str(self.value)


def fit(params,function,x, y):
    def errorfunct(params,x,y):
        return function(params,x) - y

    p = [parameter() for parameter in params]
    p1,success =  optimize.leastsq(errorfunct,p[:],args=(x,y))
    return [parameter.set(p1[i]) for i,parameter in enumerate(params)]

def reduce_xticks(factor):
    ax = gca()
    ticks = ax.get_xticks()
    ax.set_xticks(ticks[::factor])

def reduce_yticks(factor):
    ax = gca()
    ticks = ax.get_yticks()
    ax.set_yticks(ticks[::factor])
    
def unique_name(testname,existing_namelist, force_index_creation = False):
    if force_index_creation:
        uniquename = testname + "__0"
    else:
        uniquename = testname
        
    if (uniquename not in existing_namelist):
            return uniquename
    i = 0
    while(uniquename in existing_namelist):
        uniquename = testname + '__' + str(i)
        i += 1
    return uniquename
    
class Struct:
    def __init__ (self, *argv, **argd):
        if len(argd):
            # Update by dictionary
            self.__dict__.update (argd)
        else:
            # Update by position
            attrs = filter (lambda x: x[0:2] != "__", dir(self))
            for n in range(len(argv)):
                setattr(self, attrs[n], argv[n])
                
def seq_mean(seq):
    tmp = cp.copy(seq[0])
    for x in seq[1:]:
        tmp += x
    return tmp/len(seq)
        
def seq_sse(seq):
    s_mean = seq_mean(seq)
    tmp = (seq[0] - s_mean)
    tmp = tmp**2
    for x in seq[1:]:
        tmp += (x - s_mean)**2
    return tmp

def seq_std(seq):
    sse = seq_sse(seq)
    return (sse/len(seq))**0.5
    
def seq_sterr(seq):
    return seq_std/(len(seq)**0.5)
    
def seq_stats(seq):
    #get mean
    tmp = cp.copy(seq[0])
    for x in seq[1:]:
        tmp += x
    s_mean = tmp/len(seq)
    #get sse
    tmp = (seq[0] - s_mean)
    tmp = tmp**2
    for x in seq[1:]:
        tmp += (x - s_mean)**2
    s_sse = tmp
    #get std
    s_std = (s_sse/len(seq))**0.5
    #get stderr
    s_stderr = s_std/(len(seq)**0.5)
    #return as dict
    return{'mean':s_mean,'sse':s_sse,'std':s_std,'stderr':s_stderr}
    