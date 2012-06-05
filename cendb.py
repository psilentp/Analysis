#
#  cendb.py
#  
#
#  Created by Theodore Lindsay on 12/9/09.
#  Copyright (c) 2009 University of Oregon. All rights reserved.
#
""" Implements class to hold data in a cen zodb."""

import persistent
import transaction
import cenfunctions as cnf
import numpy as np
import abfloader as abl
import psilentplib as psl
import pylab as plb
import quantities as quan

class Cell(persistent.Persistent):
    def __init__(self,**kwargs):
        cdm = kwargs.pop('cdm',None)
        self.cennum = cdm['cennum']
        #code to load abf_index from cendir should sit here
        self.cdm = cdm
        #load and process synaptic dose response curve if it exists.
    
        if 'ic_o_fams' in self.cdm:
            self.presynaptic_ic = IclampDRC()
            self.presynaptic_ic.load(self.cdm)
            self.presynaptic_ic.calc_cell_stats()
        if 'post_ic_o_fams' in self.cdm:
            self.postsynaptic_ic = IclampDRC_Post()
            self.postsynaptic_ic.load(self.cdm)
            self.postsynaptic_ic.calc_cell_stats()
        if 'vc_o_fams' in self.cdm:
            self.presynaptic_vc = VclampDRC()
            self.presynaptic_vc.load(self.cdm)
            self.presynaptic_vc.calc_cell_stats()
        if 'on_cell_ifam' in self.cdm:
            self.fam = IFamCorrected(self.cdm)
            self.fam.load()
        if 'on_cell_captrans' in self.cdm:
            self.captrans = CapTrans(self.cdm)
            self.captrans.load()
            ra = self.captrans.cap.wc_ra
            ljp = quan.Quantity(15,'mV')
            self.fam.correct(ra,ljp)
            self.fam.corrected = self.fam.sub_fam.sr_ljp_correct(ra,ljp)
            self.iv = self.fam.corrected.get_iv()
            self.s_slope = cnf.StartSlope(self.fam.primary_fam)
            transaction.commit()
        if 'long_ifam' in self.cdm:
            self.long_ifam = LongIFam(self.cdm)
            self.long_ifam.load()
            try:
                ra = self.captrans.cap.wc_ra
            except AttributeError:
                ra = quan.Quantity(0.0,'ohm')
            ljp = quan.Quantity(15,'mV')
            self.long_ifam.correct(ra,ljp)
            self.long_iv = self.long_ifam.corrected.get_iv(0.3,0.4,tunits = 's')
            self.best_n = self.long_ifam.corrected.choose_iv_nslope(0.02,0.4,10)
            transaction.commit()
        if 'rev_pot' in self.cdm:
            self.rev_fam = IFamRP(self.cdm)
            self.rev_fam.load()
            self.rev_fam.correct(ra,ljp)
            #cell.rev_fam.stim_ifam.stimtrace.get_thresh_crossings(0.1,tunits = 's')
            #cell.rev_fam.stim_ifam.stims[0].get_thresh_crossings(-55,tunits = 's')
            stim_epoch = self.rev_fam.corrected.stimtrace.get_thresh_crossings(0.1,tunits = 's')
            stim_epoch[1] = stim_epoch[0]+0.05
            #base_epoch = self.rev_fam.stim_ifam.stims[0].get_thresh_crossings(-55,tunits = 's')
            base_epoch = [stim_epoch[0]-0.02,stim_epoch[0]]
            self.amp_curve = IsV(self.rev_fam.corrected.get_amp_curve(base_epoch=base_epoch,stim_epoch=stim_epoch,tunits ='s'))
            transaction.commit()
       
        if 'ic_bluepulse_long' in self.cdm:
            self.ic_bluepulse_long = BluePulseIC(self.cdm)
            self.ic_bluepulse_long.load()
            transaction.commit()
        
            
        if 'blue_pulse' in self.cdm:
            if not (self.cennum == 777):
                self.blue_pulse = BluePulse(self.cdm)
                self.blue_pulse.load()
                transaction.commit()
        if 'syn_o_fams' in self.cdm:
            self.synaptic_stim = SynVclampDRC()
            self.synaptic_stim.load(self.cdm)
            self.synaptic_stim.calc_cell_stats()
        
        if 'syn_o_fams_low' in self.cdm: ####
            self.synaptic_stim_low = SynVclampDRC_low()
            self.synaptic_stim_low.load(self.cdm)
            self.synaptic_stim_low.calc_cell_stats()
            
        if 'vc_bluepulse_long' in self.cdm:
            self.vc_bluepulse_long = BluePulseVC(self.cdm)
            self.vc_bluepulse_long.load()
            transaction.commit()
            
#base class for Experiments 
class Experiment(persistent.Persistent):
    def __init__(self,cdm):
        self.cdm = cdm

#Class def for Experiments
class CapTrans(Experiment):

    def __init__(self,cdm):
        self.cdm = cdm
        
    def __str__(self):
        return str(self.cap)
        
    def load(self):
        #parse the file names
        basedir = '/Volumes/UNTITLED/CENs/CEN' + \
                   str(self.cdm['cennum']) + '/'
        on_cell_file = basedir + self.cdm['on_cell_captrans']['files']
        in_cell_file = basedir + self.cdm['in_cell_captrans']['files']
        sig_key = self.cdm['on_cell_captrans']['sig_key']
        
        #load the raw data
        self.on_cell_cap = self.load_trans(on_cell_file,sig_key)
        self.in_cell_cap = self.load_trans(in_cell_file,sig_key)
        self.stim = self.load_stim(in_cell_file)
        
        #make the subtracted cap file
        self.sub_cap = self.in_cell_cap - self.on_cell_cap  
        #construct the CapRecord object
        self.cap = cnf.CapRecord(self.on_cell_cap,
                                 self.in_cell_cap,
                                 self.sub_cap,
                                 self.stim)
        
    def load_trans(self,filename,sig_key):
        #load the raw data
        print("loading file:" + filename)
        dat = abl.load_data(filename)
        #load the datetime info
        info = abl.load_info(filename)
        #get the time array
        tms = dat['time']
        #Hack to solve problem with old recordings
        un = 'pA'
        return cnf.TimeSeries(dt = tms[1],
                              y = dat['signals'][sig_key][0],
                              yunits = un,
                              tunits = 'us',
                              event_type = 'Captran_i',
                              event_name = 'incell_cap')
                        
    def load_stim(self,filename):
        #load the raw data
        print("loading protocol:" + filename)
        prots = abl.load_protocol(filename)
        tms = prots['time']
        return cnf.TimeSeries(dt = tms[1],
                              y = prots['channel']['Cmd 0']['analog'][0], 
                              yunits = prots['units'][0],
                              tunits = 'us')
    
    def plot(self):
        self.cap.plot()

class FamData(Experiment):
    
    def loadfam(self,filename,sig_key):
        #load the raw data
        print("loading file:" + filename)
        dat = abl.load_data(filename)
        #load the protocol
        prot = abl.load_protocol(filename)
        #load the datetime info
        info = abl.load_info(filename)
        #get the time array
        tms = dat['time']
        #create a list of sweeps
        sweeps = [cnf.TimeSeries(dt = tms[1],
                 y = d,
                 yunits = dat['units'][sig_key],
                 tunits = 'us')
                 for d in dat['signals'][sig_key]]
        #the list loads backwards so reverse it
        sweeps.reverse()
        try:
            prots = [cnf.TimeSeries(dt = tms[1],
                                y = p, 
                                yunits = prot['units'][0],
                                tunits = 'us') 
                                for p in prot['channel']['Vcmd']['analog']]
        except KeyError:
            prots = [cnf.TimeSeries(dt = tms[1],
                                y = p, 
                                yunits = prot['units'][0],
                                tunits = 'us') 
                                for p in prot['channel']['Cmd 0']['analog']]
        #make an o_fam
        if type(self) is IFamRP:
            stim = cnf.TimeSeries(dt = tms[1],
                                 y = prot['channel']['Vcmd']['digital'][0],
                                 yunits = prot['units'][0],
                                 tunits = 'us') 
            fam = cnf.StimCurrentFam(response_list = sweeps,cmnd_list = prots,stim = stim)
        else:
            fam = cnf.CurrentFamily(response_list = sweeps,cmnd_list = prots)
        #attach info attribute
        fam.info = info
        return fam
                                                
class IFamCorrected(FamData):
       
    def load(self):
        basedir = '/Volumes/UNTITLED/CENs/CEN' + \
                   str(self.cdm['cennum']) + '/'
        on_cell_file = basedir + self.cdm['on_cell_ifam']['files']
        in_cell_file = basedir + self.cdm['in_cell_ifam']['files']
        
        sig_key = self.cdm['on_cell_ifam']['sig_key']
        self.on_cell_fam = self.loadfam(on_cell_file,sig_key)
        self.in_cell_fam = self.loadfam(in_cell_file,sig_key)
        print "subtracting"
        self.sub_fam = self.in_cell_fam - self.on_cell_fam
        print "filtering"
        [x.low_filter(10000) for x in self.sub_fam]
        self.primary_fam = self.sub_fam
        transaction.commit()
        
    def correct(self,ra,ljp):
        self.corrected = self.sub_fam.sr_ljp_correct(ra,ljp)
        self.primary_fam = self.corrected
        transaction.commit()
        
    def plot(self,**kwargs):
        kwargs.update({'color':'k'})
        self.primary_fam.plot(**kwargs)
                
class LongIFam(FamData):
    
    def load(self):
        basedir = '/Volumes/UNTITLED/CENs/CEN' + \
                   str(self.cdm['cennum']) + '/'
        long_file = basedir + self.cdm['long_ifam']['files']
        sig_key = self.cdm['long_ifam']['sig_key']
        self.long_fam = self.loadfam(long_file,sig_key)
        self.primary_fam = self.long_fam
        transaction.commit()
    
    def plot(self,**kwargs):
        kwargs.update({'color':'k'})
        self.primary_fam.plot(**kwargs)
    
    def correct(self,ra,ljp):
        self.corrected = self.primary_fam.sr_ljp_correct(ra,ljp)
        transaction.commit()
        
class IFamRP(LongIFam):
    """class to hold Rev. pot. data"""
    def load(self):
        basedir = '/Volumes/UNTITLED/CENs/CEN' + \
                   str(self.cdm['cennum']) + '/'
        data_file = basedir + self.cdm['rev_pot']['files'][0]
        sig_key = self.cdm['on_cell_ifam']['sig_key']
        self.stim_ifam = self.loadfam(data_file,sig_key)
        self.primary_fam = self.stim_ifam
        try:
            data_file2 = basedir + self.cdm['rev_pot']['files'][1]
            sig_key = self.cdm['on_cell_ifam']['sig_key']
            self.stim_ifam2 = self.loadfam(data_file2,sig_key)
            self.primary_fam = self.stim_ifam2 + self.stim_ifam
            self.primary_fam /= 2.0
        except IndexError:
            print "only one rp protocol run"
            
        self.primary_fam = self.stim_ifam
        transaction.commit()
        
class IsV(Experiment):
    def __init__(self,IV):
        self.IV = IV
    
    def plot(self):
        self.IV.plot(mode = 'rp')
        
class OpticalDRC(Experiment):
    
    def __init__(self):
        Experiment.__init__(self,None)
        
    def load(self,cdm,cdmkey):
        self.cdm = cdm
        #parse the filenames
        basedir = '/Volumes/UNTITLED/CENs/CEN' + \
                   str(self.cdm['cennum']) + '/'
        filelist = [basedir + x for x in self.cdm[cdmkey]['files']]
        #make a persistent list for the trials
        self.trials = persistent.list.PersistentList()
        #load the ofams and insert into the list
        sig_key = self.cdm[cdmkey]['sig_key']
        for f in filelist:
            print("loading file:" + str(f))
            self.trials.append(self.load_trial(f,sig_key))
        transaction.commit()
        
    def plot(self):
        self.trials[0].plot()
        
class IclampDRC(OpticalDRC):
    from datafiles import powerlist
    
    def __init__(self):
        OpticalDRC.__init__(self)
        
    def load(self,cdm):
        OpticalDRC.load(self,cdm,'ic_o_fams')
        
    def load_trial(self,filename,sig_key):
        #load the raw data
        dat = abl.load_data(filename)
        #load the datetime info
        info = abl.load_info(filename)
        #get the time array
        tms = dat['time']
        #create a list of sweeps
        
        sweeps = [cnf.TimeSeries(dt = tms[1],y = x,yunits = 'mV', tunits = 'us') for x in\
                dat['signals'][sig_key]]
        #the list loads backwards so reverse it
        sweeps.reverse()
        #create a list of protocols
        prots = self.make_protocols()
        #convert the protocols to TimeSeries
        prots = [cnf.TimeSeries(dt = tms[1],y = x, yunits = 'mW/(mm*mm)',tunits = 'us') for x in prots]
        self.prots = prots
        #make an o_fam
        fam = cnf.OpticalFamily(response_list = sweeps,cmnd_list = prots)
        #attach info attribute
        fam.info = info
        return fam
        
    def make_protocols(self):
        """makes a protocol using the settings from 
        the pre-defined light powers"""
        #number of samps before the first epoch
        first_samps = 781.0
        #50ms at 50 KHz
        epoch1_samps = 2500
        #total collected in sweep
        total_samps = 50000
        #Number of samps after first sweep
        final_samps = total_samps - epoch1_samps - first_samps
        #make the prot. list and return it
        return [np.concatenate((np.zeros(first_samps)*0,
                np.ones(epoch1_samps)*x,
                np.zeros(final_samps))) for x in self.powerlist]
                
    def calc_cell_stats(self,*args):
        """calculate the mean data for the six trials - 
        average over trials first to last"""
        starttrial, stoptrial = 0,2
        try:
            starttrial = args[0]
            stoptrial = args[1]
        except IndexError:
            pass
            
        if starttrial > len(self.trials):
            starttrial = 0
        if starttrial > len(self.trials):
            stoptrial = len(self.trials)
            
        #a list to hold the stat data for each power
        statlist = list()
        #calculate the sweep stats for each power
        for i in range(len(self.powerlist)):
            #first make a sequence of the sweeps 
            seq = [x[i] for x in \
                   self.trials[starttrial:stoptrial]]
            #statlist gets with a dict of stat info for power i
            statlist.append(psl.seq_stats(seq))
        #make an o_fam for the various group stats
        self.ave_f = cnf.OpticalFamily(response_list = \
                                [x['mean'] for x in statlist],
                                cmnd_list = self.prots)
        self.sterr_f = cnf.OpticalFamily(response_list = \
                                [x['stderr'] for x in statlist],
                                cmnd_list = self.prots)                        
        transaction.commit()
        #make an object to hold the average Vm summary data
        self.calc_ave_vm()
    
    def calc_ave_vm(self):
        """create a persistent object that holds the charge transfer summary
        data"""
        self.vm_stats = persistent.mapping.PersistentMapping()
        for p in self.powerlist:
                seq = [self.get_mean_vm(i,p) for i in \
                       range(len(self.trials))]
                self.vm_stats[p] = psl.seq_stats(seq)
        transaction.commit()
                
    def get_mean_vm(self,trial,power,*args):
        """return the charge transfer for a sweep, 
        --called using trial and power"""
        #baseline is period before stim (15.6ms)
        baseline = [0,15.6]
        #allow an interval to be passed in args as start and stop. The interval
        #is defined in from the begining of the light-pulse. Defaults to 100ms.
        start, stop = baseline[1] + 0,baseline[1] + 100
        try:
            start = args[0] + baseline[1]
            stop = args[1] + baseline[1]
        except IndexError:
            pass
        tmp = self.trials[trial][self.powerlist.index(power)]
        return np.mean(tmp.ts(start,stop,tunits = 'ms'))
        
    def plot_means(self):
        mlist = [self.vm_stats[p]['mean'] for p in self.powerlist]
        plb.plot(self.powerlist,mlist,color = 'k',alpha = 0.3)
        
class IclampDRC_Post(IclampDRC):
    def __init__(self):
        IclampDRC.__init__(self)
    
    def load(self,cdm):
        OpticalDRC.load(self,cdm,'post_ic_o_fams')
        
    def plot(self,**kwargs):
        plottrial = self.trials[0]
        for i,x in enumerate(plottrial):
            plottrial[i] = x.sub_baseline(0,0.01,tunits = 's')
        plottrial.plot(**kwargs)
        
class VclampDRC(OpticalDRC):
    from datafiles import powerlist
    
    def __init__(self):
        OpticalDRC.__init__(self)
        
    def load(self,cdm):
        OpticalDRC.load(self,cdm,'vc_o_fams')
        
    def load_trial(self,filename,sig_key):
        #load the raw data
        dat = abl.load_data(filename)
        #load the datetime info
        info = abl.load_info(filename)
        #get the time array
        tms = dat['time']
        #create a list of sweeps
        sweeps = [cnf.TimeSeries(dt = tms[1],y = x,yunits = 'pA',tunits = 'us') for x in\
                dat['signals'][sig_key]]
        #the list loads backwards so reverse it
        sweeps.reverse()
        #create a list of protocols
        prots = self.make_protocols()
        #convert the protocols to TimeSeries
        prots = [cnf.TimeSeries(dt = tms[1],y = x, yunits = 'mW/(mm*mm)',tunits = 'us') for x in prots]
        self.prots = prots
        #make an o_fam
        fam = cnf.OpticalFamily(response_list = sweeps,cmnd_list = prots)
        #attach info attribute
        fam.info = info
        return fam
        
    def make_protocols(self):
        """makes a protocol using the settings from 
        the pre-defined light powers"""
        #number of samps before the first epoch
        first_samps = 781.0
        #50ms at 50 KHz
        epoch1_samps = 2500
        #total collected in sweep
        total_samps = 50000
        #Number of samps after first sweep
        final_samps = total_samps - epoch1_samps - first_samps
        #make the prot. list and return it
        return [np.concatenate((np.zeros(first_samps)*0,
                np.ones(epoch1_samps)*x,
                np.zeros(final_samps))) for x in self.powerlist]
                
    def calc_cell_stats(self,*args):
        """calculate the mean data for the six trials - 
        average over trials first to last"""
        starttrial, stoptrial = 0,2
        try:
            starttrial = args[0]
            stoptrial = args[1]
        except IndexError:
            pass
            
        if starttrial > len(self.trials):
            starttrial = 0
        if starttrial > len(self.trials):
            stoptrial = len(self.trials)
            
        #a list to hold the stat data for each power
        statlist = list()
        #calculate the sweep stats for each power
        for i in range(len(self.powerlist)):
            #first baseline subtract the sweeps 
            seq = [x[i].sub_baseline(0,15.6) for x in \
                   self.trials[starttrial:stoptrial]]
            #statlist gets with a dict of stat info for power i
            statlist.append(psl.seq_stats(seq))
        #make an o_fam for the various group stats
        self.ave_f = cnf.OpticalFamily(response_list = \
                                [x['mean'] for x in statlist],
                                cmnd_list = self.prots)
        self.sterr_f = cnf.OpticalFamily(response_list = \
                                [x['stderr'] for x in statlist],
                                cmnd_list = self.prots)                        
        transaction.commit()
        #make an object to hold the charge transfer summary data
        self.calc_cht_stats()
        
    def calc_cht_stats(self):
        """create a persistent object that holds the charge transfer summary
        data"""
        self.cht_stats = persistent.mapping.PersistentMapping()
        for p in self.powerlist:
                seq = [self.get_charge_transfer(i,p) for i in \
                       range(len(self.trials))]
                self.cht_stats[p] = psl.seq_stats(seq)
        transaction.commit()
                
    def get_charge_transfer(self,trial,power,*args):
        """return the charge transfer for a sweep, 
        --called using trial and power"""
        #baseline is period before stim (15.6ms)
        baseline = [0,15.6]
        #allow an interval to be passed in args as start and stop. The interval
        #is defined  from the begining of the light-pulse. Defaults to 100ms.
        start, stop = baseline[1] + 0,baseline[1] + 100
        try:
            start = args[0] + baseline[1]
            stop = args[1] + baseline[1]
        except IndexError:
            pass
        tmp = self.trials[trial][self.powerlist.index(power)]
        return tmp.sub_baseline(baseline[0],baseline[1],tunits = 'ms').integrate(start,stop,tunits = 'ms')
    
    def plot_means(self):
        mlist = [self.cht_stats[p]['mean'] for p in self.powerlist]
        plb.plot(self.powerlist,mlist,color = 'k', alpha = 0.3)
        
class SynVclampDRC(VclampDRC):
    def load(self,cdm):
        OpticalDRC.load(self,cdm,'syn_o_fams')
        
class SynVclampDRC_low(VclampDRC):
    def __init__(self):
        VclampDRC.__init__(self)
        self.powerlist = powerlist = [0.070000000000000021, 0.11914571614494369 , 0.20279573822416422, 0.34517490659800826, 0.58751587774119629, 1.0]
        
    def load(self,cdm):
        OpticalDRC.load(self,cdm,'syn_o_fams_low')
        
        
class BluePulseIC(Experiment):
    
    def load(self):
        basedir = '/Volumes/UNTITLED/CENs/CEN' + \
                   str(self.cdm['cennum']) + '/'
        files = [basedir + x for x in self.cdm['ic_bluepulse_long']['files']]
        self.currentlist = persistent.list.PersistentList()
        self.stimlist = persistent.list.PersistentList()
        self.tracelist = persistent.list.PersistentList()
        for fi in files:
            current,stim = self.load_stim(fi)
            stim*= 0.75
            trace = self.load_trace(fi,self.cdm['ic_bluepulse_long']['sig_key'])
            self.currentlist.append(current)
            self.stimlist.append(stim)
            self.tracelist.append(trace)
        
    def load_stim(self,filename):
        #load the raw data
        print("loading protocol:" + filename)
        prots = abl.load_protocol(filename)
        tms = prots['time']
        #return the protocl sweep
        current = cnf.TimeSeries(dt = tms[1],
                              y = prots['channel']['Icmd']['analog'][0], 
                              yunits = 'pA',
                              tunits = 'us')
                              
        stim = cnf.TimeSeries(dt = tms[1],
                              y = prots['channel']['Icmd']['digital'][0], 
                              yunits = 'mW/(mm*mm)',
                              tunits = 'us')
        return current,stim                      
        
    def load_trace(self,filename,sig_key):
        #load the raw data
        print("loading file:" + filename)
        dat = abl.load_data(filename)
        #load the datetime info
        info = abl.load_info(filename)
        #get the time array
        tms = dat['time']
        return cnf.TimeSeries(dt = tms[1],
                              y = dat['signals'][sig_key][0],
                              yunits = 'mV',
                              tunits = 'us')
    
    def plot(self):
        #self.imgplot()
        #return
        sp = cnf.SubPanel(bounds = [0.1,0.1,0.8,0.8])
        ax1_rect = [0.0,0.9,1,0.09]
        ax2_rect = [0.0,0.1,1,0.75]
        ax3_rect = [0.0,0.0,1,0.09]
        ax1 = plb.axes(sp.rect_remap(ax1_rect))
        ax2 = plb.axes(sp.rect_remap(ax2_rect),sharex = ax1)
        ax3 = plb.axes(sp.rect_remap(ax3_rect),sharex = ax1)
        plb.axes(ax1)
        for s in self.stimlist:
            s.plot()
        plb.gca().set_yscale('log')
        plb.axes(ax2)
        for t in self.tracelist:
            t.plot()
        plb.gca().axis('tight')
        plb.axes(ax3)
        for c in self.currentlist:
            c.plot()
        #plb.gca().axis('tight')
        plb.gca().set_xlim([30,36])
        
    def imgplot(self):
        #fig = plb.figure()
        #fig.subplots_adjust(top = 0.99, bottom = 0.01,left = 0.2,right = 0.99)
        nsweeps = len(self.tracelist) + 1
        plb.subplot(nsweeps,1,1)
        self.stimlist[0].plot()
        for i,s in enumerate(self.tracelist):
            ax = plb.subplot(nsweeps,1,i+2)
            s = s.get_normalized()
            s = np.vstack((s,s))
            plb.axis("off")
            plb.imshow(s,aspect = 'auto')

class BluePulseVC(Experiment):
    
    def load(self):
        basedir = '/Volumes/UNTITLED/CENs/CEN' + \
                   str(self.cdm['cennum']) + '/'
        files = [basedir + x for x in self.cdm['vc_bluepulse_long']['files']]
        self.currentlist = persistent.list.PersistentList()
        self.stimlist = persistent.list.PersistentList()
        self.tracelist = persistent.list.PersistentList()
        for fi in files:
            current,stim = self.load_stim(fi)
            stim*= 0.75
            trace = self.load_trace(fi,self.cdm['vc_bluepulse_long']['sig_key'])
            self.currentlist.append(current)
            self.stimlist.append(stim)
            self.tracelist.append(trace)
        
    def load_stim(self,filename):
        #load the raw data
        print("loading protocol:" + filename)
        prots = abl.load_protocol(filename)
        tms = prots['time']
        #return the protocl sweep
        current = cnf.TimeSeries(dt = tms[1],
                              y = prots['channel']['Vcmd']['analog'][0], 
                              yunits = 'pA',
                              tunits = 'us')
                              
        stim = cnf.TimeSeries(dt = tms[1],
                              y = prots['channel']['Vcmd']['digital'][0], 
                              yunits = 'mW/(mm*mm)',
                              tunits = 'us')
        return current,stim                      
        
    def load_trace(self,filename,sig_key):
        #load the raw data
        print("loading file:" + filename)
        dat = abl.load_data(filename)
        #load the datetime info
        info = abl.load_info(filename)
        #get the time array
        tms = dat['time']
        return cnf.TimeSeries(dt = tms[1],
                              y = dat['signals'][sig_key][0],
                              yunits = 'mV',
                              tunits = 'us')
    
    def plot(self):
        #self.imgplot()
        #return
        sp = cnf.SubPanel(bounds = [0.1,0.1,0.8,0.8])
        ax1_rect = [0.0,0.9,1,0.09]
        ax2_rect = [0.0,0.1,1,0.75]
        ax3_rect = [0.0,0.0,1,0.09]
        ax1 = plb.axes(sp.rect_remap(ax1_rect))
        ax2 = plb.axes(sp.rect_remap(ax2_rect),sharex = ax1)
        ax3 = plb.axes(sp.rect_remap(ax3_rect),sharex = ax1)
        plb.axes(ax1)
        for s in self.stimlist:
            s.plot()
        plb.gca().set_yscale('log')
        plb.axes(ax2)
        for t in self.tracelist:
            t.plot()
        plb.gca().axis('tight')
        plb.axes(ax3)
        for c in self.currentlist:
            c.plot()
        #plb.gca().axis('tight')
        plb.gca().set_xlim([30,36])
        
    def imgplot(self):
        #fig = plb.figure()
        #fig.subplots_adjust(top = 0.99, bottom = 0.01,left = 0.2,right = 0.99)
        nsweeps = len(self.tracelist) + 1
        plb.subplot(nsweeps,1,1)
        self.stimlist[0].plot()
        for i,s in enumerate(self.tracelist):
            ax = plb.subplot(nsweeps,1,i+2)
            s = s.get_normalized()
            s = np.vstack((s,s))
            plb.axis("off")
            plb.imshow(s,aspect = 'auto')


class BluePulse(Experiment):    
    def load(self):
        basedir = '/Volumes/UNTITLED/CENs/CEN' + \
                   str(self.cdm['cennum']) + '/'
        try:
            lowfiles = [basedir + x for x in self.cdm['blue_pulse']['bp_low']['files']]
            self.lowstim = self.load_stim(lowfiles[0])*0.03
            self.low = self.load_trace(lowfiles[0],self.cdm['blue_pulse']['bp_low']['sig_key'])
        except KeyError:
            self.lowstim = None
            self.low = None
        self.low_replicates = persistent.list.PersistentList()
        for x in range(1,len(lowfiles)): #just add the replicates, not the original trace
            self.low_replicates.append(self.load_trace(lowfiles[x],self.cdm['blue_pulse']['bp_low']['sig_key']))
        try:            
            medfiles = [basedir + x for x in self.cdm['blue_pulse']['bp_med']['files']]
            self.medstim = self.load_stim(medfiles[0])*0.75
            self.med = self.load_trace(medfiles[0],self.cdm['blue_pulse']['bp_med']['sig_key'])
        except KeyError:
            self.medstim = None
            self.med = None
        self.med_replicates = persistent.list.PersistentList()
        try:
            for x in range(1,len(lowfiles)): #just add the replicates, not the original trace
                self.med_replicates.append(self.load_trace(medfiles[x],self.cdm['blue_pulse']['bp_med']['sig_key']))
        except IndexError:
            print('no replicates')
        try:
            highfiles = [basedir + x for x in self.cdm['blue_pulse']['bp_high']['files']]
            self.highstim = self.load_stim(highfiles[0])*12.5
            self.high = self.load_trace(highfiles[0],self.cdm['blue_pulse']['bp_high']['sig_key'])
        except KeyError:
            self.highstim = None
            self.high = None
        self.high_replicates = persistent.list.PersistentList()
        try:
            for x in range(1,len(lowfiles)): #just add the replicates, not the original trace
                self.high_replicates.append(self.load_trace(highfiles[x],self.cdm['blue_pulse']['bp_high']['sig_key']))
        except IndexError:
            print('no replicates')
        
    def load_stim(self,filename):
        #load the raw data
        print("loading protocol:" + filename)
        prots = abl.load_protocol(filename)
        tms = prots['time']
        #return the protocl sweep
        print "HERE"
        return cnf.TimeSeries(dt = tms[1],
                              y = prots['channel']['Vcmd']['digital'][0], 
                              yunits = 'mW/(mm*mm)',
                              tunits = 'us')
                              
    def load_trace(self,filename,sig_key):
        #load the raw data
        print("loading file:" + filename)
        dat = abl.load_data(filename)
        #load the datetime info
        info = abl.load_info(filename)
        #get the time array
        tms = dat['time']
        print dat['signals'].keys()
        return cnf.TimeSeries(dt = tms[1],
                              y = dat['signals'][sig_key][0],
                              yunits = 'pA',
                              tunits = 'us')
                              
    def plot(self):
        sp = cnf.SubPanel(bounds = [0.1,0.1,0.8,0.8])
        ax1_rect = [0.0,0.9,1,0.09]
        ax2_rect = [0.0,0.0,1,0.85]
        ax1 = plb.axes(sp.rect_remap(ax1_rect))
        ax2 = plb.axes(sp.rect_remap(ax2_rect))
        plb.axes(ax1)
        for x in [self.lowstim,self.medstim,self.highstim]: 
            if not (x is None):
                x.plot()
        plb.gca().set_yscale('log')
        plb.axes(ax2)
        for x in [self.low,self.med,self.high]: 
            if not (x is None):
                x.plot()
        plb.gca().axis('tight')
                
    def animate_trial(self,trial = 'high',n=1000,speed = 0.0025,box = 14):
        import os, sys
        #speed = 0.1 # for testing low temporal rez
        files = []
        fig = plb.figure(figsize= (28,10))
        
        view_window = 14
        #box = ##########
        import time
        sp = cnf.SubPanel(bounds = [0.1,0.11,0.88,0.88])
        ax1_rect = [0.0,0.9,1,0.09]
        ax2_rect = [0.0,0.0,1,0.85]
        ax1 = plb.axes(sp.rect_remap(ax1_rect))
        ax2 = plb.axes(sp.rect_remap(ax2_rect))
        
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.set_yticks([0,12])
        #ax1.set_ylabel("Power\n(mW/mm$^2$)",rotation = 'horizontal', size = 40,multialignment = 'center')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.yaxis.set_ticks_position('left')
        ax2.xaxis.set_ticks_position('bottom')
        ax2.set_ylabel("Current(pA)",size = 40)
        ax2.set_xlabel("Time (s)",size = 40)
        
        import copy as cp
        dt = cp.copy(self.highstim.dt)
        dt.units = 's'
        tstart = time.time()
        y_stim = self.highstim.ts(0,box,tunits = 's')
        y_sig = self.high.ts(0,box,tunits = 's')
        x = y_stim.get_t_array(tunits = 's')
        pnts = len(x)
        
        plb.axes(ax1)
        padding = quan.Quantity(1,'mW/mm**2')
        stim_line, = plb.plot(x,np.ma.masked_where(x>0.1,np.array(y_stim.y)),lw = 5)
        ax1.set_ybound([self.highstim.y.min()-padding,self.highstim.y.max()+padding])
        ax1.set_xbound([0,view_window])
        [l.set_visible(False) for l in ax1.get_xticklines()]
        [l.set_visible(False) for l in ax1.get_xticklabels()]
        [lb.set_size(40) for lb in ax1.get_yticklabels()]
        
        plb.axes(ax2)
        sig_line, = plb.plot(x,np.ma.masked_where(x>0.1,np.array(y_sig.y)),color = 'k')
        ax2.set_ybound([self.high.y.min(),self.high.y.max()])
        ax2.set_xbound([0,view_window])
        [lb.set_size(40) for lb in ax2.get_yticklabels()]
        [lb.set_size(40) for lb in ax2.get_xticklabels()]

        frame = 0
        for i in np.arange(0.1,view_window,speed):
            stim_line.set_ydata(np.ma.masked_where(x>i,np.array(y_stim.y)))
            sig_line.set_ydata(np.ma.masked_where(x>i,np.array(y_sig.y)))
            plb.draw()
            fname = '_tmp%05d.png'%frame
            print 'Saving frame', fname
            fig.savefig(fname)
            files.append(fname)
            frame += 1
            
        move_distance = 0
        for i in np.arange(view_window,box,speed):
            move_distance += speed
            ax1.set_xbound([move_distance,move_distance + view_window])
            ax2.set_xbound([move_distance,move_distance + view_window])
            stim_line.set_ydata(np.ma.masked_where(x>i,np.array(y_stim.y)))
            sig_line.set_ydata(np.ma.masked_where(x>i,np.array(y_sig.y)))
            
            
            plb.draw()
            fname = '_tmp%05d.png'%frame
            print 'Saving frame', fname
            fig.savefig(fname)
            files.append(fname)
            frame += 1
        
        from animation_tools import make_movie
        make_move(files)
            
        """
        downs_factor = 8
        from PIL import Image
        im2 = Image.open(files[0])
        counter = 3
        frame = 0
        for f in files[1:]:
            counter += 1
            im = im2
            im2 = Image.open(f)
            im2 = Image.blend(im,im2,0.2)
            if counter % downs_factor == 0:
                fname = 'b_tmp%05d.png'%frame
                im2.save(fname)
                frame += 1
                print "saving" + fname
        """
        
        print 'Making movie animation.mpg - this make take a while'
        #os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800 -oac copy -o animation.mpg")
        os.system("ffmpeg -r 50 -qmax 4 -b 3200 -i b_tmp%05d.png test.mp4") #encode at 50fps
        #    plb.draw()
        #print 'FPS:' , (n+(len(x)/speed))/(time.time()-tstart)
        return files
        #plb.axes(ax2)
        #self.high.plot()
    
    