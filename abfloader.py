import numpy as np
import scipy as sp
from pylab import *
import gc
import axon
import datetime
# note bug somewhere in program - when testing library axon_get_protocol behaves differenly if axon_get_data has run already

def load_info(file_name):

    abf = axon.ABF()
    abf.Open(file_name)
    abf.ReadFileInfo()
    if abf.FileInfo.uFileSignature != 843465281:
        print 'ABF2 signature does not match!'
        abf.Close()
        return 
        
    date = abf.FileInfo.uFileStartDate
    y = int(str(date)[:4])
    m = int(str(date)[4:6])
    d = int(str(date)[6:8])
    ms = abf.FileInfo.uFileStartTimeMS
    filedate = datetime.datetime(y,m,d)
    filedate += datetime.timedelta(0,0,0,ms)
    return {'writetime':filedate}
    
def load_data(file_name):
    """
    'axon_get_data' is a helper function for Axon Library. 
    'axon_get_data(file_name)' gets (reads) raw data from ABF2 file 'file_name'
    and returns a tuple '(data, times)'. The real data is returned as an array 
    'data(x,y,z)', where 'x' is a number refering to ADC, 'y' is a number of sample and 
    'z' is a number of sweep. 'times' is an array with time samples corresponding 
    with the data. 
    """
    
    abf = axon.ABF()
    abf.Open(file_name)
    #abf.ReadFileInfo()
    abf.ReadAllSections()
    print abf.FileInfo.uFileSignature
    if abf.FileInfo.uFileSignature != 843465281:
        print 'ABF2 signature does not match!'
        abf.Close()
        return 
        
    
    # Reading raw data
    if abf.FileInfo.nDataFormat == 0:
        raw_data = np.array(abf.ReadIntData(),float);
    elif abf.FileInfo.nDataFormat == 1:
        raw_data = np.array(abf.ReadFloatData(),float);
    else:
        print 'Wrong data format!'
        abf.Close()
    
    
    #abf.ReadProtocolInfo()
    #abf.ReadADCInfo()
          
    nADC = abf.FileInfo.ADCSection.llNumEntries;       # Number of ADC channels
    ns = abf.ProtocolInfo.lNumSamplesPerEpisode;       # Number of samples per episode
    ne = abf.FileInfo.uActualEpisodes;                 # Number of episodes per run
    T  = abf.ProtocolInfo.fADCSequenceInterval;        # Sampling period in ??
    
    # Reshaping data and reading ADCInfo into list
    data = np.zeros((nADC,ns/nADC,ne))
    adc_info = []
    for i in range(nADC):
        data[i,:,:] = np.reshape(raw_data[i::nADC],(ns/nADC,ne),order='F');
        adc_info.append(abf.GetADCInfo(i));    
    
    if abf.FileInfo.nDataFormat == 0:
        for i in range(nADC):
            #Data offset
            offset = adc_info[i].fInstrumentOffset - \
                     adc_info[i].fSignalOffset;
            # Data gain
            gain =  abf.ProtocolInfo.fADCRange          / \
                    abf.ProtocolInfo.lADCResolution     / \
                    (adc_info[i].fADCProgrammableGain   * \
                     adc_info[i].fInstrumentScaleFactor * \
                     adc_info[i].fSignalGain);                 
            # Additional telegraph gain    
            if adc_info[i].nTelegraphEnable:
                gain = gain/adc_info[i].fTelegraphAdditGain;        
            # Scaling of the data
            data[i,:,:] = data[i,:,:] * gain + offset;
    
    indexes ={'ADC_name':abf.ADCInfo.lADCChannelNameIndex,'ADC_units':abf.ADCInfo.lADCUnitsIndex,'ADC_signals':abf.FileInfo.ADCSection.llNumEntries}
    signals = [abf.GetString(x) for x in range(indexes['ADC_name'],indexes['ADC_name']+indexes['ADC_signals'] * 2,2)]
    units = [abf.GetString(x) for x in range(indexes['ADC_units'],indexes['ADC_units']+indexes['ADC_signals'] * 2,2)]
    abf.Close()
    
    # Creating times array
    times = np.linspace(0, ns/nADC*T, ns/nADC);
    
    data_dict = dict()
    
    for i,s in enumerate(signals):
        data_dict[s] = [x for x in rot90(data[i,:])]
    
    unit_dict = dict()
    
    for i,s in enumerate(signals):
        unit_dict[s] = units[i]
    
    #del(abf)
    
    return {'signals':data_dict,'units':unit_dict,'time':times}


def load_protocol(file_name):
    """
    'axon_get_protocol' is a helper function for Axon Library. 
    'axon_get_protocol(file_name)' gets (reads) protocol waveform from ABF2 
    file 'file_name' and returns a tuple '(pro,times,pro_file_name)'. The times 
    series of a protocol is returned as an array 'pro(x,y)' where 'x' is 
    a number of DAC and 'y' is a number of sample and 'z' is a number of sweep. 
    Also the exact times where the samples where taken are returned in 
    an array 'times(y)'. Additionally, the name of the protocol file 
    is returned in 'pro_file_name'.
    """
    abf = axon.ABF()
    abf.Open(file_name)
    abf.ReadFileInfo()
    print abf.FileInfo.uFileSignature
    if abf.FileInfo.uFileSignature != 843465281:
        print 'ABF2 signature does not match!'
        abf.Close()
        return 
    
    # Reading raw data
    if abf.FileInfo.nDataFormat == 0:
        raw_data = np.array(abf.ReadIntData(),float);
    elif abf.FileInfo.nDataFormat == 1:
        raw_data = np.array(abf.ReadFloatData(),float);
    else:
        print 'Wrong data format!'
        abf.Close()
    print("Data in")
    abf.ReadAllSections()
    abf.ReadProtocolInfo()
    abf.ReadEpochInfo()
    abf.ReadEpochInfoPerDAC() #Call this function to load the abf.EpochInfoPerDAC matrix
         
    nADC = abf.FileInfo.ADCSection.llNumEntries       # Number of ADC channels
    nDAC = abf.FileInfo.DACSection.llNumEntries       # Number of DAC channels
    ns = abf.ProtocolInfo.lNumSamplesPerEpisode       # Number of samples per episode
    ne = abf.FileInfo.uActualEpisodes                 # Number of episodes per run (number of sweeps)
    T  = abf.ProtocolInfo.fADCSequenceInterval ##*10**-3 # Sampling period in ??

    # Reading the protocol file name from the strings
    abf.ReadStrings()
    pro_file_name = abf.GetString(abf.FileInfo.uProtocolPathIndex);
    
    # Creating times array
    times = np.linspace(0, ns/nADC*T, ns/nADC);
    
    # First and last holding (samples)
    last = np.floor(ns/nADC/64);
    
    dac_signals = dict()
    
    indexes ={'DAC_name':abf.DACInfo.lDACChannelNameIndex,'DAC_units':abf.DACInfo.lDACChannelUnitsIndex,'DAC_signals':abf.FileInfo.DACSection.llNumEntries}
    DAC_names = [abf.GetString(x) for x in range(indexes['DAC_name'],indexes['DAC_name']+indexes['DAC_signals'] * 2,2)]
    units = [abf.GetString(x) for x in range(indexes['DAC_units'],indexes['DAC_units']+indexes['DAC_signals'] * 2,2)]
    print('got startup data')
    for d,n in enumerate(DAC_names):
    # Creating a protocol array consisted of DAC's holding levels only
        pro = [np.ones(last)*abf.GetDACInfo(d).fDACHoldingLevel for x in range(ne)]
        dig = [np.ones(last)*abf.ProtocolInfo.nDigitalHolding for x in range(ne)]
    # Reading the whole epoch_per_dac array
        for j in range(abf.FileInfo.EpochPerDACSection.llNumEntries):
            ei = abf.GetEpochInfo(j)
            ei_dac = abf.GetEpochInfoPerDAC(d,j) #first arg is the Analog Output channel, second argument is the epoch
            if not(ei_dac.nEpochNum == -1): 
                pro = [sp.concatenate((epoch,np.ones(i*ei_dac.lEpochDurationInc+ei_dac.lEpochInitDuration) * float(ei_dac.fEpochInitLevel + ei_dac.fEpochLevelInc*i))) for i,epoch in enumerate(pro)] 
                dig = [sp.concatenate((epoch, np.ones(i*ei_dac.lEpochDurationInc+ei_dac.lEpochInitDuration) * float(ei.nDigitalValue))) for i,epoch in enumerate(dig)]
    
    #fill the rest of the epoch 
        total_samps = ns/nADC
        pro = [sp.concatenate((epoch,np.ones(total_samps-size(epoch)) * abf.GetDACInfo(d).fDACHoldingLevel )) for epoch in pro] 
        dig = [sp.concatenate((epoch, np.ones(total_samps-size(epoch)) * abf.ProtocolInfo.nDigitalHolding )) for epoch in dig]
        
        dac_signals[n] = {'analog':pro,'digital':dig}

    indexes ={'DAC_name':abf.DACInfo.lDACChannelNameIndex,'DAC_units':abf.DACInfo.lDACChannelUnitsIndex,'DAC_signals':abf.FileInfo.DACSection.llNumEntries}
    signals = [abf.GetString(x) for x in range(indexes['DAC_name'],indexes['DAC_name']+indexes['DAC_signals'] * 2,2)]
    units = [abf.GetString(x) for x in range(indexes['DAC_units'],indexes['DAC_units']+indexes['DAC_signals'] * 2,2)]
    
    abf.Close()
    
    #del(abf)
    
    return ({'channel':dac_signals,'time':times,'units':units})
