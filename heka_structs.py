from collections import namedtuple
import struct
#import data_structs
import copy
Tree = namedtuple('Tree','contents children')

class HEKA_type(object):
    
    def __init__(self,size):
        self.size = size
        
    def read(self,read_file,debug = False):
        rstr = read_file.read(self.size)
        if debug:
            print repr(rstr)
        self.parse(rstr)
        
    def parse(self,rstr):
        self.rstr = rstr
        self.data = struct.unpack(self.frmt,self.rstr)[0]
        
    def __getitem__(self,key):
        return self.data[key]
        
    def __int__(self):
        return int(self.data)
        
    def __repr__(self):
        return self.data.__repr__()

    def __str__(self):
        return str(self.data or 'None')
        
class INT(HEKA_type):
    def __init__(self,size):
        super(INT,self).__init__(size)
        if self.size == 4:
            self.frmt = 'L'
        elif self.size == 2:
            self.frmt = 'H'
        else:
            raise ValueError('INT size can only be 2 or 4 bytes')
            
    def __int__(self):
        return self.data
        
class StringType(HEKA_type):
    def __init__(self,size):
        super(StringType,self).__init__(size)
        self.frmt = str(self.size) +'s'
        
    def __str__(self):
        try:
            return str(self.data).strip(u'\x00') or 'None'
        except UnicodeDecodeError:
            return 'error'
        
class BYTE(HEKA_type):
    def __init__(self):
        super(BYTE,self).__init__(1)      
        self.frmt = 'B'

class BOOLEAN(HEKA_type):
    def __init__(self):
        super(BOOLEAN,self).__init__(1)  
        self.frmt = '?'
    
class CARD(HEKA_type):
    def __init__(self,size):
        super(CARD,self).__init__(size)
        if self.size == 4:
            self.frmt = 'L'
        elif self.size == 2:
            self.frmt = 'H'
        else:
            raise ValueError('CARD size can only be 2 or 4 bytes')
        
class LONGREAL(HEKA_type):
    def __init__(self):
        super(LONGREAL,self).__init__(8)
        self.frmt = 'd'
    
    def __float__(self):
        return self.data
        
class CHAR(HEKA_type):
    def __init__(self):
        super(CHAR,self).__init__(1)  
        self.frmt = 'c'

class NCHILD(HEKA_type):
    def __init__(self):
        super(NCHILD,self).__init__(4)
        self.frmt = 'L'

    def __int__(self):
        return self.data
        
class SET(HEKA_type):
    def __init__(self,size):
        super(SET,self).__init__(size)
        self.frmt = 'H'
        
class ARRAY(HEKA_type):
    def __init__(self,length,base_obj):
        self.length = length
        self.base_obj = base_obj
        self.size = length*self.base_obj.size
        self.base_class = type(base_obj)
    
    def parse(self,rstr):
        self.rstr = rstr
        self.data = [copy.copy(self.base_obj) for x in range(self.length)]
        for i in range(self.length):
            self.data[i].parse(rstr[i*self.base_obj.size:
                                    (i+1)*self.base_obj.size])
    
    def __getitem__(self,key):
        return self.data[key]
    
    def __str__(self):
        st = ''
        for c in self.data:
            st += str(c)
        return (st or 'None')

                                    
class BundleItem(HEKA_type):
    def __init__(self):
        super(BundleItem,self).__init__(16)
    
    def parse(self,rstr):
        self.data = {'oStart':INT(4),
                    'oLength' : INT(4),
                    'oExtension' : ARRAY(8,CHAR())}
        self.data['oStart'].parse(rstr[0:4])
        self.data['oLength'].parse(rstr[4:8])
        self.data['oExtension'].parse(rstr[8:16])
        
    def __getitem__(self,key):
        return self.data[key]
    
    def __getattr__(self,key):
        if not key == 'data':
            return self.data[key]
            
class UserParamDescrType(HEKA_type):
    def __init__(self):
        super(UserParamDescrType,self).__init__(40)
        
    def parse(self,rstr):
        self.data = {'name':StringType(32),
                    'unit':StringType(8)}
        self.data['name'].parse(rstr[0:32])
        self.data['unit'].parse(rstr[32:40])
        
    def __getitem__(self,key):
        return self.data[key]
    
    #def __getattr__(self,key):
    #    if not key == 'data':
    #        return self.data[key]
        
class HEKA_record(object):
    def load(self,b_file):
        #readindex = 0
        for key in self.readlist:
            self.__dict__[key].read(b_file)
            #self.__dict__[key] = self.data[key].data
            #readindex += self.data[key].size
            
    def repdict(self):
        return self.__dict__
            
    def __str__(self):
        rstr = ''
        for k in self.__dict__.keys():
            if not k == 'readlist':
                 rstr += (k + ':' + str(self.__dict__[k]) + '\n')
        return rstr

    def read_size(self):
        rsize = 0
        for key in self.readlist:
            rsize += self.__dict__[key].size
        return rsize

class BundleHeader(HEKA_record):
    def __init__(self,recordsize):
        self.readlist = ['oSignature', 'oVersion','oTime','oItems','fill',
                        'oIsLittleEndian','oBundleItems'] 
        self.oSignature             = ARRAY(8, CHAR())
        self.oVersion               = ARRAY(32,CHAR())
        self.oTime                  = LONGREAL()
        self.oItems                 = INT(4)
        self.oIsLittleEndian        = BOOLEAN()
        self.fill                   = ARRAY(11, CHAR())
        self.oBundleItems           = ARRAY(12, BundleItem())


####### Stimulus templates file (.pgf) ##########
   
class RootRecord_PGF(HEKA_record):
    def __init__(self,recordsize):
        self.recordsize = recordsize
        self.readlist = ['roVersion','roMark','roVersionName','roMaxSamples',
                        'roFiller1','roParams','roParamText','roReserved',
                        'roFiller2','roCRC','num_children']
        self.roVersion              = INT(4)
        self.roMark                 = INT(4)
        self.roVersionName          = StringType(32)
        self.roMaxSamples           = INT(4)
        self.roFiller1              = INT(4)
        self.roParams               = ARRAY(10,LONGREAL())
        self.roParamText            = ARRAY(10,StringType(32))
        self.roReserved             = ARRAY(128,BYTE())
        self.roFiller2              = INT(4)
        self.roCRC                  = CARD(4)
        self.num_children           = NCHILD()

                                   
class StimulationRecord(HEKA_record):
    def __init__(self,recordsize):
        self.recordsize = recordsize
        self.readlist = ['stMark','stEntryName','stFileName','stAnalName',
                        'stDataStartSegment','stDataStartTime',
                        'stSampleInterval','stSweepInterval','stLeakDelay',
                        'stFilterFactor','stNumberSweeps','stNumberLeaks',
                        'stNumberAverages','stActualAdcChannels',
                        'stActualDacChannels','stExtTrigger','stNoStartWait',
                        'stUseScanRates','stNoContAq','stHasSinewaves',
                        'stOldStartMacKind','stOldEndMacKind','stAutoRange',
                        'stBreakNext','stIsExpanded','stLeakCompMode',
                        'stHasChirp', 'stOldStartMacro','stOldEndMacro',
                        'sIsGapFree','sHandledExternally','stFiller1',
                        'stFiller2','stCRC','num_children']
        # removed
        #'stTag' removed
        
        self.stMark                 = INT(4)
        self.stEntryName            = StringType(32)
        self.stFileName             = StringType(32)
        self.stAnalName             = StringType(32)
        self.stDataStartSegment     = INT(4)
        self.stDataStartTime        = LONGREAL()
        self.stSampleInterval       = LONGREAL()
        self.stSweepInterval        = LONGREAL()
        self.stLeakDelay            = LONGREAL()
        self.stFilterFactor         = LONGREAL()
        self.stNumberSweeps         = INT(4)
        self.stNumberLeaks          = INT(4)
        self.stNumberAverages       = INT(4)
        self.stActualAdcChannels    = INT(4)
        self.stActualDacChannels    = INT(4)
        self.stExtTrigger           = BYTE()
        self.stNoStartWait          = BOOLEAN()
        self.stUseScanRates         = BOOLEAN()
        self.stNoContAq             = BOOLEAN()
        self.stHasSinewaves         = BOOLEAN()
        self.stOldStartMacKind      = CHAR()
        self.stOldEndMacKind        = BOOLEAN()
        self.stAutoRange            = BYTE()
        self.stBreakNext            = BOOLEAN()
        self.stIsExpanded           = BOOLEAN()
        self.stLeakCompMode         = BOOLEAN()
        self.stHasChirp             = BOOLEAN()
        self.stOldStartMacro        = StringType(32)
        self.stOldEndMacro          = StringType(32)
        self.sIsGapFree             = BOOLEAN()
        self.sHandledExternally     = BOOLEAN()
        self.stFiller1              = BOOLEAN()
        self.stFiller2              = BOOLEAN()
        self.stCRC                  = CARD(4)
        #self.stTag                  = StringType(32)
        self.num_children           = NCHILD()


        

class ChannelRecord(HEKA_record):
    def __init__(self,recordsize):
        self.recordsize = recordsize
        self.readlist = ['chMark','chLinkedChannel','chCompressionFactor',
                         'chYUnit','chAdcChannel','chAdcMode','chDoWrite',
                         'stLeakStore','chAmplMode','chOwnSegTime' ,
                         'chSetLastSegVmemb','chDacChannel','chDacMode',
                         'chFiller1','chRelevantXSegment','chRelevantYSegment',
                         'chDacUnit','chHolding','chLeakHolding','chLeakSize',
                         'chLeakHoldMode','chLeakAlternate','chAltLeakAveraging',
                         'chLeakPulseOn','chStimToDacID','chCompressionMode',
                         'chCompressionSkip','chDacBit','chHasSinewaves',
                         'chBreakMode','chZeroSeg','chFiller2','chInfoLReal',
                         'chInfoLInt','chInfoIChar','chDacOffset','chAdcOffset',
                         'chTraceMathFormat','chHasChirp','chFiller3',
                         'chCompressionOffset','chPhotoMode','chBreakLevel',
                         'chTraceMath','chOldCRC','chCRC',
                         'num_children']
        #removed 'chFiller4'
        self.chMark                 = INT(4)
        self.chLinkedChannel        = INT(4)
        self.chCompressionFactor    = INT(4)
        self.chYUnit                = StringType(8)
        self.chAdcChannel           = INT(2)
        self.chAdcMode              = BYTE()
        self.chDoWrite              = BOOLEAN()
        self.stLeakStore            = BYTE()
        self.chAmplMode             = BYTE()
        self.chOwnSegTime           = BOOLEAN()
        self.chSetLastSegVmemb      = BOOLEAN()
        self.chDacChannel           = INT(2)
        self.chDacMode              = BYTE()
        self.chFiller1              = BYTE()
        self.chRelevantXSegment     = INT(4)
        self.chRelevantYSegment     = INT(4)
        self.chDacUnit              = StringType(8)
        self.chHolding              = LONGREAL()
        self.chLeakHolding          = LONGREAL()
        self.chLeakSize             = LONGREAL()
        self.chLeakHoldMode         = BYTE()
        self.chLeakAlternate        = BOOLEAN()
        self.chAltLeakAveraging     = BOOLEAN()
        self.chLeakPulseOn          = BOOLEAN()
        self.chStimToDacID          = SET(2)
        self.chCompressionMode      = SET(2)
        self.chCompressionSkip      = INT(4)
        self.chDacBit               = INT(2)
        self.chHasSinewaves         = BOOLEAN()
        self.chBreakMode            = BYTE()
        self.chZeroSeg              = INT(4)
        self.chFiller2              = INT(4)
        self.chInfoLReal            = ARRAY(8, LONGREAL())
        self.chInfoLInt             = ARRAY(8, INT(4))
        self.chInfoIChar            = ARRAY(8, CHAR())
        self.chDacOffset            = LONGREAL()
        self.chAdcOffset            = LONGREAL()
        self.chTraceMathFormat      = BYTE()
        self.chHasChirp             = BOOLEAN()
        self.chFiller3              = ARRAY(30, CHAR())
        self.chCompressionOffset    = INT(4)
        self.chPhotoMode            = INT(4)
        self.chBreakLevel           = LONGREAL()
        self.chTraceMath            = StringType(128)
        self.chOldCRC               = CARD(4)
        #self.chFiller4              = INT(4)
        self.chCRC                  = CARD(4)
        self.num_children           = NCHILD()


class StimSegmentRecord(HEKA_record):
    def __init__(self,recordsize):
        self.recordsize = recordsize
        self.readlist = ['seMark','seClass','seDoStore','seVoltageIncMode',
                        'seDurationIncMode','seVoltage','seVoltageSource',
                        'seDeltaVFactor','seDeltaVIncrement','seDuration',
                        'seDurationSource','seDeltaTFactor','seDeltaTIncrement',
                        'seFiller1','seCRC','seScanRate','num_children']
        self.seMark                 = INT(4)
        self.seClass                = BYTE()
        self.seDoStore              = BOOLEAN()
        self.seVoltageIncMode       = BYTE()
        self.seDurationIncMode      = BYTE()
        self.seVoltage              = LONGREAL()
        self.seVoltageSource        = INT(4)
        self.seDeltaVFactor         = LONGREAL()
        self.seDeltaVIncrement      = LONGREAL()
        self.seDuration             = LONGREAL()
        self.seDurationSource       = INT(4)
        self.seDeltaTFactor         = LONGREAL()
        self.seDeltaTIncrement      = LONGREAL()
        self.seFiller1              = INT(4)
        self.seCRC                  = CARD(4)
        self.seScanRate             = LONGREAL()
        self.num_children           = NCHILD()
                                
####### PUL FILES ##########


class RootRecord_PUL(HEKA_record):
    def __init__(self,recordsize):
        self.recordsize = recordsize
        self.readlist = ['roVersion','roMark','roVersionName','roAuxFileName',
                        'roRootText','roStartTime','roMaxSamples','roCRC',
                        'roFeatures','roFiller1','roFiller2','num_children']
        self.roVersion              = INT(4)
        self.roMark                 = INT(4) 
        self.roVersionName          = StringType(32)
        self.roAuxFileName          = StringType(80)    
        self.roRootText             = StringType(400)
        self.roStartTime            = LONGREAL()
        self.roMaxSamples           = INT(4)
        self.roCRC                  = CARD(4)
        self.roFeatures             = SET(2)
        self.roFiller1              = INT(2)
        self.roFiller2              = INT(4)
        self.num_children           = NCHILD()
    
    def dictLabel(self):
        return str(self.roLabel)
                
    def repdict(self):
        return {'roVersion': self.roVersion}
        
class GroupRecord(HEKA_record):
    def __init__(self,recordsize):
        self.recordsize = recordsize
        self.readlist = ['grMark', 'grLabel','grText','grExperimentNumber',
                        'grGroupCount','grCRC','num_children']
        self.grMark                  = INT(4)
        self.grLabel                 = StringType(32)
        self.grText                  = StringType(80)
        self.grExperimentNumber      = INT(4)
        self.grGroupCount            = INT(4)
        self.grCRC                   = CARD(4)
        self.num_children            = NCHILD()
    
    def dictLabel(self):
        return str(self.grLabel)
    
    def repdict(self):
        return {'grExperimentNumber': self.grExperimentNumber}
        
    def __str__(self):
        rstr = ''
        for k in ['grLabel']:
            rstr += (k + ':' + str(self.__dict__[k]) + '\n')
        return rstr
   
class SeriesRecord(HEKA_record):
    def __init__(self,recordsize):
        self.recordsize = recordsize
        self.readlist = ['seMark', 'seLabel','seComment','seSeriesCount',
        'seNumberSweeps','seAmplStateOffset','seAmplStateSeries',
        'seSeriesType','seUseXStart','seFiller1','seFiller2','seTime',
        'sePageWidth','seSwUserParamDescr','seFiller3','seSeUserParams',
        'seLockInParams','seAmplifierState','seUsername','seUserParamDescr1',
        'seFiller4','seCRC','seSeUserParams2','seSeUserParamDescr2',
        'seScanParams','num_children']
        self.seMark                 = INT(4)
        self.seLabel                = StringType(32)
        self.seComment              = StringType(80)
        self.seSeriesCount          = INT(4)
        self.seNumberSweeps         = INT(4)
        self.seAmplStateOffset      = INT(4)
        self.seAmplStateSeries      = INT(4)
        self.seSeriesType           = CHAR()
        self.seUseXStart            = BOOLEAN()
        self.seFiller1              = CHAR()
        self.seFiller2              = CHAR()
        self.seTime                 = LONGREAL()
        self.sePageWidth            = LONGREAL()
        self.seSwUserParamDescr     = ARRAY(4, UserParamDescrType())
        self.seFiller3              = ARRAY(32,CHAR())
        self.seSeUserParams         = ARRAY(4, LONGREAL())
        self.seLockInParams         = ARRAY(96,BYTE()) # = 96, see "Pulsed.de" *)
        self.seAmplifierState       = StringType(400)
        self.seUsername             = StringType(80)
        self.seUserParamDescr1      = ARRAY(4, UserParamDescrType())
        self.seFiller4              = INT(4)
        self.seCRC                  = CARD(4)
        self.seSeUserParams2        = ARRAY(4,LONGREAL())
        self.seSeUserParamDescr2    = ARRAY(4,UserParamDescrType())
        self.seScanParams           = StringType(96)
        self.num_children           = NCHILD()
    
    def dictLabel(self):
        return str(self.seLabel)
            
    def repdict(self):
        return {'seNumberSweeps': self.seNumberSweeps}
    
    def __str__(self):
        rstr = ''
        for k in ['seLabel','seComment']:
            rstr += (k + ':' + str(self.__dict__[k]) + '\n')
        return rstr
    
class SweepRecord(HEKA_record):
    def __init__(self,recordsize):
        self.recordsize = recordsize
        self.readlist = ['swMark','swLabel','swAuxDataFileOffset',
                        'swStimCount','swSweepCount','swTime',
                        'swTimer','swSwUserParams','swTemperature',
                        'swOldIntSol','swOldExtSol','swDigitalIn',
                        'swSweepKind','swDigitalOut','swFiller1','swMarkers',
                        'swFiller2','swCRC','num_children']     
        self.swMark                 = INT(4)
        self.swLabel                = StringType(32)
        self.swAuxDataFileOffset    = INT(4)
        self.swStimCount            = INT(4)
        self.swSweepCount           = INT(4)
        self.swTime                 = LONGREAL()
        self.swTimer                = LONGREAL()
        self.swSwUserParams         = ARRAY(4,LONGREAL())
        self.swTemperature          = LONGREAL()
        self.swOldIntSol            = INT(4)
        self.swOldExtSol            = INT(4)
        self.swDigitalIn            = SET(2)
        self.swSweepKind            = SET(2)
        self.swDigitalOut           = SET(2)
        self.swFiller1              = INT(2)
        self.swMarkers              = ARRAY(4, LONGREAL())
        self.swFiller2              = INT(4)
        self.swCRC                  = CARD(4)
        self.num_children           = NCHILD()
    
    def dictLabel(self):
        return str("sweep")
        
    def repdict(self):
        return {'swStimCount': self.swStimCount}
            
    def __str__(self):
        rstr = ''
        for k in ['swLabel']:
            rstr += (k + ':' + str(self.__dict__[k]) + '\n')
        return rstr


class TraceRecord(HEKA_record):
    def __init__(self,recordsize):
        self.recordsize = recordsize
        self.readlist = ['trMark','trLabel','trTraceCount','trData',
        'trDataPoints','trInternalSolution','trAverageCount','trLeakCount',
        'trLeakTraces','trDataKind',
        'trUseXStart',
        'trFiller1','trRecordingMode','trAmplIndex',
        'trDataFormat','trDataAbscissa','trDataScaler','trTimeOffset','trZeroData',
        'trYUnit','trYMaxStr',
        'trXInterval','trXStart',
        'trXUnit','trXMaxStr',
        'trYRange','trYOffset','trBandwidth','trPipetteResistance','trCellPotential',
        'trSealResistance','trCSlow','trGSeries','trRsValue','trGLeak',
        'trMConductance','trLinkDAChannel','trValidYrange','trAdcMode',
        'trAdcChannel','trYmin','trYmax','trSourceChannel','trExternalSolution',
        'trCM','trGM','trPhase','trDataCRC','trCRC' ,'trGS',
        'trSelfChannel','trInterleaveSize','trInterleaveSkip','trImageIndex',
        'trMarkers','trSECM_X','trSECM_Y','trSECM_Z','num_children']
        
        self.trMark                  = INT(4)
        self.trLabel                 = StringType(32)
        self.trTraceCount            = INT(4)
        self.trData                  = INT(4)
        self.trDataPoints            = INT(4)
        self.trInternalSolution      = INT(4)
        self.trAverageCount          = INT(4)
        self.trLeakCount             = INT(4)
        self.trLeakTraces            = INT(4)
        self.trDataKind              = SET(2)
        self.trUseXStart             = BOOLEAN()
        self.trFiller1               = CHAR()
        self.trRecordingMode         = CHAR()
        self.trAmplIndex             = CHAR()
        self.trDataFormat            = BYTE()
        self.trDataAbscissa          = CHAR()
        self.trDataScaler            = LONGREAL()
        self.trTimeOffset            = LONGREAL()
        self.trZeroData              = LONGREAL()
        self.trYUnit                 = StringType(2) ### changed
        self.trYMaxStr               = StringType(6) ### added
        self.trXInterval             = LONGREAL()
        self.trXStart                = LONGREAL()
        self.trXUnit                 = StringType(2) ### changed
        self.trXMaxStr               = StringType(6) ### added
        self.trYRange                = LONGREAL()
        self.trYOffset               = LONGREAL()
        self.trBandwidth             = LONGREAL()
        self.trPipetteResistance     = LONGREAL()
        self.trCellPotential         = LONGREAL()
        self.trSealResistance        = LONGREAL()
        self.trCSlow                 = LONGREAL()
        self.trGSeries               = LONGREAL()
        self.trRsValue               = LONGREAL()
        self.trGLeak                 = LONGREAL()
        self.trMConductance          = LONGREAL()
        self.trLinkDAChannel         = INT(4)
        self.trValidYrange           = BOOLEAN()
        self.trAdcMode               = CHAR()
        self.trAdcChannel            = INT(2)
        self.trYmin                  = LONGREAL()
        self.trYmax                  = LONGREAL()
        self.trSourceChannel         = INT(4)
        self.trExternalSolution      = INT(4)
        self.trCM                    = LONGREAL()
        self.trGM                    = LONGREAL()
        self.trPhase                 = LONGREAL()
        self.trDataCRC               = CARD(4)
        self.trCRC                   = CARD(4)
        self.trGS                    = LONGREAL()
        self.trSelfChannel           = INT(4)
        self.trInterleaveSize        = INT(4)
        self.trInterleaveSkip        = INT(4)
        self.trImageIndex            = INT(4)
        self.trMarkers               = ARRAY(10, LONGREAL())
        self.trSECM_X                = LONGREAL()
        self.trSECM_Y                = LONGREAL()
        self.trSECM_Z                = LONGREAL()
        self.num_children            = NCHILD()
    
    def dictLabel(self):
        return str(self.trLabel)
        
    def repdict(self):
        return {'trData':self.trData,
                'trDataPoints':self.trDataPoints,
                'trDataKind':str(self.trDataKind),
                'trDataFormat':int(self.trDataFormat.data),
                'trDataScaler':str(self.trDataScaler),
                'trTimeOffset':str(self.trTimeOffset),
                'trYUnit':str(self.trYUnit),
                'trXInterval':str(self.trXInterval),
                'trXStart':str(self.trXStart),
                'trXUnit':str(self.trXUnit),
                }
        
    def __str__(self):
        rstr = ''
        for k in ['trData','trDataPoints']:
            rstr += (k + ':' + str(self.__dict__[k]) + '\n')
        return rstr