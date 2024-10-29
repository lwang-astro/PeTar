# BSE data interface
from .base import *
from .functions import *

BSE_STAR_TYPE_INDEX={'LMS':0, 'MS':1, 'HG':2, 'GB':3, 'CHeB':4, 
                     'FAGB':5, 'SAGB':6, 'HeMS':7, 'HeHG':8, 'HeGB':9, 
                     'HeWD':10, 'COWD':11, 'ONWD':12, 'NS':13, 'BH':14, 'SN':15}

BSE_STAR_TYPE_NAME={0:'LMS', 1:'MS', 2:'HG', 3:'GB', 4:'CHeB', 
                    5:'FAGB', 6:'SAGB', 7:'HeMS', 8:'HeHG', 9:'HeGB', 
                    10:'HeWD', 11:'COWD', 12:'ONWD', 13:'NS', 14:'BH', 15:'SN'}

BSE_BINARY_TYPE_INDEX={'Unset':0, 'Initial':1, 'Type_change':2, 'Start_Roche':3, 
                       'End_Roche':4, 'Contact':5, 'Start_Symbiotic':6,
                       'End_Symbiotic':7, 'Common_envelope':8, 'Giant':9,
                       'Coalescence':10, 'Blue_straggler':11, 'No_remain':12, 'Disrupt':13}

BSE_BINARY_TYPE_INDEX={0:'Unset', 1:'Initial', 2:'Type_change', 3:'Start_Roche', 
                       4:'End_Roche', 5:'Contact', 6:'Start_Symbiotic',
                       7:'End_Symbiotic', 8:'Common_envelope', 9:'Giant',
                       10:'Coalescence', 11:'Blue_straggler', 12:'No_remain', 13:'Disrupt'}


class SSEStarParameter(DictNpArrayMix):
    """ SSE star parameter class from bse_interface.h
    Keys: (class members)
        type  (1D): SSE stellar type, see help(petar.SSEType)
        mass0 (1D): initial mass at each evolution stage (Msun)
        mass  (1D): current mass (Msun)
        rad   (1D): stellar radius (Rsun)
        mcore (1D): core mass (Msun)
        rcore (1D): core radius (Rsun)
        spin  (1D): stellar rotation
        epoch (1D): time offset at each evolution stage (Myr)
        time  (1D): current physical time (Myr)
        lum   (1D): bolometric luminosity (Lsun)
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['type',np.int64],['mass0',np.float64],['mass',np.float64],['rad',np.float64],['mcore',np.float64],['rcore',np.float64],['spin',np.float64],['epoch',np.float64],['time',np.float64],['lum',np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class SSETypeChange(DictNpArrayMix):
    """ SSE type change output data from PeTar
    Keys: (class members)
        id (1D): particle id
        init (SSEStarParameter): initial status of star
        final (SSEStarParameter): final status of star after stellar evolution
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['id',np.int64],['init',SSEStarParameter],['final',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def printTable(self, column_format = 'final', print_title = True):
        """ Print merger information in a formated table

        Parameters:
        ----------
        column_format: a list of column label (class member name), format and column title, enclosed by tuple, for sub-member, use . to access
                       For exmaple: [(key1,'%s',title1), (key2,'%12.7f',title2), (key3.subkey1,'%d',title3), (key3.subkey2,'%e',title4)]
                       Some pre-defined choices (default: final):
                       final: [('id', '%10d', 'id'),
                               ('final.time','%12.4g','timef[Myr]'),('final.type','%4d','kf'),
                               ('final.mass0','%10.3f','mass0[M*]'),('final.mass','%10.3f','mass[M*]'), 
                               ('final.rad','%8.2g','r[R*]'),('final.lum','%8.2g','Lf[L*]'),
                               ('final.mcore','%9.3f','mcf[M*]'),('final.rcore','%8.2g','rcf[R*]'),
                               ('final.spin','%8.2g','spinf'),('final.epoch','%11.3g','epoch[Myr]')]
                       init-final: [('id', '%10d', 'id'),
                               ('init.time','%12.4g','timei[Myr]'),('init.type','%4d','ki'),
                               ('init.mass','%11.4f','mi[M*]'),('init.mcore','%9.3f','mci[M*]'),
                               ('final.time','%12.4g','timef[Myr]'),('final.type','%4d','kf'),
                               ('final.mass','%11.4f','mf[M*]'),('final.mcore','%9.3f','mcf[M*]')]
        print_title: print title of keys (default: True)
        """
        
        if (column_format == 'final'):
            column_format = [('id', '%10d', 'id'),
                             ('final.time','%12.4g','timef[Myr]'),('final.type','%4d','kf'),
                             ('final.mass0','%10.3f','mass0[M*]'),('final.mass','%10.3f','mass[M*]'), 
                             ('final.rad','%8.2g','r[R*]'),('final.lum','%8.2g','Lf[L*]'),
                             ('final.mcore','%9.3f','mcf[M*]'),('final.rcore','%8.2g','rcf[R*]'),
                             ('final.spin','%8.2g','spinf'),('final.epoch','%11.3g','epoch[Myr]')]

        elif (column_format =='init-final'):
            column_format = [('id', '%10d', 'id'),
                             ('init.time','%12.4g','timei[Myr]'),('init.type','%4d','ki'),
                             ('init.mass','%11.4f','mi[M*]'),('init.mcore','%9.3f','mci[M*]'),
                             ('final.time','%12.4g','timef[Myr]'),('final.type','%4d','kf'),
                             ('final.mass','%11.4f','mf[M*]'),('final.mcore','%9.3f','mcf[M*]')]
        elif (type(column_format) != list):
            raise ValueError('Error, column_format should be "final","init-final" or list of manually defined column format; given ',column_format)

        DictNpArrayMix.printTable(self, column_format, print_title)

class SSESNKick(DictNpArrayMix):
    """ SSE SN kick output data from PeTar
    Keys: (class members)
        id (1D): particle id
        vkick (1D): kick velocity value (km/s)
        star (SSEStarParameter): final status of star after stellar evolution
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['id',np.int64],['vkick',np.float64],['star',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SSEType(DictNpArrayMix):
    """ SSE stellar types
    Keys: (class members), the corresponding stellar type index used in SSE is shown in [] at the end
        LMS (1D): deeply or fully convective low mass MS star [0]
        MS  (1D): Main Sequence star [1]
        HG  (1D): Hertzsprung Gap [2]
        GB  (1D): First Giant Branch [3]
        CHeB (1D): Core Helium Burning [4]
        FAGB (1D): First Asymptotic Giant Branch [5]
        SAGB (1D): Second Asymptotic Giant Branch [6]
        HeMS (1D): Main Sequence Naked Helium star [7]
        HeHG (1D): Hertzsprung Gap Naked Helium star [8]
        HeGB (1D): Giant Branch Naked Helium star [9]
        HeWD (1D): Helium White Dwarf [10]
        COWD (1D): Carbon/Oxygen White Dwarf [11]
        ONWD (1D): Oxygen/Neon White Dwarf [12]
        NS (1D): Neutron Star [13]
        BH (1D): Black Hole [14]
        SN (1D): Massless Supernova [15]
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        key_type=np.float64
        if ('key_type' in kwargs.keys()): key_type=kwargs['key_type']
        keys = [["LMS",key_type], ["MS",key_type], ["HG",key_type], ["GB",key_type], ["CHeB",key_type], ["FABG",key_type], ["SABG",key_type], ["HeMS",key_type], ["HeHG",key_type], ["HeGB",key_type], ["HeWD",key_type], ["COWD",key_type], ["ONWD",key_type], ["NS",key_type], ["BH",key_type], ["SN",key_type]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSEBinaryEvent(DictNpArrayMix):
    """ BSE binary event data from PeTar
    Keys: (class members)
        Base:
            time (1D): evolved time (Myr)
            m1   (1D): mass component 1 (Msun)
            m2   (1D): mass component 2 (Msun)
            type1 (1D): stellar type of component 1
            type2 (1D): stellar type of component 2
            semi (1D): semi-major axis (Rsun)
            ecc  (1D): eccentricity
            radro1 (1D): stellar radius of component 1 (Roche radius 1)
            radro2 (1D): stellar radius of component 2 (Roche radius 2)
            binary_type (1D): BSE binary type
        Add:
            lum1 (1D): bolometric luminosity of component 1 (Lsun) 
            lum2 (1D): bolometric luminosity of component 2 (Lsun) 
            rad1 (1D): stellar radius of component 1 (Rsun)
            rad2 (1D): stellar radius of component 2 (Rsun)
            mcore1 (1D): core mass of component 1 (Msun)
            mcore2 (1D): core mass of component 2 (Msun)
            rcore1 (1D): core radius of component 1 (Rsun)
            rcore2 (1D): core radius of component 2 (Rsun)
            spin1  (1D): stellar rotation of component 1 
            spin2  (1D): stellar rotation of component 2 

        if keyword argument "base_output" == False (default)
           class members: Base + Add
        else 
           class members: Base
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys_base = [['time',np.float64],['m1',np.float64],['m2',np.float64],['type1',np.float64],['type2',np.float64],['semi',np.float64],['ecc',np.float64],['radro1',np.float64],['radro2',np.float64],['binary_type',np.float64]]
        keys_add = [['lum1', np.float64],['lum2', np.float64],['rad1', np.float64],['rad2', np.float64],['mcore1',np.float64],['mcore2',np.float64],['rcore1',np.float64],['rcore2',np.float64],['spin1',np.float64],['spin2',np.float64]]
        keys = keys_base + keys_add
        if ('base_output' in kwargs.keys()):
            if kwargs['base_output']:
                keys = keys_base
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSETypeChange(DictNpArrayMix):
    """ BSE type change output data from PeTar
    Keys: (class members)
        type (1D): binary type, see help(petar.BSEType)
        init (BSEBinaryEvent): initial status of binary
        final (BSEBinaryEvent): final status of binary after binary stellar evolution
        id1 (1D): particle id of component 1
        id2 (1D): particle id of component 2
        drdv (1D): two component relative position dot relative velocity (Rsun*km/s)
        dr (1D): relative distance (Rsun)

    keyword arguments:
    ----------------------
    base_output (bool): False
          False: init and final have full output (see help of petar.BSEBinaryEvent)
          True:  init and final have old output 
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['type',np.int64],['init',BSEBinaryEvent],['final',BSEBinaryEvent],['id1',np.int64],['id2',np.int64],['drdv',np.float64],['dr',np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def generateBinaryID(self):
        """ Use CantorPairing to map two components id to one binary id
            Add new member bid 
        """
        bid = cantorPairing(self.id1, self.id2)
        self.addNewMember('bid',bid)

    def printTable(self, column_format = 'final', print_title = True):
        """ Print merger information in a formated table

        Parameters:
        ----------
        column_format: a list of column label (class member name), format and column title, enclosed by tuple, for sub-member, use . to access
                       For exmaple: [(key1,'%s',title1), (key2,'%12.7f',title2), (key3.subkey1,'%d',title3), (key3.subkey2,'%e',title4)]
                       Some pre-defined choices (default: final):
                       init-final: [('type', '%3d', 'kb'), ('init.type1','%4d','k1i'),('init.type2','%4d','k2i'),
                                    ('init.time','%12.4g','timei[Myr]'),('init.m1','%11.4f','m1i[M*]'),('init.m2','%11.4f','m2i[M*]'),
                                    ('init.semi','%12.4g','ai[R*]'),('init.ecc','%13.8f','ecci'),
                                    ('final.type1','%4d','k1f'),('final.type2','%4d','k2f'), 
                                    ('final.m1','%11.4f','m1f[M*]'),('final.m2','%11.4f','m2f[M*]'),('final.semi','%12.4g','af[R*]'),('final.ecc','%13.8f','eccf')]
                       final: [('type', '%3d', 'kb'),('id1','%9d','id1'), ('id2','%9d','id2'),
                              ('final.time','%12.4g','timef[Myr]'),('final.type1','%4d','k1f'),('final.type2','%4d','k2f'), 
                              ('final.m1','%10.3f','m1f[M*]'),('final.m2','%10.3f','m2f[M*]'),
                              ('final.semi','%8.2g','af[R*]'),('final.ecc','%12.7f','eccf'),
                              ('final.rad1','%8.2g','r1[R*]'),('final.rad2','%8.2g','r2[R*]'),
                              ('final.mcore1','%9.3f','mc1f[M*]'),('final.mcore2','%9.3f','mc2f[M*]'),
                              ('final.spin1','%8.2g','spin1f'),('final.spin2','%8.2g','spin2f')]
        print_title: print title of keys (default: True)
        """
        if (column_format == 'final'):
            column_format = [('type', '%3d', 'kb'), ('id1','%9d','id1'), ('id2','%9d','id2'),
                             ('final.time','%12.4g','timef[Myr]'), ('final.type1','%4d','k1f'),('final.type2','%4d','k2f'), 
                             ('final.m1','%10.3f','m1f[M*]'),('final.m2','%10.3f','m2f[M*]'),('final.semi','%8.2g','af[R*]'),('final.ecc','%13.8f','eccf'),
                             ('final.rad1','%8.2g','r1[R*]'),('final.rad2','%8.2g','r2[R*]'),#('final.lum1','%8.2g','L1f[L*]'),('final.lum2','%8.2g','L2f[L*]'),
                             ('final.mcore1','%9.3f','mc1f[M*]'),('final.mcore2','%9.3f','mc2f[M*]'),#('final.rcore1','%9.2g','rc1f[R*]'),('final.rcore2','%9.2g','rc2f[R*]'),
                             ('final.spin1','%8.2g','spin1f'),('final.spin2','%8.2g','spin2f')]

        elif (column_format =='init-final'):
            column_format = [('type', '%3d', 'kb'), ('id1','%8d','id1'), ('id2','%8d','id2'), 
                             ('init.time','%12.4g','timei[Myr]'),('init.type1','%4d','k1i'),('init.type2','%4d','k2i'), 
                             ('init.m1','%10.3f','m1i[M*]'),('init.m2','%10.3f','m2i[M*]'),('init.semi','%8.2g','ai[R*]'),('init.ecc','%13.8f','ecci'),
                             ('final.type1','%4d','k1f'),('final.type2','%4d','k2f'), 
                             ('final.m1','%10.3f','m1f[M*]'),('final.m2','%10.3f','m2f[M*]'),('final.semi','%8.2g','af[R*]'),('final.ecc','%13.8f','eccf')]
        elif (type(column_format) != list):
            raise ValueError('Error, column_format should be "final","init-final" or list of manually defined column format; given ',column_format)

        DictNpArrayMix.printTable(self, column_format, print_title)

class BSESNKick(DictNpArrayMix):
    """ BSE SN kick output data from PeTar
    Keys: (class members)
        id1 (1D): particle id of component 1
        id2 (1D): particle id of component 2
        kindex (1D): index of component which has SN kick
        vkick (1D): kick velocity (km/s)
        star (SSEStarParameter): final status of kicked star after binary stellar evolution
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['id1',np.int64],['id2',np.int64],['kindex',np.int64],['vkick',np.float64],['star',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def generateBinaryID(self):
        """ Use CantorPairing to map two components id to one binary id
            Add new member bid 
        """
        bid = cantorPairing(self.id1, self.id2)
        self.addNewMember('bid',bid)

    
class SSEStarParameterPair(DictNpArrayMix):
    """ SSE star parameter pair
    Keys: (class members)
        p1 (SSEStarParameter): SSE status of component 1
        p2 (SSEStarParameter): SSE status of component 2
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['p1',SSEStarParameter],['p2',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSEDynamicMerge(DictNpArrayMix):
    """ BSE Dynamical merger from PeTar
    Keys: (class members)
        id1 (1D): particle id of component 1
        id2 (1D): particle id of component 2
        period (1D): period  before merge (days)
        semi (1D): semi-major axis before merge (Rsun)
        ecc (1D): eccentricity before merge
        *dr (1D): separation at merger time (Rsun)
        *t_peri (1D): if merger is binary, this is the left time to peri-center, else it is 0.0 (days)
        *sd (1D): slowdown factor for the binary
        init (SSEStarParameterPair): initial status of star
        final (SSEStarParameterPair): final status of star after stellar evolution

    keyword arguments:
    ------------------
    less_output (bool): False
        True: 
           class member (key) with the prefix '*" shown above are excluded.
           This option is for the old version of PeTar before Sep 10, 2020. 
           Notice after the version of Oct 18, 2020, it is not need to use this option anymore.
           The petar.data.gether automatically fills the three columns by zero (this means dr, t_peri and sd are not correct).
        False:
           all members exists 
        
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        less_output = False
        if 'less_output' in kwargs.keys(): less_output = kwargs['less_output']

        keys_bin = [['id1',np.int64],['id2',np.int64],['period',np.float64],['semi',np.float64],['ecc',np.float64]]
        keys_extra = [['dr',np.float64],['t_peri',np.float64],['sd',np.float64]]
        keys_p= [['init',SSEStarParameterPair],['final',SSEStarParameterPair]]
        keys = keys_bin + keys_extra + keys_p
        if (less_output): keys = keys_bin + keys_p
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


    def generateBinaryID(self):
        """ Use CantorPairing to map two components id to one binary id
            Add new member bid 
        """
        bid = cantorPairing(self.id1, self.id2)
        self.addNewMember('bid',bid)

#class BSESingleCount(DictNpArrayMix):
#    """ BSE single event record
#    """
#    def __init__(self,  _dat=None, _offset=int(0), _append=False, **kwargs):
#        keys = [["n",SSEType], ["mmax", SSEType], ["mave", SSEType]]
#        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
# 
#    def findEvents(self, data):
#        keys = self.n.keys
#        
#        for ki in range(len(keys)):
#            key = keys[ki][0]
#            sel= (data.s_type==ki)
#            self.count[key] = np.append(self.count[key], sel.sum())
#            mass = data.mass[sel]
#            self.mmax[key] = np.append(self.mmax[key], np.max(mass))
#            self.mave[key] = np.append(self.mave[key], np.average(mass))
#        self.size +=1

class BSEType(DictNpArrayMix):
    """ BSE binary types
    Keys: (class members)
        Unset   (1D): no special [0]
        Initial (1D): initialization [1]
        Type_change (1D): stellar type change [2]
        Start_Roche (1D): start Roche event [3]
        End_Roche (1D): end Roche event [4]
        Contact   (1D): contact binary [5]
        Start_Symbiotic (1D): start symbiotic evolution [6]
        End_Symbiotic (1D): end symbiotic evolution [7]
        Common_envelop (1D): common envelope [8]
        Giant (1D): giant envelope [9]
        Coalescence (1D): binary coalescence [10]
        Blue_straggler (1D): blue straggler formation [11]
        No_remain (1D): no SN remainent [12]
        Disrupt (1D): binary disrupt [13]
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        key_type=np.float64
        if ('key_type' in kwargs.keys()): key_type=kwargs['key_type']
        keys = [["Unset",key_type], ["Initial",key_type], ["Type_change",key_type], ["Start_Roche",key_type], ["End_Roche",key_type], ["Contact",key_type], ["Start_Symbiotic",key_type], ["End_Symbiotic",key_type], ["Common_envelop",key_type], ["Giant",key_type],[ "Coalescence",key_type], ["Blue_straggler",key_type], ["No_remain",key_type], ["Disrupt",key_type]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class BSENumberCount(DictNpArrayMix):
    """ BSE count single/binary numbers of different stellar types
    Keys: (class members)
        single (SSEType) number of stars for each stellar type 
        binary_one (SSEType) binary with one component being the certain stellar type
        binary_both (SSEType) binary with both components being the certain stellar type
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["single", SSEType], ["binary_one", SSEType], ["binary_both",SSEType]]
        kwargs_loc=kwargs.copy()
        kwargs_loc['key_type']=np.int64
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs_loc)


class BSEStatus(DictNpArrayMix):
    """ BSE event record
    Keys: (class members)
        time (1D): evolved time
        count (BSENumberCount): number count for different stellar types
        mmax (SSEType): maximum mass of each stellar type
        mave (SSEType): average mass of each stellar type
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["time",np.float64], ["count", BSENumberCount], ["mmax",SSEType], ["mave", SSEType]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def findEvents(self, time, single, binary):
        """ Find single and binary events from one snapshot

        Parameters
        ----------
        time: float
            current evolved time
        single: inherited SimpleParticle
            single data set
        binary: Binary
            binary data set
        """
        
        self.time = np.append(self.time, time)

        keys = self.mmax.keys
        
        for ki in range(len(keys)):
            key = keys[ki][0]
            ssel = (single.star.type==ki)
            b1sel = (binary.p1.star.type==ki)
            b2sel = (binary.p2.star.type==ki)

            bsidesel = (b1sel & np.logical_not(b2sel)) | (b2sel & np.logical_not(b1sel))
            bbothsel = b1sel & b2sel

            self.count.single[key]     = np.append(self.count.single[key],ssel.sum())
            self.count.binary_one[key] = np.append(self.count.binary_one[key],bsidesel.sum())
            self.count.binary_both[key]= np.append(self.count.binary_both[key],bbothsel.sum())
            
            smass = single.mass[ssel]
            b1mass = binary.p1.mass[b1sel]
            b2mass = binary.p2.mass[b2sel]
            
            mass = np.concatenate((smass, b1mass, b2mass))

            if (mass.size>0):
                self.mmax[key] = np.append(self.mmax[key], np.amax(mass))
                self.mave[key] = np.append(self.mave[key], np.average(mass))
            else:
                self.mmax[key] = np.append(self.mmax[key], 0.0)
                self.mave[key] = np.append(self.mave[key], 0.0)

        self.count.single.size += 1
        self.count.binary_one.size += 1
        self.count.binary_both.size += 1
        self.count.size += 1
        self.mmax.size += 1
        self.mave.size += 1
        self.size +=1
        
        #self.printSize()

class SSEStarParameterOut(DictNpArrayMix):
    """ SSE star parameter output class from bse_interface.h
    Keys: (class members)
         type0 (1D): original type before evolution, see help(petar.SSEType)
         menv (1D): mass of convective envelope (Msun)
         renv (1D): radius of convective envelope (Rsun)
         tm  (1D): Main sequence lifetime (Myr)
         vkick (2D,4): kick velocity for NS/BH formation, vx, vy, vz, |v|
         dm (1D): mass loss (Msun)
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['type0', np.int64], ['menv', np.float64], ['renv', np.float64], ['tm', np.float64], ['vkick', (np.float64, 4)], ['dm', np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SSEISO(DictNpArrayMix):
    """ Binary stellar evolution tool petar.(mo)bse output for single stars
    Keys: (class members)
        mass_init (1D): initial mass at zero age
        star (SSEStarParameter): SSE star parameter 
        out (SSEStarParameterOut): SSE star parameter output 
    """
    
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['mass_init',np.float64],['star',SSEStarParameter],['out',SSEStarParameterOut]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    

class BSEISO(DictNpArrayMix):
    """ Binary stellar evolution tool petar.(mo)bse output for binaries
    Keys: (class members)
        m1_init (1D): initial mass of component 1 (Msun)
        m2_init (1D): initial mass of component 2 (Msun)
        period_init (1D): initial period (days)
        ecc_init (1D): initial eccentricity
        period_final (1D): final period (days)
        ecc_final (1D): final eccentricity 
        star1 (SSEStarParameter): SSE star parameter of component 1
        out1 (SSEStarParameterOut): SSE star parameter output of component 1
        star2 (SSEStarParameter): SSE star parameter of component 2
        out2 (SSEStarParameterOut): SSE star parameter output of component 2
    """
    
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['m1_init',np.float64],['m2_init',np.float64],['period_init',np.float64],['ecc_init',np.float64],['period_final',np.float64],['ecc_final',np.float64],['star1',SSEStarParameter],['out1',SSEStarParameterOut],['star2',SSEStarParameter],['out2',SSEStarParameterOut]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSEMerge(DictNpArrayMix):
    """ BSE binary mergers 
    Keys: (class members)
        time  (1D): current physical time (Myr)
        bid   (1D): binary id obtained from CantorPairing of two component IDs.
        id1   (1D): particle id of component 1
        id2   (1D): particle id of component 2
        semi  (1D): semi-major axis (Rsun)
        ecc   (1D): eccentricity
        type1 (1D): stellar type of component 1 before merge
        type2 (1D): stellar type of component 2 before merge
        m1    (1D): mass component 1 (Msun) before merge
        m2    (1D): mass component 2 (Msun) before merge
        typef (1D): stellar type of merger
        mf    (1D): merger mass (Msun)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation using key list, see help(DictNpArrayMix.__init__)
        """
        keys=[['time',np.float64],['bid', np.int64],['id1',np.int64],['id2',np.int64],['semi',np.float64],['ecc',np.float64],
              ['type1',np.int64],['type2',np.int64],
              ['m1',np.float64],['m2',np.float64],['typef',np.int64],['mf',np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


    def getMergeFromTypeChange(self,type_change):
        """
        Find mergers from binary type change data, and find the first initial data and mix it with the final data as the final merger
        """
        
        def calc_ubid(_type_change):
            """
            Calculate unique id by adding bid and final time to a value of by using numpy.float128
            """
            #return cantorPairing(_type_change.bid, (np.log(_type_change.final.time)*1e6).astype(int))
            return np.float128(_type_change.bid) + np.float128(_type_change.final.time)
        
        # find mergers
        type_change.generateBinaryID()
        # select events where one component becomes massless type and the other has mass increase.
        sel_merge = ((type_change.final.m1==0) & (type_change.init.m1>0) & (type_change.final.m2 > type_change.init.m2)) | ((type_change.final.m2==0) & (type_change.init.m2>0) & (type_change.final.m1 > type_change.init.m1))
        merge_final = type_change[sel_merge]
        bid_merge = merge_final.bid
        
        # if no merger
        if (bid_merge.size==0): return merge_final
        
        # get ubid and sort merger by ubid
        ubid_merge = calc_ubid(merge_final)
        sindex = ubid_merge.argsort()
        merge_final = merge_final[sindex]
        ubid_merge = ubid_merge[sindex]
        
        # get merger progenitor history
        sel_history=np.in1d(type_change.bid, bid_merge)
        merge_history = type_change[sel_history]
        ubid_history = calc_ubid(merge_history)
        
        # find merger history with same ubid to obtain the first initial status
        temp, index = np.unique(ubid_history, return_index=True)
        ubid_history_unique = ubid_history[index]
        sel_merge_init = np.in1d(ubid_history_unique, ubid_merge)
        ubid_merge_init = ubid_history_unique[sel_merge_init]
        merge_init = merge_history[index][sel_merge_init]
        # sort merge init with ubid
        sindex = ubid_merge_init.argsort()
        merge_init = merge_init[sindex]
        ubid_merge_init = ubid_merge_init[sindex]
        
        # validate whether the ubid of init and final data match
        match_check = (ubid_merge-ubid_merge_init).sum()
        if (match_check>0):
            raise ValueError('Find merger initial status fails, ubid unmatch, the sum value is %d, should be zero!' % match_check)
        
        # get merge init status as init and final status as final
        merge_mix = merge_init
        merge_mix.final = merge_final.final
        return merge_mix

    def combine(self, type_change, dyn_merge):
        """ Find mergers from type_change and combine them with those from dynamical merge
        bid will be generated
        
        Parameters:
        -----------
        type_change: BSETypeChange data
        dyn_merge: BSEDynamicMerge data
        """
        #type_change.generateBinaryID()
        #sel = (type_change.final.m1==0) | (type_change.final.m2==0)
        #se_merge = type_change[sel]
        se_merge = self.getMergeFromTypeChange(type_change)

        dyn_merge.generateBinaryID()
        merge = self.__dict__
        merge['time']=np.concatenate((se_merge.final.time,dyn_merge.final.p1.time))
        tsort=merge['time'].argsort()
        merge['time']=merge['time'][tsort]
        merge['bid']=np.concatenate((se_merge.bid, dyn_merge.bid))[tsort]
        merge['id1']=np.concatenate((se_merge.id1,dyn_merge.id1))[tsort]
        merge['id2']=np.concatenate((se_merge.id2,dyn_merge.id2))[tsort]
        merge['semi']=np.concatenate((se_merge.init.semi,dyn_merge.semi))[tsort]
        merge['ecc']=np.concatenate((se_merge.init.ecc,dyn_merge.ecc))[tsort]
        merge['type1']=np.concatenate((se_merge.init.type1,dyn_merge.init.p1.type))[tsort]
        merge['type2']=np.concatenate((se_merge.init.type2,dyn_merge.init.p2.type))[tsort]
        merge['m1']=np.concatenate((se_merge.init.m1,dyn_merge.init.p1.mass))[tsort]
        merge['m2']=np.concatenate((se_merge.init.m2,dyn_merge.init.p2.mass))[tsort]

        se_kwf = se_merge.final.type1
        se_mf = se_merge.final.m1
        se_mf2_sel = (se_merge.final.m1 == 0)
        se_kwf[se_mf2_sel] = se_merge.final.type2[se_mf2_sel]
        se_mf[se_mf2_sel] = se_merge.final.m2[se_mf2_sel]

        dyn_kwf = dyn_merge.final.p1.type
        dyn_mf = dyn_merge.final.p1.mass
        dyn_mf2_sel = (dyn_merge.final.p1.mass == 0)
        dyn_kwf[dyn_mf2_sel] = dyn_merge.final.p2.type[dyn_mf2_sel]
        dyn_mf[dyn_mf2_sel] = dyn_merge.final.p2.mass[dyn_mf2_sel]

        merge['typef']=np.concatenate((se_kwf, dyn_kwf))[tsort]
        merge['mf']=np.concatenate((se_mf, dyn_mf))[tsort]
        self.size=merge['time'].size

    def printTable(self, 
                   column_format = [('time','%12.7f','time[Myr]'), ('bid','%14d','bid'),
                                    ('id1','%8d','id1'), ('id2','%8d','id2'), 
                                    ('semi','%10.7g','semi[R*]'), ('ecc','%15.10f','ecc'),
                                    ('type1','%4d','k1i'), ('type2','%4d','k2i'),
                                    ('m1','%12.4f','m1i[M*]'), ('m2','%12.4f','m2i[M*]'),
                                    ('typef','%4d','kf'), ('mf','%12.7f','mf[M*]')],
                   print_title = True):
        """ Print merger information in a formated table

        Parameters:
        ----------
        column_format: a list of column label (class member name), format and column title, enclosed by tuple, for sub-member, use . to access
                       For exmaple: [(key1,'%s',title1), (key2,'%12.7f',title2), (key3.subkey1,'%d',title3), (key3.subkey2,'%e',title4)]
                       Default: [('time','%12.7f','time[Myr]'), ('bid','%14d','bid'),
                                 ('id1','%8d','id1'), ('id2','%8d','id2'), 
                                 ('semi','%10.7g','semi[R*]'), ('ecc','%15.10f','ecc'),
                                 ('type1','%4d','k1i'), ('type2','%4d','k2i'),
                                 ('m1','%12.4f','m1i[M*]'), ('m2','%12.4f','m2i[M*]'),
                                 ('typef','%4d','kf'), ('mf','%12.7f','mf[M*]')]
        print_title: print title of keys (default: True)

        """
        DictNpArrayMix.printTable(self, column_format, print_title)

class GWKick(DictNpArrayMix):
    """ Gravitational wave kick data
    Keys: (class members)
        id1 (1D): particle id of component 1
        id2 (1D): particle id of component 2
        kindex (1D): index of component which has final mass and GW kick
        vkick (3D,3): kick velocity (km/s), vx, vy, vz
        time (1D): time of the merger (Myr)
        m1 (1D): particle mass of component 1 before merger
        m2 (1D): particle mass of component 2 before merger
        mf (1D): final mass after merger
        semi (1D): semi-major axis before merger (Rsun)
        ecc (1D): eccentricity before merger
        chi1 (3D,3): dimensionless spin of particle 1
        chi2 (3D,3): dimensionless spin of particle 2
        am   (2D,3): specific angular momemtum x, y, z
        dr   (2D,3): relative position vector from particle 1 to 2
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation using key list, see help(DictNpArrayMix.__init__)
        """
        keys=[['id1', np.int64], ['id2', np.int64], ['kindex', np.int64], ['vkick', (np.float64, 3)], 
              ['time', np.float64], ['m1', np.float64], ['m2', np.float64], ['mf', np.float64], 
              ['semi', np.float64], ['ecc', np.float64], ['chi1', (np.float64, 3)], ['chi2', (np.float64, 3)], 
              ['am', (np.float64, 3)], ['dr', (np.float64, 3)]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


def find_merge_tree(merger_list, merger_root):
    """ Find the merger tree for a given merger
    Parameters:
    -----------
    merger_list: BSEMerge data
    merger_root: the target merger to find tree
    Return:
    ------------
    merger_tree: numpy.ndarray(3D)
        The merger tree table for plotting, each row is one tree branch 
        The row contains 4 pair of data pointing from the leaf to the root:
           1. The position in the tree branch. The value of the leaf is based on that of the root +- 0.5^{level+1}, where the level refers to the root.
           2. indice counting from the left to the right of all leaves and roots
           3. times of components (leaves) and mergers (roots)
           4. masses of components (leaves) and mergers (roots)
    """
    def find_merge_tree_iter(merger_list, merger_root, merger_tree, binary_tree_base, binary_tree_interval, indebinary_tree_base):
        """
        
        """
        idlst=[merger_root.id1, merger_root.id2]
        btlst=[-binary_tree_interval,binary_tree_interval]
        mlst=[merger_root.m1, merger_root.m2]
        ctotlst=[indebinary_tree_base,1]
        for k in range(len(idlst)):
            if (k>0): 
                ctotlst[k] += ctotlst[k-1]
            idk = idlst[k]
            sel = ((merger_list.id1 == idk) | (merger_list.id2 == idk)) & (merger_list.time <= merger_root.time) & (merger_list.mf < merger_root.mf)
            sdat = merger_list[sel]
            if (sdat.size>0):
                c_left, c_right = find_merge_tree_iter(merger_list, sdat[-1], merger_tree, binary_tree_base+btlst[k], binary_tree_interval/2, ctotlst[k])  
                ctotlst[k] += c_right - ctotlst[k]
                merger_tree.append([[binary_tree_base+btlst[k], binary_tree_base], [c_left+1, ctotlst[0]+1], [sdat.time[-1], merger_root.time],[sdat.mf[-1], merger_root.mf]])
            else:
                ctotlst[k] += 1
                merger_tree.append([[binary_tree_base+btlst[k], binary_tree_base], [ctotlst[k], ctotlst[0]+1], [0, merger_root.time], [mlst[k], merger_root.mf]])
        return ctotlst[0], ctotlst[1]
    merger_tree=[]
    find_merge_tree_iter(merger_list, merger_root, merger_tree, 0.5, 0.25, 0)

    return np.array(merger_tree)
