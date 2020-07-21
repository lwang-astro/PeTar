# BSE data interface

from .base import *

class SSEStarParameter(DictNpArrayMix):
    """ SSE star parameter class from bse_interface.h
    Keys:
        type  (1D): SSE stellar type
        mass0 (1D): initial mass at each evolution stage (Msun)
        mass  (1D): current mass (Msun)
        rad   (1D):  stellar radius
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
        keys = [['type',1],['mass0',1],['mass',1],['rad',1],['mcore',1],['rcore',1],['spin',1],['epoch',1],['time',1],['lum',1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class SSETypeChange(DictNpArrayMix):
    """ SSE type change output data from PeTar
    Keys:
        id (1D): particle id
        init (SSEStarParameter): initial status of star
        final (SSEStarParameter): final status of star after stellar evolution
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['id',1],['init',SSEStarParameter],['final',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SSESNKick(DictNpArrayMix):
    """ SSE SN kick output data from PeTar
    Keys:
        id (1D): particle id
        vkick (1D): kick velocity value (km/s)
        star (SSEStarParameter): final status of star after stellar evolution
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['id',1],['vkick',1],['star',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SSEType(DictNpArrayMix):
    """ SSE stellar types
    Keys: 
        LMS (1D): deeply or fully convective low mass MS star
        MS  (1D): Main Sequence star
        HG  (1D): Hertzsprung Gap
        GB  (1D): First Giant Branch
        CHeB (1D): Core Helium Burning
        FAGB (1D): First Asymptotic Giant Branch
        SAGB (1D): Second Asymptotic Giant Branch
        HeMS (1D): Main Sequence Naked Helium star
        HeHG (1D): Hertzsprung Gap Naked Helium star
        HeGB (1D): Giant Branch Naked Helium star
        HeWD (1D): Helium White Dwarf
        COWD (1D): Carbon/Oxygen White Dwarf
        ONWD (1D): Oxygen/Neon White Dwarf
        NS (1D): Neutron Star
        BH (1D): Black Hole
        SN (1D): Massless Supernova
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["LMS",1], ["MS",1], ["HG",1], ["GB",1], ["CHeB",1], ["FABG",1], ["SABG",1], ["HeMS",1], ["HeHG",1], ["HeGB",1], ["HeWD",1], ["COWD",1], ["ONWD",1], ["NS",1], ["BH",1], ["SN",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSEBinaryEvent(DictNpArrayMix):
    """ BSE binary event data from PeTar
    Keys:
        time (1D): evolved time (Myr)
        m1   (1D): mass component 1 (Msun)
        m2   (1D): mass component 2 (Msun)
        type1 (1D): stellar type of component 1
        type2 (1D): stellar type of component 2
        semi (1D): semi-major axis (Rsun)
        ecc  (1D): eccentricity
        rad1 (1D): stellar radius of component 1 (Rsun)
        rad2 (1D): stellar radius of component 2 (Rsun)
        binary_type (1D): BSE binary type
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['time',1],['m1',1],['m2',1],['type1',1],['type2',1],['semi',1],['ecc',1],['rad1',1],['rad2',1],['binary_type',1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSETypeChange(DictNpArrayMix):
    """ BSE type change output data from PeTar
    Keys:
        init (BSEBinaryEvent): initial status of binary
        final (BSEBinaryEvent): final status of binary after binary stellar evolution
        id1 (1D): particle id of component 1
        id2 (1D): particle id of component 2
        drdv (1D): two component relative position dot relative velocity
        dr (1D): relative distance
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['init',BSEBinaryEvent],['final',BSEBinaryEvent],['id1',1],['id2',1],['drdv',1],['dr',1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSESNKick(DictNpArrayMix):
    """ BSE SN kick output data from PeTar
    Keys:
        id1 (1D): particle id of component 1
        id2 (1D): particle id of component 2
        kindex (1D): index of component which has SN kick
        vkick (1D): kick velocity (km/s)
        star (SSEStarParameter): final status of kicked star after binary stellar evolution
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['id1',1],['id2',1],['kindex',1],['vkick',1],['star',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
class SSEStarParameterPair(DictNpArrayMix):
    """ SSE star parameter pair
    Keys:
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
    Keys:
        id1 (1D): particle id of component 1
        id2 (1D): particle id of component 2
        period (1D): period  before merge (days)
        semi (1D): semi-major axis before merge (Rsun)
        ecc (1D): eccentricity before merge
        init (SSEStarParameter): initial status of star
        final (SSEStarParameter): final status of star after stellar evolution
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['id1',1],['id2',1],['period',1],['semi',1],['ecc',1],['init',SSEStarParameterPair],['final',SSEStarParameterPair]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

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
    Keys:
        Unset   (1D): no special
        Initial (1D): initialization
        Type_change (1D): stellar type change
        Start_Roche (1D): start Roche event
        End_Roche (1D): end Roche event
        Contact   (1D): contact binary
        Start_Symbiotic (1D): start symbiotic evolution
        End_Symbiotic (1D): end symbiotic evolution
        Common_envelop (1D): common envelope
        Giant (1D):
        Coalescence (1D): binary coalescence
        Blue_straggler (1D): blue straggler formation
        No_remain (1D): no SN remainent 
        Disrupt (1D): binary disrupt
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["Unset",1], ["Initial",1], ["Type_change",1], ["Start_Roche",1], ["End_Roche",1], ["Contact",1], ["Start_Symbiotic",1], ["End_Symbiotic",1], ["Common_envelop",1], ["Giant",1],[ "Coalescence",1], ["Blue_straggler",1], ["No_remain",1], ["Disrupt",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class BSENumberCount(DictNpArrayMix):
    """ BSE count single/binary numbers of different stellar types
    Keys:
        single (SSEType) number of stars for each stellar type 
        binary_one (SSEType) binary with one component being the certain stellar type
        binary_both (SSEType) binary with both components being the certain stellar type
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["single", SSEType], ["binary_one", SSEType], ["binary_both",SSEType]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class BSEStatus(DictNpArrayMix):
    """ BSE event record
    Keys:
        time (1D): evolved time
        count (BSENumberCount): number count for different stellar types
        mmax (SSEType): maximum mass of each stellar type
        mave (SSEType): average mass of each stellar type
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["time",1], ["count", BSENumberCount], ["mmax",SSEType], ["mave", SSEType]]
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

