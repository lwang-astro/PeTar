# BSE data interface

from .base import *

class SSEStarParameter(DictNpArrayMix):
    """ SSE star parameter class from bse_interface.h
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['type',1],['mass0',1],['mass',1],['rad',1],['mcore',1],['rcore',1],['spin',1],['epoch',1],['time',1],['lum',1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class SSETypeChange(DictNpArrayMix):
    """ SSE type change output data from PeTar
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['id',1],['init',SSEStarParameter],['final',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SSESNKick(DictNpArrayMix):
    """ SSE SN kick output data from PeTar
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['id',1],['vkick',1],['star',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SSEType(DictNpArrayMix):
    """ SSE stellar types
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["LMS",1], ["MS",1], ["HG",1], ["GB",1], ["CHeB",1], ["FABG",1], ["SABG",1], ["HeMS",1], ["HeHG",1], ["HeGB",1], ["HeWD",1], ["COWD",1], ["ONWD",1], ["NS",1], ["BH",1], ["SN",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSEBinaryEvent(DictNpArrayMix):
    """ BSE binary event data from PeTar
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['time',1],['m1',1],['m2',1],['type1',1],['type2',1],['semi',1],['ecc',1],['rad1',1],['rad2',1],['binary_type',1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSETypeChange(DictNpArrayMix):
    """ BSE type change output data from PeTar
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['init',BSEBinaryEvent],['final',BSEBinaryEvent],['id1',1],['id2',1],['drdv',1],['dr',1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSESNKick(DictNpArrayMix):
    """ BSE SN kick output data from PeTar
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['id1',1],['id2',1],['kindex',1],['vkick',1],['star',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
class SSEStarParameterPair(DictNpArrayMix):
    """ SSE star parameter pair
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['p1',SSEStarParameter],['p2',SSEStarParameter]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSEDynamicMerge(DictNpArrayMix):
    """ BSE Dynamical merger from PeTar
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
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
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["Unset",1], ["Initial",1], ["Type_change",1], ["Start_Roche",1], ["End_Roche",1], ["Contact",1], ["Start_Symbiotic",1], ["End_Symbiotic",1], ["Common_envelop",1], ["Giant",1],[ "Coalescence",1], ["Blue_straggler",1], ["No_remain",1], ["Disrupt",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class BSENumberCount(DictNpArrayMix):
    """ BSE count single/binary numbers of different stellar types
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["single", SSEType], ["binary_one", SSEType], ["binary_both",SSEType]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class BSEStatus(DictNpArrayMix):
    """ BSE event record
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["time",1], ["count", BSENumberCount], ["mmax",SSEType], ["mave", SSEType]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def findEvents(self, time, single, binary):
        """ Find single and binary events from one snapshot

        Parameters
        ----------
        single: single snapshot
        binary: binary snapshot
        """
        
        self.time = np.append(self.time, time)

        keys = self.mmax.keys
        
        for ki in range(len(keys)):
            key = keys[ki][0]
            ssel = (single.s_type==ki)
            b1sel = (binary.p1.s_type==ki)
            b2sel = (binary.p2.s_type==ki)

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

