# BSE data
from .base import *

class BSESingleType(DictNpArrayMix):
    """ BSE stellar types
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["LMS",1], ["MS",1], ["HG",1], ["GB",1], ["CHeB",1], ["FABG",1], ["SABG",1], ["HeMS",1], ["HeHG",1], ["HeGB",1], ["HeWD",1], ["COWD",1], ["ONWD",1], ["NS",1], ["BH",1], ["SN",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class BSESingleEvent(DictNpArrayMix):
    """ BSE single event record
    """
    def __init__(self,  _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["n",BSESingleType], ["mmax", BSESingleType], ["mave", BSESingleType]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
 
    def findEvents(self, data):
        keys = self.n.keys
        
        for ki in range(len(keys)):
            key = keys[ki][0]
            sel= (data.s_type==ki)
            self.count[key] = np.append(self.count[key], sel.sum())
            mass = data.mass[sel]
            self.mmax[key] = np.append(self.mmax[key], np.max(mass))
            self.mave[key] = np.append(self.mave[key], np.average(mass))
        self.size +=1

class BSEBinaryType(DictNpArrayMix):
    """ BSE binary types
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["Unset",1], ["Initial",1], ["Type_change",1], ["Start_Roche",1], ["End_Roche",1], ["Contact",1], ["Start_Symbiotic",1], ["End_Symbiotic",1], ["Common_envelop",1], ["Giant",1],[ "Coalescence",1], ["Blue_straggler",1], ["No_remain",1], ["Disrupt",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

#class BSEBinaryEvent(DictNpArrayMix):
#    """ BSE binary event record
#    """
#    def __init__(self,  _dat=None, _offset=int(0), _append=False, **kwargs):
#        keys = [["count",BSEBinaryType], ["mmax", BSEBinaryType], ["mave", BSEBinaryType]]
#        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class BSECount(DictNpArrayMix):
    """ BSE count single/binary numbers of different stellar types
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["single", BSESingleType], ["binary_one", BSESingleType], ["binary_both",BSESingleType]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class BSEEvent(DictNpArrayMix):
    """ BSE event record
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["time",1], ["count", BSECount], ["mmax",BSESingleType], ["mave", BSESingleType]]
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
