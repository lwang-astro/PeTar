import numpy as np
from .base import *
from .data import *

class SingleEscaper(Particle):
    """ Single escaper recorder output from PeTar
    Keys: (class members)
        time (1D): escaping time
        [Particle] (inherited)
    """
    
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            particle_type: string (soft)
                Basic particle type: hermite, hard, soft
                When read PeTar data, do not change this
            interrupt_mode: string (none)
                PeTar interrupt mode: base, bse, mobse, none
        """

        
        
        DictNpArrayMix.__init__(self, [['time',np.float64]], _dat, _offset, _append, **kwargs)
        Particle.__init__(self, _dat, _offset+self.ncols, True, **kwargs)

    def findEscaper(self, time, single, rcut, ecut=0.0):
        """ Find escaper from a snapshot
        Functions, calcR2, calcEkin and calcEtot, are used first for input single data set,
        then distance r<rcut and energy etot>ecut will be selected as escapers.
        The escapers will be removed from the input single
        
        Parameters
        ----------
        time: float
            evolved time of snapshot
        single: inherited SingleParticle
            single particle data set
        rcut: float
            distance criterion
        ecut: float (0.0)
            energy criterion
        """
        single.calcR2()
        single.calcEkin()
        single.calcEtot()
        
        rcut2 = rcut*rcut
        ssel = ((single.r2>rcut2) & (single.etot>ecut))
        single_esc = single[ssel]
        nssel = ssel.sum()
        single_esc.addNewMember('time',np.ones(nssel)*time)
        self.append(single_esc)
        #idsinx = self.single.time.argsort()
        #self.single=self.single[idsinx]
        if (self.size>0): self.removeDuplicate()
        return single[np.logical_not(ssel)]

    def removeDuplicate(self):
        """ removed duplicated escapers, keep the first appearing one
        """
        sindex=self.time.argsort()
        data_sort=self[sindex]
        unid, index= np.unique(data_sort.id, return_index=True)
        newdata=data_sort[index]
        self.__init__(newdata,**self.initargs)

class BinaryEscaper(Binary):
    """ Binary escaper information
    Keys: (class members)
        time (1D): escaping time
        [Binary] (inherited)
        
    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            particle_type: string (soft)
                Basic particle type: hermite, hard, soft
            interrupt_mode: string (none)
                PeTar interrupt mode: base, bse, mobse, none
            simple_mode: bool (True)
                If True, only calculate semi and ecc, save computing time significantly
            member_particle_type: type (Particle)
                Binary member particle type
            G: float (1.0) 
                gravitational constant (1.0)
        """
        if not 'member_particle_type' in kwargs.keys(): kwargs['member_particle_type']=Particle

        DictNpArrayMix.__init__(self, [['time',np.float64]], _dat, _offset, _append, **kwargs)
        Binary.__init__(self, _dat, None, _offset+self.ncols, True, **kwargs)

    def findEscaper(self, time, binary, rcut, ecut=0.0):
        """ Find escaper from a binary snapshot
        Functions, calcR2, calcEkin, calcPot and calcEtot, are used first for input binary data set,
        then c.m. distance r<rcut and c.m. energy etot>ecut will be selected as escapers.
        The escapers will be removed from the input binary
        
        Parameters
        ----------
        time: float
            evolved time of snapshot
        binary: Binary
            binary particle data set
        rcut: float
            distance criterion
        ecut: float (0.0)
            energy criterion
        """
        binary.calcR2()
        binary.calcPot()
        binary.calcEkin()
        binary.calcEtot()

        rcut2 = rcut*rcut
        bsel = (binary.r2>rcut2) & (binary.etot>ecut)
        binary_esc = binary[bsel]
        nbsel = bsel.sum()
        binary_esc.addNewMember('time',np.ones(nbsel)*time)
        self.append(binary_esc)
        #idsinx = self.binary.time.argsort()
        #self.binary=self.binary[idsinx]
        if (self.size>0): self.removeDuplicate()
        return binary[np.logical_not(bsel)]

    def removeDuplicate(self):
        """ removed duplicated escapers, keep the first appearing one
        Use the first component id to check duplicate. 
        """
        unid, index= np.unique(self.p1.id, return_index=True)
        newdata=self[index]
        self.__init__(newdata,**self.initargs)

def calcRCutIsolate(rh):
    """ For isolated star clusters, set rcut to 20 * half-mass radius

    Return
    ----------
    rcut: float 
        escaper distance criterion
    """
    return 20*rh


#def joinEscaper(*esc_list):
#    single_type = type(esc_list[0].single)
#    esc_merge = Escaper(single_type)
#    for ei in esc_list:
#        esc_merge.rcut = ei.rcut
#        esc_merge.single.append(ei.single)
#        esc_merge.binary.append(ei.binary)
#    unid, index= np.unique(esc_merge.single.id, return_index=True)
#    esc_merge.single = esc_merge.single[index]
# 
#    unid, index= np.unique(esc_merge.binary.p1.id, return_index=True)
#    esc_merge.binary = esc_merge.binary[index]
# 
#    return esc_merge
