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

    def findEscaper(self, time, single, rcut, es_cut=0.0):
        """ Find escaper from a snapshot
        Functions, calcR2, calcEkin and calcEtot, are used first for input single data set,
        then distance r<rcut and energy etot > mass*es_cut will be selected as escapers.
        The escapers will be removed from the input single
        
        Parameters
        ----------
        time: float
            evolved time of snapshot
        single: inherited SingleParticle
            single particle data set
        rcut: float
            distance criterion
        es_cut: float (0.0)
            specific energy criterion
        """
        single.calcR2()
        single.calcEkin()
        single.calcEtot()
        
        rcut2 = rcut*rcut
        ssel = ((single.r2>rcut2) & (single.etot-single.mass*es_cut>0.0))
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

    def findEscaper(self, time, binary, rcut, es_cut=0.0):
        """ Find escaper from a binary snapshot
        Functions, calcR2, calcEkin, calcPot and calcEtot, are used first for input binary data set,
        then c.m. distance r<rcut and c.m. energy etot > mass*es_cut will be selected as escapers.
        The escapers will be removed from the input binary
        
        Parameters
        ----------
        time: float
            evolved time of snapshot
        binary: Binary
            binary particle data set
        rcut: float
            distance criterion
        es_cut: float (0.0)
            specific energy criterion
        """
        binary.calcR2()
        binary.calcPot()
        binary.calcEkin()
        binary.calcEtot()

        rcut2 = rcut*rcut
        bsel = (binary.r2>rcut2) & (binary.etot-binary.mass*es_cut>0.0)
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

class Tidal(DictNpArrayMix):
    """ tidal radius and potential of the external potential of the particle system
    keys: (class members)
        time (1D): time
        rtid (1D): tidal radius 
        pot (1D): the external potential of the particle system
        mass (1D): total bound mass
        n (1D): total bound number of particles
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['time', np.float64], ['rtid', np.float64], ['pot', np.float64], ['mass', np.float64], ['n', np.int64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
    def calcTidalSphere(self, time, mass, r2, pot_ext, rc, center_pos, G):
        """ calculate tidal radius assuming R_tid = M_system/ (3*M_galaxy)^(1/3) * R_galaxy; 
        and potential of the central position by averaging the pot_ext inside core radius

        Parameters
        ----------
        time: float
            current time
        mass: 1D numpy.ndarray
            masses of particles
        r2: 1D numpy.ndarray
            the distance square to the center of the system
        pot_ext: 1D numpy.ndarray
            external potential of particles
        rc: float
            core radius
        center_pos: 1D numpy.ndarray
            the position of center in the galactic frame
        G: float
            gravitational constant

        Return
        ----------
        r_tid: float
            tidal radius
        pot_ext_cave: float
            external potential of the center
        """
        self.time = np.append(self.time, time)

        r_gal_2 = (center_pos*center_pos).sum()
        r_gal = np.sqrt(r_gal_2)
        rc2 = rc*rc
        mtot = mass.sum()
        csel = r2<rc2
        mc = mass[csel]
        mctot = mc.sum()
        pot_ext_c = pot_ext[csel]
        pot_ext_cave = (pot_ext_c*mc).sum()/mctot
        M_galaxy = - pot_ext_cave*r_gal/G
        if (M_galaxy<0):
            raise ValueError('External potential is positive! ', pot_ext_cave)
    
        r_tid = (mtot/(3*M_galaxy))**(1.0/3.0)*r_gal
        rt2 = r_tid*r_tid
        r_tid_old=r_tid*1.2
        rtsel = r2<rt2
        mcut = mass[rtsel]
        r2cut= r2[rtsel]

        while ((r_tid_old-r_tid)/r_tid_old>1e-3):
            r_tid_old = r_tid
            rt2 = r_tid*r_tid
            rtsel = r2cut<rt2
            mcut = mcut[rtsel]
            r2cut = r2cut[rtsel]
            mtot = mcut.sum()
            r_tid = (mtot/(3*M_galaxy))**(1.0/3.0)*r_gal

        self.rtid = np.append(self.rtid, r_tid)
        self.pot  = np.append(self.pot, pot_ext_cave)
        self.mass = np.append(self.mass, mtot)
        self.n = np.append(self.n, mcut.size)
        self.size += 1
    
        return r_tid, pot_ext_cave

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
