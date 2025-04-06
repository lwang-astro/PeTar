import collections
from .base import *
from .data import *

class BinaryTreeSDAR(DictNpArrayMix):
    """ Binary tree data output from SDAR 
    Keys: (class members)
        semi (1D): semi-major axis
        ecc  (1D): eccentricity
        incline (1D): inclination
        rot_horizon (1D): frame rotational angle in x-y plane (longitude of ascending node)
        rot_self (1D): frame rotational angle in orbital plane (argument of periapsis)
        t_peri (1D): time to peri-center
        period (1D): period
        ecca (1D): eccentric anomaly
        m1   (1D): component 1 mass
        m2   (1D): component 2 mass
        rrel (1D): relative distance
        am   (2D,3): specific angular momemtum x, y, z
        stab (1D): stability factor (>1: unstable)
        sd   (1D): slowdown factor
        sd_org(1D): original slowdown factor based on perturbation
        sd_max(1D): maximum slowdown factor based on timescale criterion
        p1 (member_particle_type) component one
        p2 (member_particle_type) component two
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            member_particle_type: type (HardParticle)
                Type of component particle 
        """
        member_particle_type=HardParticle
        if 'member_particle_type' in kwargs.keys(): member_particle_type=kwargs['member_particle_type']

        keys = [['semi',np.float64], ['ecc',np.float64], ['incline',np.float64],['rot_horizon',np.float64],['rot_self',np.float64],['t_peri',np.float64],['period',np.float64],['ecca',np.float64],['m1',np.float64],['m2',np.float64],['r',np.float64],['am',(np.float64,3)],['stab',np.float64],['sd',np.float64],['sd_org',np.float64],['sd_max',np.float64],['p1',member_particle_type],['p2',member_particle_type]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append,**kwargs)

    def generateBinaryID(self):
        """ Use CantorPairing to map two components id to one binary id
            Add new member bid 
        """
        bid = cantorPairing(self.p1.id, self.p2.id)
        self.addNewMember('bid',bid)
        
class GroupInfo(DictNpArrayMix):
    """ Group information output from PeTar
    Keys: (class members)
        type (1D): group type, 0: new group; 1: end group
        n    (1D): number of members in group (should be consistent with keyword argument N
        time (1D): current time
        pos  (2D,3): position of the group c.m. in the framework of the global system (without shift of global system c.m. if external_mode is on)
        vel  (2D,3): velocity of the group c.m. in the framework of the global system (without shift of global system c.m.)
        bin[X] (BinaryTreeSDAR): members of the group in a hierarchical binary tree
               Here X indicates the order. 0 represents the root (outer most) binary; 1,2,3 ... are inner binaries
               For a triple, bin0 is outer binary, bin1 is inner binary.
               p2 of bin0 is the c.m. of bin1, the id of p2 is the minimum id from the two components in bin1.
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        May receive the warning message: 
        RuntimeWarning: invalid value encountered in cast self.__dict__[key] = _dat[:,icol].astype(parameter)
        This is due to the artificial particle mode of mass_bk and status data, which are F64 instead of S64

        Parameters
        ----------
        keyword arguments:
            member_particle_type: type (HardParticle)
                Type of component particle, do not change this!
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): base, bse, mobse, none
               This option indicates whether columns of stellar evolution exist
            external_mode: string (none)
               PeTar external mode (set in configure): galpy, none 
               This option indicates whether the column of externa potential exist
            use_mpfrc: bool (False)
               If true, add three columns of pos_high indicating the high-precision parts of position
            float_type: type (np.float64)
                floating point data type
            N: int (2)
                Number of members of one group
        """
        keys=[['type',np.int64],['n',np.int64],['time',np.float64],['pos',(np.float64,3)],['vel',(np.float64,3)]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

        n=2
        if 'N' in kwargs.keys(): n = kwargs['N']
        elif (_dat!=None) & (self.size>0): n = self.N[0]

        keys_bin = [['bin'+str(i),BinaryTreeSDAR] for i in range(n-1)]
        DictNpArrayMix.__init__(self, keys_bin, _dat, _offset+self.ncols, True, **kwargs)
            
    def generateBinaryID(self, i):
        """ Use CantorPairing to map two components id to one binary id for one binary group
            Add new member bid into this group

        Parameters
        ----------
           i: int 
              binary group index, counting from 0
        """
        key='bin'+str(i)
        if key in self.__dict__.keys():
            self[key].generateBinaryID()
        else:
            raise ValueError('Error: failed to find ',key)
