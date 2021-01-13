import collections
from .base import *
from .data import *

class BinaryTree(DictNpArrayMix):
    """ Binary tree data output from SDAR and PeTar
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
            member_particle_type: type (Particle)
                Type of component particle 
        """
        member_particle_type=Particle
        if 'member_particle_type' in kwargs.keys(): member_particle_type=kwargs['member_particle_type']

        keys = [['semi',np.float64], ['ecc',np.float64], ['incline',np.float64],['rot_horizon',np.float64],['rot_self',np.float64],['t_peri',np.float64],['period',np.float64],['ecca',np.float64],['m1',np.float64],['m2',np.float64],['r',np.float64],['am',(np.float64,3)],['stab',np.float64],['sd',np.float64],['sd_org',np.float64],['sd_max',np.float64],['p1',member_particle_type],['p2',member_particle_type]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append,**kwargs)
        
class GroupInfo(DictNpArrayMix):
    """ Group information output from PeTar
    Keys: (class members)
        The keys contain n groups with the name of 'binX' where 'X' is replaced by the group index starting from 0
        Each group is the type of BinaryTree.
        n is determined in the keywords argument in initial function
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            particle_type: string (hard)
                Particle type, do not change this!
            N: int (2)
                Number of members of one group
        """
        kwargs['particle_type']='hard'
        keys=[['type',np.int64],['n',np.int64],['time',np.float64],['pos',(np.float64,3)],['vel',(np.float64,3)]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

        n=2
        if 'N' in kwargs.keys(): n = kwargs['N']
        elif (_dat!=None) & (self.size>0): n = self.N[0]

        keys_bin = [['bin'+str(i),BinaryTree] for i in range(n-1)]
        DictNpArrayMix.__init__(self, keys_bin, _dat, _offset+self.ncols, True, **kwargs)
            
