import collections
from .base import *
from .data import *

class BinaryTree(DictNpArrayMix):
    """ Binary tree
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        member_particle_type=Particle
        if 'member_particle_type' in kwargs.keys(): member_particle_type=kwargs['member_particle_type']

        keys = [['semi',1], ['ecc',1], ['incline',1],['rot_horizon',1],['rot_self',1],['t_peri',1],['period',1],['ecca',1],['m1',1],['m2',1],['r',1],['am',3],['stab',1],['sd',1],['sd_org',1],['sd_max',1],['p1',member_particle_type],['p2',member_particle_type]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append,**kwargs)
        
class GroupInfo(DictNpArrayMix):
    """ Group information
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        kwargs['particle_type']='hard'
        keys=[['type',1],['n',1],['time',1],['pos',3],['vel',3]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

        n=2
        if 'N' in kwargs.keys(): n = kwargs['N']
        elif (_dat!=None) & (self.size>0): n = self.N[0]

        keys_bin = [['bin'+str(i),BinaryTree] for i in range(n-1)]
        DictNpArrayMix.__init__(self, keys_bin, _dat, _offset+self.ncols, True, **kwargs)
            
