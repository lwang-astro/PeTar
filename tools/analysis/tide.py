# Tide data interface
from .base import *
from .bse import SSEStarParameter
from .group import BinaryTreeSDAR

class Tide(DictNpArrayMix):
    """ Tide event record
    Keys: (class members)
        id1 (1D): particle id of component 1
        id2 (1D): particle id of component 2
        pair_id1 (1D): recorded pair id of component 1 (to indicate whether this is a new two-body system)
        pair_id2 (1D): recorded pair id of component 2 (to indicate whether this is a new two-body system)
        bin_type1 (1D): recorded binary type from component 1 (to indicate what type of binary before tide)
        bin_type2 (1D): recorded binary type from component 2 (to indicate what type of binary before tide)
        poly_type1 (1D): the polynomial type for dynamical tide of component 1, if it is GW, it is zero
        poly_type2 (1D): the polynomial type for dynamical tide of component 2, if it is GW, it is zero
        drdv (1D): two component relative position dot relative velocity 
        semi0 (1D): old semi-major axis 
        ecc0 (1D): old eccentricity
        etid (1D): tide energy change
        ltid (1D): angular momentum change
        bin (BinaryTreeSDAR): binary data with two components after tide effect

    PS: 1) Here pair_id1 == id2, pair_id2 == id1, bin_type1 == bin_type2 if the two body system is not new one.
        2) The two components in bin are type of SSEStarParameter from BSE
        
    """
    
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys=[['id1', np.int64], ['id2', np.int64], ['pair_id1', np.int64], ['pair_id2', np.int64], 
              ['bin_type1', np.int64], ['bin_type2', np.int64], ['poly_type1', np.float64], ['poly_type2', np.float64], 
              ['drdv', np.float64], ['semi0', np.float64], ['ecc0', np.float64], 
              ['etid', np.float64], ['ltid', np.float64], 
              ['bin', BinaryTreeSDAR]]
        
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **{**kwargs, 'member_particle_type':SSEStarParameter})
    
    def generateBinaryID(self):
        """ Use CantorPairing to map two components id to one binary id
            Add new member bid 
        """
        bid = cantorPairing(self.p1.id, self.p2.id)
        self.addNewMember('bid',bid)

