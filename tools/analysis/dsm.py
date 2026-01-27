# DSM data interface
from sdar.base import *
from sdar.functions import *

class DSMStarParameter(DictNpArrayMix):
    """
    DSMStarParameter is a class that manages the parameters of a star in a DSM (Disk Star Merger model).
    Keys: (class members)
        type (1D): type of the star (0: Black Hole, 1: Star no growth; 2: Star with growth)
        n_merger_star (1D): time of merger with star
        n_merger_bh (1D): time of merger with black hole
        last_mass_change_time (1D): time of last mass change
        last_merger_time (1D): time of last merger
        helium (1D): helium fraction in the star
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['type',np.int64],
                ['n_merger_star',np.int64],
                ['n_merger_bh',np.int64],
                ['last_mass_change_time',np.float64],
                ['last_merger_time',np.float64],
                ['helium',np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
