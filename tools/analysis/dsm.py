# DSM data interface
from .base import *
from .functions import *

class DSMStarParameter(DictNpArrayMix):
    """
    DSMStarParameter is a class that manages the parameters of a star in a DSM (Disk Star Merger model).
    Keys: (class members)
        type (1D): type of the star (0: Black Hole, 1: Star no growth; 2: Star with growth)
        merger_star_times (1D): time of merger with star
        merger_bh_times (1D): time of merger with black hole
        growth_time_start (1D): time of start of growth
        last_merger_time (1D): time of last merger
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['type',np.int64],
                ['merger_star_times',np.int64],
                ['merger_bh_times',np.int64],
                ['growth_time_start',np.float64],
                ['last_merger_time',np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
