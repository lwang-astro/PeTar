# analysis status data
from .base import *

class Energy(DictNpArrayMix):
    """ Energy data
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["error",1],["error_cum",1],["ekin",1],["epot",1],["etot",1],["dE_modify",1],["dE_interrupt",1],["error_hard",1],["error_hard_cum",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class AngularMomentum(DictNpArrayMix):
    """ Angular momemtum
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["error",1],["error_cum",1],["Lx",1],["Ly",1],["Lz",1],["L",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class Status(DictNpArrayMix):
    """ Status data
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["time",1],["n_real_loc",1],["n_real_glb",1],["n_all_loc",1],["n_all_glb",1],["n_rm_glb",1],["n_esc_glb",1],
                ["energy",Energy], ["energy_sd",Energy], ["angular_momentum",AngularMomentum],
                ["CM_mass",1],["CM_pos",3],["CM_vel",3]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
