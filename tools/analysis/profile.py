# analysis profile data

from .base import *

class FDPSProfile(DictNpArrayMix):
    """ FDPS time profile for tree 
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["collect_sam_ptcl",1], ["decompose_domain",1], ["exchange_ptcl",1], ["make_local_tree",1], ["make_global_tree",1], ["calc_force",1], ["calc_mom_loc_tree",1], ["calc_mom_gb_tree",1], ["make_LET_1st",1], ["make_LET_2nd",1], ["exchange_LET_1st",1], ["exchange_LET_2nd",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class PeTarProfile(DictNpArrayMix):
    """ PeTar time profile 
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["total",1], ["hard_single",1], ["hard_isolated",1], ["hard_connected",1], ["hard_interrupt",1], ["tree_neighbor",1], ["tree_force",1], ["force_correct",1], ["kick",1], ["search_cluster",1], ["create_group",1], ["domain_decomp",1], ["exchange_ptcl",1], ["output",1], ["status",1],["other",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class GPUProfile(DictNpArrayMix):
    """ GPU time profile 
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["copy",1], ["send",1], ["receive",1], ["calc_force",1], ["n_walk",1], ["n_epi",1], ["n_epj",1], ["n_spj",1], ["n_call",1]]    
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class PeTarCount(DictNpArrayMix):
    """ PeTar number count
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["hard_single",1], ["hard_isolated",1], ["hard_connected",1], ["hard_interrupt",1], ["cluster_isolated",1], ["cluster_connected",1], ["AR_step_sum",1], ["AR_tsyn_step_sum",1], ["AR_group_number",1], ["Hermite_step_sum",1], ["n_neighbor_zero",1], ["Ep-Ep_interaction",1], ["Ep-Sp_interaction",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class Profile(DictNpArrayMix):
    """ Profile class
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        use_gpu=True
        if ('use_gpu' in kwargs.keys()): use_gpu=kwargs['use_gpu']
        if (use_gpu):
            keys = [['rank',1], ['time',1], ['nstep',1], ['n_loc',1], ['comp',PeTarProfile], ['comp_bar', PeTarProfile], ['tree_soft', FDPSProfile], ['tree_nb', FDPSProfile], ['gpu',GPUProfile],['count',PeTarCount]]
            DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        else:
            keys = [['rank',1], ['time',1], ['nstep',1], ['n_loc',1], ['comp',PeTarProfile], ['comp_bar', PeTarProfile], ['tree_soft', FDPSProfile], ['tree_nb', FDPSProfile], ['count',PeTarCount]]
            DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
