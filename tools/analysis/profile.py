# analysis profile data

from .base import *

class FDPSProfile(DictNpArrayMix):
    """ FDPS time profile for tree 
    """
    @InitialDictNpArrayMixMethod
    def __init__(self):
        return [["collect_sam_ptcl",1], ["decompose_domain",1], ["exchange_ptcl",1], ["make_local_tree",1], ["make_global_tree",1], ["calc_force",1], ["calc_mom_loc_tree",1], ["calc_mom_gb_tree",1], ["make_LET_1st",1], ["make_LET_2nd",1], ["exchange_LET_1st",1], ["exchange_LET_2nd",1]]

class PeTarProfile(DictNpArrayMix):
    """ PeTar time profile 
    """
    @InitialDictNpArrayMixMethod
    def __init__(self):
        return [["total",1], ["hard_single",1], ["hard_isolated",1], ["hard_connected",1], ["hard_interrupt",1], ["tree_neighbor",1], ["tree_force",1], ["force_correct",1], ["kick",1], ["search_cluster",1], ["create_group",1], ["domain_decomp",1], ["exchange_ptcl",1], ["output",1]]

class GPUProfile(DictNpArrayMix):
    """ GPU time profile 
    """
    @InitialDictNpArrayMixMethod
    def __init__(self):
        return [["copy",1], ["send",1], ["receive",1], ["calc_force",1], ["n_walk",1], ["n_epi",1], ["n_epj",1], ["n_spj",1], ["n_call",1]]    

class PeTarCount(DictNpArrayMix):
    """ PeTar number count
    """
    @InitialDictNpArrayMixMethod
    def __init__(self):
        return [["hard_single",1], ["hard_isolated",1], ["hard_connected",1], ["hard_interrupt",1], ["cluster_isolated",1], ["cluster_connected",1], ["AR_step_sum",1], ["AR_tsyn_step_sum",1], ["AR_group_number",1], ["Hermite_step_sum",1], ["Ep-Ep_interaction",1], ["Ep-Sp_interaction",1]]

class Profile(DictNpArrayMix):
    """ Profile class
    """
    def __init__ (self, _dat=0, use_gpu=True, print_tree_nb=True):
        """
        _dat: np.ndarray type data reading from profile data or Profile type data
        """
        if (isinstance(_dat, Profile)):
            self = _dat.copy()
        elif (type(_dat)==np.ndarray):
            self.rank = _dat[:,0]
            self.time = _dat[:,1]
            self.nstep= _dat[:,2]
            self.n_loc= _dat[:,3]
            icol = 4
            self.comp = PeTarProfile(_dat[:,icol:])
            icol += self.comp.ncols
            self.comp_bar= PeTarProfile(_dat[:,icol:])
            icol += self.comp_bar.ncols
            self.tree_soft = FDPSProfile(_dat[:,icol:])
            icol += self.tree_soft.ncols
            if (print_tree_nb):
                self.tree_nb = FDPSProfile(_dat[:,icol:])
                icol += self.tree_nb.ncols
            if (use_gpu):
                self.gpu = GPUProfile(_dat[:,icol:])
                icol += self.gpu.ncols
            self.count = PeTarCount(_dat[:,icol:])
            icol += self.count.ncols
            self.ncols = int(icol)
            self.size  = int(_dat.size/icol)
        else:
            raise ValueError('Initial fail, date type should be Profile or np.ndarray, given ',type(_dat))

