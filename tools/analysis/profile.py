# analysis profile data

from .base import *

class FDPSProfile(DictNpArrayMix):
    """ FDPS time profile for tree for one tree step
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["collect_sam_ptcl",1], ["decompose_domain",1], ["exchange_ptcl",1], ["make_local_tree",1], ["make_global_tree",1], ["calc_force",1], ["calc_mom_loc_tree",1], ["calc_mom_gb_tree",1], ["make_LET_1st",1], ["make_LET_2nd",1], ["exchange_LET_1st",1], ["exchange_LET_2nd",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class PeTarProfile(DictNpArrayMix):
    """ PeTar computing wallclock time profile for one tree step
    Keys:
        total (1D): total time per tree step
        hard_single (1D): short-range integration of clusters with only one particle (pure drift)
        hard_isolated (1D): short-range integration of clusters with multiple particles in local MPI process (Hermite + SDAR)
        hard_connected (1D): short-range integration of clusters with multiple particles crossing multiple MPI processes (Hermite + SDAR; MPI communication)
        hard_interrupt (1D): short-range integration of interrupted clusters
        tree_neighbor (1D): particle-tree construction of n_real and neighbor searching
        tree_force    (1D): particle-tree construction of n_all and tree forace calculattion
        force_correct (1D): force correction for changeover function
        kick (1D): kick particle velocity
        search_cluster (1D): find clusters for short-range interactions
        create_group (1D):  find particle groups and create artificial particles
        domain_decomp (1D): domain decomposition
        exchange_ptcl (1D): exchange particles between MPI processes
        output (1D): output snapshot and data
        status (1D): calculate status of system
        other (1D): other cost
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """
        Parameters
        ----------
        _dat: numpy.ndarray | same class type (None)
            If it is 2D numpy.ndarray type data, read data as readArray function; if it is the same class type, copy the data 
        _offset: int (0)
            Reading column offset of _dat if it is 2D np.ndarray
        _append: bool (False)
            If true, append keys and ncols to the current class instead of create new class members
        kwaygs: dict ()
            keyword arguments

        """
        keys = [["total",1], ["hard_single",1], ["hard_isolated",1], ["hard_connected",1], ["hard_interrupt",1], ["tree_neighbor",1], ["tree_force",1], ["force_correct",1], ["kick",1], ["search_cluster",1], ["create_group",1], ["domain_decomp",1], ["exchange_ptcl",1], ["output",1], ["status",1],["other",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class GPUProfile(DictNpArrayMix):
    """ GPU time profile for one tree force calculation
    Keys:
        copy: copy data for sending 
        send: host to GPU memory sending
        receive: GPU to host memory receiving
        calc_force: GPU force calculation
        n_walk: number of multiple walks (FDPS)
        n_epi: total number of i particles
        n_epj: total number of j particles
        n_spj: total number of super particles
        n_call: number of calls of force kernel function 
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """
        Parameters
        ----------
        _dat: numpy.ndarray | same class type (None)
            If it is 2D numpy.ndarray type data, read data as readArray function; if it is the same class type, copy the data 
        _offset: int (0)
            Reading column offset of _dat if it is 2D np.ndarray
        _append: bool (False)
            If true, append keys and ncols to the current class instead of create new class members
        kwaygs: dict ()
            keyword arguments

        """
        keys = [["copy",1], ["send",1], ["receive",1], ["calc_force",1], ["n_walk",1], ["n_epi",1], ["n_epj",1], ["n_spj",1], ["n_call",1]]    
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class PeTarCount(DictNpArrayMix):
    """ PeTar number count for one tree step
    Keys:
        hard_single: number of particles in single clusters
        hard_isolated:  number of particles in isolated clusters
        hard_connected:  number of particles in connected clusters
        hard_interrupt: number of clusters suffering interruptions
        cluster_isolated: number of clusters with multiple particles in local MPI process
        cluster_connected: number of clusters with multiple particles crosing multiple MPI processes
        AR_step_sum: total AR steps
        AR_tsyn_step_sum: total AR steps for time synchronization
        AR_group_numer: number of AR groups
        Hermite_step_sum: total Hermite steps
        n_neighbor_zero: particles have zero neighbors in Hermite 
        Ep_Ep_interaction: number of essential (active) i and j particle interactions 
        Ep_Sp_interaction: number of essential (active) i and superparticle interactions
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """
        Parameters
        ----------
        _dat: numpy.ndarray | same class type (None)
            If it is 2D numpy.ndarray type data, read data as readArray function; if it is the same class type, copy the data 
        _offset: int (0)
            Reading column offset of _dat if it is 2D np.ndarray
        _append: bool (False)
            If true, append keys and ncols to the current class instead of create new class members
        kwaygs: dict ()
            keyword arguments

        """
        keys = [["hard_single",1], ["hard_isolated",1], ["hard_connected",1], ["hard_interrupt",1], ["cluster_isolated",1], ["cluster_connected",1], ["AR_step_sum",1], ["AR_tsyn_step_sum",1], ["AR_group_number",1], ["Hermite_step_sum",1], ["n_neighbor_zero",1], ["Ep_Ep_interaction",1], ["Ep_Sp_interaction",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class Profile(DictNpArrayMix):
    """ Profile class
    Keys:
        rank (1D): MPI rank
        time (1D): evolved time
        nstep (1D): number of steps per output
        n_loc (1D): number of partiles locally
        comp (PeTarProfile): time profiling for each components of PeTar
        comp_bar (PeTarProfile): MPI barrier waiting time of each components of PeTar
        tree_soft (FDPSProfile): FDPS long-range force particle-tree profile
        tree_nb  (FDPSProfile): FDPS particle-tree for neighbor searching
        if keyword arguments "use_gpu" == True:
            gpu (GPUProfile): GPU profile for tree force calculation
        count (PeTarCount): number counts
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """
        Parameters
        ----------
        _dat: numpy.ndarray | same class type (None)
            If it is 2D numpy.ndarray type data, read data as readArray function; if it is the same class type, copy the data 
        _offset: int (0)
            Reading column offset of _dat if it is 2D np.ndarray
        _append: bool (False)
            If true, append keys and ncols to the current class instead of create new class members
        kwaygs: dict ()
            keyword arguments

        """
        use_gpu=True
        if ('use_gpu' in kwargs.keys()): use_gpu=kwargs['use_gpu']
        if (use_gpu):
            keys = [['rank',1], ['time',1], ['nstep',1], ['n_loc',1], ['comp',PeTarProfile], ['comp_bar', PeTarProfile], ['tree_soft', FDPSProfile], ['tree_nb', FDPSProfile], ['gpu',GPUProfile],['count',PeTarCount]]
            DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        else:
            keys = [['rank',1], ['time',1], ['nstep',1], ['n_loc',1], ['comp',PeTarProfile], ['comp_bar', PeTarProfile], ['tree_soft', FDPSProfile], ['tree_nb', FDPSProfile], ['count',PeTarCount]]
            DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
