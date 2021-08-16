# analysis profile data

from .base import *

class FDPSProfile(DictNpArrayMix):
    """ FDPS time profile for tree for one tree step
    Keys: (class members)
        collect_sam_ptcl (1D): collect sample
        decompose_domain (1D): decompose domains
        exchange_ptcl (1D): exchange particles
        *set_particle_local_tree (1D): set particle in local tree
        *set_particle_global_tree (1D): set particle in local tree
        make_local_tree (1D): make local tree
        make_global_tree (1D): make global tree
        *set_root_cell (1D): set root cell
        calc_force (1D): calculate force
        calc_mom_loc_tree: calculate superparticle momentum in local tree
        calc_mom_gb_tree: calcualte superparticle momentum in global tree
        make_LET_1st: make local essential tree 1st
        make_LET_2nd: make local essential tree 2nd
        exchange_LET_1st: exchange local essential tree 1st
        exchange_LET_2nd: exchange local essential tree 2nd
        *write_back (1D): write back


    PS: the prefix '*" indicates that these items do not exist for the old PeTar version before 984
        Using the keyword argument 'FDPS_version=old' in the initialization for the old version.
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [["collect_sam_ptcl",np.float64], ["decompose_domain",np.float64], ["exchange_ptcl",np.float64], 
                ["set_particle_local_tree",np.float64],["set_particle_global_tree",np.float64],
                ["make_local_tree",np.float64], ["make_global_tree",np.float64], ["set_root_cell",np.float64],
                ["calc_force",np.float64], ["calc_mom_loc_tree",np.float64], ["calc_mom_gb_tree",np.float64],
                ["make_LET_1st",np.float64], ["make_LET_2nd",np.float64], ["exchange_LET_1st",np.float64], ["exchange_LET_2nd",np.float64],
                ["write_back",np.float64]]
        if ('FPDS_version' in kwargs.keys()): 
            if (kwargs['FDPS_version']=='old'):
                keys = [["collect_sam_ptcl",np.float64], ["decompose_domain",np.float64], ["exchange_ptcl",np.float64], 
                        ["make_local_tree",np.float64], ["make_global_tree",np.float64], 
                        ["calc_force",np.float64], ["calc_mom_loc_tree",np.float64], ["calc_mom_gb_tree",np.float64], 
                        ["make_LET_1st",np.float64], ["make_LET_2nd",np.float64], ["exchange_LET_1st",np.float64], ["exchange_LET_2nd",np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class PeTarProfile(DictNpArrayMix):
    """ PeTar computing wallclock time profile for one tree step
    Keys: (class members)
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
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["total",np.float64], ["hard_single",np.float64], ["hard_isolated",np.float64], ["hard_connected",np.float64], ["hard_interrupt",np.float64], ["tree_neighbor",np.float64], ["tree_force",np.float64], ["force_correct",np.float64], ["kick",np.float64], ["search_cluster",np.float64], ["create_group",np.float64], ["domain_decomp",np.float64], ["exchange_ptcl",np.float64], ["output",np.float64], ["status",np.float64],["other",np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class GPUProfile(DictNpArrayMix):
    """ GPU time profile for one tree force calculation
    Keys: (class members)
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
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["copy",np.float64], ["send",np.float64], ["receive",np.float64], ["calc_force",np.float64], ["n_walk",np.int64], ["n_epi",np.int64], ["n_epj",np.int64], ["n_spj",np.int64], ["n_call",np.int64]]    
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class PeTarCount(DictNpArrayMix):
    """ PeTar number count for one tree step
    Keys: (class members)
        hard_single: number of particles in single clusters
        hard_isolated:  number of particles in isolated clusters
        hard_connected:  number of particles in connected clusters
        hard_interrupt: number of clusters suffering interruptions
        cluster_isolated: number of clusters with multiple particles in local MPI process
        cluster_connected: number of clusters with multiple particles crosing multiple MPI processes
        AR_step_sum: total AR steps
        AR_tsyn_step_sum: total AR steps for time synchronization
        AR_group_number: number of AR groups
        iso_group_number: number of isolated AR groups 
        Hermite_step_sum: total Hermite steps
        n_neighbor_zero: particles have zero neighbors in Hermite 
        Ep_Ep_interaction: number of essential (active) i and j particle interactions 
        Ep_Sp_interaction: number of essential (active) i and superparticle interactions
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["hard_single",np.int64], ["hard_isolated",np.int64], ["hard_connected",np.int64], ["hard_interrupt",np.int64], ["cluster_isolated",np.int64], ["cluster_connected",np.int64], ["AR_step_sum",np.int64], ["AR_tsyn_step_sum",np.int64], ["AR_group_number",np.int64], ["iso_group_number",np.int64], ["Hermite_step_sum",np.int64], ["n_neighbor_zero",np.int64], ["Ep_Ep_interaction",np.int64], ["Ep_Sp_interaction",np.int64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class Profile(DictNpArrayMix):
    """ Profile class
    Keys: (class members)
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
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            use_gpu: bool (True)
                whether cuda is used 
        """
        use_gpu=True
        if ('use_gpu' in kwargs.keys()): use_gpu=kwargs['use_gpu']
        if (use_gpu):
            keys = [['rank',np.int64], ['time',np.float64], ['nstep',np.int64], ['n_loc',np.int64], ['comp',PeTarProfile], ['comp_bar', PeTarProfile], ['tree_soft', FDPSProfile], ['tree_nb', FDPSProfile], ['gpu',GPUProfile],['count',PeTarCount]]
            DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        else:
            keys = [['rank',np.int64], ['time',np.float64], ['nstep',np.int64], ['n_loc',np.int64], ['comp',PeTarProfile], ['comp_bar', PeTarProfile], ['tree_soft', FDPSProfile], ['tree_nb', FDPSProfile], ['count',PeTarCount]]
            DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
