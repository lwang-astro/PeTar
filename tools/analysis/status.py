# analysis status data
from .base import *
from .data import *

class Energy(DictNpArrayMix):
    """ Energy data
    Keys: (class members)
        error (1D): current energy error
        error_cum (1D): cumulative energy error
        ekin (1D): kinetic energy
        epot (1D): potential energy
        etot (1D): total energy
        dE_modify (1D): modified energy due to none graivitional processes
        dE_interrupt (1D): modified energy due to interruption of AR integration (e.g. binary stellar evolution, mergers)
        error_hard (1D): current energy error in hard (short-range) integration
        error_hard_cum (1D): cumulative energy error in hard (short-range) integration
        
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["error",np.float64],["error_cum",np.float64],["ekin",np.float64],["epot",np.float64],["etot",np.float64],["dE_modify",np.float64],["dE_interrupt",np.float64],["error_hard",np.float64],["error_hard_cum",np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class AngularMomentum(DictNpArrayMix):
    """ Angular momemtum
    Keys: (class members)
        error: current error of angular momemtum value
        error_cum: cumulative error of angular momemtum value
        Lx: x component of angular momemtum
        Ly: y component of angular momemtum
        Lz: z component of angular momemtum
        L: value of angular momentum
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["error",np.float64],["error_cum",np.float64],["Lx",np.float64],["Ly",np.float64],["Lz",np.float64],["L",np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SoftParticleGroup(DictNpArrayMix):
    """ Soft particle group
    Keys: (class members)
        p[x] (Particle): particle data, [x] indicate the indice, counting from 0
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        Parameters
        ----------
        keyword arguments:
            interrupt_mode: string (none)
                Interrupt mode, base, bse, none
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, agama, none 
            N_particle int (0)
                Number of particles
        """
        n = 0
        if 'N_particle' in kwargs.keys(): n = kwargs['N_particle']
        if (n>0):
            keys_p = [['p'+str(i), Particle] for i in range(n)]
            DictNpArrayMix.__init__(self, keys_p, _dat, _offset, _append, **kwargs)
        else:
            self.ncols=0
            self.size=0

class Status(DictNpArrayMix):
    """ Status data output from PeTar
    Keys: (class members)
        time (1D): current evolved time
        n_real_loc (1D): number of real particles in the local MPI process
        n_real_glb (1D): total number of real particles in all MPI processes
        n_all_loc  (1D): number of all type of particles (real and aritifical) in the local MPI process
        n_all_glb  (1D): number of all type of particles (real and aritifical) in all MPI processes
        n_rm_glb   (1D): number of removed particles in all MPI processes
        n_esc_glb  (1D): number of escaped particles in all MPI processes
        energy    (Energy): physical energy information of the system
        energy_sd (Energy): slowdown energy information of the system
        angular_momentum (AngularMomentum): angular momentum information of the system
        CM_mass (1D): total mass of the system
        CM_pos  (2D,3): c.m. position of the system
        CM_vel  (2D,3): c.m. velocity of the system
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        Parameters:
        ----------
        keyword arguments:
            interrupt_mode: string (none)
                Interrupt mode, base, bse, none
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, agama, none 
            use_mpfrc: bool (False)
                Include three columns of high-precsion parts of particle position x, y, z
            N_particle int (0)
                Number of particles
        """
        keys = [["time",np.float64],["n_real_loc",np.int64],["n_real_glb",np.int64],["n_all_loc",np.int64],["n_all_glb",np.int64],["n_rm_glb",np.int64],["n_esc_glb",np.int64],
                ["energy",Energy], ["energy_sd",Energy], ["angular_momentum",AngularMomentum],
                ["CM_mass",np.float64],["CM_pos",(np.float64,3)],["CM_vel",(np.float64,3)]]

        n = 0
        if 'N_particle' in kwargs.keys(): 
            n = kwargs['N_particle']
            kwargs['N_column_exist'] = False
            kwargs['cm_column_exist'] = False
            kwargs['member_type'] = Particle
        if (n>0): keys += [["particles", ParticleGroup]]
            
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

