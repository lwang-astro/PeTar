from .base import *
from .data import *
from .functions import *

class SDARParticle(SimpleParticle):
    """ SDAR particle
    keys: (class members)
        mass (1D): mass
        pos (2D,3): postion x, y, z
        vel (2D,3): velocity vx, vy, vz
        radius (1D): stellar radius for interruption check
        id (1D): id of particles
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """

        SimpleParticle.__init__(self, _dat, _offset, _append, **kwargs)
        keys = [['radius', np.float64],['id', np.int64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)


class SDARProfile(DictNpArrayMix):
    """ SDAR profile
    Keys: (class members)
        n_step_sum (1D): total number of integration step
        n_step_tsyn_sum (1D): total number of steps for time synchronization 
        n_step (1D): number of integration step every output interval
        n_step_tsyn (1D): number of steps for time synchronization every output interval
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys=[["n_step_sum", np.int64], ["n_step_tsyn_sum", np.int64], ["n_step", np.int64], ["n_step_tsyn", np.int64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    

class SDARData(DictNpArrayMix):
    """ SDAR integrator print column data, used in original SDAR sample code
    Keys: (class members)
        time (1D): current evolved time (counting from zero)
        energy_err (1D): physical energy error
        etot_ref (1D): initial total energy
        ekin (1D): kinetic energy
        epot (1D): potential energy
        gt_drift (1D): time tranformation for drift step
        H (1D): extened phase space Hamiltonian
        energy_interrupt (1D): energy change due to interruption
        H_interrupt (1D): H change due to interruption
        ds (1D): integration step
        time_offset (1D): time offset to obtain the actual time (time_offset + time)
        r_break_crit (1D): distance criterion to break group (used in Hermite)
        profile (SDARProfile): SDAR profile
        particles (ParticleGroup): particle group

    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            member_type: type (SDARParticle)
                Member particle type
            N_particle: int (0)
                Number of particles, determined from file if not provided
        """
        if (not 'member_type' in kwargs.keys()):
            kwargs['member_type'] = SDARParticle
        keys=[['time', np.float64], ['energy_err', np.float64], ["etot_ref",np.float64],["ekin",np.float64],["epot",np.float64], ['gt_drift', np.float64], ['H', np.float64], ['energy_interrupt', np.float64], ['H_interrupt', np.float64], ['ds', np.float64], ['time_offset', np.float64], ['r_break_crit', np.float64], ['profile', SDARProfile], ['particles', ParticleGroup]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class HermiteEnergy(DictNpArrayMix):
    """ Hermite integrator energy data
    keys: (class members)
        error (1D): energy error 
        etot_ref (1D): initial total energy
        ekin (1D): kinetic energy
        epot (1D): potential energy
        epert (1D): perturbation energy
        dE_cum (1D): cumulative energy change
        dE_interrupt (1D): interrupt energy change (binary stellar evolution)
        dE_modify (1D): modify energy change (stellar evolution)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["error",np.float64],["etot_ref",np.float64],["ekin",np.float64],["epot",np.float64],["epert",np.float64],["dE_cum",np.float64],["dE_interrupt",np.float64],["dE_modify",np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
class SlowDown(DictNpArrayMix):
    """ SDAR slowdown data
    Keys: (class members)
        sd (1D): slowdown factor 
        sd_org (1D): slowdown original factor without limit
        sd_max (1D): maximum slowdown factor
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["sd", np.float64], ["sd_org", np.float64], ["sd_max", np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class SlowDownGroup(DictNpArrayMix):
    """ Slowdown data group
    Keys: (class members)
        sd[x] (SlowDown): slowdown data, [x] indicate the indice, counting from 0
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        Parameters
        ----------
        keyword arguments:
            N_sd int (0)
                Number of SlowDown pairs
        """
        keys=[['n', np.int64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

        n_sd = 0
        if 'N_sd' in kwargs.keys(): n_sd = kwargs['N_sd']
        elif (type(_dat)!=type(None)) & (self.size>0): n_sd = self.n[0]

        if (n_sd>0):
            keys_sd = [['sd'+str(i), SlowDown] for i in range(n_sd)]
            DictNpArrayMix.__init__(self, keys_sd, _dat, _offset+self.ncols, True, **kwargs)

class HermiteProfile(DictNpArrayMix):
    """ Hermite profile
    Keys: (class members)
        h4_step_single (1D): single particle total steps
        h4_step_group (1D): AR group total steps
        ar_step (1D): ar total steps
        ar_step_tsyn (1D): ar time synchronize steps
        break_group (1D): number of break groups
        new_group (1D): number of new groups
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys=[["h4_step_single", np.int64], ["h4_step_group", np.int64], ["ar_step", np.int64], ["ar_step_tsyn", np.int64], ["break_group", np.int64], ["new_group", np.int64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)



class HardData(DictNpArrayMix):
    """ Hermite+SDAR integrator print column data, used in petar.hard.debug
    Keys: (class members)
        time (1D): current evolved time (counting from zero)
        energy_phy (HermiteEnergy): physical energy data
        energy_sd (HermiteEnergy): slowdown energy data
        sd (SlowDownGroup): slowdown data
        time_org (1D): original physical time
        profile (HermiteProfile): hermite profile
        particles (ParticleGroup): particle group
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            N_particle: int (0)
                Number of members of one group
            N_sd: int (0)
                Number of slowdown pairs
        """
        kwargs['member_type']=Particle
        kwargs['particle_type']='hermite'
        keys=[['time', np.float64], ['energy_phy', HermiteEnergy], ['energy_sd', HermiteEnergy], ['sd', SlowDownGroup], ['time_org', np.float64], ['profile', HermiteProfile], ['particles', ParticleGroup]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
