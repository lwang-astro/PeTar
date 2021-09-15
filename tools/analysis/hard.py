from .base import *
from .data import *
from .functions import *

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

class HardParticleGroup(DictNpArrayMix):
    """ Hard particle group
    Keys: (class members)
        p[x] (Particle): particle data, [x] indicate the indice, counting from 0
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        Parameters
        ----------
        keyword arguments:
            particle_type: string (hard)
                Particle type, do not change this!
            interrupt_mode: string (none)
                Interrupt mode, base, bse, none
            N_particle int (0)
                Number of particles
        """
        kwargs['particle_type']='hard'
        keys=[['n', np.int64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

        n = 0
        if 'N_particle' in kwargs.keys(): n = kwargs['N_particle']
        elif (type(_dat)!=type(None)) & (self.size>0): n = self.n[0]

        keys_p = [['p'+str(i), Particle] for i in range(n)]
        DictNpArrayMix.__init__(self, keys_p, _dat, _offset+self.ncols, True, **kwargs)


class SDARData(DictNpArrayMix):
    """ Hermite integrator print column data, used in petar.hard.debug
    Keys: (class members)
        time (1D): current evolved time (counting from zero)
        energy_phy (HermiteEnergy): physical energy data
        energy_sd (HermiteEnergy): slowdown energy data
        sd (SlowDownGroup): slowdown data
        time_org (1D): original physical time
        profile (HermiteProfile): hermite profile
        particles (HardParticleGroup): particle group
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            particle_type: string (hard)
                Particle type, do not change this!
            N_particle: int (0)
                Number of members of one group
            N_sd: int (0)
                Number of slowdown pairs
        """
        keys=[['time', np.float64], ['energy_phy', HermiteEnergy], ['energy_sd', HermiteEnergy], ['sd', SlowDownGroup], ['time_org', np.float64], ['profile', HermiteProfile], ['particles', HardParticleGroup]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
