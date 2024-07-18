from .base import *
from .data import *
from .functions import *

class IsolatedSDARParticle(SimpleParticle):
    """ SDAR particle of isolated SDAR sample code
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
        de (1D): physical energy error
        etot_ref (1D): initial total energy
        ekin (1D): kinetic energy
        epot (1D): potential energy
        gt_drift (1D): time tranformation for drift step
        H (1D): extened phase space Hamiltonian
        de_interrupt (1D): energy change due to interruption
        dH_interrupt (1D): H change due to interruption
        if (keyword arguemnt 'data_type' == 'hard'): neighbor information defined in Hermite group (not used in isolated group case)
            r_min_index (1D):    nearest neighbor index for each ptcl  
            mass_min_index (1D): mimimum mass in neighbors         
            r_min_sq (1D):       nearest neighbor distance square      
            r_min_mass (1D):     nearest neighbor index for each ptcl 
            mass_min (1D):       mimimum mass in neighbors             
            r_neighbor_crit_sq (1D): neighbor radius criterion
            need_resolve_flag (1D): indicate whether the members need to be resolved for outside 
            initial_step_flag (1D): indicate whether the time step need to be initialized due to the change of neighbors
            n_neighbor_group (1D): number of group neighbor
            n_neighbor_single (1D): number of single neighbor
        ds (1D): integration step
        time_offset (1D): time offset to obtain the actual time (time_offset + time)
        r_break_crit (1D): distance criterion to break group (used in Hermite)
        profile (SDARProfile): SDAR profile
        if (keyword argument 'slowdown' == True):
            de_sd (1D): slowdown energy error
            etot_sd (1D): slowdown energy
            ekin_sd (1D): slowdown kinetic energy
            epot_sd (1D): slowdown potential energy
            de_sd_change (1D): slowdown energy change
            dH_sd_change (1D): slowdown H change
            de_sd_interrupt (1D): slowdown energy change due to interruption
            dH_sd_interrupt (1D): slowdown H change due to interruption
            sd (SlowDownGroup): slowdown data
        particles (ParticleGroup): particle group
            if data_type=='hard', member_type is 'petar.HardParticle'
            if data_type=='sdar', member_type is 'petar.IsolatedSDARParticle'.

    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            data_type: str (hard) 
                hard: PeTar hard particle type
                      The member_type of the particle group is 'petar.HardParticle';
                      The cm_type of the particle group is 'petar.HermiteParticle';
                      Add columns defined in Hermite neighbor, but not used in an isolated SDAR group
                sdar: isolated SDAR sample particle type
                      The member_type of the particle group is 'petar.IsolatedSDARParticle'
            N_particle: int (0)
                Number of particles, determined from file if not provided
            slowdown: bool (False)
                if True, add slowdown keys
            N_sd: int (0)
                Number of slowdown pair, used when slowdown='on'
        """
        key_add = []
        if (not 'data_type' in kwargs.keys()):
            kwargs['data_type'] = 'hard'

        if (kwargs['data_type'] == 'hard'):
            kwargs['member_type'] = HardParticle
            kwargs['cm_type'] = HermiteParticle
            key_add = [['r_min_index', np.int64], ['mass_min_index', np.int64], ['r_min_sq', np.float64], ['r_min_mass', np.float64], ['mass_min', np.float64], ['r_neighbor_crit_sq', np.float64], ['need_resolve_flag', bool], ['initial_step_flag', bool], ['n_neighbor_group', np.int64], ['n_neighbor_single', np.int64]]
        elif (kwargs['data_type'] == 'sdar'):
            kwargs['member_type'] = IsolatedSDARParticle
            kwargs['cm_type'] = IsolatedSDARParticle
        else:
            raise ValueError('data_type is not supported, should be hard or sdar, given ',kwargs['data_type'])

        keys=[['time', np.float64], ['de', np.float64], ["etot_ref",np.float64],["ekin",np.float64],["epot",np.float64], ['gt_drift', np.float64], ['H', np.float64], ['de_interrupt', np.float64], ['dH_interrupt', np.float64]]
        keys = keys + key_add + [['ds', np.float64], ['time_offset', np.float64], ['r_break_crit', np.float64], ['profile', SDARProfile]]
        if ('slowdown' in kwargs.keys()):
            if (kwargs['slowdown']):
                keys = keys + [['de_sd', np.float64], ["etot_sd",np.float64],["ekin_sd",np.float64],["epot_sd",np.float64], ['de_sd_change', np.float64], ['dH_sd_change', np.float64], ['de_sd_interrupt', np.float64], ['dH_sd_interrupt', np.float64], ['sd', (SlowDownGroup, {'with_indices':True})]]
        keys = keys + [['particles', ParticleGroup]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SDARBinary(DictNpArrayMix):
    """ SDAR binary parameter
    keys: (class members)
        semi (1D): semi-major axis
        ecc (1D): eccentricity
        incline (1D): inclination angle
        rot_horizon (1D): frame rotational angle in x-y plane (longitude of ascending node)
        rot_self (1D): frame rotational angle in orbital plane (argument of periapsis)
        t_peri (1D) : time to peri-center
        period (1D): period
        ecca (1D): eccentric anomaly (-pi, pi)
        m1   (1D): component 1 mass
        m2   (1D): component 2 mass
        rrel (1D): relative distance
        am   (2D,3): specific angular momemtum x, y, z
        stab (1D): stability factor 
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        
        keys = [['semi',np.float64],['ecc', np.float64],['incline',np.float64],['rot_horizon',np.float64],['rot_self',np.float64],['t_peri',np.float64],['period',np.float64],['ecca',np.float64],['m1',np.float64],['m2',np.float64],['rrel',np.float64],['am',(np.float64,3)],['stab',np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SDARInterruptBinary(DictNpArrayMix):
    """ SDAR interrupt binary parameter
    keys: (class members)
        time_now (1D): current time
        time_end (1D): targeted integration ending time
        status (1D): binary interruption types
        cm (particle_type): binary c.m. parameter
        bin (SDARBinary): binary orbital parameter
        p1 (particle_type): binary component 1
        p2 (particle_type): binary component 2
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        ----------
        keyword arguments:
            particle_type: type (SDARParticle)
                particle data type
        """

        if (not 'particle_type' in kwargs.keys()):
            kwargs['particle_type'] = SDARParticle
        particle_type = kwargs['particle_type']
        
        keys = [['time_now',np.float64],['time_end', np.float64],['status',np.int64],['cm',particle_type],['bin', SDARBinary],['p1',particle_type],['p2',particle_type]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    

class HermiteEnergy(DictNpArrayMix):
    """ Hermite integrator energy data
    keys: (class members)
        de (1D): energy error 
        etot_ref (1D): initial total energy
        ekin (1D): kinetic energy
        epot (1D): potential energy
        epert (1D): perturbation energy
        de_change (1D): cumulative energy change
        de_interrupt (1D): interrupt energy change (binary stellar evolution)
        de_modify (1D): modify energy change (stellar evolution)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["de",np.float64],["etot_ref",np.float64],["ekin",np.float64],["epot",np.float64],["epert",np.float64],["de_change",np.float64],["de_interrupt",np.float64],["de_modify",np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
class SlowDown(DictNpArrayMix):
    """ SDAR slowdown data
    Keys: (class members)
        if keyword argument 'with_indices' == True:
            i1 (1D): index of binary component 1
            i2 (1D): index of binary component 2
        sd (1D): slowdown factor 
        sd_org (1D): slowdown original factor without limit
        sd_max (1D): maximum slowdown factor
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        Parameters
        ----------
        keyword arguments:
            with_indices: bool (False)
                if true, first two keys (columns) are index1 and index2 of slowdown binary components
        """
        keys = []
        if ('with_indices' in kwargs.keys()):
            if kwargs['with_indices']:
                keys = [['i1',np.int64],['i2',np.int64]]
        keys = keys + [["sd", np.float64], ["sd_org", np.float64], ["sd_max", np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SlowDownGroup(DictNpArrayMix):
    """ Slowdown data group
    Keys: (class members)
        n (1D): number of slowdown pairs
        sd[x] (SlowDown): slowdown data, [x] indicate the indice, counting from 0
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        Parameters
        ----------
        keyword arguments:
            N_sd: int (0)
                Number of SlowDown pairs
            with_indices: bool (False)
                if true, first two keys in sd[x] are index1 and index2 of slowdown binary components
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

class IsolatedHermiteParticle(IsolatedSDARParticle):
    """ Hermite Particle
    keys: (class members)
        mass (1D): mass
        pos (2D,3): postion x, y, z
        vel (2D,3): velocity vx, vy, vz
        radius (1D): stellar radius for interruption check
        id (1D): id
        dt (1D): time step
        time (1D): current time
        acc (2D,3): acceleration x, y, z
        jerk (2D,3): acceleration derivative x, y, z
        pot (1D): potential
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """

        IsolatedSDARParticle.__init__(self, _dat, _offset, _append, **kwargs)
        keys_hermite_add = [['dt',np.float64],['time',np.float64],['acc',(np.float64,3)],['jerk',(np.float64,3)],['pot',np.float64]]
        #keys = [['dm', np.float64],['time_check', np.int64],['binary_state',np.int64]]
        DictNpArrayMix.__init__(self, keys_hermite_add, _dat, _offset+self.ncols, True, **kwargs)


class HermiteData(DictNpArrayMix):
    """ Hermite+SDAR integrator print column data, used in petar.hard.debug
    Keys: (class members)
        time (1D): current evolved time (counting from zero)
        time_offset (1D): time offset to calculate the global time (time+time_offset)
        energy_phy (HermiteEnergy): physical energy data
        energy_sd (HermiteEnergy): slowdown energy data
        sd (SlowDownGroup): slowdown data
        profile (HermiteProfile): hermite profile
        particles (ParticleGroup): particle group, 
            if data_type=='hard', member_type is 'petar.HardParticle', 
                                  cm_type is 'petar.HermiteParticle'
            if data_type=='hermite', member_type is 'petar.IsolatedHermiteParticle'.
                                     cm_type is 'petar.IsolatedSDARParticle'
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
            data_type: str (hard) 
                hard: PeTar hard particle type
                      the member_type of the particle group is 'petar.Particle'
                      the particle_type of 'petar.Particle' is 'hermite'
                hermite: isolated Hermite sample particle type
                      the member_type of the particle group is 'petar.HermiteParticle'
        """

        if (not 'data_type' in kwargs.keys()):
            kwargs['data_type'] = 'hard'

        if (kwargs['data_type'] == 'hard'):
            kwargs['member_type'] = HardParticle
            kwargs['cm_type'] = HermiteParticle
        elif (kwargs['data_type'] == 'hermite'):
            kwargs['member_type'] = IsolatedHermiteParticle
            kwargs['cm_type'] = IsolatedSDARParticle
        else:
            raise ValueError('data_type is not supported, should be hard or hermite, given ',kwargs['data_type'])

        keys=[['time', np.float64], ['time_offset', np.float64], ['energy_phy', HermiteEnergy], ['energy_sd', HermiteEnergy], ['sd', SlowDownGroup]]
        keys = keys + [['profile', HermiteProfile], ['particles', ParticleGroup]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
