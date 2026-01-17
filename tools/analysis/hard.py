from sdar.base import *
from sdar.hermite import *
from sdar.ar import *
from sdar.functions import *
from .data import *

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
            kwargs['member_type'] = HermiteParticle
            kwargs['cm_type'] = HermiteParticle
        elif (kwargs['data_type'] == 'hermite'):
            kwargs['member_type'] = IsolatedHermiteParticle
            kwargs['cm_type'] = IsolatedSDARParticle
        else:
            raise ValueError('data_type is not supported, should be hard or hermite, given ',kwargs['data_type'])

        keys=[['time', np.float64], ['time_offset', np.float64], ['energy_phy', HermiteEnergy], ['energy_sd', HermiteEnergy], ['sd', SlowDownGroup]]
        keys = keys + [['profile', HermiteProfile], ['particles', ParticleGroup]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
