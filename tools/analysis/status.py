# analysis status data
from .base import *

class Energy(DictNpArrayMix):
    """ Energy data
    Keys:
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
        
        keys = [["error",1],["error_cum",1],["ekin",1],["epot",1],["etot",1],["dE_modify",1],["dE_interrupt",1],["error_hard",1],["error_hard_cum",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class AngularMomentum(DictNpArrayMix):
    """ Angular momemtum
    Keys:
        error: current error of angular momemtum value
        error_cum: cumulative error of angular momemtum value
        Lx: x component of angular momemtum
        Ly: y component of angular momemtum
        Lz: z component of angular momemtum
        L: value of angular momentum
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
        keys = [["error",1],["error_cum",1],["Lx",1],["Ly",1],["Lz",1],["L",1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class Status(DictNpArrayMix):
    """ Status data output from PeTar
    Keys:
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
        keys = [["time",1],["n_real_loc",1],["n_real_glb",1],["n_all_loc",1],["n_all_glb",1],["n_rm_glb",1],["n_esc_glb",1],
                ["energy",Energy], ["energy_sd",Energy], ["angular_momentum",AngularMomentum],
                ["CM_mass",1],["CM_pos",3],["CM_vel",3]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

