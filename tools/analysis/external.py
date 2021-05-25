import numpy as np
from .base import *
from .data import *

def calcCenterPotExt(particles, rsel):
    """ Calculate the averaged external potential at the center of particles 
    Parameters:
    -------------
    rsel: float
       The distance criterion to select particles for calculating external potential 

    Return:
    -----------
    pot_ext: float
       averaged external potential within rsel of the particle system
    """
    rsel2 = rsel*rsel
    csel = particles.r2<rsel2
    msel = particles.mass[csel]
    mstot = msel.sum()
    pot_ext_c = particles.pot_ext[csel]
    pot_ext_cave = (pot_ext_c*msel).sum()/mstot
    return pot_ext_cave

def estimateGalaxyMass(pot_ext, r_gal, G):
    """ Estimate the equivalent mass of galaxy by - pot_ext*r_gal/G
    Parameters:
    -------------
    pot_ext: float
       the external potential of the center of the particle system
    r_gal: float
       the distance between the center of the particle system to the galactic center
    G: float
       gravitational constant

    Return:
    -------------
    M_galaxy: float
       the estimated mass of the galaxy
    """
    M_galaxy = - pot_ext*r_gal/G
    if (M_galaxy<0):
        raise ValueError('External potential is positive! ', pot_ext)
    return M_galaxy

def calcRCutIsolate(rh):
    """ For isolated star clusters, set rcut to 20 * half-mass radius

    Return
    ----------
    rcut: float 
        escaper distance criterion
    """
    return 20*rh

class Tidal(DictNpArrayMix):
    """ tidal radius and potential of the external potential of the particle system
    keys: (class members)
        time (1D): time
        rtid (1D): tidal radius 
        pot (1D): the external potential of the particle system
        mass (1D): total bound mass
        n (1D): total bound number of particles
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['time', np.float64], ['rtid', np.float64], ['pot', np.float64], ['mass', np.float64], ['n', np.int64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
    def calcTidalSphere(self, time, mass, r2, pot_ext, r_gal, G):
        """ calculate tidal radius assuming R_tid = M_system/ (3*M_galaxy)^(1/3) * R_galaxy; 
        and potential of the central position by averaging the pot_ext inside core radius

        Parameters
        ----------
        time: float
            current time
        mass: 1D numpy.ndarray
            masses of particles
        r2: 1D numpy.ndarray
            the distance square to the center of the system
        pot_ext: float
            the external potential of the center of the particle system
        r_gal: float
            the distance between the center of the particle system to the galactic center
        G: float
            gravitational constant

        Return
        ----------
        r_tid: float
            tidal radius
        """
        self.time = np.append(self.time, time)

        mtot = mass.sum()
    
        M_galaxy = estimateGalaxyMass(pot_ext, r_gal, G)
        r_tid = (mtot/(3*M_galaxy))**(1.0/3.0)*r_gal
        rt2 = r_tid*r_tid
        r_tid_old=r_tid*1.2
        rtsel = r2<rt2
        mcut = mass[rtsel]
        r2cut= r2[rtsel]

        while ((r_tid_old-r_tid)/r_tid_old>1e-2):
            r_tid_old = r_tid
            rt2 = r_tid*r_tid
            rtsel = r2cut<rt2
            mcut = mcut[rtsel]
            r2cut = r2cut[rtsel]
            mtot = mcut.sum()
            r_tid = (mtot/(3*M_galaxy))**(1.0/3.0)*r_gal

        self.rtid = np.append(self.rtid, r_tid)
        self.pot  = np.append(self.pot, pot_ext)
        self.mass = np.append(self.mass, mtot)
        self.n = np.append(self.n, mcut.size)
        self.size += 1
    
        return r_tid

