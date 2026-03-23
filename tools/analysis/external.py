import numpy as np
from sdar.base import *
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

def calcREscapeIsolate(rh):
    """ For isolated star clusters, set r_escape to 20 * half-mass radius

    Return
    ----------
    r_escape: float 
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
    
    def calcTidalSphere(self, time, mass, r2, M_galaxy, pot_ext, r_gal, G):
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
        M_galaxy: float
            mass of the galaxy
        pot_ext: float
            potential of the galaxy
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
    
        #M_galaxy = estimateGalaxyMass(pot_ext, r_gal, G)
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

HEADER_EXTERNAL_POT_MAP_OFFSET = 16

class ExternalPotMapConfig(DictNpArrayMix):
    """ External potential measurement parameter file
    keys: 
        time0, initial time
        dt: time step
        nstep: total evolution steps
        dt_out: output interval
        xmin, xmax, nx: x-axis mesh range and number of points
        ymin, ymax, ny: y-axis mesh range and number of points
        zmin, zmax, nz: z-axis mesh range and number of points
    """
    def __init__(self, fpar):
        fp = open(fpar, 'r')
        header = fp.readline()
        fp.close()
        t0, dt, nstep, dt_out, xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz = header.split()
        self.t0 = float(t0)
        self.dt = float(dt)
        self.nstep = int(nstep)
        self.dt_out = float(dt_out)
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.nx = int(nx)
        self.ymin = float(ymin)
        self.ymax = float(ymax)
        self.ny = int(ny)
        self.zmin = float(zmin)
        self.zmax = float(zmax)
        self.nz = int(nz)

class ExternalPotMapHeader():
    """ External potential map header information
    keys:
        time: snapshot time
        nx: number of x-axis points
        ny: number of y-axis points
    """
    def read(self, filename, snapshot_format='ascii'):
        dtype_header = np.dtype([('time', np.float64), ('nx', np.int32), ('ny', np.int32)])
        if snapshot_format == 'ascii':
            with open(filename, 'r') as f:
                header = f.readline()
                t0, nx, ny = header.split()
                self.time = float(t0)
                self.nx = int(nx)
                self.ny = int(ny)
        else:
            with open(filename, 'rb') as f:
                header_data = np.fromfile(f, dtype=dtype_header, count=1)
            self.time = header_data['time'][0]
            self.nx = header_data['nx'][0]
            self.ny = header_data['ny'][0]

    def __init__(self, filename=None, snapshot_format='ascii'):
        if filename is not None:
            self.read(filename, snapshot_format)
        else:
            self.time = 0.0
            self.nx = 0
            self.ny = 0

class ExternalPotMap(DictNpArrayMix):
    """ External potential measure point
    keys: (class members)
        mass (1D): mass
        pos (2D,3): postion x, y, z
        vel (2D,3): velocity vx, vy, vz
        acc (2D,3): acceleration ax, ay, az
        pot (1D): potential
        den (1D): density
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['mass', np.float64], ['pos', (np.float64, 3)], 
                ['vel', (np.float64, 3)], ['acc',(np.float64, 3)], 
                ['pot',np.float64], ['den',np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
    def plot(self, axes, header, plot_keys=['x','y'], log_flag=False, with_contour=False, center_offset=None, **kwargs):
        """ 
        Plot the external potential map using pcolormesh and contour
        Parameters:
        -------------
        axes: matplotlib.axes.Axes
            the axes to plot on
        header: ExternalPotMapHeader
            the header information of the potential map, used to reshape the data
        plot_keys: list of str (default: ['x','y'])
            the keys for x and y axis
        log_flag: bool (default: False)
            whether to plot log(-pot) or -pot
        with_contour: bool (default: False)
            whether to add contour lines on top of the pcolormesh
        center_offset: list or tuple of float (default: None)
            if not None, the offset to apply to the positions for plotting, used to center the data around a specific point (e.g., the center of mass of the particle system)
        kwargs: dict
            pcolormesh and contour keyword arguments, such as vmin, vmax, cmap, etc.

        Return:
        -------------
        im: matplotlib.collections.QuadMesh
            the pcolormesh object
        cset: matplotlib.contour.QuadContourSet (if with_contour is True)
            the contour set object
        """

        nx = header.nx
        ny = header.ny
        key_map = {'x': 0, 'y': 1, 'z': 2}
        x_grid = self.pos[:, key_map[plot_keys[0]]].reshape((nx, ny))
        y_grid = self.pos[:, key_map[plot_keys[1]]].reshape((nx, ny))
        pot = self['pot'].reshape((nx, ny))

        if center_offset is not None:
            x_plot = x_grid - center_offset[0]
            y_plot = y_grid - center_offset[1]
        else:
            x_plot = x_grid
            y_plot = y_grid

        count = np.log10(-pot) if log_flag else -pot        

        im = axes.pcolormesh(x_plot, y_plot, count, shading='auto', **kwargs)
        axes.set_aspect('equal', adjustable='box')
        axes.set_xlabel(plot_keys[0])
        axes.set_ylabel(plot_keys[1])

        cset = None
        if (with_contour):
            cset = axes.contour(x_plot, y_plot, count, linewidths=2, **kwargs)
            #axes[i].clabel(cset,inline=True,fmt='%1.1f',fontsize=10)

        return im, cset
