import collections
import numpy as np
from scipy import spatial as sp
from .base import *
from .data import *
from .bse import *

class Core(DictNpArrayMix):
    """ Core of a star cluster
    Keys: (class members)
        time (1D): evolved time of the star cluster
        pos  (2D,3): center of the star cluster
        vel  (2D,3): center velocity of the cluster
        rc   (1D): core radius of the star cluster
    """
    
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys  = [['time',np.float64],['pos', (np.float64,3)],['vel', (np.float64,3)], ['rc', np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def calcPotentialCenter(self, single, binary, G):
        """ Calculate potential center of the system and return the result
        r_cm = sum_i {pot_i *r_i} /sum_i pot_i (only count pot_i <0)

        Parameters
        ----------
        single: inhermited SimpleParticle
            Single particle data
        binary: Binary
            Binary particle data
        G: float
            Gravitational constant

        Return
        ----------
        cm_pos: c.m. position, numpy.ndarray with shape (*,3)
        cm_vel: c.m. velocity, numpy.ndarray with shape (*,3)

        """
        pot_s_sel = (single.pot<0)
        pot_s = single.pot[pot_s_sel]
        pos_s = single.pos[pot_s_sel]
        vel_s = single.vel[pot_s_sel]

        binary.calcPot(G)
        pot_b = binary.pot
        pos_b = binary.pos
        vel_b = binary.vel
        
        pot_sum = pot_s.sum() + pot_b.sum()
        cm_pos = np.array([(np.sum(pot_s*pos_s[:,i])+np.sum(pot_b*pos_b[:,i]))/pot_sum for i in range(3)])
        cm_vel = np.array([(np.sum(pot_s*vel_s[:,i])+np.sum(pot_b*vel_b[:,i]))/pot_sum for i in range(3)])

        self.pos = np.append(self.pos,[cm_pos],axis=0)
        self.vel = np.append(self.vel,[cm_vel],axis=0)

        return cm_pos, cm_vel
    
    def calcDensityAndCenter(self, particle, kdtree):
        """ Calculate density based on nearest six neighbors and return the result (Casertano & Hut 1985)
        
        Parameters
        ----------
        particle: inherited SimpleParticle
            Particle data set
        kdtree: scipy.spatial.cKDtree
            3D KDTree of all particle positions, can be generated from findPair

        Return
        ----------
        cm_pos: c.m. position, numpy.ndarray with shape (*,3)
        cm_vel: c.m. velocity, numpy.ndarray with shape (*,3)
        """
        # 6 nearest neighbors
        nb_r_list6, nb_index_list6 = kdtree.query(particle.pos,k=7)
        nb_mass_tot6=np.sum(particle.mass[nb_index_list6[:,:-1]],axis=1)
        
        nb_inv_r6 = 1/nb_r_list6[:,-1]
        rho = nb_mass_tot6*(nb_inv_r6*nb_inv_r6*nb_inv_r6)
        particle.addNewMember('density',rho)
        rho_tot = rho.sum()
     
        cm_pos = np.array([np.sum(rho*particle.pos[:,i])/rho_tot for i in range(3)])
        cm_vel = np.array([np.sum(rho*particle.vel[:,i])/rho_tot for i in range(3)])
        self.pos = np.append(self.pos, [cm_pos], axis=0)
        self.vel = np.append(self.vel, [cm_vel], axis=0)

        return cm_pos, cm_vel

    def calcCoreRadius(self, particle):
        """ Calculate core radius, using Casertano & Hut (1985) method with the squares of densities, also used in Nbody5 and Nbody6 (Aarseth 2003):
        rc = sqrt(sum_i rho_i^2 r_i^2 / sum_i rho_i^2)

        Parameters
        ----------
        particle: inherited SimpleParticle
            Particle data set

        Return
        ----------
        rc: core radius

        """
        rho2 = particle.density*particle.density
        rc = np.sqrt((particle.r2*rho2).sum()/(rho2.sum()))
        self.rc = np.append(self.rc, rc)

        return rc

    def addTime(self, time):
        """ Append a new time to current member 'time'
        """
        self.time = np.append(self.time, time)

    def toSkyCoord(self, **kwargs):
        """ generate astropy.coordinates.SkyCoord data in galactocentric frame

        Parameters
        -----------------
        kwargs: dict()
            pos_unit: astropy.units (units.pc)
                 position unit of the particle data
            vel_unit: astropy.units (units.pc/units.Myr)
                 velocity unit of the particle data
            galcen_distance: floating with length units (8.0*units.kpc [Galpy])
                 galactic central distance of the Sun
            z_sun: floating with length units (15.0*units.pc [Galpy])
                 z direction distance of the Sun
            galcen_v_sun: astropy.coordinates.CartesianDifferential ([10.0, 235.0, 7.0]*units.km/units.s [Galpy])
                 velocity of the Sun

        Return
        ----------------
        core_g: astropy.coordinates.SkyCoord
            core c.m. data using SkyCoord
        """
        import astropy 
        from astropy.coordinates import SkyCoord  # High-level coordinates
        from astropy.coordinates import ICRS, Galactic, Galactocentric, FK4, FK5  # Low-level frames
        from astropy.coordinates import Angle, Latitude, Longitude  # Angles
        from astropy.coordinates import CartesianDifferential
        import astropy.units as u

        pos_unit = u.pc
        if ('pos_unit' in kwargs.keys()): pos_unit = kwargs['pos_unit']
        vel_unit = u.pc/u.Myr
        if ('vel_unit' in kwargs.keys()): vel_unit = kwargs['vel_unit']

        parameters={'galcen_distance':8.0*u.kpc, 'z_sun':15.*u.pc, 'galcen_v_sun':CartesianDifferential([10.0,235.,7.]*u.km/u.s)}
        for key in parameters.keys():
            if key in kwargs.keys():
                parameter[key] = kwargs[key]

        core_g = SkyCoord(x=self.pos[:,0]*pos_unit, 
                          y=self.pos[:,1]*pos_unit, 
                          z=self.pos[:,2]*pos_unit, 
                          v_x=self.vel[:,0]*vel_unit,
                          v_y=self.vel[:,1]*vel_unit,
                          v_z=self.vel[:,2]*vel_unit,
                          frame='galactocentric', representation_type='cartesian', **parameters)
        return core_g

class LagrangianVelocity(DictNpArrayMix):
    """ Lagrangian velocity component
        Each velocity member is a 2D numpy.ndarray, the row contain the values corresponding to each Langragian radius and the core radius (last value).
        For example, in default case, mass_fraction array is np.array([0.1, 0.3, 0.5, 0.7, 0.9]). 
        The member 'abs' has the shape (*,6), each row contain the velocity referring to the 0.1, 0.3, 0.5, 0.7, 0.9 of Langragian radii and the core radius (last value).
    
    Keys: (class members)
        abs (2D,n_frac): absolute 3D velocity value
        x   (2D,n_frac): velocity in x axis
        y   (2D,n_frac): velocity in y axis
        z   (2D,n_frac): velocity in z axis
        rad (2D,n_frac): velocity in radial direction
        tan (2D,n_frac): velocity in tangential direction
        rot (2D,n_frac): rotational velocity in x-y plane

        n_frac: number of Lagrangian radii + 1, determined by keyword argument 'mass_fraction'
    """
    
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            mass_fraction: 1D numpy.ndarray (np.array([0.1, 0.3, 0.5, 0.7, 0.9]))
                An array to indicate the mass fractions to calculate lagrangian radii.
        """
        
        m_frac=np.array([0.1,0.3,0.5,0.7,0.9])
        if ('mass_fraction' in kwargs.keys()): m_frac=kwargs['mass_fraction'].copy()
        n_frac = m_frac.size + 1
        keys  = [['abs',(np.float64,n_frac)],['x',(np.float64,n_frac)],['y',(np.float64,n_frac)],['z',(np.float64,n_frac)],['rad',(np.float64,n_frac)],['tan',(np.float64,n_frac)],['rot',(np.float64,n_frac)]] # all, x, y, z, radial, tangential, rotational
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        self.initargs['mass_fraction'] = m_frac

class Lagrangian(DictNpArrayMix):
    """ Lagrangian parameters
    Keys: (class members)
        r (2D,n_frac): Lagrangian radii
        m (2D,n_frac): average mass referring to Lagrangian radii
        n (2D,n_frac): number of particles referring to Lagrangian radii
        vel   (LagrangianVelocity): average velocity referring to Lagrangian radii
        sigma (LagrangianVelocity): velocity dispersion referring to Lagrangian radii

        Additional members ff the keyword argument 'calc_energy' is true:
          epot (2D, n_frac): potential energy
          vr (2D, n_frac); virial ratio -(2*ekin/epot); notice that here the external potential energy is reduced by half.
          *epot_ext (2D, n_frac): external potential energy (exists if the keyword argument 'external_mode' is switch on)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            mass_fraction: 1D numpy.ndarray (np.array([0.1, 0.3, 0.5, 0.7, 0.9]))
                An array to indicate the mass fractions to calculate lagrangian radii.
            calc_energy: bool (False)
                If true, add kinetic, potential energies and virial ratio calculation
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, none 
                If it is not none, epot_ext will be added if calc_energy=True
            calc_multi_rc: bool (False)
                If true, calculate core radius using KDTree and correct the central position
        """
        m_frac=np.array([0.1,0.3,0.5,0.7,0.9])
        if ('mass_fraction' in kwargs.keys()): m_frac=kwargs['mass_fraction'].copy()
        calc_energy=False
        if ('calc_energy' in kwargs.keys()): calc_energy=kwargs['calc_energy']
        external_mode='none'
        if ('external_mode' in kwargs.keys()): external_mode=kwargs['external_mode']
        calc_multi_rc=False
        if ('calc_multi_rc' in kwargs.keys()): calc_multi_rc=kwargs['calc_multi_rc']
        
        n_frac = m_frac.size + 1
        keys = [['r', (np.float64, n_frac)],['m', (np.float64, n_frac)],['n', (np.int64, n_frac)], ['vel', LagrangianVelocity], ['sigma', LagrangianVelocity]]
        if (calc_energy):
            keys = keys + [['epot',(np.float64, n_frac)], ['vr',(np.float64, n_frac)]]
            if (external_mode!='none'):
                keys = keys + [['epot_ext',(np.float64, n_frac)]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        self.initargs['mass_fraction'] = m_frac
        self.initargs['calc_energy'] = calc_energy
        self.initargs['external_mode'] = external_mode
        self.initargs['calc_multi_rc'] = calc_multi_rc

        #keys = [['r', (np.float64, n_frac)],['m', (np.float64, n_frac)],['n', (np.float64, n_frac)]] # radius, mass, number
        #DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        #self.vel  = LagrangianVelocity(_dat, _offset+self.ncols, False, **kwargs)
        #self.ncols += self.vel.ncols
        #self.keys.append(['vel',LagrangianVelocity])
        #self.sigma= LagrangianVelocity(_dat, _offset+self.ncols, False, **kwargs)
        #self.ncols += self.sigma.ncols
        #self.keys.append(['sigma',LagrangianVelocity])
        #self.initargs['mass_fraction'] = m_frac

    def calcTrh(self, G, gamma=0.02, mode='sphere'):
        """ Calculate Spitzer one-component half-mass relaxation time
            Trh = 0.138 N^0.5 Rh^1.5 /( G^0.5 m^0.5 ln(gamma N))
            Then add the result as a new member 'trh'.
            Notice mass fraction must contain 0.5 to use this function

        Parameters
        ----------
        G: float
           Gravitational constant
        gamma: float (0.02 # Giersz M., Heggie D. C., 1996, MNRAS, 279, 1037)
           The coefficient for Coulomb logarithm
        mode: string (sphere)
            sphere: calculate averaged properties from center to Lagrangian radii
                   m use the average mass within half-mass radius
                   N use the number of objects within half-mass radius
            shell: calculate properties between two neighbor Lagrangian radii

        """
        rhindex=np.where(self.initargs['mass_fraction']==0.5)[0][0]
        if (mode=='shell'): 
            n=self.n[:,0:(rhindex+1)].sum(axis=1)
            m=(self.m[:,0:(rhindex+1)]*self.n[:,0:(rhindex+1)]).sum(axis=1)/n
            trh=calcTrh(n, self.r[:,rhindex], m, G, gamma)
            self.addNewMember('trh',trh)
        else:
            trh=calcTrh(self.n[:,rhindex], self.r[:,rhindex], self.m[:,rhindex], G, gamma)
            self.addNewMember('trh',trh)

    def calcTcr(self, G, mode='sphere'):
        """ Calculate half-mass crossing time
            Tcr = Rh^1.5/sqrt(G M)
            Then add the result as a new member 'tcr'.

            Notice mass fraction must contain 0.5 to use this function

        Parameters
        ----------
        G: float
           Gravitational constant
        mode: string (sphere)
            sphere: calculate averaged properties from center to Lagrangian radii
            shell: calculate properties between two neighbor Lagrangian radii

        """
        rhindex=np.where(self.initargs['mass_fraction']==0.5)[0][0]
        if (mode=='shell'): 
            M=(self.m[:,0:(rhindex+1)]*self.n[:,0:(rhindex+1)]).sum(axis=1)*2
            tcr=calcTcr(M, self.r[:,rhindex], G)
            self.addNewMember('tcr',tcr)
        else:
            tcr=calcTcr(self.n[:,rhindex]*self.m[:,rhindex]*2, self.r[:,rhindex], G)
            self.addNewMember('tcr',tcr)
            

    def calcOneSnapshot(self, _particle, _rc, _mode='sphere', **kwargs):
        """ Calculate one snapshot lagrangian parameters

        Parameters
        ----------
        _particle: inherited SimpleParticle
            particles data set sorted by distance to the coordinate center
        _rc: float
            core radius
        _mode: string (sphere)
            sphere: calculate averaged properties from center to Lagrangian radii
            shell: calculate properties between two neighbor Lagrangian radii
        ----------
        keyword arguments:
            read_rlagr: 1D np.array([])
                 When this option is used, Lagragian radii are not calculated,
                 but are given by the argument (an array that contain Lagragian radii).
                 Other properties (e.g. average mass, velocity, disperion) are calculated based on these.
        ----------
        return
            rlagr: 1D np.array
                Lagrangian radii of current data
        """
        shell_mode = True if (_mode == 'shell') else False

        mass_fraction = self.initargs['mass_fraction']
        calc_energy = self.initargs['calc_energy']
        external_mode = self.initargs['external_mode']
        calc_multi_rc = self.initargs['calc_multi_rc']
        n_frac = mass_fraction.size+1

        rc = _rc
        if calc_multi_rc:
            if (_particle.size>6):
                core = Core()
                kdtree = sp.cKDTree(_particle.pos)
                cm_pos, cm_vel = core.calcDensityAndCenter(_particle, kdtree)
                _particle.correctCenter(cm_pos, cm_vel)
                _particle.calcR2()
                rc = core.calcCoreRadius(_particle)
            else:
                rc = 0
                

        def find_mass_index(_mass_cum,mass_fraction):
            mass_bins=np.append(0,mass_fraction*_mass_cum[-1])
            index,count=np.histogram(_mass_cum,bins=mass_bins)
            index=index.cumsum()
            index[index>=_mass_cum.size]=_mass_cum.size-1
            return index

        if (_particle.size<=1):
            size = self.size
            empty = Lagrangian(np.zeros([1,n_frac*self.ncols]),**self.initargs)
            self.append(empty)
            if (size+1!=self.size):
                raise ValueError('Size should increase one, but increase', self.size-size)
            return np.zeros(mass_fraction.size)
        else:
            self.size += 1
            self.vel.size += 1
            self.sigma.size += 1
            r = np.sqrt(_particle.r2)
            rindex=np.array([])
            mcum=_particle.mass.cumsum()
            if ('read_rlagr' in kwargs.keys()):
                rbin = np.append(0,kwargs['read_rlagr'])
                rindex, e = np.histogram(r, bins=rbin)
                rindex=rindex.cumsum()
                rindex[rindex>=r.size]=r.size-1
            else:
                rindex= find_mass_index(mcum, mass_fraction)
            rlagr = r[rindex]
            if(len(self.r.shape)!=2):
                raise ValueError('r shape is wrong',self.r.shape)
            self.r = np.append(self.r, [np.append(rlagr,rc)], axis=0)
            nlagr = rindex+1
            if (shell_mode): nlagr[1:] -= nlagr[:-1]
            nc = (_particle.r2<(rc*rc)).sum()
            self.n = np.append(self.n, [np.append(nlagr,nc)], axis=0)
            mlagr = mcum[rindex]
            if (shell_mode): mlagr[1:] -= mlagr[:-1]
            sel = nlagr>0
            mlagr[sel] /= nlagr[sel]
            if (nc>0): mc = mcum[nc-1]/nc
            else:      mc = 0.0
            self.m = np.append(self.m, [np.append(mlagr,mc)], axis=0)

            m   =_particle.mass 
            vel = _particle.vel
            pos = _particle.pos
            rx = pos[:,0]
            ry = pos[:,1]
            rz = pos[:,2]
            vx = vel[:,0]
            vy = vel[:,1]
            vz = vel[:,2]

            # x-y plane radial velocity * rxy
            rvxy = rx*vx + ry*vy
            # radial velocity value
            vr   = (rvxy + rz*vz)/r
            # tangential velocity vector
            vt = [None]*3
            vt[0] = vx - vr*rx/r
            vt[1] = vy - vr*ry/r
            vt[2] = vz - vr*rz/r
            # x-y plane radial position square
            rxy2 = rx*rx + ry*ry
            # rotational velocity
            vrotx = vx - rvxy*rx/rxy2
            vroty = vy - rvxy*ry/rxy2
            vrot = np.sqrt(vrotx*vrotx+vroty*vroty)
            # rotational direction sign
            vrotd = vrotx*ry - vroty*rx
            ng_sig = (vrotd<0.0)
            vrot[ng_sig] = -vrot[ng_sig]
            
            n_offset = np.append(np.zeros(1), rindex+1).astype(int)
            vlst = [vx, vy, vz, vr, vt[0], vt[1], vt[2], vrot]
            vave = [None]*len(vlst)
            for k in range(len(vlst)):
                vlagr = []
                if (shell_mode):
                    vlagr = [np.average(m[n_offset[i]:n_offset[i+1]]*vlst[k][n_offset[i]:n_offset[i+1]])/mlagr[i] if (n_offset[i]<n_offset[i+1]) else 0.0 for i in range(mlagr.size)]
                else:
                    vlagr = [np.average(m[:n_offset[i+1]]*vlst[k][:n_offset[i+1]])/mlagr[i] if (n_offset[i+1]>0) else 0.0 for i in range(mlagr.size)]
                # core radius
                vlagr.append(np.average(m[0:nc]*vlst[k][0:nc])/mc if (nc>0) else 0.0)
                vave[k] = np.array(vlagr)
            
            self.vel.x = np.append(self.vel.x, [vave[0]], axis=0)
            self.vel.y = np.append(self.vel.y, [vave[1]], axis=0)
            self.vel.z = np.append(self.vel.z, [vave[2]], axis=0)
            self.vel.abs= np.append(self.vel.abs, [np.sqrt(vave[0]*vave[0]+vave[1]*vave[1]+vave[2]*vave[2])], axis=0)
            self.vel.rad = np.append(self.vel.rad, [vave[3]], axis=0)
            self.vel.tan = np.append(self.vel.tan, [np.sqrt(vave[4]*vave[4]+vave[5]*vave[5]+vave[6]*vave[6])], axis=0)
            self.vel.rot = np.append(self.vel.rot, [vave[7]], axis=0)
            
            sigma = [None]*len(vlst)
            for k in range(len(vlst)):
                slagr = []
                if (shell_mode):
                    slagr = [np.average(m[n_offset[i]:n_offset[i+1]] * (vlst[k][n_offset[i]:n_offset[i+1]] - vave[k][i])**2) / mlagr[i] if (n_offset[i]<n_offset[i+1]) else 0.0 for i in range(mass_fraction.size)]
                else:
                    slagr = [np.average(m[:n_offset[i+1]] * (vlst[k][:n_offset[i+1]] - vave[k][i])**2) / mlagr[i] if (n_offset[i+1]>0) else 0.0 for i in range(mass_fraction.size)]
                # core radius
                slagr.append(np.average(m[0:nc] * (vlst[k][0:nc] - vave[k][-1])**2) / mc if (nc>0) else 0.0)
                sigma[k] = np.array(slagr)

            sigma_all = sigma[0]+sigma[1]+sigma[2]
            self.sigma.x = np.append(self.sigma.x, [np.sqrt(sigma[0])], axis=0)
            self.sigma.y = np.append(self.sigma.y, [np.sqrt(sigma[1])], axis=0)
            self.sigma.z = np.append(self.sigma.z, [np.sqrt(sigma[2])], axis=0)
            self.sigma.abs= np.append(self.sigma.abs, [np.sqrt(sigma_all)], axis=0)
            self.sigma.rad = np.append(self.sigma.rad, [np.sqrt(sigma[3])], axis=0)
            self.sigma.tan = np.append(self.sigma.tan, [np.sqrt(sigma[4]+sigma[5]+sigma[6])], axis=0)
            self.sigma.rot = np.append(self.sigma.rot, [np.sqrt(sigma[7])], axis=0)

            if (calc_energy):
                epot = []
                if (shell_mode):
                    epot = [np.sum(_particle.pot[n_offset[i]:n_offset[i+1]] * _particle.mass[n_offset[i]:n_offset[i+1]]) if (n_offset[i]<n_offset[i+1]) else 0.0 for i in range(mass_fraction.size)]

                else:
                    epot = [np.sum(_particle.pot[:n_offset[i+1]] * _particle.mass[:n_offset[i+1]]) if (n_offset[i+1]>0) else 0.0 for i in range(mass_fraction.size)]

                epot.append(np.sum(_particle.pot[:nc]*_particle.mass[:nc]  if (nc>0) else 0.0))
                ekin = sigma_all*np.append(mlagr,mc)*np.append(nlagr,nc)
                epot = np.array(epot)
                self.epot = np.append(self.epot, [epot], axis=0)
                self.vr = np.append(self.vr, [-2.0*ekin/epot], axis=0)
                if (external_mode != 'none'):
                    epot_ext = []
                    if (shell_mode):
                        epot_ext = [np.sum(_particle.pot_ext[n_offset[i]:n_offset[i+1]] * _particle.mass[n_offset[i]:n_offset[i+1]]) if (n_offset[i]<n_offset[i+1]) else 0.0 for i in range(mass_fraction.size)]
                    else:
                        epot_ext = [np.sum(_particle.pot_ext[:n_offset[i+1]] * _particle.mass[:n_offset[i+1]]) if (n_offset[i+1]>0) else 0.0 for i in range(mass_fraction.size)]
                    epot_ext.append(np.sum(_particle.pot_ext[:nc]*_particle.mass[:nc]  if (nc>0) else 0.0))
                    self.epot_ext = np.append(self.epot_ext, [np.array(epot_ext)], axis=0)
            return rlagr

class LagrangianMultiple(DictNpArrayMix):
    """ Lagrangian for single, binaries and all
    Keys: (class members)
        single (Lagrangian): Lagrangian data for single particles
        binary (Lagrangian): Lagrangian data for binaries
        all (Lagrangian): Lagrangian data for all data (binary is treated as c.m., count once)
        Additional members that depend on the keyword argument 'add_star_type':
            [add_star_type] (Lagrangian): Lagrangian data for specific type of stars, see the following content for details
            [add_mass_range] (Lagrangian): Lagrangian data for specific mass ranges of stars, see the following content for details

    Keyword argument:
    --------------
    'add_star_type' (list of string):       
        In case when interrupt_mode='(mo)bse', a list of star types can be given,
        so that the Lagrangian properties for these specific types will be calculated.
        There are four styles of star types:
           (1) a single SSE type name (see help(petar.SSEType))
               For example, if 'BH' is given, one additional class member (type is class Lagrangian), BH, is added.
               The Lagrangian radii are determined by only counting BHs, 
               then the properties (average mass, velocity...) of BHs are calculated using these Lagrangian radii.
           (2) a combination of different SSE types connected by '_'.
               This will include mutliple SSE types in the calculation.
               For example, if 'BH_NS_WD' is given, the additional class member, BH_NS_WD,
               count BHs, NSs and WDs together to calculate the Lagrangian properties.
           (3) a single SSE type name with the prefix 'no'
               This will exclude the given SSE type name in the calculation.
               For exmaple, if 'noBH' is given, the additional class member, noBH, 
               count all stars except BHs to calculate the Lagrangian properties.
           (4) two types (can be any case of the style 1-3) are given by '[type 1]__in__[type 2]'
               Then, the properties of type 1 stars within the spheres or shells of 
               Lagragian radii of type 2 stars are calculated. The Lagrangian radii of type 1 are not calculated.
               For example, if 'BH__in__all' is given, Lagrangian properties of BHs are calculated 
               within the shell or sphere of Lagrangian radii of all stars (instead of Lagrangian radii of BHs).
               Notice that the two type names (except 'all') should also be added separately in the list.
               In the case of 'BH__in__all', add_star_type must contain 'BH', i.e. add_star_type=['BH','BH__in__all', ...].
               In another example, 'BH__in__MS', add_star_type=['BH','MS','BH__in__MS',...].

    'add_mass_range' (list of string): 
         This argument contains a list of mass ranges to select stars for calculating the Lagrangian properties.
         The format of mass range have two styles:
             (1) [minimum mass]_[maximum mass]: the minimum mass and the maximum mass to select stars.
                 For example, if '0.08_1' is given, Lagrangian properties are calculated by selecting stars with masses from 0.08 to 1.0.
                 A new class member 'mass_0.08_1' is added. 
                 The minimum mass must > 0 to ensure the result is correct.
             (2) [mass range 1]__in__[mass range 2]: the mass range 1 and mass range 2 have the format of style (1): [minimum mass]_[maxmimum mass]
                 This type calculates the Lagrangian properties by using the mass range 1 within the Lagrangian radii of the mass range 2.
                 The Lagrangian radii of mass range 1 are not calculated. 
                 For example, if '1_150__in__0.08_1' is given, 'mass_1_150__in__mass_0.08_1' is added, where Lagrangian properties with masses 
                 between 1 and 150 are calculated within the shell or sphere of Lagrangian radii of stars with masses between 0.08 and 1.0.

    For binaries:
        One binary is counted once in class member 'binary', its c.m. data is used.

        In case of [add_star_type] and [add_mass_range], one binary is counted once using its c.m. position and velocity.
        If two components are both this type of star or within the mass range, the total mass is used.
        If one component is this type of star or within the mass range, its mass and c.m. position and velocity are used.

    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            mass_fraction: 1D numpy.ndarray (np.array([0.1, 0.3, 0.5, 0.7, 0.9]))
                An array to indicate the mass fractions to calculate lagrangian radii.
            calc_energy: bool (False)
                If true, add kinetic, potential energies and virial ratio calculation
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, none 
                If it is not none, epot_ext will be added if calc_energy=True
            calc_multi_rc: bool (False)
                If true, calculate core radius using KDTree and correct the central position
            add_star_type: list of string ([])
                A list containing the star type names to calculate additional Lagrangian properties for specific star types.
            add_mass_range: list of mass bins ([])
                A list containing the mass bins for selecting stars to calculate additional Lagrangian properties, mass boundaries should >0
        """
        m_frac=np.array([0.1,0.3,0.5,0.7,0.9])
        if ('mass_fraction' in kwargs.keys()): m_frac=kwargs['mass_fraction'].copy()
        calc_energy=False
        if ('calc_energy' in kwargs.keys()): calc_energy=kwargs['calc_energy']
        external_mode='none'
        if ('external_mode' in kwargs.keys()): external_mode=kwargs['external_mode']
        calc_multi_rc=False
        if ('calc_multi_rc' in kwargs.keys()): calc_multi_rc=kwargs['calc_multi_rc']

        keys=[['time',np.float64], ['single',Lagrangian], ['binary', Lagrangian], ['all', Lagrangian]]
        add_star_type=[]
        if ('add_star_type' in kwargs.keys()): 
            add_star_type=kwargs['add_star_type'].copy()
            for name in add_star_type:
                keys += [[name, Lagrangian]]
        add_mass_range=[]
        if ('add_mass_range' in kwargs.keys()): 
            add_mass_range=kwargs['add_mass_range'].copy()
            for name in add_mass_range:
                name_cross= name.split('__in__')
                if (len(name_cross)>1):
                    keys += [['mass_'+name_cross[0]+'__in__mass_'+name_cross[1], Lagrangian]]
                else:
                    keys += [['mass_'+name, Lagrangian]]

        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        self.initargs['mass_fraction'] = m_frac
        self.initargs['calc_energy'] = calc_energy
        self.initargs['external_mode'] = external_mode
        self.initargs['calc_multi_rc'] = calc_multi_rc
        self.initargs['add_star_type'] = add_star_type
        self.initargs['add_mass_range'] = add_mass_range

        #DictNpArrayMix.__init__(self, [['time',np.float64]], _dat, _offset, _append, **kwargs)
        #self.single = Lagrangian(_dat, _offset+self.ncols, False, **kwargs)
        #self.ncols += self.single.ncols
        #self.keys.append(['single',Lagrangian])
        #self.binary = Lagrangian(_dat, _offset+self.ncols, False, **kwargs)
        #self.ncols += self.binary.ncols
        #self.keys.append(['binary',Lagrangian])
        #self.all    = Lagrangian(_dat, _offset+self.ncols, False, **kwargs)
        #self.ncols += self.all.ncols
        #self.keys.append(['all',Lagrangian])
        #self.size = self.all.size

    def calcTrh(self, G, gamma=0.02, mode='sphere'):
        """ Calculate Spitzer one-component half-mass relaxation time
            the branch "all" is used, where binaries are treated as unresolved single objests.
            Trh = 0.138 N^0.5 Rh^1.5 /( G^0.5 m^0.5 ln(gamma N))
            Then add the result as a new member in the branch all: 'all.trh'

            Notice mass fraction must contain 0.5 to use this function

        Parameters
        ----------
        G: float
           Gravitational constant
        gamma: float (0.02 # Giersz M., Heggie D. C., 1996, MNRAS, 279, 1037)
           The coefficient for Coulomb logarithm
        mode: string (sphere)
            sphere: calculate averaged properties from center to Lagrangian radii
                   m use the average mass within half-mass radius
                   N use the number of objects within half-mass radius
            shell: calculate properties between two neighbor Lagrangian radii

        """
        ncols_old = self.all.ncols
        self.all.calcTrh(G, gamma, mode)
        self.ncols += self.all.ncols - ncols_old
        
    def calcTcr(self, G, mode='sphere'):
        """ Calculate half-mass crossing time
            Tcr = Rh^1.5/sqrt(G M)
            Then add the result as a new member in the branch all: 'all.tcr'

            Notice mass fraction must contain 0.5 to use this function

        Parameters
        ----------
        G: float
           Gravitational constant
        mode: string (sphere)
            sphere: calculate averaged properties from center to Lagrangian radii
            shell: calculate properties between two neighbor Lagrangian radii

        """
        ncols_old = self.all.ncols
        self.all.calcTcr(G,mode)
        self.ncols += self.all.ncols - ncols_old

    def calcPsi(self, rindex=2, star_type1='noBH__in__all', star_type2='BH__in__all'):
        """ Calculate psi factor for coverting one-component relaxation time (trh1) to two-component one (trh2): trh2 = trh1/psi
            Used when two mass components exists, such as black holes and light objects. 
            Add the result as a new member in the branch all: 'all.psi'

            See the theoretical description in Wang, 2020, MNRAS, 491, 2413 (https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2413W/abstract)
            
            s1 = <n_1>*<m_1>^2/|sigma_1|
            s2 = <n_2>*<m_2>^2/|sigma_2|
            s = <n_all>*<m_all>^2/|sigma_all|
            psi = (s1+s2)/s

        Parameters
        ----------
        rindex: int (2)
            the index of mass fraction to select data of n, m and sigma
            default is index for half-mass (50%), assuming the mass fraction used in petar.data.process to calculate lagrangian radii is also the defaulted case 
        star_type1: string ('noBH__in__all')
            name of first component, need to pick up names from 'single, binary' or names defined by 'add_star_type' used in petar.data.process
            default is all non BH components.
        star_type2: string ('BH__in__all')
            name of second component, default is all black holes 
        """

        lagr1=self[star_type1]
        lagr2=self[star_type2]
        lagra=self.all
        s1 = lagr1.n[:,rindex]*lagr1.m[:,rindex]**2/lagr1.sigma.abs[:,rindex]
        s2 = lagr2.n[:,rindex]*lagr2.m[:,rindex]**2/lagr2.sigma.abs[:,rindex]
        s  = lagra.n[:,rindex]*lagra.m[:,rindex]**2/lagra.sigma.abs[:,rindex]
        psi = (s1+s2)/s
        psi[np.isnan(psi)] = 1
        ncols_diff = self.all.addNewMember('psi',psi)
        self.ncols += ncols_diff

    def calcOneSnapshot(self, time, single, binary, rc, mode):
        """ Calculate Lagrangian radii and related properties for one snapshot

        Parameters
        ----------
        time: float
            current evolved time of the system
        single: inherited SimpleParticle
            single partilces (center is corrected and r2 exist)
        binary: Binary
            binaries (center is corrected and r2 exist)
        rc: float
            Core radius
        mode: string
            sphere: calculate averaged properties from center to Lagrangian radii
            shell: calculate properties between two neighbor Lagrangian radii
        """    
        self.time = np.append(self.time, time)
        single_sim = SimpleParticle(single)
        single_sim.calcR2()
        single_sim.addNewMember('pot',single.pot)
        if (self.initargs['external_mode']!='none'):
            single_sim.addNewMember('pot_ext',single.pot_ext)
        binary_sim = SimpleParticle(binary)
        binary_sim.calcR2()
        binary.calcPot()
        binary_sim.addNewMember('pot',binary.pot)
        if (self.initargs['external_mode']!='none'):
            binary.calcPotExt()
            binary_sim.addNewMember('pot_ext',binary.pot_ext)
        if (len(self.initargs['add_star_type'])>0):
            single_sim.addNewMember('type1',single.star.type)
            single_sim.addNewMember('m1',single.mass)
            single_sim.addNewMember('type2',single.star.type)
            single_sim.addNewMember('m2',np.zeros(single.size))
            binary_sim.addNewMember('type1',binary.p1.star.type)
            binary_sim.addNewMember('m1',binary.p1.mass)
            binary_sim.addNewMember('type2',binary.p2.star.type)
            binary_sim.addNewMember('m2',binary.p2.mass)
        elif (len(self.initargs['add_mass_range'])>0):
            single_sim.addNewMember('m1',single.mass)
            single_sim.addNewMember('m2',np.zeros(single.size))
            binary_sim.addNewMember('m1',binary.p1.mass)
            binary_sim.addNewMember('m2',binary.p2.mass)
        all_sim = join(single_sim, binary_sim)
        n_single = single.size
        n_binary = binary.size

        idx = all_sim.r2.argsort()
        idx_single = idx[idx<n_single]
        idx_binary = idx[idx>=n_single] - n_single
        single_sort = single_sim[idx_single]
        binary_sort = binary_sim[idx_binary]
        all_sort = all_sim[idx]
    
        self.single.calcOneSnapshot(single_sort, rc, mode)
        self.binary.calcOneSnapshot(binary_sort, rc, mode)
        rlagr=dict()
        rlagr['all'] = self.all.calcOneSnapshot(all_sort, rc, mode)
        if (self.binary.size != self.single.size):
            raise ValueError('Size inconsistence: single.size:', self.single.size, ' binary.size:', self.binary.size)

        name_list=[]
        cross_list=[]
        obj_sort=dict()
        for name_org in self.initargs['add_star_type']:
            name_cross= name_org.split('__in__')
            if (len(name_cross)>1):
                cross_list.append(name_cross)
                obj_sort[name_cross[0]] = None
            else:
                name_list.append(name_cross[0])

        for name in name_list:
            sel = np.zeros(all_sort.size).astype(bool)
            all_sort.mass = np.zeros(all_sort.size) # clear up mass, final mass is the sum of matched component masses
            for subname in name.split('_'):
                sel1 = np.array([]).astype(bool)
                sel2 = np.array([]).astype(bool)
                if (subname[:2]=='no'):
                    type_index = BSE_STAR_TYPE_INDEX[subname[2:]]
                    sel1 = (all_sort.type1!=type_index) 
                    sel2 = (all_sort.type2!=type_index)
                else:
                    type_index = BSE_STAR_TYPE_INDEX[subname]
                    sel1 = (all_sort.type1==type_index)
                    sel2 = (all_sort.type2==type_index)
                sel = sel | (sel1 | sel2)
                all_sort.mass[sel1] += all_sort.m1[sel1]
                all_sort.mass[sel2] += all_sort.m2[sel2]
            all_sel_sort = all_sort[sel]
            rlagr[name]=self.__dict__[name].calcOneSnapshot(all_sel_sort, rc, mode)
            # save particle list for cross check
            if (name in obj_sort.keys()): obj_sort[name] = all_sel_sort
            #self.__dict__[name+'_cross'].calcOneSnapshot(all_sel_sort, rc, mode, read_rlagr=rlagr)

        name_list=[]
        for name_org in self.initargs['add_mass_range']:
            name_cross= name_org.split('__in__')
            if (len(name_cross)>1):
                for i in range(len(name_cross)):
                    name_cross[i] = 'mass_'+name_cross[i]
                cross_list.append(name_cross)
                obj_sort[name_cross[0]] = None
            else:
                name_list.append('mass_' + name_cross[0])

        for name in name_list:
            label, mmin, mmax = name.split('_')
            mmin = float(mmin)
            mmax = float(mmax)
            all_sort.mass = np.zeros(all_sort.size) # clear up mass, final mass is the sum of matched component masses
            sel1 = (all_sort.m1 >= mmin) & (all_sort.m2 < mmax)
            sel2 = (all_sort.m1 >= mmin) & (all_sort.m2 < mmax)
            sel = (sel1 | sel2)            

            all_sort.mass = np.zeros(all_sort.size) # clear up mass, final mass is the sum of matched component masses
            all_sort.mass[sel1] += all_sort.m1[sel1]
            all_sort.mass[sel2] += all_sort.m2[sel2]
            all_sel_sort = all_sort[sel]

            rlagr[name]=self.__dict__[name].calcOneSnapshot(all_sel_sort, rc, mode)
            # save particle list for cross check
            if (name in obj_sort.keys()): obj_sort[name] = all_sel_sort
            #self.__dict__[name+'_cross'].calcOneSnapshot(all_sel_sort, rc, mode, read_rlagr=rlagr)

        # cross check
        for name_pair in cross_list:
            self.__dict__[name_pair[0]+'__in__'+name_pair[1]].calcOneSnapshot(obj_sort[name_pair[0]], rc, mode, read_rlagr=rlagr[name_pair[1]])

        self.size += 1
        
