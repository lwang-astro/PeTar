# useful functions
import numpy as np

# vector dot of x, y 
vecDot = lambda x,y: np.sum(x*y,axis=1)

# rotate 3D vector using Euler angles
def vecRot(data, seq, euler_angles, **kwargs):
    """ Rotate the 3D vector array by three rotational angles 
       Using scipy.spatial.transform.Rotation

    Parameters
    ----------
    data: numpy.ndarray (floating with shape of (3,))
        3D vectors to rotate
    seq: string 
        sequencing of x,y,z coordinates (e.g. 'xyz'), see the reference of Rotation. 
    euler_angles: numpy.ndarray (floating of 3) or list
        three euler angles (yaw, pitch, roll) for rotation
    kwargs: dict()
        keyword arguments passed to scipy.spatial.transform.Rotation.from_euler

    Return
    data_rot: numpy.ndarray (floating with shape of (3,))
        rotated vector array
    """
    from scipy.spatial.transform import Rotation as R
    r = R.from_euler(seq, euler_angles, **kwargs)
    return r.apply(data)

def cantorPairing(id1, id2):
    """ Use CantorPairing to map two components id to one binary id

    Parameters
    ----------
    id1: 1D numpy.ndarray or int
        ID for component 1
    id2: 1D numpy.ndarray or int
        ID for component 2
    
    return: binary id
    """ 
    i1=np.minimum(id1,id2).astype('int64')
    i2=np.maximum(id1,id2).astype('int64')
    return ((i1+i2+1)*(i1+i2)/2+i2).astype('int64')

def calcTrh(N, rh, m, G, gamma=0.02): 
    """ Calculate Spitzer one-component half-mass relaxation time
        Trh = 0.138 N^0.5 Rh^1.5 /( G^0.5 m^0.5 ln(gamma N))

    Parameters
    ----------
    N: 1D numpy.ndarray or int
       Total number of particles
    rh: 1D numpy.ndarray or float
       Half-mass radius
    m: 1D numpy.ndarray or float
       mass of one particle
    G: float
       Gravitational constant
    gamma: float (0.02 # Giersz M., Heggie D. C., 1996, MNRAS, 279, 1037)
       The coefficient for Coulomb logarithm

    return: half-mass relaxation time
    """
    return 0.138*N**0.5*rh**1.5/(m**0.5*np.log(gamma*N)*G**0.5)

def calcTcr(M, rh, G): 
    """ Calculate half-mass crossing time
        Tcr = Rh^1.5/sqrt(G M)

    Parameters
    ----------
    M: 1D numpy.ndarray or float
       total mass of the system
    rh: 1D numpy.ndarray or float
       Half-mass radius
    G: float
       Gravitational constant

    return: half-mass crossing time
    """
    return rh**1.5/np.sqrt(G*M)

def calcGWMyr(m1, m2, semi, ecc):
    """ Calculate GW merge timescale in Myr using Peters (1964) formula
    If ecc >1.0, return np.NaN

    Parameters
    ----------
    m1: 1D numpy.ndarray or float
        mass 1 (Msun)
    m2: 1D numpy.ndarray or float
        mass 1 (Msun)
    semi: 1D numpy.ndarray or float
        semi-major axies (pc)
    ecc: 1D numpy.ndarray or float
        eccentricity
    """    
 
    pc_to_au = 206264.81

    # Merging time in myr for one, hyperbolic orbit returns nan
    def time_gw_myr_one(_m1_msun, _m2_msun, _semi_au, _ecc):
     
        G=4*np.pi**2 # gravitational constant in Msun, AU, year
        c=173*365.25 # speed of light
        beta = (64.0/5.0)*G**3*_m1_msun*_m2_msun*(_m1_msun+_m2_msun)/c**5
        
        if (_ecc==0.0):
            return _semi_au**4/(4*beta)*1e-6
        elif (_ecc>=1.0):
            return np.NaN
        else:
            c0 = _semi_au/(_ecc**(12.0/19.0)/(1-_ecc**2)*(1+(121.0/304.0)*_ecc**2)**(870.0/2299.0))
            def e_fun(ecc_):
                return ecc_**(29.0/19.0)*(1+(121.0/304.0)*ecc_**2)**(1181.0/2299.0)/(1-ecc_**2)**1.5
            eint=integrate.quad(e_fun,0,_ecc)
            return (12.0/19.0)*c0**4/beta*(eint[0])*1e-6

    semi_au = semi*pc_to_au
    if (type(m1) == np.ndarray | type(m1) == list):
        return np.array(list(map(time_gw_myr_one,m1,m2,semi_au,ecc)))
    else: 
        return time_gw_yr(m1, m2, semi_au, ecc)

def convergentPointCheck(data, velocity):
    """ calculate proper motions in the frame of convergent point based on the given velocity and calculate the residuals 
        The method is described in e.g., van Leeuwen F., 2009, A\&A, 497, 209. doi:10.1051/0004-6361/200811382; 
        and Jerabkova T., Boffin H.~M.~J., Beccari G., de Marchi G., de Bruijne J.~H.~J., Prusti T., 2021, A\&A, 647, A137. doi:10.1051/0004-6361/202039949

        The algorithm can be described as follows
       1. Assume that the stars have the given velocity and rotate it by RA (alpha) and DEC (delta) of stars. Thus this predicted velocity is in the same frame as the observation.
       2. Calculate the angle between two components of the predicted proper motions.
       3. Rotate both the predicted and observed proper motions to the regular frame where the x component is proprotional to v*sin(lambda), 
          where lambda is the angle between the radial velocity and the full velocity, and the y component is zero.
       4. Calculate the difference of the predicted and observed proper motions in the new frame in x, y directions.

    Parameters
    ----------
    data: astropy.coordinates.SkyCoord
        postions and velocities of stars
    velocity: numpy.array (floating with astropy.units and a size of 3)
        The reference co-moving velocity in the Cartesian frame

    Return
    ----------
    dmu_parallel: numpy.ndarray(floating), size of input data
       The difference of the predicted proper motion and the observed proper motion along the parallel direction (x-axis) in the regular frame
    dmu_perpendicular: numpy.ndarray(floating), size of input data
       The difference of the predicted proper motion and the observed proper motion along the perperndicular direction (x-axis) in the regular frame
    angle: numpy.ndarray(floating), size of input data
       The angle for rotating the proper motions to the regular frame
    """
    import astropy.units as u
    
    parallax = 1000/data.icrs.distance
    ra = data.icrs.ra.to(u.rad)
    dec= data.icrs.dec.to(u.rad)
    pm_ra_cosdec = data.icrs.pm_ra_cosdec*4.74047/parallax
    pm_dec = data.icrs.pm_dec*4.74047/parallax

    vpred = vecRot(velocity,'zyx', np.transpose(np.array([-ra, dec, np.zeros(ra.size)])))  # return vr, pm_ra_cosdec, pm_dec
    psi = np.arctan2(vpred[:,2],vpred[:,1]) # angle 
    vpred_rot = vecRot(vpred, 'x', -psi)
    vobs = np.transpose([data.icrs.radial_velocity, pm_ra_cosdec, pm_dec])
    vobs_rot = vecRot(vobs, 'x', -psi)
    
    dmu_parallel = vobs_rot[:,1] - vpred_rot[:,1]
    dmu_perpendicular = vobs_rot[:,2] - vpred_rot[:,2]
    
    return dmu_parallel, dmu_perpendicular, psi

def coordinateCorrection(data, snap_center, obs_center, **kwargs):
    """ Correct c.m. coordinate based on the difference between snapshot center and observational center in the galactocentric frame
    First the snap and obs centers are transferred to the galactocentric frame
    The difference of coordinates in the sphercial representational type are calculated and added to particle data as the correction
    
    Parameters:
    -----------
    data: astropy.coordinates.SkyCoord
         particle data
    snap_center: astropy.coordinates.SkyCoord (one particle)
         snapshot center coordinate
    obs_center:  astropy.coordinates.SkyCoord (one particle)
         observational center coordinate

    keyword arguments:
       galcen_distance: floating with length units (8.0*units.kpc [Galpy])
            galactic central distance of the Sun
       z_sun: floating with length units (15.0*units.pc [Galpy])
            z direction distance of the Sun
       galcen_v_sun: astropy.coordinates.CartesianDifferential ([10.0, 235.0, 7.0]*units.km/units.s [Galpy])
            velocity of the Sun    

    Return:
    ----------
    data_c: astropy.coordinates.SkyCoord
       corrected particle data
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

    obs_cg = obs_center.transform_to(Galactocentric(**parameters))
    obs_cg.representation_type = 'spherical'

    snap_cg = snap_center.transform_to(Galactocentric(**parameters))
    snap_cg.representation_type = 'spherical'

    data_g = data.transform_to(Galactocentric(**parameters))
    data_g.representation_type = 'spherical'
    
    dlon = obs_cg.lon - snap_cg.lon
    dlat = obs_cg.lat - snap_cg.lat
    lon_new = data_g.lon + dlon
    lat_new = data_g.lat + dlat

    sel = lat_new > 90*u.deg
    lat_new[sel] = 180*u.deg - lat_new[sel]
    lon_new[sel] = 180*u.deg + lon_new[sel]
    
    sel = lat_new < -90*u.deg
    lat_new[sel] = -180*u.deg - lat_new[sel]
    lon_new[sel] = 180*u.deg + lon_new[sel]

    sel = (lon_new > 360*u.deg) | (lon_new < -360*u.deg)
    lon_new[sel] = lon_new[sel] - (lon_new[sel].to(u.deg).value/360.0).astype(int)*360*u.deg

    ddis = obs_cg.distance - snap_cg.distance
    dpm_lon = obs_cg.pm_lon - snap_cg.pm_lon
    dpm_lat = obs_cg.pm_lat - snap_cg.pm_lat
    drv  = obs_cg.radial_velocity - snap_cg.radial_velocity
    
    
    data_c = SkyCoord(lon = lon_new,
                      lat = lat_new,
                      distance = data_g.distance + ddis,
                      pm_lon = data_g.pm_lon + dpm_lon,
                      pm_lat = data_g.pm_lat + dpm_lat,
                      radial_velocity = data_g.radial_velocity + drv,
                      frame='galactocentric', representation_type='spherical', **parameters)
    
    return data_c

