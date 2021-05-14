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

