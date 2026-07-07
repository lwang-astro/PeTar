# read snapshot and obtain multiple systems
import collections
from scipy import spatial as sp
from sdar.base import *
from sdar.functions import *
from sdar.ar import SDARInterruptBinary
from sdar.ar import SDARData as sdar_SDARData
from sdar.particle import SimpleParticle, ParticleGroup
from sdar.particle import Binary as sdar_Binary
from sdar.hermite import HermiteData as hermite_HermiteData
import sdar.group as hermite_group
from .bse import *
from .dsm import *

G_MSUN_PC_MYR=0.00449830997959438 # Msun, pc, myr
G_HENON=1 # Henon unit
HEADER_OFFSET=24 # header offset in bytes for snapshots with the BINARY format
HEADER_OFFSET_F128=32 # header offset in bytes for snapshots with the BINARY format
HEADER_OFFSET_WITH_CM=72 # header offset with center-of-the-mass data in bytes for snapshots with the BINARY format
HEADER_OFFSET_WITH_CM_F128=128 # header offset with center-of-the-mass data in bytes for snapshots with the BINARY format

class PeTarDataHeader():
    """ Petar snapshot data header
    members:
        file_id: int 
           file id
        n: int 
           number of particles
        time: float 
           time of snapshot
        *pos_offset: list of float (length of 3)
           position offset of particle system
        *vel_offset: list of float (length of 3)
           velocity offset of particle system
  
        pos_offset and vel_offset only exist when keyword argument 'external_mode' is not none
    """

    def __init__(self, _filename=None, **kwargs):
        """ Initial data header
        
        Parameters:
        -----------
        _filename: string
            PeTar snapshot file name to read the header, if not provide, all members are initialized to zero (None)
        kwargs: dict
            Keyword arguments:
            float_type: type (np.float64)
                floating point data type
            snapshot_format: string (binary)
                Data format of snapshot files: binary or ascii
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, agama, none 
                If not none, this option indicates the pos_offset and vel_offset exists 
        """
        self.file_id = int(0)
        self.n = int(0)
        self.time = 0.0
        self.pos_offset=[0.0,0.0,0.0]
        self.vel_offset=[0.0,0.0,0.0]
        self.offset_flag=False
        
        if (_filename!=None): self.read(_filename,**kwargs)

    def read(self, _filename, **kwargs):
        """ Read snapshot file to obtain the header information

        Parameters:
        -----------
        _filename: string
            PeTar snapshot file name to read the header
        kwargs: dict
            Keyword arguments:
            float_type: type (np.float64)
                floating point data type
            snapshot_format: string (binary)
                Data format of snapshot files: binary or ascii
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, agama, none 
                If not none, this option indicates the pos_offset and vel_offset exists 
        """
        snapshot_format='binary'
        if ('snapshot_format' in kwargs.keys()): snapshot_format=kwargs['snapshot_format']
        if ('external_mode' in kwargs.keys()):
            if (kwargs['external_mode']!='none'): self.offset_flag=True
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        if (snapshot_format=='ascii'):
            fp = open(_filename, 'r')
            header=fp.readline()
            header_items=header.split()
            if (self.offset_flag):
                if (len(header_items)!=9):
                    raise ValueError('Snapshot header item number mismatch! Need 9 (file_id, N, time, xcm, ycm, zcm, vxcm, vycm, vzcm), got %d. Make sure the external_mode keyword set correctly.' % len(header_items))

                file_id, n_glb, t, x, y, z, vx, vy, vz = header_items
                fp.close()

                self.file_id = int(file_id)
                self.n = int(n_glb)
                self.time = float_type(t)
                self.pos_offset = [float_type(x),float_type(y),float_type(z)]
                self.vel_offset = [float_type(vx),float_type(vy),float_type(vz)]
            else:
                if (len(header_items)!=3):
                    raise ValueError('Snapshot header item number mismatch! Need 3 (file_id, N, time), got %d. Make sure the external_mode keyword set correctly.' % len(header_items))

                file_id, n_glb, t = header_items
                fp.close()

                self.file_id = int(file_id)
                self.n = int(n_glb)
                self.time = float_type(t)

        else:
            if (self.offset_flag):
                fp = np.fromfile(_filename, dtype=np.dtype([('file_id',np.int64),('n_glb',np.int64),('time',float_type),('x',float_type),('y',float_type),('z',float_type),('vx',float_type),('vy',float_type),('vz',float_type)]),count=1)
                self.file_id = fp['file_id'][0]
                self.n = fp['n_glb'][0]
                self.time = fp['time'][0]
                self.pos_offset = np.array([fp['x'][0], fp['y'][0], fp['z'][0]])
                self.vel_offset = np.array([fp['vx'][0], fp['vy'][0], fp['vz'][0]])
            else:
                fp = np.fromfile(_filename, dtype=np.dtype([('file_id',np.int64),('n_glb',np.int64),('time',float_type)]),count=1)
                self.file_id = fp['file_id'][0]
                self.n = fp['n_glb'][0]
                self.time = fp['time'][0]

    def savetxt(self, fname, **kwargs):
        """ Save class member data to a file
        Use the getherDataToArray and then numpy.savetxt

        Parameters
        ----------
        fname: string of filename or file handler
        kwargs: dict
            keyword arguments for numpy.savetxt
        """
        offset_flag=False
        close_flag = False
        if (type(fname)==str):
            f = open(fname, 'w')
            close_flag = True
        elif (hasattr(fname, 'write')):
            f = fname
        if (self.offset_flag):
            f.write("%d %d %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n" % 
                    (self.file_id, self.n, self.time, 
                     *self.pos_offset, *self.vel_offset))
        else:
            f.write("%d %d %.20g\n" % (self.file_id, self.n, self.time))

        if close_flag:
            f.close()

    def tofile(self, fname):
        """ Write class member data to a file using numpy.save
        Use numpy.save to write data, the dtype is defined by keys (members)

        Parameters
        ----------
        fname: string of filename or file header
        kwargs: dict
            keyword arguments for numpy.save, notice dtype is already defined, do not provide that
        """

        close_flag = False
        if (type(fname)==str):
            f = open(fname, 'wb')
            close_flag = True
        elif (hasattr(fname, 'write')):
            f = fname

        import struct
        if (self.offset_flag):
            header_buffer = struct.pack('qqddddddd', self.file_id, self.n, self.time, *self.pos_offset, *self.vel_offset)
        else:
            header_buffer = struct.pack('qqd', self.file_id, self.n, self.time)
        f.write(header_buffer)

        if close_flag:
            f.close()
        

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
                parameters[key] = kwargs[key]

        sky = SkyCoord(x=self.pos_offset[0]*pos_unit, 
                       y=self.pos_offset[1]*pos_unit, 
                       z=self.pos_offset[2]*pos_unit, 
                       v_x=self.vel_offset[0]*vel_unit,
                       v_y=self.vel_offset[1]*vel_unit,
                       v_z=self.vel_offset[2]*vel_unit,
                       frame='galactocentric', representation_type='cartesian', **parameters)
        return sky
    
class BaseParticle(SimpleParticle):
    """ Base particle type of PeTar
        The members include simple particle information, binary status and stellar evolution data
       
    keys: (class members)
        Members inherited from SimpleParticle: mass (1D), pos (2D,3), *pos_high (2D,3) vel (2D,3) 
            see help(petar.SimpleParticle)
        binary_state (1D): binary interruption state
        if (keyword argument 'interrupt_mode' == 'base', 'bse', 'bseEmp', 'mobse', 'dsm'):
            radius:        (1D): radius for merger checker
            dm:            (1D): mass loss
            time_record    (1D): last time of interruption check
            time_interrupt (1D): next interruption time
        if (keyword argument 'interrupt_mode' == 'bse', 'bseEmp', 'mobse'):
            star  (SSEStarParameter): BSE based stellar evolution parameters
        if (keyword argument 'interrupt_mode' == 'dsm'):
            star  (DSMStarParameter): Disk star merger parameters
    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): merger, base, bse, mobse, none, dsm
               This option indicates whether columns of stellar evolution exist
            use_mpfrc: bool (False)
                if true, add three columns of pos_high indicating the high-precision parts of position
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys_bstat = [['binary_state',np.int64]]
        keys_se  = [['radius',float_type],['dm',float_type],['time_record',float_type],['time_interrupt',float_type]]    
        
        keys = keys_bstat
        if ('interrupt_mode' in kwargs.keys()):
            if (kwargs['interrupt_mode'] in ('base', 'merger')):
                keys = keys_bstat+keys_se
            elif ('bse' in kwargs['interrupt_mode']):
                keys = keys_bstat+keys_se+[['star',SSEStarParameter]]
            elif (kwargs['interrupt_mode']=='dsm'):
                keys = keys_bstat+keys_se+[['star',DSMStarParameter]]
            
        SimpleParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

class HardParticle(BaseParticle):
    """ Hard particle type of PeTar
        The member include BaseParticle and searching radius, id, artificial particle data and changeover radii
        
    keys: (class members)
        Members inherited from BaseParticle: see help(petar.BaseParticle)
             Please set the keyword argument 'interrupt_mode' to determine the members of stellar evolution data
        r_search (1D): searching radius
        id       (1D): identification
        mass_bk  (1D): artificial particle parameter 1 
        status   (1D): artificial particle parameter 2
        r_in     (1D): changeover function inner boundary
        r_out    (1D): changeover function outer boundary
    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        keyword arguments:
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): merger, base, bse, mobse, none, dsm
               This option indicates whether columns of stellar evolution exist
            use_mpfrc: bool (False)
                if true, add three columns of pos_high indicating the high-precision parts of position
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['r_search',float_type], ['id',np.int64], ['mass_bk',np.float64], ['status',np.float64], ['r_in',float_type], ['r_out',float_type]]

        BaseParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

class HermiteParticle(HardParticle):
    """ Hermite particle type of PeTar
        The member include HardParticle and dt, time, acc, jerk and pot
        
    keys: (class members)
        Members inherited from HardParticle: see help(petar.HardParticle)
             Please set the keyword argument 'interrupt_mode' to determine the members of stellar evolution data
        dt    (1D): time step
        time  (1D): current time
        acc   (2D,3): acceleration x, y, z
        jerk  (2D,3): acceleration derivative x, y, z
        pot   (1D): potential
    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        keyword arguments:
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): merger, base, bse, mobse, none, dsm
               This option indicates whether columns of stellar evolution exist
            use_mpfrc: bool (False)
                if true, add three columns of pos_high indicating the high-precision parts of position
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['dt',float_type],['time',float_type],['acc',(float_type,3)],['jerk',(float_type,3)],['pot',float_type]]

        HardParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

class Particle(HardParticle):
    """ (Soft) Particle type of PeTar, also used in snapshot
        The particle data of PeTar. Depending on the compile configuration of PeTar, 
        The data structures (columns) of the particle snapshots are different.
        Using the correct keyword arguments in the initialization to control the member definition (Keys)

    keys: (class members)
        Members inherited from HardParticle: see help(petar.HardParticle)
             Please set the keyword argument 'interrupt_mode' to determine the members of stellar evolution data
             Please set the keyword argument 'use_mpfrc' to determine whether high-precision parts of particle position are included
        acc_soft (2D,3): long-range interaction acceleration (particle-tree) x, y, z
        if (keyword argument 'collect_sp_acc' == True'):
            acc_sp (2D,3): superparticle acceleration (only used when superparticle is enabled) x, y, z
        pot      (1D): total potential
        pot_soft (1D): long-range interaction potential
        if (keyword argument 'external_mode' != 'none'):
             pot_ext  (1D): external potential (only exist when keyword argument 'external_mode' is not 'none')
        n_nb:    (1D): number of neighbors (short-interaction)
    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): merger, base, bse, mobse, none, dsm
               This option indicates whether columns of stellar evolution exist
            external_mode: string (none)
               PeTar external mode (set in configure): galpy, agama, none 
               This option indicates whether the column of externa potential exist
            use_mpfrc: bool (False)
               If true, add three columns of pos_high indicating the high-precision parts of position
            collect_sp_acc: bool (False)
               If true, the superparticle acceleration is collected and the column acc_sp exists
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['acc_soft',(float_type,3)]]
        if ('collect_sp_acc' in kwargs.keys()):
            if (kwargs['collect_sp_acc']):
                keys += [['acc_sp',(float_type,3)]]
        keys += [['pot',float_type], ['pot_soft',float_type]]
        if ('external_mode' in kwargs.keys()):
            if (kwargs['external_mode']!='none'):
                keys += [['pot_ext',float_type]]
        keys += [['n_nb',np.int64]]

        HardParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

    def calcEtot(self):
        """ Calculate total energy and add it as the member, etot
        """
        etot = self.ekin + self.mass*self.pot
        self.addNewMember('etot',etot)

class Binary(sdar_Binary):
    """ Binary class
        The binary (tree) data. Depending on the definition of two members 
        (keyword argument member_particle_type(|_one|_two), 
        The binary can refer to any type of multiple system.

    Keys:
        The final keys depends on kwargs of initial function
  
        kwargs['simple_mode'] (bool)
            True: (default)
                mass (1D): total mass of two components
                pos  (2D,3): c.m. position x, y, z
                vel  (2D,3): c.m. velocity vx, vy, vz
                rrel (1D): relative distance
                semi (1D): semi-major axis
                ecc  (1D): eccentricity
                p1   (member_particle_type_one) component one
                p2   (member_particle_type_two) component two
            False:
                mass (1D): total mass of two components
                pos  (2D,3): c.m. position x, y, z
                vel  (2D,3): c.m. velocity vx, vy, vz
                m1   (1D): component 1 mass
                m2   (1D): component 2 mass
                rrel (1D): relative distance
                semi (1D): semi-major axis
                am   (2D,3): specific angular momemtum x, y, z
                L    (2D,3): angular momemtum x, y, z
                eccvec  (2D,3): eccentric vector
                incline (1D): inclination
                rot_horizon (1D): frame rotational angle in x-y plane (longitude of ascending node)
                ecc  (1D): eccentricity
                rot_self (1D): frame rotational angle in orbital plane (argument of periapsis)
                ecca (1D): eccentric anomaly
                period (1D): period
                t_peri (1D): time to peri-center
                p1 (member_particle_type_one) component one
                p2 (member_particle_type_two) component two

        The member_particle_type(|_one|_two) is given by keyword arguments:
           'member_particle_type' (for both members),'member_particle_type_one','member_particle_type_two'.
        In default, it is petar.Particle.
        If a type (e.g., petar.Particle) is given, the member is a single star.
        If a list with two members (e.g., [petar.Particle, petar.Particle]) is given, 
        the member is a binary with two single stars.
        A hierarchical list can be provided, e.g., [petar.Particle, [petar.Particle, petar.Particle]]
        to indicate a triple system.
    """

    def __init__ (self, _p1=None, _p2=None, _offset=int(0), _append=False, **kwargs):
        """
        Parameters
        ----------
        _p1: particle data | 2D numpy.ndarray | petar.Binary | None
            If _p1 is a particle type data or a petar.Binary type data, it is treated as the first component of binary
            If _p1 is a petar.Binary type data and _p2 is None, the class instance is initialized by copy the data of _p1.
            If _p1 is None, initialize class with empty data
        _p2: particle data | None
            If _p2 is a particle type data or a petar.Binary type data, it is treated as the second component of binary
            If _p2 is None, _p1 should be petar.Binary data or None
        _offset: int (0)
            Reading column offset of _dat if it is 2D np.ndarray
        _append: bool (False)
            If true, append keys and ncols to the current class instead of create new class members

        keyword arguments:
            simple_mode: bool (True)
                If True, only calculate semi and ecc, save computing time significantly
            G: float (1.0)
                Gravitational constant
            interrupt_mode: string (none)
                PeTar interrupt mode (set in configure): merger, base, bse, mobse, none, dsm
                This option indicates whether columns of stellar evolution exist
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, agama, none 
                This option indicates whether the column of external potential exist
            use_mpfrc: bool (False)
               If true, add three columns of pos_high indicating the high-precision parts of position
            member_particle_type: type or list (SimpleParticle)
                Type of component particle (both)
            member_particle_type_one: type or list (SimpleParticle)
                Type of 1st component
            member_particle_type_two: type or list (SimpleParticle)
                Type of 2nd component 
            float_type: type (np.float64)
                floating point data type
        """
        kwargs_local = dict(kwargs)
        kwargs_local['member_particle_type'] = Particle
        super().__init__(_p1, _p2, _offset, _append, **kwargs_local)


class InterruptBinary(SDARInterruptBinary):
    """ Data of stellar evolution interrupted binary in base mode
        Inherit from sdar.ar.SDARInterruptBinary
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ Initial InterruptBinary class
        Parameters
        ----------
        ----------
        keyword arguments:
            particle_type: type (HardParticle)
                particle data type
        """

        if (not 'particle_type' in kwargs.keys()):
            kwargs['particle_type'] = HardParticle
        particle_type = kwargs['particle_type']

        SDARInterruptBinary.__init__(self, _dat, _offset, _append, **kwargs)

def calculateParticleCMDict(pcm, _p1, _p2):
    """ Calculate the center-of-the-mass of two particle sets
    
    Parameters
    ----------
    _p1: inherited SimpleParticle
        particle set 1
    _p2: inherited SimpleParticle 
        particle set 2, should have the same size as _p1
    pcm: dict 
        particle center-of-the-mass, should include keys: 'mass','pos','vel'.
    """
    if (issubclass(type(_p1), SimpleParticle)) & (issubclass(type(_p2),SimpleParticle)):
        pcm['mass'] = _p1.mass + _p2.mass
        pcm['pos']  = (_p1.mass[:,None]*_p1.pos + _p2.mass[:,None]*_p2.pos)/pcm['mass'][:,None]
        pcm['vel']  = (_p1.mass[:,None]*_p1.vel + _p2.mass[:,None]*_p2.vel)/pcm['mass'][:,None]
    elif (isinstance(_p1, collections.OrderedDict)) & (isinstance(_p2,collections.OrderedDict)) | (isinstance(_p1, dict)) & (isinstance(_p2, dict)):
        pcm['mass'] = _p1['mass'] + _p2['mass']
        pcm['pos']  = (_p1['mass'][:,None]*_p1['pos'] + _p2['mass'][:,None]*_p2['pos'])/pcm['mass'][:,None]
        pcm['vel']  = (_p1['mass'][:,None]*_p1['vel'] + _p2['mass'][:,None]*_p2['vel'])/pcm['mass'][:,None]
    else:
        raise ValueError('Initial fail, date type should be Particle or collections.OrderDict, given',type(_p1))

class GroupInfo(hermite_group.GroupInfo):
    """ Group information output from PeTar
    Keys: (class members)
        type (1D): group type, 0: new group; 1: end group
        n    (1D): number of members in group (should be consistent with keyword argument N
        time (1D): current time
        pos  (2D,3): position of the group c.m. in the framework of the global system (without shift of global system c.m. if external_mode is on)
        vel  (2D,3): velocity of the group c.m. in the framework of the global system (without shift of global system c.m.)
        bin[X] (BinaryTreeSDAR): members of the group in a hierarchical binary tree
               Here X indicates the order. 0 represents the root (outer most) binary; 1,2,3 ... are inner binaries
               For a triple, bin0 is outer binary, bin1 is inner binary.
               p2 of bin0 is the c.m. of bin1, the id of p2 is the minimum id from the two components in bin1.
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ GroupInfo class for petar, inherit from hermite_group.GroupInfo

            Parameters
            ----------
            keyword arguments:
                member_particle_type: type (HardParticle)
                    Type of component particle, do not change this!
                interrupt_mode: string (none)
                    PeTar interrupt mode (set in configure): merger, base, bse, mobse, none, dsm
                    This option indicates whether columns of stellar evolution exist
                external_mode: string (none)
                    PeTar external mode (set in configure): galpy, agama, none 
                    This option indicates whether the column of external potential exist
                use_mpfrc: bool (False)
                    If true, add three columns of pos_high indicating the high-precision parts of position
                float_type: type (np.float64)
                    floating point data type
                N: int (2)
                    Number of members of one group
        """
        kwargs_local = dict(kwargs)
        kwargs_local['member_particle_type'] = HardParticle
        super().__init__(_dat, _offset, _append, **kwargs_local)
    
class SDARData(sdar_SDARData):
    """ SDARData class for petar, inherit from sdar.SDARData

        Keys: (class members)
            time (1D): current evolved time (counting from zero)
            de (1D): physical energy error
            etot_ref (1D): initial total energy
            ekin (1D): kinetic energy
            epot (1D): potential energy
            gt_drift (1D): time tranformation for drift step
            H (1D): extened phase space Hamiltonian
            if (keyword argument 'include_H_approx' == True):
                H_approx (1D): approximated phase space Hamiltonian
            de_interrupt (1D): energy change due to interruption
            dH_interrupt (1D): H change due to interruption
            perturber (perturber_type): perturber data, depending on the keyword argument 'perturber_type'
                                        if not given, this key is not included
            info (SDARInfo): SDAR information shown as follows:
                ds (1D): integration step
                time_offset (1D): time offset to obtain the actual time (time_offset + time)
                r_break_crit (1D): distance criterion to break group (used in Hermite)
            if (keyword argument 'hybrid' == True):
                hybrid_flag (1D): if 1, hybrid method is used, else, normal method
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
            particles (ParticleGroup): particle group, depending on the keyword argument 'member_type' and 'cm_type'

    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ SDARData class for petar, inherit from sdar.SDARData

            Parameters
            ----------
            keyword arguments:
                member_type: member particle type (HardParticle)
                cm_type: c.m. particle type (HermiteParticle)
                perturber_type: perturber particle type (None)
                N_particle: int (0)
                    Number of particles, determined from file if not provided
                slowdown: bool (False)
                    if True, add slowdown keys
                N_sd: int (0)
                    Number of slowdown pair, used when slowdown='on'
                time_measure: bool (False)
                    if True, add time measure keys in profile
                include_H_approx: bool (False)
                    if True, add H_approx key
                interrupt_mode: string (none)
                    PeTar interrupt mode (set in configure): merger, base, bse, mobse, none, dsm
                    This option indicates whether columns of stellar evolution exist
                external_mode: string (none)
                    PeTar external mode (set in configure): galpy, agama, none 
                    This option indicates whether the column of external potential exist
                use_mpfrc: bool (False)
                    If true, add three columns of pos_high indicating the high-precision parts of position
                float_type: type (np.float64)
                    floating point data type
        """
        kwargs_local = dict(kwargs)
        kwargs_local['member_type'] = HardParticle
        kwargs_local['cm_type'] = HermiteParticle

        super().__init__(_dat, _offset, _append, **kwargs_local)

class HermiteData(hermite_HermiteData):
    """ HermiteData class for petar, inherit from sdar.HermiteData
    Keys: (class members)
        time (1D): current evolved time (counting from zero)
        time_offset (1D): time offset to calculate the global time (time+time_offset)
        energy_phy (HermiteEnergy): physical energy data
        energy_sd (HermiteEnergy): slowdown energy data
        sd (SlowDownGroup): slowdown data
        profile (HermiteProfile): hermite profile
        particles (ParticleGroup): particle group, depending on the keyword argument 'member_type' and 'cm_type'
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ HermiteData class for petar, inherit from sdar.HermiteData

            Parameters
            ----------
            keyword arguments:
                member_type: member particle type (HermiteParticle)
                cm_type: c.m. particle type (HermiteParticle)
                N_particle: int (0)
                    Number of members of one group
                N_sd: int (0)
                    Number of slowdown pairs
                time_measure: bool (False)
                    if True, add time measure keys in profile
                interrupt_mode: string (none)
                    PeTar interrupt mode (set in configure): merger, base, bse, mobse, none, dsm
                    This option indicates whether columns of stellar evolution exist
                external_mode: string (none)
                    PeTar external mode (set in configure): galpy, agama, none 
                    This option indicates whether the column of external potential exist
                use_mpfrc: bool (False)
                    If true, add three columns of pos_high indicating the high-precision parts of position
                float_type: type (np.float64)
                    floating point data type
        """
        kwargs_local = dict(kwargs)
        kwargs_local['member_type'] = HermiteParticle
        kwargs_local['cm_type'] = HermiteParticle

        super().__init__(_dat, _offset, _append, **kwargs_local)