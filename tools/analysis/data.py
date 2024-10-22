# read snapshot and obtain multiple systems
import collections
from scipy import spatial as sp
from .base import *
from .bse import *
from .functions import *

G_MSUN_PC_MYR=0.00449830997959438 # Msun, pc, myr
G_HENON=1 # Henon unit
HEADER_OFFSET=24 # header offset in bytes for snapshots with the BINARY format
HEADER_OFFSET_WITH_CM=72 # header offset with center-of-the-mass data in bytes for snapshots with the BINARY format

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
            snapshot_format: string (ascii)
                Data format of snapshot files: binary or ascii
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, none 
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
            snapshot_format: string (ascii)
                Data format of snapshot files: binary or ascii
            external_mode: string (none)
                PeTar external mode (set in configure): galpy, none 
                If not none, this option indicates the pos_offset and vel_offset exists 
        """
        snapshot_format='ascii'
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
                self.pos_offset = [fp['x'][0], fp['y'][0], fp['z'][0]]
                self.vel_offset = [fp['vx'][0], fp['vy'][0], fp['vz'][0]]
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
                parameter[key] = kwargs[key]

        sky = SkyCoord(x=self.pos_offset[0]*pos_unit, 
                       y=self.pos_offset[1]*pos_unit, 
                       z=self.pos_offset[2]*pos_unit, 
                       v_x=self.vel_offset[0]*vel_unit,
                       v_y=self.vel_offset[1]*vel_unit,
                       v_z=self.vel_offset[2]*vel_unit,
                       frame='galactocentric', representation_type='cartesian', **parameters)
        return sky
    

class SimpleParticle(DictNpArrayMix):
    """ Simple particle class with only mass, postion, velocity
    keys: (class members)
        mass (1D): mass
        pos (2D,3): postion x, y, z
        *pos_high (2D,3): high-precision parts of position x, y, z, only exist when use_mpfrc is True
        vel (2D,3): velocity vx, vy, vz
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters:
        -----------
        Keyword arguments:
            float_type: type (np.float64)
                floating point data type
            use_mpfrc: bool (False)
                if true, add three columns of pos_high indicating the high-precision parts of position
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64
        keys = [['mass', float_type], ['pos', (float_type, 3)]]
        if ('use_mpfrc' in kwargs.keys()): use_mpfrc = kwargs['use_mpfrc']
        else: use_mpfrc = False
        if (use_mpfrc):
            keys += [['pos_high', (float_type, 3)]]
        keys += [['vel', (float_type, 3)]]
        
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def calcR2(self):
        """ calculate distance square, r2, and add/update it as a class member
        """
        if (self.size>0):
            r2 = vecDot(self.pos,self.pos)
        else:
            r2 = np.array([])
        self.addNewMember('r2',r2)

    def calcEkin(self):
        """ calculate kinetic energy, ekin, and add/update it as a class member
        """
        if (self.size>0):
            ekin = 0.5*vecDot(self.vel,self.vel)*self.mass
        else:
            ekin = np.array([])
        self.addNewMember('ekin',ekin)

    def correctCenter(self, cm_pos, cm_vel):
        self.pos -= cm_pos
        self.vel -= cm_vel


    def toSkyCoord(self, **kwargs):
        """ generate astropy.coordinates.SkyCoord data
            Be careful when external_mode is used, remember to use the keyword arguments pos_offset and vel_offset to add the center shift in the header of PeTar snapshot.
        Parameters
        -----------------
        kwargs: dict()
            pos_offset: numpy.ndarray ([0.0,0.0,0.0])
                 position offset to add
            vel_offset: numpy.ndarray ([0.0,0.0,0.0]) 
                 velocity offset to add
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
        snap: astropy.coordinates.SkyCoord
            snapshot data using SkyCoord
        """
        import astropy 
        from astropy.coordinates import SkyCoord  # High-level coordinates
        from astropy.coordinates import ICRS, Galactic, Galactocentric, FK4, FK5  # Low-level frames
        from astropy.coordinates import Angle, Latitude, Longitude  # Angles
        from astropy.coordinates import CartesianDifferential
        import astropy.units as u

        cm_cor=np.zeros(6)
        if ('pos_offset' in kwargs.keys()):
            pos_offset = kwargs['pos_offset']
            if(type(pos_offset)==np.ndarray) | (type(pos_offset)==list):
                cm_cor[0] = pos_offset[0]
                cm_cor[1] = pos_offset[1]
                cm_cor[2] = pos_offset[2]
            else:
                raise ValueError('pos_offset should be an array or a list with size of 3, given ', pos_offset)
        if ('vel_offset' in kwargs.keys()):
            vel_offset = kwargs['vel_offset']
            if(type(vel_offset)==np.ndarray) | (type(vel_offset)==list):
                cm_cor[3] = vel_offset[0]
                cm_cor[4] = vel_offset[1]
                cm_cor[5] = vel_offset[2]
            else:
                raise ValueError('vel_offset should be an array or a list with size of 3, given ', vel_offset)

        pos_unit = u.pc
        if ('pos_unit' in kwargs.keys()): pos_unit = kwargs['pos_unit']
        vel_unit = u.pc/u.Myr
        if ('vel_unit' in kwargs.keys()): vel_unit = kwargs['vel_unit']

        parameters={'galcen_distance':8.0*u.kpc, 'z_sun':15.*u.pc, 'galcen_v_sun':CartesianDifferential([10.0,235.,7.]*u.km/u.s)}
        for key in parameters.keys():
            if key in kwargs.keys():
                parameters[key] = kwargs[key]

        snap = SkyCoord(x=(self.pos[:,0]+cm_cor[0])*pos_unit, 
                        y=(self.pos[:,1]+cm_cor[1])*pos_unit, 
                        z=(self.pos[:,2]+cm_cor[2])*pos_unit, 
                        v_x=(self.vel[:,0]+cm_cor[3])*vel_unit,
                        v_y=(self.vel[:,1]+cm_cor[4])*vel_unit,
                        v_z=(self.vel[:,2]+cm_cor[5])*vel_unit,
                        frame='galactocentric', representation_type='cartesian', **parameters)
        return snap
        
class BaseParticle(SimpleParticle):
    """ Base particle type of PeTar
        The members include simple particle information, binary status and stellar evolution data
       
    keys: (class members)
        Members inherited from SimpleParticle: mass (1D), pos (2D,3), vel (2D,3) 
            see help(petar.SimpleParticle)
        binary_state (1D): binary interruption state
        if (keyword argument 'interrupt_mode' == 'base', 'bse', 'bseEmp', 'mobse'):
            radius:        (1D): radius for merger checker
            dm:            (1D): mass loss
            time_record    (1D): last time of interruption check
            time_interrupt (1D): next interruption time
        if (keyword argument 'interrupt_mode' == 'bse', 'bseEmp', 'mobse'):
            star  (SSEStarParameter): BSE based stellar evolution parameters
    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): base, bse, mobse, none
               This option indicates whether columns of stellar evolution exist
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys_bstat = [['binary_state',np.int64]]
        keys_se  = [['radius',float_type],['dm',float_type],['time_record',float_type],['time_interrupt',float_type]]    
        
        keys = keys_bstat
        if ('interrupt_mode' in kwargs.keys()):
            if (kwargs['interrupt_mode']=='base'):
                keys = keys_bstat+keys_se
            elif ('bse' in kwargs['interrupt_mode']):
                keys = keys_bstat+keys_se+[['star',SSEStarParameter]]
            
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
               PeTar interrupt mode (set in configure): base, bse, mobse, none
               This option indicates whether columns of stellar evolution exist
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['r_search',float_type], ['id',np.int64], ['mass_bk',np.int64], ['status',np.int64], ['r_in',float_type], ['r_out',float_type]]

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
               PeTar interrupt mode (set in configure): base, bse, mobse, none
               This option indicates whether columns of stellar evolution exist
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
               PeTar interrupt mode (set in configure): base, bse, mobse, none
               This option indicates whether columns of stellar evolution exist
            external_mode: string (none)
               PeTar external mode (set in configure): galpy, none 
               This option indicates whether the column of externa potential exist
            use_mpfrc: bool (False)
               If true, add three columns of pos_high indicating the high-precision parts of position
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['acc_soft',(float_type,3)], ['pot',float_type], ['pot_soft',float_type], ['n_nb',np.int64]]
        if ('external_mode' in kwargs.keys()):
            if (kwargs['external_mode']!='none'):
                keys = [['acc_soft',(float_type,3)], ['pot',float_type], ['pot_soft',float_type], ['pot_ext',float_type], ['n_nb',np.int64]]

        HardParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

    def calcEtot(self):
        """ Calculate total energy and add it as the member, etot
        """
        etot = self.ekin + self.mass*self.pot
        self.addNewMember('etot',etot)

class ParticleGroup(DictNpArrayMix):
    """ A group of particles
    Keys: (class members)
        n (1D): number of particles, when keyword argument N_column_exist=False, this member does not exist
        cm (cm_type): center-of-the-mass particle data 
        p[x] (member_type): particle data, [x] indicate the indice, counting from 0
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        Parameters
        ----------
        keyword arguments:
            member_type: type (SimpleParticle)
                Member particle type
            cm_type: type (SimpleParticle)
                Center-of-the-mass particle type
            N_particle: int (0)
                Number of particles (ignore the value in the column N)
                If data path is provided in the initialization; the column N exists and this argument is not provided, 
                the first value in column N is used to determine N_particle
            N_column_exist: bool (True)
                if True, the class member n exists, otherwise not.
            cm_column_exist: bool (True)
                if True, the center-of-the-mass particle data exist, otherwise not
            float_type: type (np.float64)
                floating point data type
        """

        N_column_exist = True
        if 'N_column_exist' in kwargs.keys():
            N_column_exist = kwargs['N_column_exist']

        n = 0
        if (N_column_exist):
            keys = [['n', np.int64]]
            DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
            if (type(_dat)!=type(None)) & (self.size>0): n = self.n[0]  
        else:
            self.keys = []
            self.ncols = 0

        if 'N_particle' in kwargs.keys(): n = kwargs['N_particle']

        cm_column_exist = True
        if 'cm_column_exist' in kwargs.keys():
            cm_column_exist = kwargs['cm_column_exist']

        cm_type = SimpleParticle
        if ('cm_type' in kwargs.keys()): cm_type = kwargs['cm_type']

        if (cm_column_exist):
            keys_cm = [['cm', cm_type]]
            DictNpArrayMix.__init__(self, keys_cm, _dat, _offset+self.ncols, True, **kwargs)

        member_type = SimpleParticle
        if ('member_type' in kwargs.keys()): member_type = kwargs['member_type']
        
        if (n>0):
            keys_p = [['p'+str(i), member_type] for i in range(n)]
            DictNpArrayMix.__init__(self, keys_p, _dat, _offset+self.ncols, True, **kwargs)
        elif (not N_column_exist):
            self.ncols = 0
            self.size = 0

    def get(self, index):
        """ Get particle by using index
        Parameters
        ----------
        index: int
            particle index
        """

        key = 'p'+str(index)
        return self.__dict__[key]

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
        pcm['pos']  = np.array(list(map(lambda m1,x1,m2,x2:(m1*x1+m2*x2)/(m1+m2), _p1.mass, _p1.pos, _p2.mass, _p2.pos)))
        pcm['vel']  = np.array(list(map(lambda m1,x1,m2,x2:(m1*x1+m2*x2)/(m1+m2), _p1.mass, _p1.vel, _p2.mass, _p2.vel)))
    elif (isinstance(_p1, collections.OrderedDict)) & (isinstance(_p2,collections.OrderedDict)) | (isinstance(_p1, dict)) & (isinstance(_p2, dict)):
        pcm['mass'] = _p1['mass'] + _p2['mass']
        pcm['pos']  = np.array(list(map(lambda m1,x1,m2,x2:(m1*x1+m2*x2)/(m1+m2), _p1['mass'], _p1['pos'], _p2['mass'], _p2['pos'])))
        pcm['vel']  = np.array(list(map(lambda m1,x1,m2,x2:(m1*x1+m2*x2)/(m1+m2), _p1['mass'], _p1['vel'], _p2['mass'], _p2['vel'])))
    else:
        raise ValueError('Initial fail, date type should be Particle or collections.OrderDict, given',type(_p1))

class Binary(SimpleParticle):
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
            member_particle_type: type or list (Particle)
                Type of component particle (both)
            member_particle_type_one: type or list (Particle)
                Type of 1st component
            member_particle_type_two: type or list (Particle)
                Type of 2nd component 
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        G=1
        simple_mode=True
        member_particle_type=Particle
        member_particle_type_one=member_particle_type
        member_particle_type_two=member_particle_type
        
        if 'G' in kwargs.keys(): G=kwargs['G']
        if 'simple_mode' in kwargs.keys(): simple_mode=kwargs['simple_mode']
        if 'member_particle_type' in kwargs.keys(): 
            member_particle_type=kwargs['member_particle_type']
            member_particle_type_one=member_particle_type
            member_particle_type_two=member_particle_type
        if 'member_particle_type_one' in kwargs.keys(): member_particle_type_one=kwargs['member_particle_type_one']
        if 'member_particle_type_two' in kwargs.keys(): member_particle_type_two=kwargs['member_particle_type_two']

        if (issubclass(type(_p1), SimpleParticle)) & (issubclass(type(_p2),SimpleParticle)):
            if (simple_mode): 
                self.keys = [['mass',float_type],['pos',(float_type,3)],['vel',(float_type,3)],['rrel',float_type],['semi',float_type],['ecc',float_type],['p1',(type(_p1),_p1.initargs)], ['p2', (type(_p2),_p2.initargs)]]
                self.particleToSemiEcc(_p1, _p2, G)
                self.ncols= int(10)
            else:
                self.keys = [['mass',float_type],['pos',(float_type,3)],['vel',(float_type,3)],['m1',float_type],['m2',float_type],['rrel',float_type],['semi',float_type],['am',(float_type,3)],['L',(float_type,3)],['eccvec',(float_type,3)],['incline',float_type],['rot_horizon',float_type],['ecc',float_type],['rot_self',float_type],['ecca',float_type],['period',float_type],['t_peri',float_type],['p1',(type(_p1),_p1.initargs)], ['p2', (type(_p2),_p2.initargs)]]
                self.particleToBinary(_p1, _p2, G)
                self.ncols= int(27)
            self.p1 = _p1
            self.p1.setHost(self)
            self.p2 = _p2
            self.p2.setHost(self)
            if (not 'host' in self.__dict__.keys()):
                self.host = None
            self.size = _p1.size
            self.ncols += self.p1.ncols + self.p2.ncols
            self.initargs = kwargs.copy()
            binary_tree = self.createMemberParticleTypeTree()
            self.initargs['member_particle_type_one']=binary_tree[0]
            self.initargs['member_particle_type_two']=binary_tree[1]
        elif (_p2==None):
            type_one = member_particle_type_one
            if (type(member_particle_type_one) == list):
                type_one = (Binary, {'member_particle_type_one':member_particle_type_one[0],'member_particle_type_two':member_particle_type_one[1]})
            type_two = member_particle_type_two
            if (type(member_particle_type_two) == list):
                type_two = (Binary, {'member_particle_type_one':member_particle_type_two[0],'member_particle_type_two':member_particle_type_two[1]})
            if (simple_mode):
                keys = [['rrel',float_type],['semi',float_type],['ecc',float_type],['p1',type_one], ['p2', type_two]]
                SimpleParticle.__init__(self, _p1, _offset, _append, **kwargs)
                DictNpArrayMix.__init__(self, keys, _p1, _offset+self.ncols, True, **kwargs)
            else:
                keys=[['m1',float_type],['m2',float_type],['rrel',float_type],['semi',float_type],['am',(float_type,3)],['L',(float_type,3)],['eccvec',(float_type,3)],['incline',float_type],['rot_horizon',float_type],['ecc',float_type],['rot_self',float_type],['ecca',float_type],['period',float_type],['t_peri',float_type],['p1', type_one],['p2', type_two]]
                SimpleParticle.__init__(self, _p1, _offset, _append, **kwargs)
                DictNpArrayMix.__init__(self, keys, _p1, _offset+self.ncols, True, **kwargs)
            self.initargs = kwargs.copy()
        else:
            raise ValueError('Initial fail, date type should be Particle (2), Binary (1) or no argument (0)')

    def calcEkin(self, member_also=False):
        """ Calculate c.m. kinetic energy, ekin, and add it as a member
        """
        if (self.size>0):
            ekin = 0.5*vecDot(self.vel,self.vel)*self.mass
        else:
            ekin = np.array([])
        self.addNewMember('ekin',ekin)
        if (member_also):
            self.p1.calcEkin()
            self.p2.calcEkin()

    def calcEtot(self, member_also=False):
        """ Calculate c.m. total energy (binary energy is excluded) , etot, and add it as a member
        """
        etot = self.ekin + self.mass*self.pot
        self.addNewMember('etot',etot)
        if (member_also):
            self.p1.calcEtot()
            self.p2.calcEtot()

    def calcR2(self, member_also=False):
        """ Calculate c.m. distance square, r2, and add it as a member
        """
        if (self.size>0):
            r2 = vecDot(self.pos,self.pos)
        else:
            r2 = np.array([])
        self.addNewMember('r2',r2)
        if (member_also):
            self.p1.calcR2()
            self.p2.calcR2()

    def calcEbin(self):
        """ Calculate binding energy, ebin, and add it as a member 
            Notice G should be given the correct value in initialization (keyword argument 'G')
        """
        ebin = self.initargs['G']*self.p1.mass*self.p2.mass/(2*self.semi)
        self.addNewMember('ebin',ebin)

    def calcPot(self):
        """ Calculate potential of c.m., pot, and add it as a member
            Notice G should be given the correct value in initialization (keyword argument 'G')
        """
        G = self.initargs['G']
        pos_b1 = self.p1.pos
        pos_b2 = self.p2.pos
        m_b1 = self.p1.mass
        m_b2 = self.p2.mass
        dr = pos_b1-pos_b2
        dr2 = vecDot(dr,dr)
        invr = 1/np.sqrt(dr2)
        pot_b1 = self.p1.pot + G*m_b2*invr
        pot_b2 = self.p2.pot + G*m_b1*invr
        pot = (m_b2*pot_b1 + m_b1*pot_b2)/self.mass
        self.addNewMember('pot',pot)

    def calcPotExt(self):
        """ Calculate external potential of c.m., pot_ext, and add it as a member
        """
        pot_ext = (self.p1.mass*self.p1.pot_ext + self.p2.mass*self.p2.pot_ext)/self.mass
        self.addNewMember('pot_ext',pot_ext)


    def generateBinaryID(self):
        """ Use CantorPairing to map two components id to one binary id
            Add new member bid 
        """
        bid = cantorPairing(self.p1.id, self.p2.id)
        self.addNewMember('bid',bid)
            
    def correctCenter(self, cm_pos, cm_vel):
        """ Corrent c.m and component position and velocity by subtracting cm_pos and cm_vel
        """
        self.pos -= cm_pos
        self.vel -= cm_vel
        self.p1.correctCenter(cm_pos, cm_vel)
        self.p2.correctCenter(cm_pos, cm_vel)

    def particleToSemiEcc(self, _p1, _p2, _G):
        """ Calculate relative distance, semi-major axis and eccentricity from particle pairs

        Parameters
        ----------
        _p1, _p2: inherited SimpleParticle
            Particle pair data set
        _G: float
            Gravitational constant

        """
        calculateParticleCMDict(self.__dict__, _p1, _p2)

        if ('use_mpfrc' in _p1.initargs.keys()):
            if (_p1.initargs['use_mpfrc']):
                pos1_mp = np.zeros((_p1.size,3),dtype=np.float128)
                pos1_mp += _p1.pos
                pos1_mp += _p1.pos_high
                pos2_mp = np.zeros((_p2.size,3),dtype=np.float128)
                pos2_mp += _p2.pos
                pos2_mp += _p2.pos_high
                dr = (pos1_mp - pos2_mp)
            else:
                dr = (_p1.pos - _p2.pos)
        else:
            dr = (_p1.pos - _p2.pos)
        dv = (_p1.vel - _p2.vel)
        
        dr2  = (dr*dr).sum(axis=1)
        dv2  = (dv*dv).sum(axis=1)
        rvdot= (dr*dv).sum(axis=1)
    
        dr   = np.sqrt(dr2)
        m    = (_p1.mass+_p2.mass)
        semi = 1.0/(2.0/dr - dv2/(_G*m))

        dr_semi = 1.0 - dr/semi
        ecc = np.sqrt(dr_semi*dr_semi + rvdot*rvdot/(_G*m*semi))

        self.rrel = dr
        self.semi = semi
        self.ecc  = ecc

    def particleToBinary(self, _p1, _p2, _G):
        """ Calculate binary orbit from particle pairs

        Parameters
        ----------
        _p1, _p2: inherited SimpleParticle
            Particle pair data set
        _G: float
            Gravitational constant

        """
        binary=self.__dict__
     
        def regular_sign(_a,_a_err):
            _a[(_a<0) & (_a>-_a_err)] *= -1
     
        f_err = 1e-2
        calculateParticleCMDict(binary, _p1, _p2)

        binary['m1'] = _p1.mass
        binary['m2'] = _p2.mass
        m_tot = binary['mass']
        Gm_tot = _G*m_tot
        
        dx = _p1.pos-_p2.pos
        dv = _p1.vel-_p2.vel
        dr2  = vecDot(dx,dx)
        dv2  = vecDot(dv,dv)
        rvdot= vecDot(dx,dv)
        dr   = np.sqrt(dr2)
        binary['rrel'] = np.sqrt(dr2)
     
        inv_dr = 1.0 / binary['rrel']
        binary['semi'] = 1.0 / (2.0*inv_dr - dv2 / Gm_tot)
        binary['am'] = np.cross(dx,dv)
        dp = _p1.vel*_p1.mass[:,None] - _p2.vel*_p2.mass[:,None]
        binary['L'] = np.cross(dx,dp)
        binary['eccvec'] = np.cross(dv,binary['am'])/Gm_tot[:,None]-dx/dr[:,None]
     
        binary['incline'] = np.arctan2(np.sqrt(binary['am'][:,0]*binary['am'][:,0]+binary['am'][:,1]*binary['am'][:,1]),binary['am'][:,2])
        binary['rot_horizon'] = np.arctan2(binary['am'][:,0],-binary['am'][:,1])
        regular_sign(binary['am'][:,0],f_err)
        regular_sign(binary['am'][:,1],f_err)
        #binary['rot_horizon'][binary['rot_horizon']<0] += np.pi
        binary['rot_horizon'][binary['am'][:,1]==0.0]=0.0
     
        cosOMG = np.cos(binary['rot_horizon'])
        sinOMG = np.sin(binary['rot_horizon'])
        cosinc = np.cos(binary['incline'])
        sininc = np.sin(binary['incline'])
     
        pos_bar_x =   dx[:,0]*cosOMG + dx[:,1]*sinOMG
        pos_bar_y = (-dx[:,0]*sinOMG + dx[:,1]*cosOMG)*cosinc + dx[:,2]*sininc
        pos_bar_z = 0.0
        vel_bar_x =   dv[:,0]*cosOMG + dv[:,1]*sinOMG
        vel_bar_y = (-dv[:,0]*sinOMG + dv[:,1]*cosOMG)*cosinc + dv[:,2]*sininc
        vel_bar_z = 0.0
     
        h = np.sqrt(np.sum(binary['am']*binary['am'],axis=1))
        ecccosomg =  h/Gm_tot*vel_bar_y - pos_bar_x*inv_dr
        eccsinomg = -h/Gm_tot*vel_bar_x - pos_bar_y*inv_dr
        binary['ecc'] = np.sqrt( ecccosomg*ecccosomg + eccsinomg*eccsinomg )
        regular_sign(ecccosomg,f_err)
        regular_sign(eccsinomg,f_err)
        binary['rot_self'] = np.arctan2(eccsinomg,ecccosomg)
        #binary['rot_self'][binary['rot_self']<-np.pi+1e-5] += 2*np.pi 
        #binary['rot_self'][binary['rot_self']>=np.pi-1e-5] -= 2*np.pi
     
        regular_sign(pos_bar_y,f_err)
        regular_sign(pos_bar_x,f_err)
        phi = np.arctan2(pos_bar_y, pos_bar_x)
        #phi[phi<-np.pi+1e-5] += 2*np.pi
        #phi[phi>=np.pi-1e-5] -= 2*np.pi
     
        f = phi - binary['rot_self']
        binary['ecca'] = np.arctan(np.sin(f)*np.sqrt(np.abs(binary['ecc']*binary['ecc'] - 1.0))/(binary['ecc']+np.cos(f)))
        n = np.sqrt(Gm_tot/np.abs(binary['semi']*binary['semi']*binary['semi']))
        binary['period'] = 8.0*np.arctan(1.0)/n
        l = binary['ecca'] - binary['ecc']*np.sin(binary['ecca'])
        binary['t_peri'] = l / n

    def createMemberParticleTypeTree(self):
        """ scan the members to create the member particle type tree list
            For example, if the binary structure is a triple: p1: single, p2: binary.
            Then the returned tree list is [particle_typename, [particle_typename, particle_typename]]
        """
        binary_tree=[None,None]
        if (type(self.p1) == Binary):
            binary_tree[0] = self.p1.createMemberParticleTypeTree()
        else:
            binary_tree[0] = type(self.p1)
        if (type(self.p2) == Binary):
            binary_tree[1] = self.p2.createMemberParticleTypeTree()
        else:
            binary_tree[1] = type(self.p2)
        return binary_tree

    def generateBSEInput(self, fpath, time=0.0):
        """
        Generate input for petar.bse
        line: m1, m2, type1, type2, period, ecc, time
        
        Parameters:
        ------------
        fpath: file path to save the input data
        time: time of the snapshot (0.0)
        """
        period = periodToSemi(self.p1.mass, self.p2.mass, self.semi, self.initargs['G'])
        out_data = np.transpose((self.p1.mass, self.p2.mass, self.p1.star.type, self.p2.star.type, period, self.ecc, np.ones(self.size)*time))
        f = open(fpath, 'w')
        f.write('%d\n' % self.size)
        np.savetxt(f, out_data, fmt='%.24g %.24g %d %d %.24g %.24g %.24g')
        f.close()

        
def findPair(_dat, _G, _rmax, use_kdtree=False, simple_binary=True):
    """  Find binaries in a particle data set
    The scipy.spatial.cKDTree is used to find pairs

    Parameters
    ----------
    _dat: inhermited SimpleParticle
        Particle data set
    _G: float
        Gravitational constant
    _rmax: float
        Maximum binary separation
    use_kdtree: bool (False)
        If True, use KDtree to find all binaries (slow); otherwise use information from PeTar, only hard binaries are detected (fast)
    simple_binary: bool (True)
        If True, only calculate semi and ecc (fast); otherwise calculating all binary parameters (slow)

    Return
    ----------
    kdt: KDtree structure if use_kdtree=True
    single: type of _dat
        single particle data set
    binary: Binary(simple_mode=simple_binary, member_particle_type=type(single), G=_G)
        binary data set
    """
    if (not issubclass(type(_dat), SimpleParticle)):
        raise ValueError("Data type wrong",type(_dat)," should be subclass of ", SimpleParticle)

    if (use_kdtree):
        # create KDTree
        #print('create KDTree')
        kdt=sp.cKDTree(_dat.pos)
     
        # find all close pairs
        #pairs=kdt.query_pairs(_rmax*AU2PC)
            
        # only check nearest index
        #pair_index=np.unique(np.transpose(np.array([np.array([x[0],x[1]]) for x in pairs])),axis=0)
         
        # find pair index and distance
        #print('Get index')
        r,index=kdt.query(_dat.pos,k=2)
        pair_index=np.transpose(np.unique(np.sort(index,axis=1),axis=0))
        #pair_index = np.transpose(index)

        #index = kdt.query_pairs(_rmax,output_type='ndarray')
        #pair_index = np.transpose(index)
     
        # two members
        p1 = _dat[pair_index[0]]
        p2 = _dat[pair_index[1]]
     
        # check orbits
        #print('Create binary')
        binary = Binary(p1, p2, G=_G, simple_mode=simple_binary)
        apo =binary.semi*(binary.ecc+1.0)
     
        bsel= ((binary.semi>0) & (apo<_rmax))
        binary = binary[bsel]
        
        single_mask = np.ones(_dat.size).astype(bool)
        single_mask[pair_index[0][bsel]]=False
        single_mask[pair_index[1][bsel]]=False
        single = _dat[single_mask]
        return kdt, single, binary
    else:
        idx = _dat.status.argsort()
        dat_sort = _dat[idx]
        status, index, inverse, counts = np.unique(dat_sort.status, return_index=True, return_inverse=True, return_counts=True)
        binary_i1 = index[counts==2]
        binary_i2 = binary_i1+1
        binary = Binary(dat_sort[binary_i1], dat_sort[binary_i2], _G)
        single = dat_sort[index[-1]:]

        return single, binary

def findMultiple(_single, _binary, _G, _rmax, simple_binary=True):
    """  Find triples and quadruples from single and binary data
    The scipy.spatial.cKDTree is used to find pairs

    Parameters
    ----------
    _single: inhermited SimpleParticle
        Single particle data set
    _binary: Binary
        Binary data set
    _G: float
        Gravitational constant
    _rmax: float
        Maximum binary separation
    simple_binary: bool (True)
        If True, only calculate semi and ecc (fast); otherwise calculating all binary parameters (slow)

    Return
    ----------
    kdt: KDtree structure if use_kdtree=True
    single: type of _dat
        single particle data set
    binary: Binary(simple_mode=simple_binary, member_particle_type=type(single), G=_G)
        binary data set
    triple: Binary(p1: type(single), p2: type(binary), G=_G)
        triple data set
    quadruple: Binary(p1: type(binary), p2: type(binary), G=_G)
        quadruple (binary-binary) data set
    """
    if (not issubclass(type(_single), SimpleParticle)):
        raise ValueError("Data type wrong",type(_single)," should be subclass of ", SimpleParticle)

    single_sin = SimpleParticle(_single)
    binary_sin = SimpleParticle(_binary)
    all_sin = join(single_sin, binary_sin)

    # create KDTree
    kdt=sp.cKDTree(all_sin.pos)
     
    # find pair index and distance
    r,index=kdt.query(all_sin.pos,k=2)
    pair_index=np.transpose(np.unique(np.sort(index,axis=1),axis=0))

    bout_i1 = pair_index[0]
    bout_i2 = pair_index[1]

    Ns = _single.size
    Nb = _binary.size
    quad_pre_sel= (bout_i1>=Ns) & (bout_i2>=Ns)
    tri_pre_sel = (bout_i1<Ns) & (bout_i2>=Ns)
    bin_pre_sel = (bout_i1<Ns) & (bout_i2<Ns)

    n_quad_pre = quad_pre_sel.sum()
    n_tri_pre = tri_pre_sel.sum()
    n_bin_pre = bin_pre_sel.sum()
    if (bout_i1.size != n_quad_pre + n_tri_pre + n_bin_pre):
        raise ValueError('Error: multiple index selection size miss match: dat:',bout_i1.size,'quad:',n_quad_pre,'tri:',n_tri_pre,'bin:',n_bin_pre)

    s_del_index=np.array([]).astype(int)
    b_del_index=np.array([]).astype(int)

    quadruple = Binary(member_particle_type = [type(_single), type(_single)], **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
    if (quad_pre_sel.sum()):
        q1_index = bout_i1[quad_pre_sel]-Ns
        q2_index = bout_i2[quad_pre_sel]-Ns
        quad_pre = Binary(_binary[q1_index], _binary[q2_index], **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
        apo = quad_pre.semi*(quad_pre.ecc+1.0)
        quad_sel = (quad_pre.semi>0) & (apo<_rmax)
        quadruple = quad_pre[quad_sel]
        b_del_index=np.append(q1_index[quad_sel],q2_index[quad_sel])

    triple = Binary(member_particle_type_one = type(_single), 
                    member_particle_type_two = [type(_single), type(_single)], 
                    **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
    if (tri_pre_sel.sum()):
        s_index = bout_i1[tri_pre_sel]
        b_index = bout_i2[tri_pre_sel]-Ns
        tri_pre = Binary(_single[s_index], _binary[b_index], **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
        apo = tri_pre.semi*(tri_pre.ecc+1.0)
        tri_sel = (tri_pre.semi>0) & (apo<_rmax)
        triple = tri_pre[tri_sel]
        b_del_index=np.append(b_del_index,b_index[tri_sel])
        s_del_index=s_index[tri_sel]
        
    bmask=np.ones(Nb).astype(bool)
    if (b_del_index.size>0): bmask[b_del_index]=False;
    binary = _binary[bmask]

    if (bin_pre_sel.sum()):
        s1_index = bout_i1[bin_pre_sel]
        s2_index = bout_i2[bin_pre_sel]
        bin_pre = Binary(_single[s1_index], _single[s2_index], **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
        apo = bin_pre.semi*(bin_pre.ecc+1.0)
        bin_sel = (bin_pre.semi>0) & (apo<_rmax)
        binary.append(bin_pre[bin_sel])
        s_del_index = np.concatenate((s_del_index, s1_index[bin_sel], s2_index[bin_sel]))

    smask=np.ones(Ns).astype(bool)
    smask[s_del_index]=False
    single = _single[smask]

    return single, binary, triple, quadruple

