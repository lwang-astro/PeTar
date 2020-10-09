# read snapshot and obtain multiple systems
import collections
from scipy import spatial as sp
from .base import *
from .bse import *

G_MSUN_PC_MYR=0.00449830997959438 # Msun, pc, myr
G_HENON=1 # Henon unit
HEADER_OFFSET=24 # header offset in bytes for snapshots with the BINARY format

class PeTarDataHeader():
    """ Petar snapshot data header
    members:
        fid: file id
        n: number of particles
        time: time of snapshot
    """

    def __init__(self, _filename=None, **kwargs):
        """ Initial data header
        
        Parameters:
        -----------
        _filename: string
            PeTar snapshot file name to read the header, if not provide, all members are initialized to zero (None)
        kwargs: dict
            Keyword arguments:
            snapshot_format: string (ascii)
                Data format of snapshot files: binary or ascii
        """
        self.fid = int(0)
        self.n = int(0)
        self.time = 0.0
        
        if (_filename!=None): self.read(_filename,**kwargs)

    def read(self, _filename, **kwargs):
        """ Read snapshot file to obtain the header information

        Parameters:
        -----------
        _filename: string
            PeTar snapshot file name to read the header
        kwargs: dict
            Keyword arguments:
            snapshot_format: string (ascii)
                Data format of snapshot files: binary or ascii
        """
        snapshot_format='ascii'
        if ('snapshot_format' in kwargs.keys()): snapshot_format=kwargs['snapshot_format']

        if (snapshot_format=='ascii'):
            fp = open(_filename, 'r')
            header=fp.readline()
            file_id, n_glb, t = header.split()
            fp.close()

            self.fid = int(file_id)
            self.n = int(n_glb)
            self.time = float(t)

        elif (snapshot_format=='binary'):
            fp = np.fromfile(_filename, dtype=np.dtype([('file_id',np.int64),('n_glb',np.int64),('time',np.float64)]),count=1)
            self.file_id = fp['file_id'][0]
            self.n_glb = fp['n_glb'][0]
            self.time = fp['time'][0]
        else: 
            raise ValueError('Snapshot format unknown, should be binary or ascii, given', snapshot_format)

class SimpleParticle(DictNpArrayMix):
    """ Simple particle class with only mass, postion, velocity
    keys: (class members)
        mass (1D): mass
        pos (2D,3): postion x, y, z
        vel (2D,3): velocity vx, vy, vz
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['mass', np.float64], ['pos', (np.float64, 3)], ['vel', (np.float64, 3)]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def calcR2(self):
        """ calculate distance square, r2, and add it as a class member
        """
        if (not 'r2' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['r2',1])
        self.r2 = vecDot(self.pos,self.pos)

    def calcEkin(self):
        """ calculate kinetic energy
        """
        if (not 'ekin' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['ekin',1])
        self.ekin = 0.5*vecDot(self.vel,self.vel)*self.mass

    def correctCenter(self, cm_pos, cm_vel):
        self.pos -= cm_pos
        self.vel -= cm_vel

class Particle(SimpleParticle):
    """ Particle class 
    keys: (class members)
        The final keys are a combination of sub keys depending on kwargs of initial function

        Sub key list:
        basic: [inherit SimpleParticle]
        add: binary_state: binary interruption state 
        se: radius:        (1D): radius for merger checker
            dm:            (1D): mass loss
            time_record    (1D): last time of interruption check
            time_interrupt (1D): next interruption time
        ptcl: r_search (1D): searching radius
              id       (1D): identification
              mass_bk  (1D): artificial particle parameter 1 
              status   (1D): artificial particle parameter 2
              r_in     (1D): changeover function inner boundary
              r_out    (1D): changeover function outer boundary
        hermite: dt    (1D): time step
                 time  (1D): current time
                 acc   (2D,3): acceleration x, y, z
                 jerk  (2D,3): acceleration derivative x, y, z
                 pot   (1D): potential
        soft: acc_soft (2D,3): long-range interaction acceleration (particle-tree) x, y, z
              pot      (1D): total potential
              pot_soft (1D): long-range interaction potential
              n_nb:    (1D): number of neighbors (short-interaction)

        Combination: 
        ends:
            kwargs['particle_type']:
                hermite:   ptcl + hermite
                hard:      ptcl
                soft (default): ptcl + soft
        Final:
        keys:
            kwargs['interrupt_mode']:
                base:      basic + add + se + ends
                bse:       basic + add + se + ['star',SSEStarParameter] + ends
                none (default): basic + add + ends

    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            particle_type: basic particle type: hermite, hard, soft (soft)
            interrupt_mode: PeTar interrupt mode: base, bse, none (none)
        """

        keys_add = [['binary_state',np.int64]]
        keys_se  = [['radius',np.float64],['dm',np.float64],['time_record',np.float64],['time_interrupt',np.float64]]
        keys_ptcl_add = [['r_search',np.float64], ['id',np.int64], ['mass_bk',np.int64], ['status',np.int64], ['r_in',np.float64], ['r_out',np.float64]]
        keys_hermite_add = [['dt',np.float64],['time',np.float64],['acc',(np.float64,3)],['jerk',(np.float64,3)],['pot',np.float64]]
        keys_soft_add = [['acc_soft',(np.float64,3)], ['pot',np.float64], ['pot_soft',np.float64], ['n_nb',np.int64]]
        keys_end =  keys_ptcl_add + keys_soft_add
        if ('particle_type' in kwargs.keys()):
            if (kwargs['particle_type']=='hermite'):
                keys_end = keys_ptcl_add + keys_hermite_add
            elif (kwargs['particle_type']=='hard'):
                keys_end = keys_ptcl_add
        keys=keys_add+keys_end
        if ('interrupt_mode' in kwargs.keys()):
            if (kwargs['interrupt_mode']=='base'):
                keys = keys_add+keys_se+keys_end
            elif (kwargs['interrupt_mode']=='bse'):
                keys = keys_add+keys_se+[['star',SSEStarParameter]]+keys_end
            
        SimpleParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

    def calcEtot(self):
        """ Calculate total energy and add it as the member, etot
        """
        if (not 'etot' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['etot',np.float64])
        self.etot = self.ekin + self.mass*self.pot

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
                p1   (member_particle_type) component one
                p2   (member_particle_type) component two
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
                p1 (member_particle_type) component one
                p2 (member_particle_type) component two

        member_particle_type is given by kwargs['member_particle_type'], in default it is SimpleParticle
               
    """
    def __init__ (self, _p1=None, _p2=None, _offset=int(0), _append=False, **kwargs):
        """
        Parameters
        ----------
        _p1: inherited SimpleParticle | 2D numpy.ndarray | Binary | None
            If the type is inherited SimpleParticle, it is the first component of binary (_p2 should be the same type).
            If the type is Binary, the class instance is initialized by copy the data of _p1.
            If it is None, initialize class with empty data
        _p2: inherited SimpleParticle | None
            If the type is inherited SimpleParticle, it is the second component of binary 
            If it is None, _p1 should be either 2D numpy.ndarray or Bina
        _offset: int (0)
            Reading column offset of _dat if it is 2D np.ndarray
        _append: bool (False)
            If true, append keys and ncols to the current class instead of create new class members
        kwaygs: dict ()
            keyword arguments:
                simple_mode: only calculate semi and ecc, save computing time significantly (True)
                G: gravitational constant (1.0)
                member_particle_type: type of component particle (SimpleParticle)
        """
        G=1
        simple_mode=True
        member_particle_type=SimpleParticle
        
        if 'G' in kwargs.keys(): G=kwargs['G']
        if 'simple_mode' in kwargs.keys(): simple_mode=kwargs['simple_mode']
        if 'member_particle_type' in kwargs.keys(): member_particle_type=kwargs['member_particle_type']

        if (issubclass(type(_p1), SimpleParticle)) & (issubclass(type(_p2),SimpleParticle)):
            member_particle_type = type(_p1)
            if (simple_mode): 
                self.keys = [['mass',np.float64],['pos',(np.float64,3)],['vel',(np.float64,3)],['rrel',np.float64],['semi',np.float64],['ecc',np.float64],['p1',member_particle_type], ['p2', member_particle_type]]
                self.particleToSemiEcc(_p1, _p2, G)
                self.ncols= int(10)
            else:
                self.keys = [['mass',np.float64],['pos',(np.float64,3)],['vel',(np.float64,3)],['m1',np.float64],['m2',np.float64],['rrel',np.float64],['semi',np.float64],['am',(np.float64,3)],['L',(np.float64,3)],['eccvec',(np.float64,3)],['incline',np.float64],['rot_horizon',np.float64],['ecc',np.float64],['rot_self',np.float64],['ecca',np.float64],['period',np.float64],['t_peri',np.float64],['p1', member_particle_type],['p2', member_particle_type]]
                self.particleToBinary(_p1, _p2, G)
                self.ncols= int(27)
            self.p1 = _p1
            self.p2 = _p2
            self.size = _p1.size
            self.ncols += self.p1.ncols + self.p2.ncols
        elif (_p2==None):
            if (simple_mode):
                keys = [['rrel',np.float64],['semi',np.float64],['ecc',np.float64],['p1',member_particle_type], ['p2', member_particle_type]]
                SimpleParticle.__init__(self, _p1, _offset, _append, **kwargs)
                DictNpArrayMix.__init__(self, keys, _p1, _offset+self.ncols, True, **kwargs)
            else:
                keys=[['m1',np.float64],['m2',np.float64],['rrel',np.float64],['semi',np.float64],['am',(np.float64,3)],['L',(np.float64,3)],['eccvec',(np.float64,3)],['incline',np.float64],['rot_horizon',np.float64],['ecc',np.float64],['rot_self',np.float64],['ecca',np.float64],['period',np.float64],['t_peri',np.float64],['p1', member_particle_type],['p2', member_particle_type]]
                SimpleParticle.__init__(self, _p1, _offset, _append, **kwargs)
                DictNpArrayMix.__init__(self, keys, _p1, _offset+self.ncols, True, **kwargs)
        else:
            raise ValueError('Initial fail, date type should be Particle (2), Binary (1) or no argument (0)')
        self.initargs = kwargs.copy()
        self.initargs['G'] = G

    def calcEkin(self):
        """ Calculate c.m. kinetic energy, ekin, and add it as a member
        """
        if (not 'ekin' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['ekin',np.float64])
        self.ekin = 0.5*vecDot(self.vel,self.vel)*self.mass

    def calcEtot(self):
        """ Calculate c.m. total energy (binary energy is excluded) , etot, and add it as a member
        """
        if (not 'etot' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['etot',np.float64])
        self.etot = self.ekin + self.mass*self.pot

    def calcR2(self, member_also=False):
        """ Calculate c.m. distance square, r2, and add it as a member
        """
        if (not 'r2' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['r2',np.float64])
        self.r2 = vecDot(self.pos,self.pos)
        if (member_also):
            ncols = self.p1.ncols + self.p2.ncols
            self.p1.calcR2()
            self.p2.calcR2()
            ncols = self.p1.ncols + self.p2.ncols - ncols
            self.ncols += ncols

    def calcEbin(self):
        """ Calculate binding energy, ebin, and add it as a member 
            Notice G should be given the correct value in initialization (keyword argument 'G')
        """
        if (not 'ebin' in self.__dict__.keys()):
            self.ncols += 1
            self.keys.append(['ebin',np.float64])
        self.ebin = self.initargs['G']*self.p1.mass*self.p2.mass/(2*self.semi)

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
        if (not 'pot' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['pot',np.float64])
        self.pot = (m_b2*pot_b1 + m_b1*pot_b2)/self.mass
            
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
        binary['am'] = np.array(list(map(lambda x,y:np.cross(x,y),dx,dv)))
        dp = np.array(list(map(lambda m1,x1,m2,x2:m1*x1-m2*x2,_p1.mass,_p1.vel,_p2.mass,_p2.vel)))
        binary['L'] = np.array(list(map(lambda x,y:np.cross(x,y),dx,dp)))
        binary['eccvec'] = np.array(list(map(lambda v,am,gm,dx,dr:np.cross(v,am)/gm-dx/dr,dv,binary['am'],Gm_tot,dx,dr)))
     
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
     
        h = np.array(list(map(lambda x:np.sqrt(np.inner(x,x)),binary['am'])))
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

    # two members
    p1 = all_sin[pair_index[0]]
    p2 = all_sin[pair_index[1]]
     
    # check orbits
    binary_out = Binary(p1, p2, G=_G, simple_mode=simple_binary)
    apo =binary_out.semi*(binary_out.ecc+1.0)
     
    bsel= ((binary_out.semi>0) & (apo<_rmax))
    bout_i1=pair_index[0][bsel]
    bout_i2=pair_index[1][bsel]
    # check single or cm
    Ns = _single.size
    Nb = _binary.size
    quad_sel= (bout_i1>=Ns) & (bout_i2>=Ns)
    tri_sel = (bout_i1<Ns) & (bout_i2>=Ns)
    bin_sel = (bout_i1<Ns) & (bout_i2<Ns)

    if (bout_i1.size != quad_sel.sum()+tri_sel.sum()+bin_sel.sum()):
        raise ValueError('Error: select size miss match: dat:',bout_i1.size,'quad:',quad_sel.sum(),'tri:',tri_sel.sum(),'bin:',bin_sel.sum())

    multiple = binary_out[bsel]

    quadruple = multiple[quad_sel]
    ncol_diff_bin = _binary.ncols - binary_sin.ncols
    ncol_diff_sin = _single.ncols - single_sin.ncols
    if (quad_sel.sum()):
        q1_index = bout_i1[quad_sel]-Ns
        quadruple.p1 = _binary[q1_index]
        q2_index = bout_i2[quad_sel]-Ns
        quadruple.p2 = _binary[q2_index]
        quadruple.initargs['member_particle_type'] = Binary
        quadruple.ncols += 2*ncol_diff_bin
        if (quadruple.size != quadruple.p1.size):
            raise ValueError('Error: quadruple size', quadruple.size,' mismatch member size ',quadruple.p1.size)

    triple = multiple[tri_sel]
    if (tri_sel.sum()):
        s_index = bout_i1[tri_sel]
        triple.p1 = _single[s_index]
        b_index = bout_i2[tri_sel]-Ns
        triple.p2 = _binary[b_index]
        triple.ncols += ncol_diff_bin + ncol_diff_sin
        if (triple.size != triple.p1.size):
            raise ValueError('Error: triple size', triple.size,' mismatch member size ',triple.p1.size)

    binary_new = multiple[bin_sel]
    if (bin_sel.sum()):
        s1_index = bout_i1[bin_sel]
        binary_new.p1 = _single[s1_index]
        s2_index = bout_i2[bin_sel]
        binary_new.p2 = _single[s2_index]
        binary_new.ncols += 2*ncol_diff_sin
        binary_new.initargs['member_particle_type'] = type(_single)
        if (binary_new.size != binary_new.p1.size):
            raise ValueError('Error: binary_new size', binary_new.size,' mismatch member size ',binary_new.p1.size)

        if (binary_new.ncols != _binary.ncols): 
            raise ValueError('Error: old binary ncols', _binary.ncols,' mismatch new binary ncols', binary_new.ncols)

    # delete used singles
    s_del_index = np.append(bout_i1[bout_i1<Ns],bout_i2[bout_i2<Ns])
    smask=np.ones(Ns).astype(bool)
    smask[s_del_index]=False
    single = _single[smask]

    # delete used binaries
    b_del_index = np.append(bout_i1[bout_i1>=Ns],bout_i2[bout_i2>=Ns])-Ns
    bmask=np.ones(Nb).astype(bool)
    bmask[b_del_index]=False
    binary = _binary[bmask]
    if (bmask.sum()+b_del_index.size != Nb):
        raise ValueError('Error: binary size not match: origin:',Nb,'remain:',bmask.sum(),'del:',b_del_index.size,'del(unique):',np.unique(b_del_index).size)

    if (binary_new.size>0):
        binary = join(binary,binary_new)

    return single, binary, triple, quadruple
