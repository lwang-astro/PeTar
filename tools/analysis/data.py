# read snapshot and obtain multiple systems
import collections
from scipy import spatial as sp
from .base import *

class SimpleParticle(DictNpArrayMix):
    """ Simple particle class with only mass, postion, velocity and r2
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['mass',1], ['pos',3], ['vel',3]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def calcR2(self):
        """ calculate distance square
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
    """
    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys_se  = [['radius',1],['dm',1],['time_record',1],['time_interrupt',1],['binary_state',1]]
        keys_bse = [['s_type',1],['s_mass0',1],['s_mass',1],['s_rad',1],['s_mcore',1],['s_rcore',1],['s_spin',1],['s_epoch',1],['s_time',1],['s_lum',1]]
        keys_std = [['r_search',1], ['id',1], ['mass_bk',1], ['status',1], ['r_in',1], ['r_out',1], ['acc_soft',3], ['pot',1], ['pot_soft',1], ['n_nb',1]]
        keys=keys_std
        if ('interrupt_mode' in kwargs.keys()):
            if (kwargs['interrupt_mode']=='base'):
                keys = keys_se+keys_std
            elif (kwargs['interrupt_mode']=='bse'):
                keys = keys_se+keys_bse+keys_std
            
        SimpleParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

    def calcEtot(self):
        if (not 'etot' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['etot',1])
        self.etot = self.ekin + self.mass*self.pot

def calculateParticleCMDict(pcm, _p1, _p2):
    """ calculate cm of particle pair"""
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

class Binary(DictNpArrayMix):
    """ Binary class
    """
    def __init__ (self, _p1=None, _p2=None, _offset=int(0), _append=False, **kwargs):
        """
        simple_mode: only calculate semi and ecc
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
                self.keys = [['mass',1],['pos',3],['vel',3],['rrel',1],['semi',1],['ecc',1],['p1',member_particle_type], ['p2', member_particle_type]]
                self.particleToSemiEcc(_p1, _p2, G)
                self.ncols= int(10)
            else:
                self.keys = [['mass',1],['pos',3],['vel',3],['m1',1],['m2',1],['rrel',1],['semi',1],['am',3],['L',3],['eccvec',3],['incline',1],['rot_horizon',1],['ecc',1],['rot_self',1],['ecca',1],['period',1],['t_peri',1],['p1', member_particle_type],['p2', member_particle_type]]
                self.particleToBinary(_p1, _p2, G)
                self.ncols= int(27)
            self.p1 = _p1
            self.p2 = _p2
            self.size = _p1.size
            self.ncols += self.p1.ncols + self.p2.ncols
        elif (_p2==None):
            if (simple_mode):
                keys = [['mass',1],['pos',3],['vel',3],['rrel',1],['semi',1],['ecc',1],['p1',member_particle_type], ['p2', member_particle_type]]
                DictNpArrayMix.__init__(self, keys, _p1, _offset, _append, **kwargs)
            else:
                keys=[['mass',1],['pos',3],['vel',3],['m1',1],['m2',1],['rrel',1],['semi',1],['am',3],['L',3],['eccvec',3],['incline',1],['rot_horizon',1],['ecc',1],['rot_self',1],['ecca',1],['period',1],['t_peri',1],['p1', member_particle_type],['p2', member_particle_type]]
                DictNpArrayMix.__init__(self, keys, _p1, _offset, _append, **kwargs)
        else:
            raise ValueError('Initial fail, date type should be Particle (2), Binary (1) or no argument (0)')
        self.initargs = kwargs.copy()

    def calcEkin(self):
        """ calculate kinetic energy
        """
        if (not 'ekin' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['ekin',1])
        self.ekin = 0.5*vecDot(self.vel,self.vel)*self.mass

    def calcEtot(self):
        if (not 'etot' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['etot',1])
        self.etot = self.ekin + self.mass*self.pot

    def calcR2(self, member_also=False):
        """ calculate distance square
        """
        if (not 'r2' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['r2',1])
        self.r2 = vecDot(self.pos,self.pos)
        if (member_also):
            ncols = self.p1.ncols + self.p2.ncols
            self.p1.calcR2()
            self.p2.calcR2()
            ncols = self.p1.ncols + self.p2.ncols - ncols
            self.ncols += ncols

    def calcPot(self, G):
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
            self.keys.append(['pot',1])
        self.pot = (m_b2*pot_b1 + m_b1*pot_b2)/self.mass
            
    def correctCenter(self, cm_pos, cm_vel):
        self.pos -= cm_pos
        self.vel -= cm_vel
        self.p1.correctCenter(cm_pos, cm_vel)
        self.p2.correctCenter(cm_pos, cm_vel)

    def particleToSemiEcc(self, _p1,_p2, _G):
        """
        calculate binary semi-major axis and eccentricity from particle pairs
        _p1, _p2: data class
        _G: gravitational constant
        return: semi, ecc
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
        """ 
        Calculate binary orbit from particle pairs
        _p1, _p2: particle array of member 1 and 2
        _G: gravitaitonal constant
        return: binary dicto
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
    """
    Find paris
    _dat: Particle type data 
    _G: gravitational constant
    _rmax: maximum binary separation
    use_kdtree: use KDtree to find all binaries (slow); otherwise use information from PeTar, only hard binaries are detected (fast)
    simple_binary: only calculate semi and ecc (fast); otherwise calculating all binary parameters (slow)
    return: [KDtree], single, binary
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
        binary = Binary(p1, p2, G=_G)
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

