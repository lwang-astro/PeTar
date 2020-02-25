# read snapshot and obtain multiple systems
import numpy as np
from scipy import spatial as sp

class ParticleArray:
    """ Particle class 
    """
    def __init__ (self, _dat=np.empty([0,18]), _sel=True):
        """
        _dat: np.ndarray type data reading from snapshot or ParticleArray type
        _sel: if the value is true(bool), copy the data, else give the address (if _dat is ParticleArray type); if the value is bool array (same size as the particle number)CCCCCC, select the data based on array
        """
        # if _dat is same type
        if (isinstance(_dat,ParticleArray)):
            if (type(_sel)==bool):
                if (_sel==True):
                    self = _dat.copy()
            elif (type(_sel)==np.ndarray):
                for key in _dat.__dict__.keys():
                    self.__dict__[key] = _dat.__dict__[key][_sel]
            else:
                raise ValueError('Initial fail: _sel type should be bool or np.ndarray, given ',type(_sel))
        # if _dat is array
        elif (type(_dat)==np.ndarray):
            dat = _dat
            if (type(_sel)==bool):
                if (_sel==True):    
                    dat = _dat
            elif (type(_sel)==np.ndarray):
                dat = _dat[_sel]
            else:
                raise ValueError('Initial fail: _sel type should be bool or np.ndarray, given ',type(_sel))
            self.mass     = dat[:,0]
            self.pos      = dat[:,1:4]
            self.vel      = dat[:,4:7]
            self.r_search = dat[:,7]
            self.id       = dat[:,8]
            self.mass_bk  = dat[:,9]
            self.status   = dat[:,10]
            self.r_in     = dat[:,11]
            self.r_out    = dat[:,12]
            self.acc      = dat[:,13:16]
            self.pot      = dat[:,16]
            self.n_nb     = dat[:,17]
        else:
            raise ValueError('Initial fail, date type should be ParticleArray or np.ndarray, given ',type(_dat))

class Binary:
    """ Binary class
    """
    def __init__ (self, _p1, _p2, _G):
        if (isinstance(_p1,ParticleArray)) & (isinstance(_p2,ParticleArray)):
            self.__dict__ = particleToBinary(_p1.__dict__, _p2.__dict__, _G)
        elif (isinstance(_p1,collections.OrderedDict)) & (isinstance(_p2,collections.OrderedDict)):
            self.__dict__ = particleToBinary(_p1, _p2, _G)

    def particleToBinary(_p1, _p2, _G):
        """ 
        Calculate binary orbit from particle pairs
        _p1, _p2: particle array of member 1 and 2
        _G: gravitaitonal constant
        return: binary dicto
        """

        def regular_sign(_a,_a_err):
            _a[(_a<0) & (_a>-_a_err)] *= -1

        f_err = 1e-2
        binary=collections.OrderedDict()
        binary['m1'] = _p1['mass']
        binary['m2'] = _p2['mass']
        m_tot = _p1['mass'] + _p2['mass']
        Gm_tot = _G*m_tot
        
        dx = _p1['pos']-_p2['pos']
        dv = _p1['vel']-_p2['vel']
        dr2  = np.array(list(map(lambda x:np.inner(x,x),dx)))
        dv2  = np.array(list(map(lambda x:np.inner(x,x),dv)))
        rvdot= np.array(list(map(lambda x,y:np.inner(x,y),dx,dv)))
        dr   = np.sqrt(dr2)
        binary['r'] = np.sqrt(dr2)
    
        inv_dr = 1.0 / binary['r']
        binary['semi'] = 1.0 / (2.0*inv_dr - dv2 / Gm_tot)
        binary['am'] = np.array(list(map(lambda x,y:np.cross(x,y),dx,dv)))
        dp = np.array(list(map(lambda m1,x1,m2,x2:m1*x1-m2*x2,_p1['mass'],_p1['vel'],p2['mass'],_p2['vel'])))
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
        sinu = binary['r']*np.sin(f) / (binary['semi']*np.sqrt(1.0 - binary['ecc']*binary['ecc']))
        cosu = (binary['r']*np.cos(f) / binary['semi']) + binary['ecc']
        binary['ecca'] = np.arctan2(sinu,cosu)
        n = np.sqrt(Gm_tot/(binary['semi']*binary['semi']*binary['semi']))
        binary['period'] = 8.0*np.arctan(1.0)/n
        l = binary['ecca'] - binary['ecc']*np.sin(binary['ecca'])
        binary['t_peri'] = l / n

        return binary

def JoinParticleArray(*_dat):
    """
    Join multiple particle array to one
    """
    for idat in _dat:
        if (not isinstance(idat,ParticleArray)):
            raise ValueError('Initial fail, date type should be ParticleArray or np.ndarray, given ',type(idat))
    keys = _dat[0].__dict__.keys()
    new_dat = ParticleArray()
    for key in keys:
        new_dat.__dict__[key] = np.concatenate(tuple(map(lambda x:x.__dict__[key], _dat)))
    return new_dat

class Units:
    """
    Unit class for scaling
    """
    G = 1
    r = 1
    v = 1
    m = 1

def particleToSemiEcc(_p1,_p2, _G):
    """
    calculate binary semi-major axis and eccentricity from particle pairs
    _p1, _p2: data class
    _G: gravitational constant
    return: semi, ecc
    """
    dr = (_p1.pos - _p2.pos)*_units.r
    dv = (_p1.vel - _p2.vel)*_units.v
    
    dr2  = (dr*dr).sum(axis=1)
    dv2  = (dv*dv).sum(axis=1)
    rvdot= (dr*dv).sum(axis=1)
    
    dr   = np.sqrt(dr2)
    m    = (_p1.mass+_p2.mass)
    semi = 1.0/(2.0/dr - dv2/(_G*m))

    dr_semi = 1.0 - dr/semi
    ecc = np.sqrt(dr_semi*dr_semi + rvdot*rvdot/(_G*m*semi))
    return semi, ecc


def findPair(_dat):
    """
    Find paris
    _dat: snapshot array
    return: member 1, member 2, binary
    """
    if (not isinstance(_dat,ParticleArray)):
        raise ValueError("Data type wrong",type(_dat))
    
    # create KDTree
    kdt=sp.KDTree(_dat.pos)

    # find all close pairs
    #pairs=kdt.query_pairs(_rmax*AU2PC)
        
    # only check nearest index
    #pair_index=np.unique(np.transpose(np.array([np.array([x[0],x[1]]) for x in pairs])),axis=0)
     
    # find pair index and distance
    r,index=kdt.query(_dat.pos,k=2)

    pair_index=np.transpose(np.unique(np.sort(index,axis=1),axis=0))

    # two members
    p1 = ParticleArray(_dat,pair_index[0])
    p2 = ParticleArray(_dat,pair_index[1])
            
    # check orbits
    unit = Units()
    binary = Binary(p1, p2, unit)
    
    return p1, p2, binary
