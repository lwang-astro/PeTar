# read snapshot and obtain multiple systems
import numpy as np
import collections
from scipy import spatial as sp
from tools.base import DictNpArrayMix, InitialDictNpArrayMixMethod

class Particle(DictNpArrayMix):
    """ Particle class 
    """
    @InitialDictNpArrayMixMethod
    def __init__ (self):
        return [['mass',1], ['pos',3], ['vel',3], ['r_search',1], ['id',1], ['mass_bk',1], ['status',1], ['r_in',1], ['r_out',1], ['acc',3], ['pot',1], ['n_nb',1]]

class Binary(DictNpArrayMix):
    """ Binary class
    """
    def __init__ (self, _p1=0, _p2=0, _G=0):
        if (isinstance(_p1, Particle)) & (isinstance(_p2,Particle)):
            self.particleToBinary(_p1.__dict__, _p2.__dict__, _G)
        elif (isinstance(_p1, collections.OrderedDict)) & (isinstance(_p2,collections.OrderedDict)):
            self.particleToBinary(_p1, _p2, _G)
        elif (_p2==0) & (_G==0) & (type(_p1)==Binary):
            self = _p1.copy()
        elif (_p1==0) & (_p2==0) & (_G==0):
            keys=[['m1',1],['m2',1],['r',1],['semi',1],['am',3],['L',3],['eccvec',3],['incline',1],['rot_horizon',1],['ecc',1],['rot_self',1],['ecca',1],['period',1],['t_peri',1]]
            for key,dimension in keys:
                self.__dict__[key] = np.empty([0,dimension])
        else:
            raise ValueError('Initial fail, date type should be Particle (2), Binary (1) or no argument (0)')
            

    def particleToBinary(self, _p1, _p2, _G):
        """ 
        Calculate binary orbit from particle pairs
        _p1, _p2: particle array of member 1 and 2
        _G: gravitaitonal constant
        return: binary dicto
        """
     
        def regular_sign(_a,_a_err):
            _a[(_a<0) & (_a>-_a_err)] *= -1
     
        f_err = 1e-2
        binary=self.__dict__
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
        dp = np.array(list(map(lambda m1,x1,m2,x2:m1*x1-m2*x2,_p1['mass'],_p1['vel'],_p2['mass'],_p2['vel'])))
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

def findPair(_dat, _G):
    """
    Find paris
    _dat: Particle type data 
    return: member 1, member 2, binary
    """
    if (not isinstance(_dat,Particle)):
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
    p1 = _dat[pair_index[0]]
    p2 = _dat[pair_index[1]]
            
    # check orbits
    binary = Binary(p1, p2, _G)
    
    return p1, p2, binary

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

