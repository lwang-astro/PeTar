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
            self.m        = dat[:,0]
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

def particleToBinary(_p1,_p2, _units):
    """
    calculate binary orbit from particles
    _p1, _p2: data class
    _units: unit class
    return: semi, ecc
    """
    dr = (_p1.pos - _p2.pos)*_units.r
    dv = (_p1.vel - _p2.vel)*_units.v
    
    dr2  = (dr*dr).sum(axis=1)
    dv2  = (dv*dv).sum(axis=1)
    rvdot= (dr*dv).sum(axis=1)
    
    dr   = np.sqrt(dr2)
    m    = (_p1.m+_p2.m)*_units.m
    semi = 1.0/(2.0/dr - dv2/(_units.G*m))

    dr_semi = 1.0 - dr/semi
    ecc = np.sqrt(dr_semi*dr_semi + rvdot*rvdot/(_units.G*m*semi))
    return semi,ecc

def findPair(_dat):
    """
    Find paris
    _dat: snapshot array
    return: members, semi, ecc, separation
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
    semi, ecc = particleToBinary(p1, p2, unit)
    
    return p1, p2, semi, ecc,r
