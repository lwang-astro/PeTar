# read snapshot and obtain multiple systems
import numpy as np
from scipy import spatial as sp

# data class
class ParticleArray:
    def __init__ (self, _dat, _sel=True):
        if (isinstance(_dat,ParticleArray)):
            if (type(_sel)==bool):
                if (_sel==True):
                    self = _dat
            elif (type(_sel)==np.ndarray):
                self.m        = _dat.m       [_sel] 
                self.pos      = _dat.pos     [_sel] 
                self.vel      = _dat.vel     [_sel] 
                self.r_search = _dat.r_search[_sel] 
                self.mass_bk  = _dat.mass_bk [_sel] 
                self.id       = _dat.id      [_sel] 
                self.status   = _dat.status  [_sel] 
                self.r_in     = _dat.r_in    [_sel] 
                self.r_out    = _dat.r_out   [_sel] 
                self.acc      = _dat.acc     [_sel] 
                self.pot      = _dat.pot     [_sel] 
                self.n_nb     = _dat.n_nb    [_sel] 
            else:
                raise ValueError('Initial fail: _sel type error',type(_sel))
        elif (type(_dat)==np.ndarray):
            dat = _dat
            if (type(_sel)==bool):
                if (_sel==True):    
                    dat = _dat
            elif (type(_sel)==np.ndarray):
                dat = _dat[_sel]
            else:
                raise ValueError('Initial fail: _sel type error',type(_sel))
            self.m        = dat[:,0]
            self.pos      = dat[:,1:4]
            self.vel      = dat[:,4:7]
            self.r_search = dat[:,7]
            self.mass_bk  = dat[:,8]
            self.id       = dat[:,9]
            self.status   = dat[:,10]
            self.r_in     = dat[:,11]
            self.r_out    = dat[:,12]
            self.acc      = dat[:,13:16]
            self.pot      = dat[:,16]
            self.n_nb     = dat[:,17]
        else:
            raise ValueError("Initial fail, date type error",type(_dat))

# Unit class for scaling
class Units:
    G = 1
    r = 1
    v = 1
    m = 1

# calculate binary orbit from particles
# _p1, _p2: data class
# _units: unit class
# return: semi, ecc
def particleToBinary(_p1,_p2, _units):
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

# _dat: snapshot array
# return: members, semi, ecc, separation
def findPair(_dat):
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
