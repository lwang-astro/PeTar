import collections
import numpy as np
from .base import *
from .data import *

class Core(DictNpArrayMix):
    """ center of mass 
    """
    
    def __init__(self, _dat=None, _offset=int(0)):
        keys  = [['time',1],['pos', 3],['vel', 3], ['rc', 1]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset)

    def calcPotentialCenter(self, single, binary, G):
        """
        get potential center of the system
        r_cm = \sum_i pot_i *r_i /\sum_i pot_i (only count pot_i <0)
        """
        pot_s_sel = (single.pot<0)
        pot_s = single.pot[pot_s_sel]
        pos_s = single.pos[pot_s_sel]
        vel_s = single.vel[pot_s_sel]
        
        pos_b1 = binary.p1.pos
        pos_b2 = binary.p2.pos
        m_b1 = binary.p1.mass
        m_b2 = binary.p2.mass
        dr = pos_b1-pos_b2
        dr2 = vec_dot(dr,dr)
        invr = 1/np.sqrt(dr2)
        pot_b1 = binary.p1.pot + G*m_b2*invr
        pot_b2 = binary.p2.pot + G*m_b1*invr
        pot_b = (m_b2*pot_b1 + m_b1*pot_b2)/binary.mass
        pos_b = binary.pos
        vel_b = binary.vel
        
        pot_sum = pot_s.sum() + pot_b.sum()
        cm_pos = np.array([(np.sum(pot_s*pos_s[:,i])+np.sum(pot_b*pos_b[:,i]))/pot_sum for i in range(3)])
        cm_vel = np.array([(np.sum(pot_s*vel_s[:,i])+np.sum(pot_b*vel_b[:,i]))/pot_sum for i in range(3)])

        self.pos = np.append(self.pos,[cm_pos],axis=0)
        self.vel = np.append(self.vel,[cm_vel],axis=0)

        return cm_pos, cm_vel
    
    def calcDensityAndCenter(self, particle, kdtree):
        """ calculate density based on six neighbors and return density center
        particle: all particles
        kdtree: 3D KDTree of all particle positions
        return: density_center_position
        """
        # 6 nearest neighbors
        nb_r_list6, nb_index_list6 = kdtree.query(particle.pos,k=6)
        nb_mass_tot6=np.sum(particle[nb_index_list6].mass,axis=1) + particle.mass
        
        nb_inv_r6 = 1/nb_r_list6[:,5]
        rho = nb_mass_tot6*(nb_inv_r6*nb_inv_r6*nb_inv_r6)
        particle.addNewMember('density',rho)
        rho_tot = rho.sum()
     
        cm_pos = np.array([np.sum(rho*particle.pos[:,i])/rho_tot for i in range(3)])
        cm_vel = np.array([np.sum(rho*particle.vel[:,i])/rho_tot for i in range(3)])
        self.pos = np.append(self.pos, [cm_pos], axis=0)
        self.vel = np.append(self.vel, [cm_vel], axis=0)

        return cm_pos, cm_vel

    def calcCoreRadius(self, particle):
        """ 
        calculate core radius, using Casertano & Hut (1985) method:
        rc = sqrt(\sum_i rho_i^2 r_i^2 / \sum_i rho_i^2)
        """
        rho2 = particle.density*particle.density
        rc = np.sqrt((particle.r2*rho2).sum()/(rho2.sum()))
        self.rc = np.append(self.rc, rc)

        return rc

    def addTime(self, time):
        self.time = np.append(self.time, time)

def joinLagrangian(*_dat):
    """
    Join multiple data to one
    """
    type0 = type(_dat[0])
    for idat in _dat:
        if (type(idat) != type0):
            raise ValueError('Initial fail, date type not consistent, type [0] is ',type0,' given ',type(idat))
    new_dat = type0(_dat[0].n_frac)
    for key, item in _dat[0].__dict__.items():
        if (type(item) == np.ndarray):
            new_dat.__dict__[key] = np.concatenate(tuple(map(lambda x:x.__dict__[key], _dat)),axis=0)
        elif(issubclass(type(item), DictNpArrayMix)):
            new_dat.__dict__[key] = joinLagrangian(*tuple(map(lambda x:x.__dict__[key], _dat)))
        else:
            new_dat.__dict__[key] = _dat[0].__dict__[key]
    new_dat.size = np.sum(tuple(map(lambda x:x.size, _dat)))
    return new_dat

class LagrangianVelocity(DictNpArrayMix):
    """ Lagrangian velocity component
    """
    
    def __init__(self, _n_frac, _dat=None, _offset=int(0)):
        keys  = [['abs',_n_frac],['x',_n_frac],['y',_n_frac],['z',_n_frac],['rad',_n_frac],['tan',_n_frac],['rot',_n_frac]] # all, x, y, z, radial, tangential, rotational
        DictNpArrayMix.__init__(self, keys, _dat, _offset)
        self.n_frac = _n_frac

class Lagrangian(DictNpArrayMix):
    """ Lagrangian parameters
    """

    def __init__(self, _n_frac, _dat=None, _offset=int(0)):
        keys = [['r', _n_frac],['m', _n_frac],['n', _n_frac]] # radius, mass, number
        DictNpArrayMix.__init__(self, keys, _dat, _offset)
        self.vel  = LagrangianVelocity(_n_frac,_dat, self.ncols+_offset)
        self.ncols += self.vel.ncols
        self.sigma= LagrangianVelocity(_n_frac,_dat, self.ncols+_offset)
        self.ncols += self.sigma.ncols
        self.n_frac = _n_frac

    def calcOneSnapshot(self, _particle, _mass_fraction, _rc, _mode='sphere'):
        """ calculate one snapshot lagrangian parameters
        _particle: sorted particles
        _mass_fraction: fraction of mass 
        _rc: core radius
        _mode: sphere: calculate averaged properties from center to Lagrangian radii; shell: calculate properties between two neighbor radii
        """
        shell_mode = True if (_mode == 'shell') else False

        def find_mass_index(_mass_cum,_mass_fraction):
            mass_bins=np.append(0,_mass_fraction*_mass_cum[-1])
            index,count=np.histogram(_mass_cum,bins=mass_bins)
            index=index.cumsum()
            index[index>=_mass_cum.size]=_mass_cum.size-1
            return index

        if (self.n_frac != int(1+_mass_fraction.size)):
            raise ValueError('Mass fraction size ',mass_fraction.size,' is not consistent with Lagrangian array size ',self.n_frac)

        if (_particle.size<=1):
            size = self.size
            empty = Lagrangian(self.n_frac,np.zeros([1,self.n_frac*self.ncols]))
            self.append(empty)
            if (size+1!=self.size):
                raise ValueError('Size should increase one, but increase', self.size-size)
        else:
            self.size += 1
            self.vel.size += 1
            self.sigma.size += 1
            mcum=_particle.mass.cumsum()
            r = np.sqrt(_particle.r2)
            rindex= find_mass_index(mcum, _mass_fraction)
            rlagr = r[rindex]
            if(len(self.r.shape)!=2):
                raise ValueError('r shape is wrong',self.r.shape)
            self.r = np.append(self.r, [np.append(rlagr,_rc)], axis=0)
            nlagr = rindex+1
            if (shell_mode): nlagr[1:] -= nlagr[:-1]
            nc = (_particle.r2<(_rc*_rc)).sum()
            self.n = np.append(self.n, [np.append(nlagr,nc)], axis=0)
            mlagr = mcum[rindex]
            if (shell_mode): mlagr[1:] -= mlagr[:-1]
            sel = nlagr>0
            mlagr[sel] /= nlagr[sel]
            if (nc>0): mc = mcum[nc-1]/nc
            else:      mc = 0.0
            self.m = np.append(self.m, [np.append(mlagr,mc)], axis=0)

            m   =_particle.mass 
            vel = _particle.vel
            pos = _particle.pos
            rx = pos[:,0]
            ry = pos[:,1]
            rz = pos[:,2]
            vx = vel[:,0]
            vy = vel[:,1]
            vz = vel[:,2]

            # x-y plane radial velocity * rxy
            rvxy = rx*vx + ry*vy
            # radial velocity value
            vr   = (rvxy + rz*vz)/r
            # tangential velocity vector
            vt = [None]*3
            vt[0] = vx - vr*rx/r
            vt[1] = vy - vr*ry/r
            vt[2] = vz - vr*rz/r
            # x-y plane radial position square
            rxy2 = rx*rx + ry*ry
            # rotational velocity
            vrotx = vx - rvxy*rx/rxy2
            vroty = vy - rvxy*ry/rxy2
            vrot = np.sqrt(vrotx*vrotx+vroty*vroty)
            # rotational direction sign
            vrotd = vrotx*ry - vroty*rx
            ng_sig = (vrotd<0.0)
            vrot[ng_sig] = -vrot[ng_sig]
            
            n_offset = np.append(np.zeros(1),nlagr).astype(int)
            vlst = [vx, vy, vz, vr, vt[0], vt[1], vt[2], vrot]
            vave = [None]*len(vlst)
            for k in range(len(vlst)):
                vlagr = []
                if (shell_mode):
                    vlagr = [np.average(m[n_offset[i]:n_offset[i+1]]*vlst[k][n_offset[i]:n_offset[i+1]])/mlagr[i] if (n_offset[i]<n_offset[i+1]) else 0.0 for i in range(mlagr.size)]
                else:
                    vlagr = [np.average(m[:n_offset[i+1]]*vlst[k][:n_offset[i+1]])/mlagr[i] if (n_offset[i+1]>0) else 0.0 for i in range(mlagr.size)]
                # core radius
                vlagr.append(np.average(m[0:nc]*vlst[k][0:nc])/mc if (nc>0) else 0.0)
                vave[k] = np.array(vlagr)
            
            self.vel.x = np.append(self.vel.x, [vave[0]], axis=0)
            self.vel.y = np.append(self.vel.y, [vave[1]], axis=0)
            self.vel.z = np.append(self.vel.z, [vave[2]], axis=0)
            self.vel.abs= np.append(self.vel.abs, [np.sqrt(vave[0]*vave[0]+vave[1]*vave[1]+vave[2]*vave[2])], axis=0)
            self.vel.rad = np.append(self.vel.rad, [vave[3]], axis=0)
            self.vel.tan = np.append(self.vel.tan, [np.sqrt(vave[4]*vave[4]+vave[5]*vave[5]+vave[6]*vave[6])], axis=0)
            self.vel.rot = np.append(self.vel.rot, [vave[7]], axis=0)
            
            sigma = [None]*len(vlst)
            for k in range(len(vlst)):
                slagr = []
                if (shell_mode):
                    slagr = [np.average(m[n_offset[i]:n_offset[i+1]] * (vlst[k][n_offset[i]:n_offset[i+1]] - vave[k][i])**2) / mlagr[i] if (n_offset[i]<n_offset[i+1]) else 0.0 for i in range(_mass_fraction.size)]
                else:
                    slagr = [np.average(m[:n_offset[i+1]] * (vlst[k][:n_offset[i+1]] - vave[k][i])**2) / mlagr[i] if (n_offset[i+1]>0) else 0.0 for i in range(_mass_fraction.size)]
                # core radius
                slagr.append(np.average(m[0:nc] * (vlst[k][0:nc] - vave[k][-1])**2) / mc if (nc>0) else 0.0)
                sigma[k] = np.array(slagr)

            self.sigma.x = np.append(self.sigma.x, [np.sqrt(sigma[0])], axis=0)
            self.sigma.y = np.append(self.sigma.y, [np.sqrt(sigma[1])], axis=0)
            self.sigma.z = np.append(self.sigma.z, [np.sqrt(sigma[2])], axis=0)
            self.sigma.abs= np.append(self.sigma.abs, [np.sqrt(sigma[0]+sigma[1]+sigma[2])], axis=0)
            self.sigma.rad = np.append(self.sigma.rad, [np.sqrt(sigma[3])], axis=0)
            self.sigma.tan = np.append(self.sigma.tan, [np.sqrt(sigma[4]+sigma[5]+sigma[6])], axis=0)
            self.sigma.rot = np.append(self.sigma.rot, [np.sqrt(sigma[7])], axis=0)

class LagrangianMultiple(DictNpArrayMix):
    """ Lagrangian for single, binaries and all
    """
    def __init__ (self, _n_frac, _dat=None, _offset=int(0)):
        DictNpArrayMix.__init__(self, [['time',1]], _dat, _offset)
        ncols = self.ncols
        self.single = Lagrangian(_n_frac, _dat, ncols+_offset)
        ncols += self.single.ncols
        self.binary = Lagrangian(_n_frac, _dat, ncols+_offset)
        ncols += self.binary.ncols
        self.all    = Lagrangian(_n_frac, _dat, ncols+_offset)
        ncols += self.all.ncols
        self.ncols = ncols
        self.size = self.all.size
        self.n_frac= _n_frac

    def calcOneSnapshot(self, time, single, binary, mass_fraction, rc, mode):
        """ Calculate Lagrangian radii and related properties
        single: single partilces (cm corrected and r2 exist)
        binary: binaries (cm corrected and r2 exist)
        mass_fraction: Lagragian radius mass fraction
        rc: core radius
        """    
        self.time = np.append(self.time, time)
        single_sim = SimpleParticle(single)
        single_sim.calc_r2()
        binary_sim = SimpleParticle(binary)
        binary_sim.calc_r2()
        all_sim = join(single_sim, binary_sim)
        n_single = single.size
        n_binary = binary.size

        idx = all_sim.r2.argsort()
        idx_single = idx[idx<n_single]
        idx_binary = idx[idx>=n_single] - n_single
        single_sort = single_sim[idx_single]
        binary_sort = binary_sim[idx_binary]
        all_sort = all_sim[idx]
    
        self.single.calcOneSnapshot(single_sort, mass_fraction, rc, mode)
        self.binary.calcOneSnapshot(binary_sort, mass_fraction, rc, mode)
        self.all.calcOneSnapshot(all_sort, mass_fraction, rc)
        if (self.binary.size != self.single.size):
            raise ValueError('Size inconsistence: single.size:', self.single.size, ' binary.size:', self.binary.size)

        self.size += 1
        
