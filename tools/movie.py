#!/usr/bin/env python3
import sys
import collections
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
import multiprocessing as mp
import getopt
import petar
import imageio
import os
#from pygifsicle import optimize

plt.style.use('dark_background')

def create_movie(filenames, fps, output_file):
    with imageio.get_writer(output_file, fps=fps) as writer:
        for filename in filenames:
            writer.append_data(imageio.imread(filename))
    writer.close()


class PlotXY:
    def __init__(self):
        self.plot_mode = 'x-y'
        self.boxsize = 2
        self.cm_mode = 'core'
        self.lum_min = 1e-5
        self.lum_max = 1e6
        self.x_min = -2
        self.x_max = 2
        self.y_min = -2
        self.y_max = 2
        self.cm_boxsize = 2
        self.framescale = 1
        self.nlayer_cross = 5
        self.nlayer_point = 10
        self.alpha_amplifier = 2.5
        self.marker_scale = 1.0
        self.mass_power = 1.0
        self.size_mode = 'mass'
        self.ptcls=[]

    def init(self, axe, **kwargs):
        for key in self.__dict__.keys():
            if (key in kwargs.keys()): self.__dict__[key] = kwargs[key]
        if ('boxsize' in kwargs.keys()):
            for key in ['x_min','y_min']:
                if (not key in kwargs.keys()): self.__dict__[key] = -self.boxsize
            for key in ['x_max','y_max']:
                if (not key in kwargs.keys()): self.__dict__[key] = self.boxsize

        #if ('boxsize' in kwargs.keys()): boxsize = kwargs['boxsize']
        #if ('nlayer_cross' in kwargs.keys()): nlayer_cross = kwargs['nlayer_cross']
        #if ('alpha_amplifier' in kwargs.keys()): alpha_amplifier = kwargs['alpha_amplifier']
        if (not 'cm_boxsize' in kwargs.keys()):
            self.cm_boxsize = self.boxsize

        nlayer = self.nlayer_cross + self.nlayer_point
        self.nlayer = nlayer
        alphascale = np.linspace(1,nlayer,nlayer)*2.0/(nlayer*(nlayer+1))*self.alpha_amplifier
        #print('Alpha layer sequence:',alphascale,' sum:',alphascale.sum())
        self.sizescale = np.logspace(0,3,nlayer)[::-1]
        #print('Size layer sequence:',self.sizescale)

        axe.set_xlim(self.x_min, self.x_max)
        axe.set_ylim(self.y_min, self.y_max)
        axe.set_aspect(1.0)

        unit_length=''
        unit_vel=''
        if ('unit_length' in kwargs.keys()):
            unit_length = '['+kwargs['unit_length']+']'
        if ('unit_vel' in kwargs.keys()):
            unit_vel = '['+kwargs['unit_vel']+']'

        plot_mode = self.plot_mode
        axes_name = plot_mode.split('-')
        if (len(axes_name)!=2):
            raise ValueError('Plot mode %s is not supported, the format should be [x-axis name]-[y-axis name], check petar.movie -h for the options of -m' % plot_mode)
        labels=[None,None]
        for i in range(2):
            name = axes_name[i]
            if name in ['x','y','z','r']:
                labels[i] = r'%s %s' % (name, unit_length)
            elif name == 'rxy':
                labels[i] = r'$r_{xy}$ %s' % unit_length
            elif name in ['vx','vy','vz','vxy','vr','v']:
                labels[i] = r'$v_{%s}$ %s' % (name[1:], unit_vel)
            elif name == 'ra': 
                labels[i] = r'RA [$\degree$]'
            elif name == 'dec': 
                labels[i] = r'Dec [$\degree$]'
            elif name == 'pmracosdec':
                labels[i] = r'$\mu_{RA,cosDec}$ [mas/yr]'
            elif name == 'pmdec':
                labels[i] = r'$\mu_{Dec}$ [mas/yr]'
            elif name in ['lon','lat']:
                labels[i] = r'%s [$\degree$]' % name
            elif name in ['pmlon','pmlat']:
                labels[i] = r'$\mu_{%s}$ [mas/yr]' % name[2:]
            else:
                raise ValueError('Plot mode axis name %s is not supported, check petar.movie -h for the options of -m' % name)
        axe.set_xlabel(labels[0])
        axe.set_ylabel(labels[1])
        self.xy_labels=labels

        for i in range(self.nlayer_cross):
            pt =axe.scatter([],[],marker='+',alpha=alphascale[i])
            self.ptcls.append(pt)
        for i in range(self.nlayer_point):
            pt =axe.scatter([],[],alpha=alphascale[i],linewidth=0)
            self.ptcls.append(pt)
        return self.ptcls

    def correctCM(self, x, y):
        xcm = 0
        ycm = 0
        boxsize = self.cm_boxsize
        if (self.cm_mode=='density'):
            nbins=1000
            xmid = np.average(np.abs(x))
            ymid = np.average(np.abs(y))
            #print(xmid,ymid)
            xmin = -5*xmid
            xmax =  5*xmid
            ymin = -5*ymid
            ymax =  5*ymid
            xbins=np.linspace(xmin,xmax,nbins)
            ybins=np.linspace(ymin,ymax,nbins)
            lbin=[xbins,ybins]
            counts, _, _ = np.histogram2d(x, y, bins=lbin)
            xp,yp=np.where(counts>0.2*counts.max())
            m = counts[xp,yp]
            mtot = m.sum()
            xcm=((xp*m).sum()/mtot+0.5)*(xmax-xmin)/nbins+xmin
            ycm=((yp*m).sum()/mtot+0.5)*(ymax-ymin)/nbins+ymin
            x = x - xcm
            y = y - ycm
     
            nbins=500
            xmid = np.average(np.abs(x))
            ymid = np.average(np.abs(y))
            xmin = -boxsize
            xmax =  boxsize
            ymin = -boxsize
            ymax =  boxsize
            xbins=np.linspace(xmin,xmax,nbins)
            ybins=np.linspace(ymin,ymax,nbins)
            lbin=[xbins,ybins]
            counts, _, _ = np.histogram2d(x,y, bins=lbin)
            xp,yp=np.where(counts>0.2*counts.max())
            m = counts[xp,yp]
            mtot = m.sum()
            xcm2=((xp*m).sum()/mtot+0.5)*(xmax-xmin)/nbins+xmin
            ycm2=((yp*m).sum()/mtot+0.5)*(ymax-ymin)/nbins+ymin
            x = x - xcm2
            y = y - ycm2
            xcm += xcm2
            ycm += ycm2
        elif (self.cm_mode=='average'):
            xcm = x.sum()/x.size
            ycm = y.sum()/y.size
            x = x - xcm
            y = y - ycm
            sel=(x>-boxsize) & (x<boxsize) & (y>-boxsize) & (y<boxsize)
            nsel=sel.sum()
            xcm2 = (x[sel]).sum()/float(nsel)
            ycm2 = (y[sel]).sum()/float(nsel)
            x = x - xcm2 
            y = y - ycm2
            xcm += xcm2
            ycm += ycm2
        return xcm, ycm

    def plot(self, data, xcm_text, ycm_text):

        plot_mode = self.plot_mode
        core_correct = (self.cm_mode=='core') & (data.generate_binary != 2) 
        cm_is_core = (data.generate_binary == 2)
        read_core = cm_is_core | core_correct
        origin_mode = (self.cm_mode=='none')
        axes_name = plot_mode.split('-')
        if (len(axes_name)!=2):
            raise ValueError('Plot mode %s is not supported, the format should be [x-axis name]-[y-axis name], check petar.movie -h for the options of -m' % plot_mode)
        labels=[None,None]
        cm_text=[None,None]
        xy = [[],[]]
        xycm = [0,0]
        xyz_index={'x':0,'y':1,'z':2}
        for i in range(2):
            name = axes_name[i]
            if name in ['x','y','z','vx','vy','vz']:
                v_flag = (name[0]=='v')
                if v_flag:
                    x = data.data.vel[:,xyz_index[name[1]]]
                    xycm[i] = data.header.vel_offset[xyz_index[name[1]]]
                else:
                    x = data.data.pos[:,xyz_index[name]]
                    xycm[i] = data.header.pos_offset[xyz_index[name]]
                xyc = 0
                if (read_core):
                    if v_flag:
                        xyc = data.core.vel[0,xyz_index[name[1]]]
                    else:
                        xyc = data.core.pos[0,xyz_index[name]]
                if (cm_is_core):
                    xycm[i] = xyc
                if (core_correct):
                    xy[i] = x + xycm[i] - xyc
                    xycm[i] = xyc
                elif (origin_mode):
                    xy[i] = x + xycm[i]
                else:
                    xy[i] = x
                if v_flag: 
                    cm_text[i] = r'$v_{%s,cm}=$' % name[1:]
                else:
                    cm_text[i] = r'$%s_{cm}=$' % name
            elif name in ['rxy','vxy']:
                v_flag = (name[0]=='v')
                if v_flag:
                    rx = data.data.vel[:,0]
                    ry = data.data.vel[:,1]
                    rxcm = data.header.vel_offset[0]
                    rycm = data.header.vel_offset[1]
                else:
                    rx = data.data.pos[:,0]
                    ry = data.data.pos[:,1]
                    rxcm = data.header.pos_offset[0]
                    rycm = data.header.pos_offset[1]
                rxycm = np.sqrt(rxcm*rxcm + rycm*rycm)
                rxc = 0
                ryc = 0
                if (read_core):
                    if v_flag:
                        rxc = data.core.vel[0,0]
                        ryc = data.core.vel[0,1]
                    else:
                        rxc = data.core.pos[0,0]
                        ryc = data.core.pos[0,1]
                if (cm_is_core):
                    rxcm = rxc
                    rycm = ryc
                    xycm[i] = np.sqrt(rxc*rxc + ryc*ryc)
                if (core_correct):
                    xy[i] = np.sqrt((rx + rxcm - rxc)**2 + (ry + rycm - ryc)**2)
                    xycm[i] = np.sqrt(rxc*rxc + ryc*ryc)
                elif (origin_mode):
                    xy[i] = np.sqrt((rx + rxcm)**2 + (ry + rycm)**2)
                else:
                    xy[i] = np.sqrt(rx*rx + ry*ry)
                if v_flag:
                    cm_text[i] = r'$v_{xy,cm}=$'
                else:
                    cm_text[i] = r'$r_{xy,cm}=$'
            elif name in ['r','v']:
                v_flag = (name=='v')
                if v_flag: 
                    rx = data.data.vel[:,0]
                    ry = data.data.vel[:,1]
                    rz = data.data.vel[:,2]
                    rxcm = data.header.vel_offset[0]
                    rycm = data.header.vel_offset[1]
                    rzcm = data.header.vel_offset[2]
                else:
                    rx = data.data.pos[:,0]
                    ry = data.data.pos[:,1]
                    rz = data.data.pos[:,2]
                    rxcm = data.header.pos_offset[0]
                    rycm = data.header.pos_offset[1]
                    rzcm = data.header.pos_offset[2]
                xycm[i] = np.sqrt(rxcm*rxcm + rycm*rycm + rzcm*rzcm)
                rxc = 0
                ryc = 0
                rzc = 0
                if (read_core):
                    if v_flag:
                        rxc = data.core.vel[0,0]
                        ryc = data.core.vel[0,1]
                        rzc = data.core.vel[0,2]
                    else:
                        rxc = data.core.pos[0,0]
                        ryc = data.core.pos[0,1]
                        rzc = data.core.pos[0,2]
                if (cm_is_core):
                    rxcm = rxc
                    rycm = ryc
                    rzcm = rzc 
                    xycm[i] = np.sqrt(rxc*rxc + ryc*ryc + rzc*rzc)
                if (core_correct):
                    xy[i] = np.sqrt((rx + rxcm - rxc)**2 + (ry + rycm - ryc)**2 + (rz + rzcm - rzc)**2)
                    xycm[i] = np.sqrt(rxc*rxc + ryc*ryc + rzc*rzc)
                elif (origin_mode):
                    xy[i] = np.sqrt((rx + rxcm)**2 + (ry + rycm)**2 + (rz + rzcm)**2)
                else:
                    xy[i] = np.sqrt(rx*rx + ry*ry + rz*rz)
                if v_flag:
                    cm_text[i] = r'$v_{cm}=$' 
                else:
                    cm_text[i] = r'$r_{cm}=$' 
            elif name == 'vr':
                vx = data.data.vel[:,0]
                vy = data.data.vel[:,1]
                vz = data.data.vel[:,2]
                vxcm = data.header.vel_offset[0]
                vycm = data.header.vel_offset[1]
                vzcm = data.header.vel_offset[2]
                rx = data.data.pos[:,0]
                ry = data.data.pos[:,1]
                rz = data.data.pos[:,2]
                rxcm = data.header.pos_offset[0]
                rycm = data.header.pos_offset[1]
                rzcm = data.header.pos_offset[2]
                rxc = 0
                ryc = 0
                rzc = 0
                vxc = 0
                vyc = 0
                vzc = 0
                rcm = np.sqrt(rxcm*rxcm + rycm*rycm + rzcm*rzcm)                
                drdvcm = rxcm*vxcm + rycm*vycm + rzcm*vzcm
                xycm[i] = drdvcm/rcm
                if (read_core):
                    vxc = data.core.vel[0,0]
                    vyc = data.core.vel[0,1]
                    vzc = data.core.vel[0,2]
                    rxc = data.core.pos[0,0]
                    ryc = data.core.pos[0,1]
                    rzc = data.core.pos[0,2]
                if (cm_is_core):
                    rxcm = rxc
                    rycm = ryc
                    rzcm = rzc 
                    vxcm = vxc
                    vycm = vyc
                    vzcm = vzc 
                    rcm = np.sqrt(rxcm*rxcm + rycm*rycm + rzcm*rzcm)                
                    drdvcm = rxcm*vxcm + rycm*vycm + rzcm*vzcm
                    xycm[i] = drdvcm/rcm
                if (core_correct):
                    rxf = rx + rxcm - rxc
                    ryf = ry + rycm - ryc
                    rzf = rz + rzcm - rzc
                    vxf = vx + vxcm - vxc
                    vyf = vy + vycm - vyc
                    vzf = vz + vzcm - vzc
                    rcm = np.sqrt(rxc*rxc + ryc*ryc + rzc*rzc)                
                    drdvcm = rxc*vxc + ryc*vyc + rzc*vzc
                    xycm[i] = drdvcm/rcm
                    r = np.sqrt(rxf*rxf + ryf*ryf + rzf*rzf)
                    drdv = rxf*vxf + ryf*vyf + rzf*vzf
                    xy[i] = drdv/r
                elif (origin_mode):
                    rxf = rx + rxcm
                    ryf = ry + rycm
                    rzf = rz + rzcm
                    vxf = vx + vxcm
                    vyf = vy + vycm
                    vzf = vz + vzcm
                    r = np.sqrt(rxf*rxf + ryf*ryf + rzf*rzf)
                    drdv = rxf*vxf + ryf*vyf + rzf*vzf
                    xy[i] = drdv/r
                else:
                    r = np.sqrt(rx*rx + ry*ry + rz*rz)
                    drdv = rx*vx + ry*vy + rz*vz
                    xy[i] = drdv/r
                cm_text[i] = r'$v_{r}=$'
            elif name == 'ra':
                xycm[i] = data.skycm.icrs.ra.value
                xy[i] = data.sky.icrs.ra.value - xycm[i]
                cm_text[i] = r'$RA_{cm}=$'
            elif name == 'dec':
                xycm[i] = data.skycm.icrs.dec.value
                xy[i] = data.sky.icrs.dec.value - xycm[i]
                cm_text[i] = r'$Dec_{cm}=$'
            elif name == 'pmracosdec':
                xycm[i] = data.skycm.icrs.pm_ra_cosdec.value
                xy[i] = data.sky.icrs.pm_ra_cosdec.value - xycm[i]
                cm_text[i] = r'$\mu_{RA,cosdec,cm}=$'
            elif name == 'pmdec':
                xycm[i] = data.skycm.icrs.pm_dec.value
                xy[i] = data.sky.icrs.pm_dec.value - xycm[i]
                cm_text[i] = r'$\mu_{Dec,cm}=$'
            elif name == 'lon':
                data.skycm.representation_type = 'spherical'
                xycm[i] = data.skycm.lon.value
                data.sky.representation_type = 'spherical'
                xy[i] = data.sky.lon.value - xycm[i]
                cm_text[i] = r'$lon_{cm}=$'
            elif name == 'lat':
                data.skycm.representation_type = 'spherical'
                xycm[i] = data.skycm.lat.value
                data.sky.representation_type = 'spherical'
                xy[i] = data.sky.lat.value - xycm[i]
                cm_text[i] = r'$lat_{cm}=$'
            elif name == 'pmlon':
                data.skycm.representation_type = 'spherical'
                xycm[i] = data.skycm.pm_lon.value
                data.sky.representation_type = 'spherical'
                xy[i] = data.sky.pm_lon.value - xycm[i]
                cm_text[i] = r'$\mu_{lon,cm}=$'
            elif name == 'pmlat':
                data.skycm.representation_type = 'spherical'
                xycm[i] = data.skycm.pm_lat.value
                data.sky.representation_type = 'spherical'
                xy[i] = data.sky.pm_lat.value - xycm[i]
                cm_text[i] = r'$\mu_{lat,cm}=$'
            else:
                raise ValueError('Plot mode axis name %s is not supported, check petar.movie -h for the options of -m' % name)
        
        # not for cm_mode=core
        dxcm, dycm = self.correctCM(xy[0], xy[1])
        xycm[0] += dxcm
        xycm[1] += dycm
        xcm_text.set_text(cm_text[0]+('%f' % xycm[0]))
        ycm_text.set_text(cm_text[1]+('%f' % xycm[1]))

        mass = data.data.mass
        lum = data.lum
        colors=data.getColor()
        for i in range(self.nlayer):
            if (self.size_mode == 'loglum'):
                sizes = (np.log10(lum)-np.log10(self.lum_min))
            else:
                sizes = (mass**self.mass_power)*self.sizescale[i]*self.framescale*self.marker_scale
            #sizes = luminosity
            #print(self.plot_mode,xy[0].shape,xy[1].shape)
            self.ptcls[i].set_offsets(np.array([xy[0],xy[1]]).transpose())
            self.ptcls[i].set_sizes(sizes)
            self.ptcls[i].set_color(colors)
        return self.ptcls

class PlotHR:
    def __init__(self):
        self.lum_min = 1e-5
        self.lum_max = 1e6
        self.temp_min = 1000
        self.temp_max = 50000
        self.ptcls = []

    def init(self, axe, **kwargs):
        for key in self.__dict__.keys():
            if (key in kwargs.keys()): self.__dict__[key] = kwargs[key]
        self.ptcls.append(axe.scatter([],[],marker='.',linewidth=0))
        axe.set_xlim(self.temp_max,self.temp_min);
        axe.set_ylim(self.lum_min,self.lum_max);
        axe.set_yscale('log');
        axe.set_xscale('log');
        axe.set_xlabel(r'$T_{eff}[K]$')
        axe.set_ylabel(r'$L[L_\odot]$')

        return self.ptcls


    def plot(self, lum, temp, types):
        colors = cm.rainbow(types/15.0)
        self.ptcls[0].set_offsets(np.array([temp, lum]).transpose())
        self.ptcls[0].set_color(colors)
        return self.ptcls

class PlotSemiEcc:
    def __init__(self):
        self.semi_min = 1e-7
        self.semi_max = 0.1
        self.ecc_min = 0.0
        self.ecc_max = 1.0
        self.mass_power = 1.0
        self.marker_scale = 1.0
        self.cm_mode = 'core'
        self.bin_color = 'white'
        self.bin_rmax = 2.0
        self.bin_rmin = 1e-5
        self.ptcls = []

    def init(self, axe, **kwargs):
        for key in self.__dict__.keys():
            if (key in kwargs.keys()): self.__dict__[key] = kwargs[key]
        self.ptcls.append(axe.scatter([],[],marker='.',linewidth=0))
        axe.set_xscale('log')
        axe.set_xlim(self.semi_min,self.semi_max)
        axe.set_ylim(self.ecc_min,self.ecc_max)
        if ('unit_length' in kwargs.keys()):
            unit_label = '['+kwargs['unit_length']+']'
            axe.set_xlabel('semi'+unit_label)
        else:
            axe.set_xlabel('semi')
        axe.set_ylabel('ecc')
        return self.ptcls

    def plot(self, data):
        #colors = cm.rainbow(types/13.0)
        sizes = data.binary.mass**self.mass_power*self.marker_scale
        core_correct = (self.cm_mode=='core') & (data.generate_binary != 2)
        origin_mode = (self.cm_mode=='none')
        xcm = data.header.pos_offset[0]
        ycm = data.header.pos_offset[1]
        zcm = data.header.pos_offset[2]
        x = data.binary.pos[:,0]
        y = data.binary.pos[:,1]
        z = data.binary.pos[:,2]
        if (core_correct):
            xc = data.core.pos[0,0]
            yc = data.core.pos[0,1]
            zc = data.core.pos[0,2]
            x += xcm - xc
            y += ycm - yc
            z += zcm - xc
            xcm = xc
            ycm = yc
            zcm = zc
        elif (origin_mode):
            x += xcm
            y += ycm
            z += zcm
        r = np.sqrt(x*x+y*y+z*z)
        self.ptcls[0].set_offsets(np.array([data.binary.semi, data.binary.ecc]).transpose())
        if (self.bin_color == 'white'):
            self.ptcls[0].set_color('white')
        elif (self.bin_color == 'distance'):
            colors = cm.hot_r(np.log(r/self.bin_rmin)/np.log(self.bin_rmax/self.bin_rmin))
            self.ptcls[0].set_color(colors)
        else:
            raise ValueError("Color mode %s is not supported for semi-ecc plot " % (self.bin_color))
        self.ptcls[0].set_sizes(sizes)
        return self.ptcls

class PlotLagr:
    def __init__(self):
        self.rlagr_min = 0
        self.rlagr_max = 1
        self.time_min = 0
        self.time_max = 1
        self.rlagr_scale= 'linear'
        self.ptcls = []

    def init(self, axe, lagr, **kwargs):
        for key in self.__dict__.keys():
            if (key in kwargs.keys()): self.__dict__[key] = kwargs[key]
        mfrac=list(lagr.initargs['mass_fraction']) + ['Rc']
        for mi in range(len(mfrac)):
            pt, = axe.plot(lagr.time, lagr.all.r[:,mi],'-', color='grey')
            self.ptcls.append(pt)
        for mi in range(len(mfrac)):
            mlabel=mfrac[mi]
            pt, = axe.plot([],[], '*', markersize=8, markeredgewidth=0, markeredgecolor='none', label=mlabel)
            self.ptcls.append(pt)
        axe.set_yscale(self.rlagr_scale)
        axe.set_xlim(self.time_min, self.time_max)
        axe.set_ylim(self.rlagr_min, self.rlagr_max)
        if ('unit_length' in kwargs.keys()):
            unit_label = '['+kwargs['unit_length']+']'
            axe.set_ylabel(r'$R_{lagr}$'+unit_label)
        else:
            axe.set_ylabel(r'$R_{lagr}$')
        if ('unit_time' in kwargs.keys()):
            unit_label = ' '+kwargs['unit_time']
            axe.set_xlabel('Time'+unit_label)
        else:
            axe.set_xlabel('Time')
        axe.legend(loc='upper right')


    def plot(self, lagr, tnow):
        nfrac=lagr.initargs['mass_fraction'].size + 1
        sel = (lagr.time==tnow)
        for mi in range(nfrac):
            self.ptcls[mi+nfrac].set_data(lagr.time[sel], lagr.all.r[sel,mi])
        return self.ptcls

class Data:
    def __init__(self, **kwargs):
        self.skiprows = 0
        self.generate_binary=2
        self.interrupt_mode = 'bse'
        self.external_mode = 'none'
        self.G = 0.00449830997959438 # pc^3/(Msun*Myr^2)
        self.semi_max = 0.1
        self.cm_mode = 'core'
        self.snapshot_format = 'ascii'
        self.lum_min = 1e-5
        self.lum_max = 1e6
        self.temp_min = 1000
        self.temp_max = 50000
        self.ekin_min = 0
        self.ekin_max = 100
        self.pot_min = -200
        self.pot_max = 0
        self.etot_min = -100
        self.etot_max = 100
        self.get_skycoord=False
        self.color_mode = 'white'
        self.galev_filter = ['Johnson', 'SDSS', 'HST', 'CSST', 'Gaia']
        self.galev_mode = 'abs_mag'
        self.galev_color = ['HST.F555W','HST.F814W']
        self.galev_mag_range = np.array([0,2])

        for key in self.__dict__.keys():
            if (key in kwargs.keys()): self.__dict__[key] = kwargs[key]

    def __getitem__(self, k):
        return self.__dict__[k]

    def keys(self):
        return self.__dict__.keys()

    def read(self, file_path, core):
        data=self.__dict__
        skiprows = self.skiprows

        header=petar.PeTarDataHeader(file_path, snapshot_format=self.snapshot_format, external_mode=self.external_mode)
        header_offset = petar.HEADER_OFFSET
        if (self.external_mode!='none'): header_offset=petar.HEADER_OFFSET_WITH_CM

        data['header'] = header
        data['t'] = header.time
        if (self.cm_mode=='core') | (self.generate_binary==2):
            tsel=(core.time==data['t'])
            data['core']=core[tsel]

        if (self.generate_binary == 2):
            single = petar.Particle(interrupt_mode=self.interrupt_mode, external_mode=self.external_mode)
            p1 = petar.Particle(interrupt_mode=self.interrupt_mode, external_mode=self.external_mode)
            p2 = petar.Particle(interrupt_mode=self.interrupt_mode, external_mode=self.external_mode)
            binary = petar.Binary(p1,p2, G=self.G)
            if os.path.getsize(file_path+'.single')>0:
                if (self.snapshot_format=='ascii'):
                    single.loadtxt(file_path+'.single')
                elif (self.snapshot_format=='binary'):
                    single.fromfile(file_path+'.single')
                elif (self.snapshot_format=='npy'):
                    single.load(file_path+'.single.npy')
                else:
                    raise ValueError('Snapshot format %s unknown, should be ascii, binary or npy.' % self.snapshot_format)
                if (self.color_mode == 'galev'):
                    mag = petar.GalevMag(filter = self.galev_filter, mode = self.galev_mode)
                    fpath = file_path+'.single.galev.mag'
                    if os.path.exists(fpath):
                        if (self.snapshot_format=='ascii'):
                            mag.loadtxt(fpath)
                        elif (self.snapshot_format=='binary'):
                            mag.fromfile(fpath)
                        elif (self.snapshot_format=='npy'):
                            mag.load(fpath)
                        else:
                            raise ValueError('Snapshot format %s unknown, should be ascii, binary or npy.' % self.snapshot_format)      
                        single.addNewMember('galev',mag)
                    else:
                        raise ValueError('Galev data %s is not found.' % fpath)
            single.calcEkin()
            single.calcEtot()
            if os.path.getsize(file_path+'.binary')>0:
                if (self.snapshot_format=='ascii'):
                    binary.loadtxt(file_path+'.binary')
                elif (self.snapshot_format=='binary'):
                    binary.fromfile(file_path+'.binary')
                elif (self.snapshot_format=='npy'):
                    binary.load(file_path+'.binary.npy')
                else:
                    raise ValueError('Snapshot format %s unknown, should be ascii, binary or npy.' % self.snapshot_format)
                if (self.color_mode == 'galev'):
                    mag = petar.GalevMag(filter = self.galev_filter, mode = self.galev_mode)
                    fpath = file_path+'.binary.galev.mag'
                    if os.path.exists(fpath):
                        if (self.snapshot_format=='ascii'):
                            mag.loadtxt(fpath)
                        elif (self.snapshot_format=='binary'):
                            mag.fromfile(fpath)
                        elif (self.snapshot_format=='npy'):
                            mag.load(fpath)
                        else:
                            raise ValueError('Snapshot format %s unknown, should be ascii, binary or npy.' % self.snapshot_format)
                        mag1 = mag[:binary.size]
                        mag2 = mag[binary.size:]
                        binary.p1.addNewMember('galev',mag1)
                        binary.p2.addNewMember('galev',mag2)
                        # get binary flux
                        mag.convertToFlux()
                        bmag = mag[:binary.size] + mag[binary.size:]
                        bmag.convertToMag()
                        binary.addNewMember('galev',bmag)
                    else:
                        raise ValueError('Galev data %s is not found.' % fpath)
            binary.calcEkin(True)
            binary.calcPot()
            binary.calcEtot(True)

            single_sim = petar.SimpleParticle(single)
            p1_sim = petar.SimpleParticle(binary.p1)
            p2_sim = petar.SimpleParticle(binary.p2)
            all_sim = petar.join(single_sim, p1_sim, p2_sim)
            data['data'] = all_sim
            data['ekin'] = np.concatenate((single.ekin, binary.p1.ekin, binary.p2.ekin))
            data['pot'] = np.concatenate((single.pot, binary.p1.pot, binary.p2.pot))
            data['etot'] = np.concatenate((single.etot, binary.p1.etot, binary.p2.etot))

            if (self.get_skycoord):
                data['sky'] = all_sim.toSkyCoord(pos_offset=data['core'].pos[0], vel_offset=data['core'].vel[0])
                data['skycm'] = data['core'].toSkyCoord()[0]
                #print(data['core'].pos,data['core'].vel,data['skycm'].icrs)

            data['binary'] = binary
            if ('bse' in self.interrupt_mode):
                data['lum'] = np.concatenate((single.star.lum, binary.p1.star.lum, binary.p2.star.lum))
                data['rad'] = np.concatenate((single.star.rad, binary.p1.star.rad, binary.p2.star.rad))
                data['type']= np.concatenate((single.star.type,binary.p1.star.type,binary.p2.star.type))
                temp_single = 5778*(single.star.lum/(single.star.rad*single.star.rad))**0.25
                temp_b1 = 5778*(binary.p1.star.lum/(binary.p1.star.rad*binary.p1.star.rad))**0.25
                temp_b2 = 5778*(binary.p2.star.lum/(binary.p2.star.rad*binary.p2.star.rad))**0.25
                temp_binary = (temp_b1*binary.p1.star.lum+temp_b2*binary.p2.star.lum)/(binary.p1.star.lum+binary.p2.star.lum)
                data['temp']= np.concatenate((temp_single, temp_b1, temp_b2))
                data['lum_cm'] = np.append(single.star.lum, binary.p1.star.lum+binary.p2.star.lum)
                data['temp_cm']= np.append(temp_single,temp_binary)
                data['type_cm']= np.append(single.star.type,np.max([binary.p1.star.type,binary.p2.star.type],axis=0))
            if (self.color_mode == 'galev'):
                mag_list = []
                mag_cm_list = []
                for key in self.galev_color:
                    mag_list.append(np.concatenate((single.galev[key], binary.p1.galev[key], binary.p2.galev[key])))
                    mag_cm_list.append(np.concatenate((single.galev[key], binary.galev[key])))
                data['galev_mag'] = np.array(mag_list)
                data['galev_mag_cm'] = np.array(mag_cm_list)
        else:
            particles=petar.Particle(interrupt_mode=self.interrupt_mode, external_mode=self.external_mode)
            if (self.snapshot_format=='ascii'): 
                particles.loadtxt(file_path, skiprows=1)
            else: 
                particles.fromfile(file_path, offset=header_offset)
            if (self.color_mode == 'galev'):
                mag = petar.GalevMag(filter = self.galev_filter, mode = self.galev_mode)
                fpath = file_path+'.galev.mag'
                if os.path.exists(fpath):
                    if (self.snapshot_format=='ascii'):
                        mag.loadtxt(fpath)
                    elif (self.snapshot_format=='binary'):
                        mag.fromfile(fpath)
                    elif (self.snapshot_format=='npy'):
                        mag.load(fpath)
                    else:
                        raise ValueError('Snapshot format %s unknown, should be ascii, binary or npy.' % self.snapshot_format)
                particles.addNewMember('galev',mag)
            data['data'] = particles
            if (self.get_skycoord):
                data['sky'] = particles.toSkyCoord(pos_offset=header.pos_offset, vel_offset=header.vel_offset)
                if (self.cm_mode == 'core'):
                    data['skycm'] = data['core'].toSkyCoord()[0]
                else:
                    data['skycm'] = header.toSkyCoord()
            
            if ('bse' in self.interrupt_mode):
                data['lum'] = particles.star.lum
                data['rad'] = particles.star.rad
                data['type']= particles.star.type
                data['temp']= 5778*(data['lum']/(data['rad']*data['rad']))**0.25

            if (self.color_mode == 'galev'):
                mag_list = []
                for key in self.galev_color:
                    mag_list.append(np.concatenate((single.galev[key], binary.p1.galev[key], binary.p2.galev[key])))
                data['galev_mag'] = np.array(mag_list)

            if (self.generate_binary>0):
                kdtree,single,binary = petar.findPair(particles, self.G, self.semi_max*2.0, True)
                data['binary'] = binary
                if ('bse' in self.interrupt_mode):
                    data['lum_cm'] = np.append(single.star.lum, binary.p1.star.lum+binary.p2.star.lum)
                    temp_single = 5778*(single.star.lum/(single.star.rad*single.star.rad))**0.25
                    temp_b1 = 5778*(binary.p1.star.lum/(binary.p1.star.rad*binary.p1.star.rad))**0.25
                    temp_b2 = 5778*(binary.p2.star.lum/(binary.p2.star.rad*binary.p2.star.rad))**0.25
                    temp_binary = (temp_b1*binary.p1.star.lum+temp_b2*binary.p2.star.lum)/(binary.p1.star.lum+binary.p2.star.lum)
                    data['temp_cm']= np.append(temp_single,temp_binary)
                    data['type_cm']= np.append(single.star.type,np.max([binary.p1.star.type,binary.p2.star.type],axis=0))
                    
                if (self.color_mode == 'galev'):
                    binary.p1.galev.convertToFlux()
                    binary.p2.galev.convertToFlux()
                    bmag = binary.p1.galev + binary.p2.galev
                    binary.p1.galev.convertToMag()
                    binary.p2.galev.convertToMag()
                    mag_cm_list = []
                    for key in self.galev_color:
                        mag_cm_list.append(np.concatenate((single.galev[key], bmag[key])))
                    data['galev_mag_cm'] = np.array(mag_cm_list)
            

    def getColor(self):
        colors='w'
        if (self.color_mode=='logtemp'):
            log_temp=(np.log10(self.temp)-np.log10(self.temp_min))/(np.log10(self.temp_max)-np.log10(self.temp_min))
            colors=cm.rainbow(1.0-log_temp)
        elif (self.color_mode=='loglum'):
            log_lum=(np.log10(self.lum)-np.log10(self.lum_min))/(np.log10(self.lum_max)-np.log10(self.lum_min))
            colors=cm.rainbow(1.0-log_lum)
        elif (self.color_mode=='ekin'):
            colors=cm.hot((self.ekin - self.ekin_min)/(self.ekin_max-self.ekin_min))
        elif (self.color_mode=='pot'):
            colors=cm.hot((self.pot - self.pot_min)/(self.pot_max-self.pot_min))
        elif (self.color_mode=='etot'):
            colors=cm.hot((self.etot - self.etot_min)/(self.etot_max-self.etot_min))
        elif (self.color_mode=='galev'):
            mag_min = self.galev_mag_range[0]
            mag_max = self.galev_mag_range[1]

            #colors= (self.galev_mag[0] - mag_min[None,:])/(mag_max-mag_min)[None,:]
            dmag = (self.galev_mag[0] - self.galev_mag[1] - mag_min)/(mag_max-mag_min)
            colors=cm.rainbow(dmag)
            #print(colors.max(axis=0),colors.min(axis=0))
        return colors


def plotOne(file_path, axe, plots, core, lagr, **kwargs):
    data = Data(**kwargs)
    data.read(file_path, core)

    if ('format_time' in kwargs.keys()):
        axe[0].set_title('T = '+kwargs['format_time'] % data['t'])
    else:
        axe[0].set_title('T = '+str(data['t']))
    
    ptcls=[]
    plot_item = plots['label']
    iaxe=0
    for pi in plot_item:
        if pi[0] == 'main':
            ptcls += plots['plot'][iaxe].plot(data, plots['xcm'][iaxe], plots['ycm'][iaxe])
            ptcls += [plots['xcm'][iaxe], plots['ycm'][iaxe]]

        if pi[0] == 'plot_HRdiagram':
            if ('lum_cm' in data.keys()):
                ptcls = ptcls + plots['plot'][iaxe].plot(data.lum_cm, data.temp_cm, data.type_cm)
            else:
                ptcls = ptcls + plots['plot'][iaxe].plot(data.lum, data.temp, data.type)

        if pi[0] == 'plot_semi_ecc':
            ptcls = ptcls + plots['plot'][iaxe].plot(data)

        if pi[0] == 'plot_lagr':
            ptcls = ptcls + plots['plot'][iaxe].plot(lagr, data['t'])
        iaxe +=1
                              
    return ptcls

def initPlot(axe, model_title, plot_item, lagr, **kwargs):
    plots=dict()

    if ('compare_in_column' in kwargs.keys()):
        plots['title'] = axe[0].text(0.5, 1.2, model_title, verticalalignment='center', horizontalalignment='center',  transform = axe[0].transAxes, color='white')
    else:
        plots['title'] = axe[0].text(-0.2, 0.5, model_title, verticalalignment='center', horizontalalignment='center',  transform = axe[0].transAxes, color='white', rotation=90)
    plots['plot'] = dict()
    plots['xcm'] = dict()
    plots['ycm'] = dict()
    plots['label'] = plot_item

    iaxe=0
    imain=0
    for pi in plot_item:
        if pi[0] == 'main':
            plots['plot'][iaxe]=PlotXY()
            kwargs_sub=dict()
            kwargs_sub['plot_mode'] = pi[1]
            for key in ['boxsize','x_min','x_max','y_min','y_max']:
                if (key+'_list' in kwargs.keys()):
                    kwargs_sub[key] = kwargs[key+'_list'][imain]
            plots['plot'][iaxe].init(axe[iaxe], **kwargs_sub , **kwargs)
            plots['xcm'][iaxe] = axe[iaxe].text(.05, 0.95, '', transform = axe[iaxe].transAxes, color='white')
            plots['ycm'][iaxe] = axe[iaxe].text(.05, 0.9, '', transform = axe[iaxe].transAxes, color='white')
            imain += 1
        
        if pi[0] == 'plot_HRdiagram':
            plots['plot'][iaxe]=PlotHR()
            plots['plot'][iaxe].init(axe[iaxe],**kwargs)

        if pi[0] == 'plot_semi_ecc':
            plots['plot'][iaxe]=PlotSemiEcc()
            plots['plot'][iaxe].init(axe[iaxe],**kwargs)

        if pi[0] == 'plot_lagr':
            if not 'time_min' in kwargs.keys():
                kwargs['time_min'] = lagr.time[0]
            if not 'time_max' in kwargs.keys():
                kwargs['time_max'] = lagr.time[-1]
            plots['plot'][iaxe]=PlotLagr()
            plots['plot'][iaxe].init(axe[iaxe], lagr, **kwargs)
        iaxe +=1 

    return plots

def initFig(model_list, frame_xsize, frame_ysize, ncol, nplots, compare_in_column):
    xsize=frame_xsize
    ysize=frame_ysize
    nrow= 1
    if (len(model_list)>1): 
        nrow = len(model_list)
        if (compare_in_column):
            nrow = ncol
            ncol = len(model_list)
        ysize = frame_ysize*nrow
        xsize = frame_xsize*ncol
    elif (ncol<0):
        ncol = nplots
        xsize = frame_xsize*nplots
    else:
        xsize = ncol*frame_xsize
        nrow = int(nplots/ncol)
        if (nrow*ncol<nplots): nrow+=1
        ysize = nrow*frame_ysize

    fig, axe = plt.subplots(nrow,ncol,figsize=(xsize, ysize))
    # axe should be two dimension, rows contain different models and columns contain different types of plots
    if (len(model_list)>1):
        if compare_in_column: 
            if (nrow>1): # swap rows and columns
                axe=[[axe[i][j] for i in range(nrow)] for j in range(ncol)]
            else: 
                axe=[[axe[i]] for i in range(ncol)] # enclose each model in additional []
        else:
            if (ncol==1): axe=[[axe[i]] for i in range(nrow)] # enclose each model in additional []
    else:
        # in only one model case, only one row exist, cols of plots should be flatten into one row and enclosed by []
        if (nrow>1) & (ncol>1): axe=[axe.flatten()]  
        if (nrow==1) | (ncol==1): axe=[axe] # one row or one column case, enclosed by []
        if (nrow==1) & (ncol==1): axe=[axe] # only one plot case, one more [] is needed: axe should be 2D [[axe]]
    return fig, axe


def createImage(_path_list, model_list, frame_xsize, frame_ysize, ncol, plot_item, core, lagr, dpi, **kwargs):
    n_frame = len(_path_list)
    use_previous = False
    if ('use_previous' in kwargs.keys()): use_previous = kwargs['use_previous']

    if (n_frame>0):
        nplots = len(plot_item)

        if (len(model_list)>1): ncol = nplots
        compare_in_column = False
        if ('compare_in_column' in kwargs.keys()): compare_in_column = kwargs['compare_in_column']

        fig, axe = initFig(model_list, frame_xsize, frame_ysize, ncol, nplots, compare_in_column)

        plots=dict()
        for i in range(len(model_list)):
            plots[i] = initPlot(axe[i], model_list[i][2],plot_item, lagr[i], **kwargs)

        for k in range(n_frame):
            file_path = _path_list[k]
            if (use_previous):
                if (os.path.exists(file_path+'.png')):
                    print('find existing %s' % file_path)
                    continue

            for i in range(len(model_list)):
                mi = model_list[i][0]
                prefix = model_list[i][1]+'.' if model_list[i][1] != '' else ''
                plotOne(mi+'/'+prefix+file_path, axe[i], plots[i], core[i], lagr[i], **kwargs)
            print('processing ',file_path)
            fig.savefig(file_path+'.png',bbox_inches = "tight", dpi=dpi)

    return n_frame

if __name__ == '__main__':

    filename='dat.lst'
    model_path=''
    core_file='data.core'
    lagr_file='data.lagr'
    read_lagr_data = False
    fps = 30
    dpi = None
    output_file = 'movie'
    plot_item=[]
    ncol= -1
    frame_xsize= 4
    frame_ysize= 4
    plot_format='mp4'
    n_cpu = 0
    plot_images=True
    format_time = '%g'

    pxy = PlotXY()
    phr = PlotHR()
    pse = PlotSemiEcc()
    plagr=PlotLagr()
    data = Data()
    
    def usage():
        print("A tool to generate a movie for the evolution of positions, velocities, HR diagram, ")
        print("binary semi-major axis-eccentricity and Lagrangian radii by reading the snapshots.")
        print("Different plots can be combined together based on the provided options.")
        print("Multiple simulations can also be compared if the numbers and times of snapshots are the same.")
        print("Usage: petar.movie [options] data_list_filename")
        print("data_list_filename: 1) In one model mode, this file contains a list of snapshot data pathes, each line for one snapshot. For example:")
        print("                    ./data.0")
        print("                    ./data.1")
        print("                     ...")
        print("        Or (generated by petar.data.gether):")
        print("                    data.0")
        print("                    data.1")
        print("                    ...")
        print("    2) In the comparison mode (-l), this file only contains the indices of snapshots without path and prefix. For example:")
        print("                    0")
        print("                    1")
        print("                    ...")
        print("option: default values are shown at last")
        print("  -h(--help): help")
        print("  -f [F]: output frame FPS: ",fps)
        print("  -m [S]: add one panel for ploting distribution of individual particles: ",pxy.plot_mode)
        print("          The data types of x- and y-axes are combined by '-'.")
        print("          For exmaple, -m x-y and -m vx-vy plot particle positions x vs. y and velocities v_x vs. v_y, respectively")
        print("          The support data type names: ")
        print("              x, y, z: position x, y and z")
        print("              vx, vy, vz: velocity v_x, v_y and v_z")
        print("              rxy: projected distance at x-y plane, sqrt(x*x+y*y)")
        print("              vxy: projected velocity at x-y plane, sqrt(vx*vx+vy*vy)")
        print("              r: distance to the center, sqrt(x*x+y*y+z*z)")
        print("              v: velocity value, sqrt(vx*vx+vy*vy+vz*vz)")
        print("              vr: radial velocity referring to the center, (x*vx+y*vy+z*vz)/r")
        print("              ra: RA in the ICRS frame")
        print("              dec: Dec in the ICRS frame")
        print("              pmracosdec: proper motion of RA*cos(Dec) in the ICRS frame")
        print("              pmdec: proper motion of Dec in the ICRS frame")
        print("              lon: Longtitude in the Galactocentric frame")
        print("              lat: Latitude in the Galactocentric frame")
        print("              pmlon: proper motion of Longtitude in the Galactocentric frame")
        print("              pmlat: proper motion of Latitude in the Galactocentric frame")
        print("          For ra, dec, lon, lat and proper motion, the snapshot data must use the astronomical units (pc, pc/Myr)")
        print("          The mode can be combined by ','. For exmaple, '-m x-y,y-z' provide two plots: x-y and y-z")
        print("          In this case, -R or --x-min/max, --y-min/max should also set the range for each plot,")
        print("          such as -R 10,10, --x-min -10,-10")
        print("  -R [F]: x- and y-axis length of -m; suppressed when --x-min/max, --y-min/max are used: ",pxy.boxsize)
        print("  -H    : add one panel of HR-diagram")
        print("  -b    : add one panel of semi-ecc diagram for binaries")
        print("          sizes indicate the mass based on scaling option --markser-scale and --mass-power")
        print("  -L [S]: add one panel of Lagrangian radii evolution, argument is filename of lagrangian data", lagr_file)
        print("          Here the filename is not used in the comparison mode")
        print("  -G [F]: gravitational constant for calculating binary orbit: ",data.G)
        print("  -o [S]: output movie filename: ",output_file)
        print("  -p    : Use previous generated png images to speed up the movie generation")
        print("  -s [S]: snapshot format: ascii, binary or npy : ", data.snapshot_format)
        print("              ascii:  all snapshots are in ascii format")
        print("              binary: all snapshots are in binary format")
        print("              npy: snapshots from petar (data.*) are in binary format;")
        print("                   snapshots from petar.data.process (data.*.[single/binary]) are in npy format")
        print("  -l [S]: the filename of a list of pathes to different models, this will switch on the comparison mode.")
        print("          Each line contains three values: directory path of model, prefix of filename of snapshots, name of model.")
        print("          When data is reading, the path and the prefix will be added in front of the filenames. ")
        print("          For example, path:'./' ; prefix:'data' ; reading file: './data.[0-9*]'.")
        print("          The number of snapshots should also be the same for all models.")
        print("  -i  [S]: interrupt mode used in petar: no, base, bse, mobse: ",data.interrupt_mode)
        print("  -t  [S]: external mode used in petar: no, galpy: ",data.external_mode)
        print("  -c  [S]: color type for particles plot of -m option: loglum, logtemp, ekin, pot, etot, galev, white: ",data.color_mode)
        print("              loglum: log(luminosity)")
        print("              logtemp: log(temperature)")
        print("              ekin: kinetic energy")
        print("              pot: potential")
        print("              etot: total energy")
        print("              galev: use galev photometry color, need to read [prefix].galev.mag files; also need to properly set galev related options")
        print("              white: pure white color")
        print("  --compare-in-column: in comparison mode, models are compared in columns instead of rows.")
        print("  --generate-binary [I]: 0: no binary, 1: detect binary by using KDtree (slow), 2: read single and binary data generated by petar.data.process: ", data.generate_binary)
        print("  --n-cpu       [I]: number of CPU processors to use: all CPU cores")
        print("  --x-min       [F]: minimum of the main plot x-axis range",pxy.x_min)
        print("  --x-max       [F]: maximum of the main plot x-axis range",pxy.x_max)
        print("  --y-min       [F]: minimum of the main plot y-axis range",pxy.y_min)
        print("  --y-max       [F]: maximum of the main plot y-axis range",pxy.y_max)
        print("  --lum-min     [F]: minimum lumonisity: ",phr.lum_min)
        print("  --lum-max     [F]: maximum lumonisity: ",phr.lum_max)
        print("  --temp-min    [F]: minimum temperature: ",phr.temp_min)
        print("  --temp-max    [F]: maximum temperature: ",phr.temp_max)
        print("  --semi-min    [F]: minimum semi-major axis: ",pse.semi_min)
        print("  --semi-max    [F]: minimum semi-major axis: ",pse.semi_max)
        print("  --ecc-min     [F]: minimum eccentricity: ",pse.ecc_min)
        print("  --ecc-max     [F]: maximum eccentricity: ",pse.ecc_max)
        print("  --bin-rmax    [F]: maximum distance for color scaling in semi-ecc plot: ",pse.bin_rmax)
        print("  --bin-rmin    [F]: minimum distance for color scaling in semi-ecc plot: ",pse.bin_rmin)
        print("  --bin-color   [I]: color mode for binary (semi-ecc) plot: ",pse.bin_color)
        print("                         distance: scale with the log(distance) of binaries to the center, normalized by --bin-rmin and --bin-rmax")
        print("                                   color map 'hot_r' is used, redder color represent binaries more distant to the center")
        print("                         white: pure white")
        print("  --time-min    [F]: minimum time in evolution plot (x-axis): auto determine from Lagrangian data")
        print("  --time-max    [F]: maximum time in evolution plot (x-axis)): auto determine from Lagrangian data")
        print("  --ekin-min    [F]: minimum kinetic energy: ",data.ekin_min)
        print("  --ekin-max    [F]: maximum kinetic energy: ",data.ekin_max)
        print("  --pot-min     [F]: minimum potential: ",data.pot_min)
        print("  --pot-max     [F]: maximum potential: ",data.pot_max)
        print("  --etot-min    [F]: minimum total energy: ",data.etot_min)
        print("  --etot-max    [F]: maximum total energy: ",data.etot_max)
        print("  --rlagr-min   [F]: minimum radius in Lagrangian plot: ", plagr.rlagr_min)
        print("  --rlagr-max   [F]: maximum radius in Lagrangian plot: ", plagr.rlagr_max)
        print("  --rlagr-scale [S]: scaling of Lagrangian radii in the plot (y-axis): ",plagr.rlagr_scale)
        print("  --lagr-energy    : option calc_energy for reading lagr data")
        print("  --lagr-type   [S]: option add_star_type for reading lagr data")
        print("  --lagr-mfrac  [S]: option for mass fraction for reading lagr data")
        print("  --unit-length [S]: set label of length unit for x, y, z and semi: no print")
        print("  --unit-vel    [S]: set label of velocity unit: no print")
        print("  --format-time [S]: set print format of time: ",format_time)
        print("  --skiprows    [I]: number of rows to escape when read snapshot: Unset")
        print("  --plot-ncols  [I]: column number of panels: same as panels")
        print("  --plot-xsize  [F]: x size of panel: ",frame_xsize)
        print("  --plot-ysize  [F]: y size of panel: ",frame_ysize)
        print("  --cm-mode     [S]: plot origin position determination: ",data.cm_mode)
        print("                       density: density center;")
        print("                       average: average of x,y;")
        print("                       core: use core data file generated from petar.data.process;")
        print("                       none: use origin of snapshots.")
        print("  --core-file   [S]: core data file name, not used in the comparison mode: ",core_file)
        print("  --cm-boxsize  [F]: boxsize to search the coordinate center for the x-y plot: 5.0 times ploting size (-R)")
        print("  --n-layer-cross [I]: number of layers of crosses for particles in the x-y plot: 5")
        print("  --n-layer-point [I]: number of layers of points for particles in the x-y plot: 10")
        print("  --layer-alpha   [F]: transparency factor of layers in the x-y plot: 2.5")
        print("  --marker-scale  [F]: amplify the size of markers in x-y plot: 1.0")
        print("  --mass-power    [F]: the power index of mass to obtain sizes of markers: 1.0")
        print("  --size-mode    [S]: point size mode: ", pxy.size_mode)
        print("                      mass: scale with mass with --mass-power")
        print("                      loglum: scale with loglum")
        print("  --galev-filter [S]: filter list, seperate by comma:", ','.join(data.galev_filter))
        print("  --galev-mode   [S]: galev mode of values with choices: abs_mag, app_mag, abs_flux, abs_flux", data.galev_mode)
        print("  --galev-color  [S]: two filters for mapping rainbow colors:", ','.join(data.galev_color))
        print("  --galev-mag-range[S]: magnitude range (min,max) to normalize colores:", ','.join([str(x) for x in data.galev_mag_range]))
        print("  --suppress-images: do not plot snapshot images (png files) and use matplotlib.animation instead of imageio, this cannot use multi-processing, much slower")
        print("  --format        [S]: video format, require imageio installed, for some formats (e.g. avi, mp4) may require ffmpeg and imageio-ffmpeg installed: ", plot_format)
        print("  --dpi           [F]: dpi of image: ", dpi)
        print("PS:: Each panel of plots can be added mutliple times (the order is recored)")

    try:
        shortargs = 'm:s:f:R:z:o:c:G:l:L:i:t:psHbh'
        longargs = ['help','n-cpu=','lum-min=','lum-max=','temp-min=','temp-max=',
                    'semi-min=','semi-max=','ecc-min=','ecc-max=','bin-rmax=','bin-rmin=','bin-color=',
                    'ekin-min=','ekin-max=','pot-min=','pot-max=','etot-min=','etot-max=',
                    'rlagr-min=','rlagr-max=','rlagr-scale=',
                    'lagr-energy','lagr-type=','lagr-mfrac=',
                    'time-min=','time-max=','x-min=','x-max=','y-min=','y-max=',
                    'unit-length=','unit-vel=','format-time=',
                    'skiprows=','generate-binary=',
                    'plot-ncols=','plot-xsize=','plot-ysize=',
                    'suppress-images','format=','cm-mode=','core-file=',
                    'n-layer-cross=','n-layer-point=','layer-alpha=','marker-scale=','mass-power=','size-mode=',
                    'galev-filter=','galev-mode=','galev-color=','galev-mag-range=',
                    'cm-boxsize=','compare-in-column','dpi=']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        kwargs=dict()
        sky_modes=['ra','dec','pmracosdec','pmdec','lon','lat','pmlon','pmlat']

        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-f'):
                fps = float(arg)
            elif opt in ('-m'):
                plot_mode_str = arg
                for i in plot_mode_str.split(','):
                    plot_item.append(['main', i])
                    xname,yname = i.split('-')
                    if (xname in sky_modes) | (yname in sky_modes):
                        kwargs['get_skycoord'] = True
            elif opt in ('-R'):
                boxsize_str = arg
                kwargs['boxsize_list'] = []
                for i in boxsize_str.split(','):
                    kwargs['boxsize_list'].append(float(i))
            elif opt in ('-H'):
                plot_item.append(['plot_HRdiagram'])
            elif opt in ('-b'):
                plot_item.append(['plot_semi_ecc'])
            elif opt in ('-o'):
                output_file = arg
            elif opt in ('-G'):
                kwargs['G'] = float(arg)
            elif opt in ('-p'):
                kwargs['use_previous'] = True
            elif opt in ('-l'):
                model_path = arg
            elif opt in ('-L'):
                plot_item.append(['plot_lagr'])
                lagr_file = arg
                read_lagr_data=True
            elif opt in ('-i'):
                kwargs['interrupt_mode'] = arg
            elif opt in ('-t'):
                kwargs['external_mode'] = arg
            elif opt in ('-s'):
                kwargs['snapshot_format'] = arg
            elif opt in ('-c'):
                kwargs['color_mode'] = arg
            elif opt in ('--n-cpu'):
                n_cpu = int(arg)
            elif opt in ('--x-min'):
                x_min_str = arg
                kwargs['x_min_list'] = []
                for i in x_min_str.split(','):
                    kwargs['x_min_list'].append(float(i))
            elif opt in ('--x-max'):
                x_max_str = arg
                kwargs['x_max_list'] = []
                for i in x_max_str.split(','):
                    kwargs['x_max_list'].append(float(i))
            elif opt in ('--y-min'):
                y_min_str = arg
                kwargs['y_min_list'] = []
                for i in y_min_str.split(','):
                    kwargs['y_min_list'].append(float(i))
            elif opt in ('--y-max'):
                y_max_str = arg
                kwargs['y_max_list'] = []
                for i in y_max_str.split(','):
                    kwargs['y_max_list'].append(float(i))
            elif opt in ('--lum-min'):
                kwargs['lum_min'] = float(arg)
            elif opt in ('--lum-max'):
                kwargs['lum_max'] = float(arg)
            elif opt in ('--temp-min'):
                kwargs['temp_min'] = float(arg)
            elif opt in ('--temp-max'):
                kwargs['temp_max'] = float(arg)
            elif opt in ('--semi-min'):
                kwargs['semi_min'] = float(arg)
            elif opt in ('--semi-max'):
                kwargs['semi_max'] = float(arg)
            elif opt in ('--ecc-min'):
                kwargs['ecc_min'] = float(arg)
            elif opt in ('--ecc-max'):
                kwargs['ecc_max'] = float(arg)
            elif opt in ('--bin-rmax'):
                kwargs['bin_rmax'] = float(arg)
            elif opt in ('--bin-rmin'):
                kwargs['bin_rmin'] = float(arg)
            elif opt in ('--bin-color'):
                kwargs['bin_color'] = arg
            elif opt in ('--ekin-min'):
                kwargs['ekin_min'] = float(arg)
            elif opt in ('--ekin-max'):
                kwargs['ekin_max'] = float(arg)
            elif opt in ('--pot-min'):
                kwargs['pot_min'] = float(arg)
            elif opt in ('--pot-max'):
                kwargs['pot_max'] = float(arg)
            elif opt in ('--etot-min'):
                kwargs['etot_min'] = float(arg)
            elif opt in ('--etot-max'):
                kwargs['etot_max'] = float(arg)
            elif opt in ('--rlagr-min'):
                kwargs['rlagr_min'] = float(arg)
            elif opt in ('--rlagr-max'):
                kwargs['rlagr_max'] = float(arg)
            elif opt in ('--time-min'):
                kwargs['time_min'] = float(arg)
            elif opt in ('--time-max'):
                kwargs['time_max'] = float(arg)
            elif opt in ('--rlagr-scale'):
                kwargs['rlagr_scale'] = arg
            elif opt in ('--lagr-energy'):
                kwargs['calc_energy'] = True
            elif opt in ('--lagr-type'):
                kwargs['add_star_type'] = [x for x in arg.split(',')]
            elif opt in ('--lagr-mfrac'):
                kwargs['mass_fraction'] = np.array([float(x) for x in arg.split(',')])
            elif opt in ('--cm-mode'):
                kwargs['cm_mode']= arg
            elif opt in ('--cm-boxsize'):
                kwargs['cm_boxsize'] = float(arg)
            elif opt in ('--core-file'):
                core_file = arg
            elif opt in ('--unit-length'):
                kwargs['unit_length'] = arg
            elif opt in ('--unit-vel'):
                kwargs['unit_vel'] = arg
            elif opt in ('--format-time'):
                kwargs['format_time'] = arg
            elif opt in ('--skiprows'):
                kwargs['skiprows'] = int(arg)
            elif opt in ('--generate-binary'):
                kwargs['generate_binary']=int(arg)
            elif opt in ('--plot-ncols'):
                ncol = int(arg)
            elif opt in ('--plot-xsize'):
                frame_xsize = float(arg)
            elif opt in ('--plot-ysize'):
                frame_ysize = float(arg)
            elif opt in ('--n-layer-cross'):
                kwargs['nlayer_cross'] = int(arg)
            elif opt in ('--n-layer-point'):
                kwargs['nlayer_point'] = int(arg)
            elif opt in ('--layer-alpha'):
                kwargs['alpha_amplifier'] = float(arg)
            elif opt in ('--marker-scale'):
                kwargs['marker_scale'] = float(arg)
            elif opt in ('--mass-power'):
                kwargs['mass_power'] = float(arg)
            elif opt in ('--size-mode'):
                kwargs['size_mode'] = arg
            elif opt in ('--galev-filter'):
                kwargs['galev_filter'] =  [x for x in arg.split(',')]
            elif opt in ('--galev-color'):
                kwargs['galev_color'] = [x for x in arg.split(',')]
            elif opt in ('--galev-mode'):
                kwargs['galev_mode'] = arg
            elif opt in ('--galev-mag-range'):
                kwargs['galev_mag_range'] = np.array([float(x) for x in arg.split(',')])
            elif opt in ('--suppress-images'):
                plot_images = False
            elif opt in ('--format'):
                plot_format = arg
            elif opt in ('--dpi'):
                dpi = float(arg)
            elif opt in ('--compare-in-column'):
                kwargs['compare_in_column'] = True
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    for key in kwargs.keys():
        print(key,kwargs[key])

    filename = remainder[0]

    fl = open(filename,'r')
    file_list = fl.read()
    path_list = file_list.splitlines()
    fl.close()

    core=dict()
    read_core=False
    if ('cm_mode' in kwargs.keys()):
        if (kwargs['cm_mode']=='core'): read_core=True
    if (not 'generate_binary' in kwargs.keys()): 
        read_core = True
    elif (kwargs['generate_binary']==2): read_core=True

    lagr=dict()

    model_list= [['.','','']]
    if (model_path!=''):
        fl = open(model_path,'r')
        lines = fl.read().splitlines()
        model_list = [lines[i].split() for i in range(len(lines))]
        fl.close()

        for i in range(len(model_list)):
            core[i] = petar.Core()
            if (read_core):
                core[i].loadtxt(model_list[i][0]+'/'+model_list[i][1]+'.core')

            lagr[i] = petar.LagrangianMultiple(**kwargs)
            if (read_lagr_data):
                lagr[i].loadtxt(model_list[i][0]+'/'+model_list[i][1]+'.lagr')
    else:
        core[0] = petar.Core()
        if (read_core):
            core[0].loadtxt(core_file)
        lagr[0] = petar.LagrangianMultiple(**kwargs)
        if (read_lagr_data):
            lagr[0].loadtxt(lagr_file)

    if (len(plot_item)==0): plot_item=[['main','x-y']]

    if (plot_images):
        if (n_cpu==int(0)):
            n_cpu = mp.cpu_count()
        pool = mp.Pool(n_cpu)
        
        n_files=len(path_list)
        n_pieces = np.ones(n_cpu)*int(n_files/n_cpu)
        n_left = n_files%n_cpu
        n_pieces[:n_left]+=1
        n_offset=np.append([0],n_pieces.cumsum()).astype(int)
        
        file_part = [path_list[n_offset[i]:n_offset[i+1]] for i in range(n_cpu)]
        results=[None]*n_cpu
        if (n_cpu==1):
            createImage(file_part[0], model_list, frame_xsize, frame_ysize, ncol, plot_item, core, lagr, dpi, **kwargs)
        else:
            for rank in range(n_cpu):
                results[rank]=pool.apply_async(createImage, (file_part[rank], model_list, frame_xsize, frame_ysize, ncol, plot_item, core, lagr, dpi), kwargs)

        # Step 3: Don't forget to close
        pool.close()
        pool.join()

        png_list = [path_list[i]+'.png' for i in range(n_files)]
        create_movie(png_list, fps, output_file+'.'+plot_format)

    else:
        nplots = len(plot_item)

        if (len(model_list)>1): ncol = nplots
        compare_in_column = False
        if ('compare_in_column' in kwargs.keys()): compare_in_column = kwargs['compare_in_column']

        fig, axe = initFig(model_list, frame_xsize, frame_ysize, ncol, nplots, compare_in_column)
        plots=dict()
        for i in range(len(model_list)):
            mi = model_list[i][0]
            plots[i]=initPlot(axe[i], model_list[i][2], plot_item, lagr[i], **kwargs)

        def init():
            iaxe=0
            ptcls=[]
            for i in range(len(model_list)):
                if ('format_time' in kwargs.keys()):
                    axe[i][0].set_title('T = '+kwargs['format_time'] % 0)
                else:
                    axe[i][0].set_title('T = 0')
                
            return ptcls
     
        def animate(k):
     
            file_path = path_list[k]
            print('process ',file_path)
            ptcls=[]
            for i in range(len(model_list)):
                mi = model_list[i][0]
                prefix = model_list[i][1]+'.' if model_list[i][1] != '' else ''
                pi=plotOne(mi+'/'+prefix+file_path, axe[i], plots[i], core[i], lagr[i], **kwargs)
                ptcls=ptcls+pi
     
            return ptcls

        n_frame = len(path_list)
        anime = animation.FuncAnimation(fig, animate, init_func=init,
                                        frames=n_frame, interval=200, blit=True)
    
        anime.save(output_file+'.'+plot_format, fps=fps, extra_args=['-vcodec', 'libx264'])    

