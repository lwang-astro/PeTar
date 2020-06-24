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
        self.boxsize = 2
        self.framescale = 1
        self.nlayer_cross = 5
        self.nlayer_point = 10
        self.alpha_amplifier = 2.5
        self.temp_min = 1000
        self.temp_max = 50000
        self.ptcls=[]

    def init(self, axe, **kwargs):
        for key in self.__dict__.keys():
            if (key in kwargs.keys()): self.__dict__[key] = kwargs[key]
        #if ('boxsize' in kwargs.keys()): boxsize = kwargs['boxsize']
        #if ('nlayer_cross' in kwargs.keys()): nlayer_cross = kwargs['nlayer_cross']
        #if ('alpha_amplifier' in kwargs.keys()): alpha_amplifier = kwargs['alpha_amplifier']

        nlayer = self.nlayer_cross + self.nlayer_point
        self.nlayer = nlayer
        alphascale = np.linspace(1,nlayer,nlayer)*2.0/(nlayer*(nlayer+1))*self.alpha_amplifier
        #print('Alpha layer sequence:',alphascale,' sum:',alphascale.sum())
        self.sizescale = np.logspace(0,3,nlayer)[::-1]
        #print('Size layer sequence:',self.sizescale)

        axe.set_xlim(-self.boxsize,self.boxsize)
        axe.set_ylim(-self.boxsize,self.boxsize)
        axe.set_aspect(1.0)
        if ('unit_length' in kwargs.keys()):
            unit_label = '['+kwargs['unit_length']+']'
            axe.set_xlabel('x'+unit_label)
            axe.set_ylabel('y'+unit_label)
        else:
            axe.set_xlabel('x')
            axe.set_ylabel('y')

        for i in range(self.nlayer_cross):
            pt =axe.scatter([],[],marker='+',alpha=alphascale[i],edgecolors='none')
            self.ptcls.append(pt)
        for i in range(self.nlayer_point):
            pt =axe.scatter([],[],alpha=alphascale[i],edgecolors='none')
            self.ptcls.append(pt)
        return self.ptcls

    def getColor(self, temperature_eff):
        log_temp=(np.log10(temperature_eff)-np.log10(self.temp_min))/(np.log10(self.temp_max)-np.log10(self.temp_min))
        colors=cm.rainbow(1.0-log_temp)
        return colors

    def plot(self, x, y, mass, colors):
        for i in range(self.nlayer):
            sizes = mass*self.sizescale[i]*self.framescale
            #sizes = (np.log10(luminosity)-np.log10(lum_min)+1)
            #sizes = luminosity
            self.ptcls[i].set_offsets(np.array([x,y]).transpose())
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
        self.ptcls.append(axe.scatter([],[],marker='.',edgecolors='none'))
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
        self.ptcls = []

    def init(self, axe, **kwargs):
        for key in self.__dict__.keys():
            if (key in kwargs.keys()): self.__dict__[key] = kwargs[key]
        self.ptcls.append(axe.scatter([],[],marker='+',edgecolors='none'))
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

    def plot(self, semi, ecc, types):
        colors = cm.rainbow(types/13.0)
        self.ptcls[0].set_offsets(np.array([semi, ecc]).transpose())
        self.ptcls[0].set_color(colors)
        return self.ptcls

class Data:
    def __init__(self, **kwargs):
        self.xcol = -1
        self.ycol = -1
        self.mcol = -1
        self.skiprows = 0
        self.generate_binary=2
        self.interrupt_mode = 'bse'
        self.G = 0.00449830997959438 # pc^3/(Msun*Myr^2)
        self.semi_max = 0.1
        self.cm_mode = 'density'

        for key in self.__dict__.keys():
            if (key in kwargs.keys()): self.__dict__[key] = kwargs[key]

    def __getitem__(self, k):
        return self.__dict__[k]

    def keys(self):
        return self.__dict__.keys()

    def read(self, file_path):
        data=self.__dict__
        xcol = self.xcol
        ycol = self.ycol
        mcol = self.mcol
        skiprows = self.skiprows

        if (xcol>=0) & (ycol>=0):
            if (mcol>=0):
                data['x'], data['y'] ,data['mass'] = np.loadtxt(file_path, unpack=True, usecols=(xcol,ycol,mcol), ndmin=2, skiprows=skiprows)
            else:
                data['x'], data['y'] = np.loadtxt(file_path, unpack=True, usecols=(xcol,ycol), ndmin=2, skiprows=skiprows)
                data['mass'] = np.ones(data['x'].size)
            data['t'] = file_path
        else:
            fp = open(file_path, 'r')
            header=fp.readline()
            file_id, n_glb, data['t'] = header.split()
            fp.close()

            if (self.generate_binary>0):
                if (self.generate_binary==2):
                    single = petar.Particle(interrupt_mode=self.interrupt_mode)
                    p1 = petar.Particle(interrupt_mode=self.interrupt_mode)
                    p2 = petar.Particle(interrupt_mode=self.interrupt_mode)
                    binary = petar.Binary(p1,p2)
                    if os.path.getsize(file_path+'.single')>0:
                        single.loadtxt(file_path+'.single')
                    if os.path.getsize(file_path+'.binary')>0:
                        binary.loadtxt(file_path+'.binary')
                    data['x'] = np.concatenate((single.pos[:,0], binary.p1.pos[:,0], binary.p2.pos[:,0])) 
                    data['y'] = np.concatenate((single.pos[:,1], binary.p1.pos[:,1], binary.p2.pos[:,1])) 
                    data['mass'] = np.concatenate((single.mass, binary.p1.mass, binary.p2.mass))
                    data['semi'] = binary.semi
                    data['ecc'] = binary.ecc
                    data['state'] = binary.p1.binary_state
                    if (self.interrupt_mode=='bse'):
                        data['lum'] = np.concatenate((single.s_lum, binary.p1.s_lum, binary.p2.s_lum))
                        data['rad'] = np.concatenate((single.s_rad, binary.p1.s_rad, binary.p2.s_rad))
                        data['type']= np.concatenate((single.s_type,binary.p1.s_type,binary.p2.s_type))
                        temp_single = 5778*(single.s_lum/(single.s_rad*single.s_rad))**0.25
                        temp_b1 = 5778*(binary.p1.s_lum/(binary.p1.s_rad*binary.p1.s_rad))**0.25
                        temp_b2 = 5778*(binary.p2.s_lum/(binary.p2.s_rad*binary.p2.s_rad))**0.25
                        temp_binary = (temp_b1*binary.p1.s_lum+temp_b2*binary.p2.s_lum)/(binary.p1.s_lum+binary.p2.s_lum)
                        data['temp']= np.concatenate((temp_single, temp_b1, temp_b2))
                        data['lum_cm'] = np.append(single.s_lum, binary.p1.s_lum+binary.p2.s_lum)
                        data['temp_cm']= np.append(temp_single,temp_binary)
                        data['type_cm']= np.append(single.s_type,np.max([binary.p1.s_type,binary.p2.s_type],axis=0))
                else:
                    particles=petar.Particle(interrupt_mode=self.interrupt_mode)
                    particles.loadtxt(file_path,skiprows=1)
                    kdtree,single,binary = petar.findPair(particles, self.G, self.semi_max*2.0, True)
                    data['x'] = particles.pos[:,0]
                    data['y'] = particles.pos[:,1]
                    data['mass'] =particles.mass
                    data['semi'] = binary.semi
                    data['ecc'] = binary.ecc
                    data['state'] = binary.p1.binary_state
                    if (self.interrupt_mode=='bse'):
                        data['lum'] = particles.s_lum
                        data['rad'] = particles.s_rad
                        data['type']= particles.s_type
                        data['lum_cm'] = np.append(single.s_lum, binary.p1.s_lum+binary.p2.s_lum)
                        data['temp']= 5778*(data['lum']/(data['rad']*data['rad']))**0.25
                        temp_single = 5778*(single.s_lum/(single.s_rad*single.s_rad))**0.25
                        temp_b1 = 5778*(binary.p1.s_lum/(binary.p1.s_rad*binary.p1.s_rad))**0.25
                        temp_b2 = 5778*(binary.p2.s_lum/(binary.p2.s_rad*binary.p2.s_rad))**0.25
                        temp_binary = (temp_b1*binary.p1.s_lum+temp_b2*binary.p2.s_lum)/(binary.p1.s_lum+binary.p2.s_lum)
                        #temp_binary = temp_b1
                        #sel = (binary.p1.s_lum<binary.p2.s_lum)
                        #temp_binary[sel] = temp_b2[sel]
                        #print(temp_single.size,temp_binary.size,data['lum_cm'].size)
                        data['temp_cm']= np.append(temp_single,temp_binary)
                        data['type_cm']= np.append(single.s_type,np.max([binary.p1.s_type,binary.p2.s_type],axis=0))
            else:
                particles=petar.Particle(interrupt_mode=interrupt_mode)
                particles.loadtxt(file_path,skiprows=1)
                data['x'] = particles.pos[:,0]
                data['y'] = particles.pos[:,1]
                data['mass'] =particles.mass
                if (interrupt_mode=='bse'):
                    data['lum'] = particles.s_lum
                    data['rad'] = particles.s_rad
                    data['type']= particles.s_type
        
    def correctCM(self, boxsize):
        if (self.cm_mode=='density'):
            nbins=100
            xmid = np.average(np.abs(self.x))
            ymid = np.average(np.abs(self.y))
            #print(xmid,ymid)
            xmin = -5*xmid
            xmax =  5*xmid
            ymin = -5*ymid
            ymax =  5*ymid
            xbins=np.linspace(xmin,xmax,nbins)
            ybins=np.linspace(ymin,ymax,nbins)
            lbin=[xbins,ybins]
            counts, _, _ = np.histogram2d(self.x,self.y, bins=lbin)
            xp,yp=np.where(counts>0.2*counts.max())
            m = counts[xp,yp]
            mtot = m.sum()
            xcm=(xp*m).sum()/mtot/nbins*(xmax-xmin)+xmin
            ycm=(yp*m).sum()/mtot/nbins*(ymax-ymin)+ymin
            self.x = self.x - xcm
            self.y = self.y - ycm

            nbins=200
            xmid = np.average(np.abs(self.x))
            ymid = np.average(np.abs(self.y))
            xmin = -boxsize
            xmax =  boxsize
            ymin = -boxsize
            ymax =  boxsize
            xbins=np.linspace(xmin,xmax,nbins)
            ybins=np.linspace(ymin,ymax,nbins)
            lbin=[xbins,ybins]
            counts, _, _ = np.histogram2d(self.x,self.y, bins=lbin)
            xp,yp=np.where(counts>0.2*counts.max())
            m = counts[xp,yp]
            mtot = m.sum()
            xcm2=(xp*m).sum()/mtot/nbins*(xmax-xmin)+xmin
            ycm2=(yp*m).sum()/mtot/nbins*(ymax-ymin)+ymin
            self.x = self.x - xcm2
            self.y = self.y - ycm2
            xcm += xcm2
            ycm += ycm2
        else:
            xcm = self.x.sum()/self.x.size
            ycm = self.y.sum()/self.y.size
            self.x = self.x - xcm
            self.y = self.y - ycm
            sel=(self.x>-boxsize) & (self.x<boxsize) & (self.y>-boxsize) & (self.y<boxsize)
            nsel=sel.sum()
            xcm2 = (self.x[sel]).sum()/float(nsel)
            ycm2 = (self.y[sel]).sum()/float(nsel)
            self.x = self.x - xcm2 
            self.y = self.y - ycm2
            xcm += xcm2
            ycm += ycm2
        return xcm, ycm

def plotOne(file_path, axe, plots, **kwargs):
    data = Data(**kwargs)
    data.read(file_path)
    xcm, ycm=data.correctCM(plots['xy'].boxsize)

    if ('unit_time' in kwargs.keys()):
        unit_label = ' '+kwargs['unit_time']
        axe[0].set_title('T = '+data['t']+unit_label)
    else:
        axe[0].set_title('T = '+data['t'])
    
    colors='w'
    if ('lum' in data.keys()): colors = plots['xy'].getColor(data['temp'])

    ptcls=[]
    ptcls = ptcls + plots['xy'].plot(data.x, data.y, data.mass, colors)
    plots['xcm'].set_text(r'$x_{cm}=%f$' % xcm)
    plots['ycm'].set_text(r'$y_{cm}=%f$' % ycm)
    ptcls = ptcls + [plots['xcm'], plots['ycm']]
                                     
    if ('xyzoom' in plots.keys()): 
        ptcls = ptcls + plots['xyzoom'].plot(data.x, data.y, data.mass, colors)

    if ('hr' in plots.keys()): 
        if ('lum_cm' in data.keys()):
            ptcls = ptcls + plots['hr'].plot(data.lum_cm, data.temp_cm, data.type_cm)
        else:
            ptcls = ptcls + plots['hr'].plot(data.lum, data.temp. data.type)

    if ('se' in plots.keys()):
        ptcls = ptcls + plots['se'].plot(data.semi, data.ecc, data.state)
    return ptcls

def initPlot(axe, plot_zoom, plot_HRdiagram, plot_semi_ecc, **kwargs):
    plots=dict()

    iaxe=0
    framescale = 1
    if 'framescale' in kwargs.keys():
        framescale = kwargs['framescale']
        kwargs['framescale'] = 1
    plots['xy'] = PlotXY()
    plots['xy'].init(axe[iaxe], **kwargs)
    plots['xcm'] = axe[0].text(.05, 0.95, '', transform = axe[0].transAxes, color='white')
    plots['ycm'] = axe[0].text(.35, 0.95, '', transform = axe[0].transAxes, color='white')
    
    if (plot_zoom): 
        iaxe +=1 
        kwargs['boxsize'] = plots['xy'].boxsize / framescale
        kwargs['framescale'] = framescale
        plots['xyzoom']=PlotXY()
        plots['xyzoom'].init(axe[iaxe], **kwargs)
        
    if (plot_HRdiagram): 
        iaxe +=1
        plots['hr']=PlotHR()
        plots['hr'].init(axe[iaxe],**kwargs)

    if (plot_semi_ecc): 
        iaxe +=1 
        plots['se']=PlotSemiEcc()
        plots['se'].init(axe[iaxe],**kwargs)

    return plots

def initFig(frame_xsize, frame_ysize, ncol, nplots):
    xsize=frame_xsize
    ysize=frame_ysize
    nrow= 1
    if (ncol<0):
        ncol = nplots
        xsize = frame_xsize*nplots
    else:
        xsize = ncol*frame_xsize
        nrow = nplots/ncol
        if (nrow*ncol<nplots): nrow+=1
        ysize = nrow*frame_ysize

    fig, axe = plt.subplots(nrow,ncol,figsize=(xsize, ysize))
    if (nrow>1) & (ncol>1): axe=axe.flatten()
    if (nrow==1) & (ncol==1): axe=[axe]
    return fig, axe


def createImage(_path_list, frame_xsize, frame_ysize, ncol, plot_zoom, plot_HRdiagram, plot_semi_ecc, **kwargs):
    n_frame = len(_path_list)
    if (n_frame>0):
        nplots = 1
        if (plot_zoom): nplots+=1
        if (plot_HRdiagram): nplots+=1
        if (plot_semi_ecc): nplots+=1

        fig, axe = initFig(frame_xsize, frame_ysize, ncol, nplots)

        plots=initPlot(axe, plot_zoom, plot_HRdiagram, plot_semi_ecc, **kwargs)

        for k in range(n_frame):
            file_path = _path_list[k]
            plotOne(file_path, axe, plots,  **kwargs)
            fig.savefig(file_path+'.png',bbox_inches = "tight")
            print(file_path)
    return n_frame

if __name__ == '__main__':

    filename='dat.lst'
    fps = 30
    output_file = 'movie'
    plot_HRdiagram=False
    plot_zoom=False
    plot_semi_ecc=False
    ncol= -1
    frame_xsize= 8
    frame_ysize= 8
    plot_format='mp4'
    n_cpu = 0
    plot_images=True

    pxy = PlotXY()
    pxy_zoom = PlotXY()
    phr = PlotHR()
    pse = PlotSemiEcc()
    data = Data()
    
    def usage():
        print("A tool for processing a list of snapshot data to detect binaries, calculate Langragian radii and properties, get the density center and core radius")
        print("Usage: petar.motion.movie [options] data_list_filename")
        print("data_list_filename: A list of snapshot data path, each line for one snapshot")
        print("option: default values are shown at last")
        print("  -h(--help): help")
        print("  -f: output frame FPS: ",fps)
        print("  -R: movie box length: ",pxy.boxsize)
        print("  -z: plot zoom frame with emplification factor: not set")
        print("  -H: plot HR-diagram")
        print("  -b: plot semi-ecc diagram for binaries")
        print("  -G: gravitational constant for calculating binary orbit: ",data.G)
        print("  -o: output movie filename: ",output_file)
        print("  --n-cpu: number of CPU processors to use, only work for gif format: all CPU cores")
        print("  --lum-min: minimum lumonisity: ",phr.lum_min)
        print("  --lum-max: maximum lumonisity: ",phr.lum_max)
        print("  --temp-min: minimum temperature: ",phr.temp_min)
        print("  --temp-max: maximum temperature: ",phr.temp_max)
        print("  --semi-min: minimum semi-major axis: ",pse.semi_min)
        print("  --semi-max: minimum semi-major axis: ",pse.semi_max)
        print("  --ecc-min: minimum eccentricity: ",pse.ecc_min)
        print("  --ecc-max: maximum eccentricity: ",pse.ecc_max)
        print("  --interrupt-mode: no, base, bse: ",data.interrupt_mode)
        print("  --xcol: column index for x-axis: Unset")
        print("  --ycol: column index for x-axis: Unset")
        print("  --mcol: column index for mass, if not set and not PeTar output, assume equal mass: Unset")
        print("  --unit-length: set label of length unit for x, y and semi: no print")
        print("  --unit-time: set label of time unit: no print")
        print("  --skiprows: number of rows to escape when read snapshot: Unset")
        print("  --generate-binary: 0: no binary, 1: detect binary by using KDtree (slow), 2: read single and binary data generated by petar.data.process: ", data.generate_binary)
        print("  --plot-ncols: column number of panels: same as panels")
        print("  --plot-xsize: x size of panel: ",frame_xsize)
        print("  --plot-ysize: y size of panel: ",frame_ysize)
        print("  --cm-mode: plot origin position determination: density: density center; average: average of x,y: ", data.cm_mode)
        print("  --suppress-images: do not plot snapshot images (png files) and use matplotlib.animation instead of imageio, this cannot use multi-processing, much slower")
        print("  --format: video format, require imageio installed, for some formats (e.g. avi, mp4) may require ffmpeg and imageio-ffmpeg installed: ", plot_format)
        print("PS:: when xcol, ycol, skiprows are not provided, the snapshot files are assumed to be the output of PeTar")

    try:
        shortargs = 's:f:R:z:o:G:Hbh'
        longargs = ['help','lum-min=','lum-max=','temp-min=','temp-max=','semi-min=','semi-max=','ecc-min=','ecc-max=','interrupt-mode=','xcol=','ycol=','mcol=','unit-length=','unit-time=','skiprows=','generate-binary=','plot-ncols=','plot-xsize=','plot-ysize=','suppress-images','format=','cm-mode=']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        kwargs=dict()
        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-s'):
                fps = float(arg)
            elif opt in ('-R'):
                kwargs['boxsize'] = float(arg)
            elif opt in ('-z'):
                kwargs['framescale'] = float(arg)
                plot_zoom = True
            elif opt in ('-H'):
                plot_HRdiagram = True
            elif opt in ('-b'):
                plot_semi_ecc = True
            elif opt in ('-o'):
                output_file = arg
            elif opt in ('-G'):
                kwargs['G'] = float(arg)
            elif opt in ('-f'):
                fps = int(arg)
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
            elif opt in ('--interrupt-mode'):
                kwargs['interrupt_mode'] = arg
            elif opt in ('--cm-mode'):
                kwargs['cm_mode']= arg
            elif opt in ('--xcol'):
                kwargs['xcol'] = int(arg)
            elif opt in ('--ycol'):
                kwargs['ycol'] = int(arg)
            elif opt in ('--mcol'):
                kwargs['mcol'] = int(arg)
            elif opt in ('--unit-length'):
                kwargs['unit_length'] = arg
            elif opt in ('--unit-time'):
                kwargs['unit_time'] = arg
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
            elif opt in ('--suppress-images'):
                plot_images = False
            elif opt in ('--format'):
                plot_format = arg
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)


    filename = remainder[0]

    fl = open(filename,'r')
    file_list = fl.read()
    path_list = file_list.splitlines()

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
        for rank in range(n_cpu):
#            createImage(file_part[rank], frame_xsize, frame_ysize, ncol, plot_zoom, plot_HRdiagram, plot_semi_ecc, **kwargs)
            results[rank]=pool.apply_async(createImage, (file_part[rank], frame_xsize, frame_ysize, ncol, plot_zoom, plot_HRdiagram, plot_semi_ecc, ), kwargs)

        # Step 3: Don't forget to close
        pool.close()
        pool.join()

        png_list = [path_list[i]+'.png' for i in range(n_files)]
        create_movie(png_list, fps, output_file+'.'+plot_format)

    else:
        nplots = 1
        if (plot_zoom): nplots+=1
        if (plot_HRdiagram): nplots+=1
        if (plot_semi_ecc): nplots+=1

        fig, axe = initFig(frame_xsize, frame_ysize, ncol, nplots)
        plots=initPlot(axe, plot_zoom, plot_HRdiagram, plot_semi_ecc, **kwargs)

        def init():
            iaxe=0
            ptcls=[]
            if ('unit_time' in kwargs.keys()):
                unit_label = ' '+kwargs['unit_time']
                axe[0].set_title('T = ' + str(0) + unit_label)
            else:
                axe[0].set_title('T = %f' % 0)
                
            return ptcls
     
        def animate(k):
     
            file_path = path_list[k]
            print('process ',file_path)
            ptcls=plotOne(file_path, axe, plots,  **kwargs)
     
            return ptcls

        n_frame = len(path_list)
        anime = animation.FuncAnimation(fig, animate, init_func=init,
                                        frames=n_frame, interval=200, blit=True)
    
        anime.save(output_file+'.'+plot_format, fps=fps, extra_args=['-vcodec', 'libx264'])    

