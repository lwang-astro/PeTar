#!/usr/bin/env python3
import sys
import collections
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation, rc
from matplotlib import cm
from IPython.display import HTML
import getopt
import petar

plt.style.use('dark_background')

if __name__ == '__main__':

    filename='dat.lst'
    dt = 60
    R = 1
    framescale=20
    lum_min = 5e-4
    lum_max = 1e5
    temp_min = 2000
    temp_max = 50000
    output_file = 'movie'
    skiprows = 0
    xcol = -1
    ycol = -1
    mcol = -1
    interrupt_mode='bse'
    plot_HRdiagram=True

    def usage():
        print("A tool for processing a list of snapshot data to detect binaries, calculate Langragian radii and properties, get the density center and core radius")
        print("Usage: petar.motion.movie [options] data_list_filename")
        print("data_list_filename: A list of snapshot data path, each line for one snapshot")
        print("option: default values are shown at last")
        print("  -h(--help): help")
        print("  -s: output frame speed: ",dt)
        print("  -R: movie box length: ",R)
        print("  -z: plot zoom frame instead of HR digram, framescale to zoom in: ",framescale)
        print("  -o: output movie filename: ",output_file)
        print("  --lum-min: mimimum lumonisity: ",lum_min)
        print("  --lum-max: maximum lumonisity: ",lum_max)
        print("  --temp-min: mimimum temperature: ",temp_min)
        print("  --temp-max: maximum temperature: ",temp_max)
        print("  --interrupt-mode: no, base, bse: ",interrupt_mode)
        print("  --xcol: column index for x-axis: Unset")
        print("  --ycol: column index for x-axis: Unset")
        print("  --mcol: column index for mass, if not set and not PeTar output, assume equal mass: Unset")
        print("  --skiprows: number of rows to escape when read snapshot: Unset")
        print("PS:: when xcol, ycol, skiprows are not provided, the snapshot files are assumed to be the output of PeTar")

    try:
        shortargs = 's:R:z:o:h'
        longargs = ['help','lum-min=','lum-max=','temp-min=','temp-max=','interrupt-mode=','xcol=','ycol=','mcol=','skiprows=']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-s'):
                dt = float(arg)
            elif opt in ('-R'):
                R = float(arg)
            elif opt in ('-z'):
                framescale = float(arg)
                plot_HRdiagram = False
            elif opt in ('-o'):
                output_file = arg
            elif opt in ('--lum-min'):
                lum_min = float(arg)
            elif opt in ('--lum-max'):
                lum_max = float(arg)
            elif opt in ('--temp-min'):
                temp_min = float(arg)
            elif opt in ('--temp-max'):
                temp_max = float(arg)
            elif opt in ('--interrupt-mode'):
                interrupt_mode = arg
                if (interrupt_mode!='bse'): plot_HRdiagram = False
            elif opt in ('--xcol'):
                xcol = int(arg)
            elif opt in ('--ycol'):
                ycol = int(arg)
            elif opt in ('--mcol'):
                mcol = int(arg)
            elif opt in ('--skiprows'):
                skiprows = int(arg)
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)


    ncol= 2
    xsize=16
    fig, axe = plt.subplots(1,ncol,figsize=(xsize, 8))

    nlayer_cross=5
    nlayer_point=10
    nlayer = nlayer_cross + nlayer_point
    alpha_amplifier = 2.5
    alphascale=np.linspace(1,nlayer,nlayer)*2.0/(nlayer*(nlayer+1))*alpha_amplifier
    print('Alpha layer sequence:',alphascale,' sum:',alphascale.sum())
    sizescale=np.logspace(0,3,nlayer)[::-1]
    print('Size layer sequence:',sizescale)
    ptcls=[]

    filename = remainder[0]

    fl = open(filename,'r')
    file_list = fl.read()
    path_list = file_list.splitlines()
    pos_cm=[0,0]

    def init():
        boxsize=2.0*R
        axe[0].set_xlim(-boxsize,boxsize)
        axe[0].set_ylim(-boxsize,boxsize)
        axe[0].set_aspect(1.0)
        axe[0].set_xlabel('x')
        axe[0].set_ylabel('y')

        if plot_HRdiagram:
            axe[1].set_xlim(30e3,1e3);
            axe[1].set_ylim(1e-5,1e5);
            axe[1].set_yscale('log');
            axe[1].set_xscale('log');
            axe[1].set_xlabel(r'$T_{eff}[K]$')
            axe[1].set_ylabel(r'$L[L_\odot]$')
        else:
            boxsize /= framescale
            axe[1].set_xlim(-boxsize,boxsize)
            axe[1].set_ylim(-boxsize,boxsize)
            axe[1].set_aspect(1.0)
            axe[1].set_xlabel('x')
            axe[1].set_ylabel('y')
     
        for i in range(nlayer_cross):
            pt =axe[0].scatter([],[],marker='+',alpha=alphascale[i],edgecolors='none')
            ptcls.append(pt)
        for i in range(nlayer_point):
            pt =axe[0].scatter([],[],alpha=alphascale[i],edgecolors='none')
            ptcls.append(pt)
     
        for i in range(nlayer_cross):
            pt =axe[1].scatter([],[],marker='+',alpha=alphascale[i],edgecolors='none')
            ptcls.append(pt)
        for i in range(nlayer_point):
            pt =axe[1].scatter([],[],alpha=alphascale[i],edgecolors='none')
            ptcls.append(pt)
     
        axe[0].set_title('T = %f' % 0)

        return ptcls

    def animate(k):

        file_path = path_list[k]
        print('process ',file_path)
        data=dict()
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

            particles=petar.Particle(use_BSE=(interrupt_mode=='bse'))
            particles.loadtxt(file_path,skiprows=1)
            data['x'] = particles.pos[:,0]
            data['y'] = particles.pos[:,1]
            data['mass'] =particles.mass
            if (interrupt_mode=='bse'):
                data['lum'] = particles.s_lum
                data['rad'] = particles.s_rad
                data['type']= particles.s_type

        x = data['x']
        y = data['y']
        mass = data['mass']
        boxsize = 2.0*R
        sel=(x-pos_cm[0]>-boxsize) & (x-pos_cm[1]<boxsize) & (y-pos_cm[0]>-boxsize) & (y-pos_cm[1]<boxsize)
        nsel=sel.sum()
        pos_cm[0] = (x[sel]).sum()/float(nsel)
        pos_cm[1] = (y[sel]).sum()/float(nsel)
        x = x - pos_cm[0]
        y = y - pos_cm[1]
        axe[0].set_title('T = '+data['t'])

        colors='w'
        if ('lum' in data.keys()):
            luminosity = data['lum']
            radius = data['rad']
            temperature_eff = 5778*(luminosity/(radius*radius))**0.25
            LogTeff=(np.log10(temperature_eff)-np.log10(temp_min))/(np.log10(temp_max)-np.log10(temp_min))
            colors=cm.rainbow(1.0-LogTeff)
        #sizes = (np.log10(luminosity)-np.log10(lum_min)+1)*100
        #sizes = luminosity

        for i in range(nlayer):
            sizes = mass*sizescale[i]
            #sizes = (np.log10(luminosity)-np.log10(lum_min)+1)
            #sizes = luminosity
            ptcls[i].set_offsets(np.array([x,y]).transpose())
            ptcls[i].set_sizes(sizes)
            ptcls[i].set_color(colors)

        if plot_HRdiagram:
            types = data['type']
            colors = cm.rainbow(types/15.0)
            ptcls[2*nlayer-1].set_offsets(np.array([temperature_eff, luminosity]).transpose())
            ptcls[2*nlayer-1].set_color(colors)
        else:
            for i in range(nlayer):
                sizes = mass*framescale*sizescale[i]
                ptcls[nlayer+i].set_offsets(np.array([x,y]).transpose())
                ptcls[nlayer+i].set_sizes(sizes)
                ptcls[nlayer+i].set_color(colors)

        return ptcls

    n_frame = len(path_list)

    anime = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=n_frame, interval=dt, blit=True)
    
    #gravity.stop()
    #stellar_evolution.stop()

    anime.save(output_file+'.mp4', fps=30, extra_args=['-vcodec', 'libx264'])    

