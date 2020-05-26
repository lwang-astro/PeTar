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

    def usage():
        print("A tool for processing a list of snapshot data to detect binaries, calculate Langragian radii and properties, get the density center and core radius")
        print("Usage: petar.data.process [options] data_filename")
        print("data_filename: A list of snapshot data path, each line for one snapshot")
        print("option:")
        print("  -h(--help): help")
        print("  -s: output frame speed: ",dt)
        print("  -R: movie box length: ",R)
        print("  -z: framescale to zoom in: ",framescale)
        print("  -o: output movie filename: ",output_file)
        print("  --lum-min: mimimum lumonisity: ",lum_min)
        print("  --lum-max: maximum lumonisity: ",lum_max)
        print("  --temp-min: mimimum temperature: ",temp_min)
        print("  --temp-max: maximum temperature: ",temp_max)

    try:
        shortargs = 's:R:z:o:h'
        longargs = ['help','lum-min=','lum-max=','temp-min=','temp-max=']
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
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError:
        print('getopt error!')
        usage()
        sys.exit(1)


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

        boxsize /= framescale
        axe[1].set_xlim(-boxsize,boxsize)
        axe[1].set_ylim(-boxsize,boxsize)
        axe[1].set_aspect(1.0)
        axe[1].set_xlabel('x[pc]')
        axe[1].set_ylabel('y[pc]')
     
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
        fp = open(file_path, 'r')
        header=fp.readline()
        file_id, n_glb, t = header.split()
        fp.close()

        time = float(t)
        particles=petar.Particle(use_BSE=True)
        particles.loadtxt(file_path,skiprows=1)
        x = particles.pos[:,0]
        y = particles.pos[:,1]
        mass =particles.mass
        boxsize = 2.0*R
        sel=(x-pos_cm[0]>-boxsize) & (x-pos_cm[1]<boxsize) & (y-pos_cm[0]>-boxsize) & (y-pos_cm[1]<boxsize)
        nsel=sel.sum()
        pos_cm[0] = (x[sel]).sum()/float(nsel)
        pos_cm[1] = (y[sel]).sum()/float(nsel)
        x = x - pos_cm[0]
        y = y - pos_cm[1]
        axe[0].set_title('T = %f' % time)

        luminosity = particles.s_lum
        radius = particles.s_rad
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

