#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp
import getopt
import imageio
import os
import pylab as pb
import petar

def create_movie(filenames, fps, output_file):
    with imageio.get_writer(output_file, fps=fps) as writer:
        for filename in filenames:
            writer.append_data(imageio.imread(filename))
    writer.close()

def createImage(index_list, xyscale, plot_format, log_flag, with_contour, only_xy, **kwargs):
    """
    index_list: file index list
    xyscale: plot x, y scale
    plot_format: plot format
    log_flag: if true, color is in logscale
    with_contour: if true, add contour plot
    only_xy: if true, only plot xy plane
    kwargs: vmin, vmax
    """
    
    if (only_xy):    
        fig, axes=plt.subplots(1,2,figsize=(4,3), gridspec_kw={'width_ratios': [20, 1]})
    else:
        fig, axes=plt.subplots(1,3,figsize=(8,3), gridspec_kw={'width_ratios': [20, 20, 1]})
    plt.subplots_adjust(wspace=0.4)

    for k in index_list:    
        print('process index ',k)

        if (only_xy):
            labels=['xy']
        else:
            labels=['xy','xz']

        for i,label in enumerate(labels):
            axes[i].clear()
            fname = label+str(k)

            if (read_binary):
                header = petar.ExternalPotMapHeader(fname, snapshot_format='binary')
                data= petar.ExternalPotMap()
                data.fromfile(fname, offset=petar.HEADER_EXTERNAL_POT_MAP_OFFSET)
            else:
                header = petar.ExternalPotMapHeader(fname, snapshot_format='ascii')
                data = petar.ExternalPotMap()
                data.loadtxt(fname, skiprows=1)

            pdata = data.plot(axes[i], header=header, plot_keys=list(label), with_contour=with_contour, log_flag=log_flag, **kwargs)
        
        if (with_contour):
            im = pdata[0]
        else:
            im = pdata

        axes[0].set_title('Time = %s' % header.time)
        if (only_xy):
            cbar = plt.colorbar(im, cax = axes[1]) 
        else:
            cbar = plt.colorbar(im, cax = axes[2]) 
        if log_flag: cbar.set_label(r'$\log{(-P)}$')
        else: cbar.set_label(r'-P')

        fig.savefig('pot'+str(k)+'.png', dpi=300, bbox_inches = "tight")

if __name__ == '__main__':

    fps = 30
    output_file = 'pot_movie'
    plot_format='mp4'
    n_cpu = 0
    log_flag = False
    read_binary = True
    with_contour = False
    only_xy = False
    kwargs=dict()

    def usage():
        print("A tool to generate a movie for the evolution of external potential")
        print("Need to use petar.external first to generate snapshots of potential map")
        print("Usage: petar.external.pot.movie [options] [petar.external parameter file]")
        print("Options:")
        print("  -h(--help): help")
        print("  -f [F]: output frame FPS: ",fps)
        print("  -o [S]: output movie filename: ",output_file)
        print("  -A    : read snapshot in ASCII format, default is BINARY format")
        print("  --with-contour: add contour plot, default is pure imshow")
        print("  --only-xy: only plot xy plane, default is both xy and xz planes")
        print("  --vmin [F]: (positive) potential minimum for color map, if not provided, use first snapshot for reference")
        print("  --vmax [F]: (positive) potential maximum for color map")
        print("  --log: color map is in logscale")
        print("  --n-cpu   [I]: number of CPU processors to use: all CPU cores")
        print("  --format  [S]: video format, require imageio installed, for some formats (e.g. avi, mp4) may require ffmpeg and imageio-ffmpeg installed: ", plot_format)
    try:
        shortargs = 'f:o:hA'
        longargs = ['help','vmax=','vmin=','format=','n-cpu=','only-xy','with-contour','log']
        opts, remainder = getopt.getopt(sys.argv[1:], shortargs, longargs)

        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-f'):
                fps = float(arg)
            elif opt in ('-o'):
                output_file = arg
            elif opt in ('-A'):
                read_binary = False
            elif opt in ('--vmin'):
                kwargs['vmin'] = float(arg)
            elif opt in ('--vmax'):
                kwargs['vmax'] = float(arg)
            elif opt in ('--format'):
                plot_format = arg
            elif opt in ('--with-contour'):
                with_contour = True
            elif opt in ('--only-xy'):
                only_xy = True
            elif opt in ('--n-cpu'):
                n_cpu = int(arg)
            elif opt in ('--log'):
                log_flag = True
            else:
                assert False, "unhandeld option"
            
    
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)


    fpar = remainder[0]

    fp = open(fpar, 'r')
    header = fp.readline()
    fp.close()
    t0, dt, nstep, dt_out, xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz = header.split()
    xmin=float(xmin)*1e-3
    xmax=float(xmax)*1e-3
    ymin=float(ymin)*1e-3
    ymax=float(ymax)*1e-3
    zmin=float(zmin)*1e-3
    zmax=float(zmax)*1e-3
    nstep= int(nstep)
    
    xyscale=[[xmin,xmax,ymin,ymax],[xmin,xmax,zmin,zmax]]

    if (not 'vmin' in kwargs.keys()) | (not 'vmax' in kwargs.keys()) :
        if (read_binary):
            data = petar.ExternalPotMap()
            data.fromfile('xy0', offset=petar.HEADER_EXTERNAL_POT_MAP_OFFSET)
            pot = data.pot
            pot_min = pot.min()
            pot_max = pot.max()
            data = petar.ExternalPotMap()
            data.fromfile('xz0', offset=petar.HEADER_EXTERNAL_POT_MAP_OFFSET)
            pot = data.pot
            pot_min = np.minimum(pot_min, pot.min())
            pot_max = np.maximum(pot_max, pot.max())
        else:
            # read the first snapshot to get the potential range
            # the first snapshot is xy0, the second is xz0
            # the potential is in the 10th column
            # the first row is the header
            data = petar.ExternalPotMap()
            data.loadtxt('xy0', skiprows=1)
            pot = data.pot
            pot_min = pot.min()
            pot_max = pot.max()
            data = petar.ExternalPotMap()
            data.loadtxt('xz0', skiprows=1)
            pot = data.pot
            pot_min = np.minimum(pot_min, pot.min())
            pot_max = np.maximum(pot_max, pot.max())
        if (not 'vmin' in kwargs.keys()):
            if (log_flag): kwargs['vmin'] = np.log10(-pot_max)
            else: kwargs['vmin'] = -pot_max
        if (not 'vmax' in kwargs.keys()):
            if (log_flag): kwargs['vmax'] = np.log10(-pot_min)
            else: kwargs['vmax'] = -pot_min
        

    if (n_cpu==int(0)):
        n_cpu = mp.cpu_count()
    pool = mp.Pool(n_cpu)
    
    n_files=int(nstep)+1
    file_list = range(0,nstep+1,int(float(dt_out)/float(dt)))
    n_pieces = np.ones(n_cpu)*int(len(file_list)/n_cpu)
    n_left = n_files%n_cpu
    n_pieces[:n_left]+=1
    n_offset=np.append([0],n_pieces.cumsum()).astype(int)
    
    file_part = [file_list[n_offset[i]:n_offset[i+1]] for i in range(n_cpu)]
    results=[None]*n_cpu
    if (n_cpu==1):
        createImage(file_part[0], xyscale, plot_format, log_flag, with_contour, only_xy, **kwargs)
    else:
        for rank in range(n_cpu):
            results[rank]=pool.apply_async(createImage, (file_part[rank], xyscale, plot_format, log_flag, with_contour, only_xy), kwargs)

    # Step 3: Don't forget to close
    pool.close()
    pool.join()

    png_list = ['pot'+str(file_list[i])+'.png' for i in range(len(file_list))]
    create_movie(png_list, fps, output_file+'.'+plot_format)
    
