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

def create_movie(filenames, fps, output_file):
    with imageio.get_writer(output_file, fps=fps) as writer:
        for filename in filenames:
            writer.append_data(imageio.imread(filename))
    writer.close()

def createImage(index_list, xyscale, plot_format, log_flag, dtype_data, dtype_header, offset_header, with_countour, only_xy, **kwargs):
    """
    index_list: file index list
    xyscale: plot x, y scale
    plot_format: plot format
    log_flag: if true, color is in logscale
    dtype_data: data types for reading binary files
    dtype_header: data types for reading header of binary files
    offset_header: offset for reading binary files
    with_countour: if true, add contour plot
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

        data=dict()
        if (only_xy):
            labels=['xy']
        else:
            labels=['xy','xz']
        xylabels=[['X','Y'],['X','Z']]

        for i in range(len(labels)):
            axes[i].clear()
            key=labels[i]
            data[key]=dict()
            fname = key+str(k)

            if (read_binary):
                header = np.fromfile(fname, dtype=dtype_header, count=1)
                time = header['time'][0]
                nx = int(header['nx'][0])
                ny = int(header['ny'][0])
                data[key] = np.fromfile(fname, dtype=dtype_data, offset=offset_header, count=nx*ny)
                data[key] = data[key].reshape(nx, ny)
            else:
                fp = open(fname, 'r')
                header = fp.readline()
                fp.close()
                time, nx, ny = header.split()
                x, y, z, ax, ay, az, pot = np.loadtxt(fname, unpack=True, usecols=(1,2,3,7,8,9,10),skiprows=1)

                nx=int(nx)
                ny=int(ny)
                data[key]['x']=x.reshape(nx,ny)
                data[key]['y']=y.reshape(nx,ny)
                data[key]['z']=z.reshape(nx,ny)
                data[key]['pot']=pot.reshape(nx,ny)

            count = np.log10(-data[key]['pot']) if log_flag else -data[key]['pot']
            if (with_countour):
                cset = axes[i].contour(count,linewidths=2,extent=xyscale[i], **kwargs)
                #axes[i].clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
            axes[i].set_xlabel(xylabels[i][0])
            axes[i].set_ylabel(xylabels[i][1])            

            im = axes[i].imshow(count,cmap=pb.cm.RdBu,
                                aspect=(xyscale[i][1]-xyscale[i][0])/(xyscale[i][3]-xyscale[i][2]),
                                interpolation='bilinear', origin='lower', 
                                extent=xyscale[i], **kwargs)

        axes[0].set_title('Time = %s' % time)
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
    with_countour = False
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
        print("  --with-countour: add contour plot, default is pure imshow")
        print("  --only-xy: only plot xy plane, default is both xy and xz planes")
        print("  --vmin [F]: (positive) potential minimum for color map, if not provided, use first snapshot for reference")
        print("  --vmax [F]: (positive) potential maximum for color map")
        print("  --log: color map is in logscale")
        print("  --n-cpu   [I]: number of CPU processors to use: all CPU cores")
        print("  --format  [S]: video format, require imageio installed, for some formats (e.g. avi, mp4) may require ffmpeg and imageio-ffmpeg installed: ", plot_format)
    try:
        shortargs = 'f:o:hA'
        longargs = ['help','vmax=','vmin=','format=','n-cpu=','only-xy','with-countour','log']
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
            elif opt in ('--with-countour'):
                with_countour = True
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

    dtype_data=np.dtype([('mass',np.float64),
                        ('x',np.float64), ('y',np.float64), ('z',np.float64),
                        ('vx',np.float64), ('vy',np.float64), ('vz',np.float64),
                        ('ax',np.float64), ('ay',np.float64), ('az',np.float64),
                        ('pot',np.float64),('den',np.float64)])
    dtype_header = np.dtype([('time', np.float64), ('nx', np.int32), ('ny', np.int32)])
    offset_header=8+4*2 # 8 bytes for time, 4 bytes for nx, 4 bytes for ny

    if (not 'vmin' in kwargs.keys()) | (not 'vmax' in kwargs.keys()) :
        if (read_binary):
            data = np.fromfile('xy0', dtype=dtype_data, offset=offset_header)
            pot = data['pot']
            pot_min = pot.min()
            pot_max = pot.max()
            data = np.fromfile('xz0', dtype=dtype_data, offset=offset_header)
            pot = data['pot']
            pot_min = np.minimum(pot_min, pot.min())
            pot_max = np.maximum(pot_max, pot.max())
        else:
            # read the first snapshot to get the potential range
            # the first snapshot is xy0, the second is xz0
            # the potential is in the 10th column
            # the first row is the header
            pot = np.loadtxt('xy0', unpack=True, usecols=(10),skiprows=1)
            pot_min = pot.min()
            pot_max = pot.max()
            pot = np.loadtxt('xz0', unpack=True, usecols=(10),skiprows=1)
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
    file_list = range(nstep+1)
    n_pieces = np.ones(n_cpu)*int(n_files/n_cpu)
    n_left = n_files%n_cpu
    n_pieces[:n_left]+=1
    n_offset=np.append([0],n_pieces.cumsum()).astype(int)
    
    file_part = [file_list[n_offset[i]:n_offset[i+1]] for i in range(n_cpu)]
    results=[None]*n_cpu
    if (n_cpu==1):
        createImage(file_part[0], xyscale, plot_format, log_flag, dtype_data, dtype_header, offset_header, with_countour, only_xy, **kwargs)
    else:
        for rank in range(n_cpu):
            results[rank]=pool.apply_async(createImage, (file_part[rank], xyscale, plot_format, log_flag, dtype_data, dtype_header, offset_header, with_countour, only_xy), kwargs)

    # Step 3: Don't forget to close
    pool.close()
    pool.join()

    png_list = ['pot'+str(file_list[i])+'.png' for i in range(n_files)]
    create_movie(png_list, fps, output_file+'.'+plot_format)
    
