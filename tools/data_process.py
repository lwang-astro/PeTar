#!/usr/bin/env python3

import numpy as np
import multiprocessing as mp
import sys
import petar
import getopt

if __name__ == '__main__':

    filename='dat.lst'
    lagr_filename='lagr.dat'
    core_filename='core.dat'
    esc_prefix='esc'
    average_mode='sphere'
    n_cpu=0

    def usage():
        print("A tool for processing a list of snapshot data to detect binaries, calculate Langragian radii and properties, get the density center and core radius")
        print("Usage: petar.data.process [options] data_filename")
        print("data_filename: A list of snapshot data path, each line for one snapshot")
        print("option:")
        print("  -h(--help): help")
        print("  -l(--lagr-filename): Lagrangian radii data filename (lagr.dat)")
        print("  -m(--mass-fraction): Lagrangian radii mass fraction (0.1,0.3,0.5,0.7,0.9)")
        print("  -G(--gravitational-constant): Gravitational constant (1.0)")
        print("  -r(--r-max-binary): maximum sepration for detecting binaries (0.1)")
        print("  -a(--average-mode): Lagrangian properity average mode: sphere: average from center to Lagragian radii; shell: average between two neighbor radii (sphere)")
        print("  -c(--core-filename): core data (time, density center, core radius) filename (core.dat)")
        print("  -e(--esc-filename): esc data filename prefix (esc)")
        print("  -n(--n-cpu): number of CPU threads for parallel processing (all threads)")

    try:
        shortargs = 'l:m:G:r:a:c:e:n:h'
        longargs = ['lagr-filename=','mass-fraction=','gravitational-constant=','r-max-binary=','average-mode=', 'core-filename=','esc-filename=','n-cpu=','help']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        kwargs=dict()
        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-l','--lagr-filename'):
                lagr_filename = arg
            elif opt in ('-m','--mass-fraction'):
                kwargs['mass_fraction'] = np.array([float(x) for x in arg.split(',')])
            elif opt in ('-G','--gravitational-constant'):
                kwargs['G'] = float(arg)
            elif opt in ('-r','--r-max-binary'):
                kwargs['r_max_binary'] = float(arg)
            elif opt in ('-a','--average-mode'):
                kwargs['average_mode'] = arg
            elif opt in ('-c','--core-filename'):
                core_filename = arg
            elif opt in ('-e','--esc-filename'):
                esc_filename = arg
            elif opt in ('-n','--n-cpu'):
                n_cpu = int(arg)
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError:
        print('getopt error!')
        usage()
        sys.exit(1)

    filename = remainder[0]

    for key, item in kwargs.items(): print(key,':',item)

    fl = open(filename,'r')
    file_list = fl.read()
    path_list = file_list.splitlines()
     
    lagr,core,esc,time_profile = petar.parallelDataProcessList(path_list,n_cpu=n_cpu,**kwargs)
     
    lagr.savetxt(lagr_filename)
    core.savetxt(core_filename)
    esc.single.savetxt(esc_prefix+'.single')
    esc.binary.savetxt(esc_prefix+'.binary')
     
    print ('Time profile:')
    for key, item in time_profile.items():
        print (key,item,)
     
    print ("Lagr data is saved in file: ",lagr_filename, ' core data is saved in file: ', core_filename, ' esc data are saved in file: ', esc_prefix,'.[single/binary]')

