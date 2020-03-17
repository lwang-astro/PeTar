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
    average_mode='sphere'
    n_cpu=0

    def usage():
        print("A tool for processing a list of snapshot data to detect binaries, calculate Langragian radii and properties, get the density center and core radius")
        print("Usage: petar.data.process [options] data_filename")
        print("data_filename: A list of snapshot data path, each line for one snapshot")
        print("option:")
        print("  -h(--help): help")
        print("  -l(--lagr-filename): Lagrangian radii data filename (lagr.dat)")
        print("  -m(--average-mode): Lagrangian properity average mode: sphere: average from center to Lagragian radii; shell: average between two neighbor radii (sphere)")
        print("  -c(--core-filename): core data (time, density center, core radius) filename (core.dat)")
        print("  -n(--n-cpu): number of CPU threads for parallel processing (all threads)")

    try:
        shortargs = 'N:M:R:t:o:s:h'
        longargs = ['lagr-filename=','average-mode=', 'core-filename=','help','n-cpu=']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-l','--lagr-filename'):
                lagr_filename = arg
            elif opt in ('-m','--average-mode'):
                average_mode = arg
            elif opt in ('-c','--core-filename'):
                core_filename = arg
            elif opt in ('-n','--n-cpu'):
                n_cpu = int(arg)
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError:
        print('getopt error!')
        usage()
        sys.exit(1)

    filename = remainder[0]

    fl = open(filename,'r')
    file_list = fl.read()
    file_list = file_list.splitlines()
    path_list = [file for file in file_list]
     
    lagr,core,time_profile = petar.parallel_data_process_list(path_list,n_cpu=n_cpu)
     
    lagr.savetxt(lagr_filename)
    core.savetxt(core_filename)
     
    print ('Time profile:')
    for key, item in time_profile.items():
        print (key,item,)
     
    print ("Lagr data is saved in file: ",lagr_filename, ' core data is saved in file: ', core_filename)

