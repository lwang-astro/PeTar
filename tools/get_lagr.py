import numpy as np
import multiprocessing as mp
import sys
import petar

filename='dat.lst'
lagr_filename='lagr.dat'
n_cpu=0

narg=len(sys.argv)-1
if (narg>=1):
    filename = sys.argv[1]
    if (filename=='-h'): 
        print("arguments: snapshot path list [dat.lst], Lagr output dat file [lagr.dat], number of parallel cores [Use all threads]")
        quit()
if (narg>=2):    lagr_filename=sys.argv[2]
if (narg>=3):    n_cpu=int(sys.argv[3])

fl = open(filename,'r')
file_list = fl.read()
file_list = file_list.splitlines()
path_list = [file for file in file_list]

lagr,time_profile = petar.parallel_data_process_list(path_list,n_cpu=n_cpu)

lagr.savetxt(lagr_filename)

print ('Time profile:')
for key, item in time_profile.items():
    print (key,item,)

print ("Lagr data is saved in file: ",lagr_filename)

