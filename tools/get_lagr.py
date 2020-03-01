import numpy as np
import multiprocessing as mp
import sys
sys.path.append('/home/lwang/code/PeTar/tools')
import analysis

path='./'
fl = open(path+'dat.lst','r')
file_list = fl.read()
file_list = file_list.splitlines()
path_list = [path+file for file in file_list]

lagr = analysis.parllel_data_process_list(path_list)
#result = analysis.data_process_list(path_list)

lagr.savetxt('lagr.dat')
#print(lagr.time)
#print(lagr.all.r[:,0])
