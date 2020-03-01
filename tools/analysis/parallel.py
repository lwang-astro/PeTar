import numpy as np
import multiprocessing as mp
from .base import *
from .data import *
from .lagrangian import *

def data_process_one(file_path, lagr, m_frac=np.array([0.1,0.3,0.5,0.7,0.95]), G=1.0, r_bin=0.1):
    fp = open(file_path, 'r')
    header=fp.readline()
    file_id, n_glb, time = header.split()
    
    #print('Loadfile')
    snap=np.loadtxt(file_path, skiprows=1)
    particle=Particle(snap)
    # find binary
    #print('Find pair')
    kdtree,single,binary=findPair(particle,G,r_bin)
    
    # get cm, density
    #print('Get density')
    cm_pos, cm_vel=density_six_nb(particle,kdtree)
    #print('cm pos:',cm_pos,' vel:',cm_vel)
    #print('Correct center')
    particle.correct_center(cm_pos, cm_vel)

    # r2
    particle.calc_r2()
    # rc

    #print('Core radius')
    rc = core_radius(particle)
    #print('rc: ',rc)

    n_frac=m_frac.size+1
    single.correct_center(cm_pos, cm_vel)
    binary.correct_center(cm_pos, cm_vel)

    #print('Lagrangian radius')
    lagr.calcOneSnapshot(float(time), single,binary,m_frac,rc)


def data_process_list(file_list, m_frac=np.array([0.1,0.3,0.5,0.7,0.95]), G=1.0, r_bin=0.1):
    """ process lagragian calculation for a list of file snapshots
    file_list: file path list
    """
    lagr=LagrangianMultiple(m_frac.size+1)
    for path in file_list:
        print(' data:',path)
        data_process_one(path, lagr, m_frac, G, r_bin)
    return lagr


def parllel_data_process_list(file_list, m_frac=np.array([0.1,0.3,0.5,0.7,0.95]), G=1.0, r_bin=0.1):
    """ parellel process lagragian calculation for a list of file snapshots
    file_list: file path list
    """
    n_cpu = mp.cpu_count()
    print('n_cpu:',n_cpu)
    pool = mp.Pool(n_cpu)

    n_files=len(file_list)
    n_pieces = np.ones(n_cpu)*int(n_files/n_cpu)
    n_left = n_files%n_cpu
    n_pieces[:n_left]+=1
    n_offset=np.append([0],n_pieces.cumsum()).astype(int)
    print('work_pieces',n_pieces)

    file_part = [file_list[n_offset[i]:n_offset[i+1]] for i in range(n_cpu)]

    result=[None]*n_cpu
    for rank in range(n_cpu):
        result[rank] = pool.apply_async(data_process_list, args=(file_part[rank], m_frac, G, r_bin))

    # Step 3: Don't forget to close
    pool.close()

    lagri=[]
    for i in range(n_cpu):
        lagri.append(result[i].get())
    lagr = joinLagrangian(*lagri)
    return lagr

