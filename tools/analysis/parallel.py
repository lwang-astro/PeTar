import numpy as np
import multiprocessing as mp
from .base import *
from .data import *
from .lagrangian import *
import time


def dataProcessOne(file_path, lagr, core, time_profile, m_frac=np.array([0.1,0.3,0.5,0.7,0.95]), G=1.0, r_bin=0.1, average_mode='sphere'):
    fp = open(file_path, 'r')
    header=fp.readline()
    file_id, n_glb, t = header.split()
    fp.close()
    
    start_time = time.time()
    #print('Loadfile')
    snap=np.loadtxt(file_path, skiprows=1)
    particle=Particle(snap)
    read_time = time.time()

    # find binary
    #print('Find pair')
    kdtree,single,binary=findPair(particle,G,r_bin,use_kdtree=True)
    find_pair_time = time.time()
    
    # get cm, density
    #print('Get density')
    cm_pos, cm_vel=core.calcDensityAndCenter(particle,kdtree)
    #print('cm pos:',cm_pos,' vel:',cm_vel)
    get_density_time = time.time()

    #print('Correct center')
    particle.correct_center(cm_pos, cm_vel)

    # r2
    particle.calc_r2()
    # rc

    #print('Core radius')
    rc = core.calcCoreRadius(particle)
    #print('rc: ',rc)

    core.addTime(float(t))
    core.size+=1

    n_frac=m_frac.size+1
    single.correct_center(cm_pos, cm_vel)
    binary.correct_center(cm_pos, cm_vel)
    center_and_r2_time = time.time()

    single.savetxt(file_path+'.single')
    binary.savetxt(file_path+'.binary')

    #print('Lagrangian radius')
    lagr.calcOneSnapshot(float(t), single, binary, m_frac, rc, average_mode)
    lagr_time = time.time()

    time_profile['read'] += read_time-start_time
    time_profile['find_pair'] += find_pair_time-read_time
    time_profile['density'] += get_density_time-find_pair_time
    time_profile['center_core'] += center_and_r2_time-get_density_time
    time_profile['lagr'] += lagr_time-center_and_r2_time

    return time_profile

def dataProcessList(file_list, m_frac=np.array([0.1,0.3,0.5,0.7,0.95]), G=1.0, r_bin=0.1, average_mode='sphere'):
    """ process lagragian calculation for a list of file snapshots
    file_list: file path list
    """
    lagr=LagrangianMultiple(m_frac.size+1)
    time_profile=dict()
    time_profile['read'] = 0.0
    time_profile['find_pair'] = 0.0
    time_profile['density'] = 0.0   
    time_profile['center_core'] = 0.0   
    time_profile['lagr'] = 0.0
    core=Core()
    for path in file_list:
        #print(' data:',path)
        dataProcessOne(path, lagr, core, time_profile, m_frac, G, r_bin, average_mode)
    for key, item in time_profile.items():
        item /= len(file_list)
    return lagr, core, time_profile


def parallelDataProcessList(file_list, m_frac=np.array([0.1,0.3,0.5,0.7,0.95]), G=1.0, r_bin=0.1, average_mode='sphere', n_cpu=int(0)):
    """ parellel process lagragian calculation for a list of file snapshots
    file_list: file path list
    """
    if (n_cpu==int(0)):
        n_cpu = mp.cpu_count()
        #print('n_cpu:',n_cpu)
    pool = mp.Pool(n_cpu)

    n_files=len(file_list)
    n_pieces = np.ones(n_cpu)*int(n_files/n_cpu)
    n_left = n_files%n_cpu
    n_pieces[:n_left]+=1
    n_offset=np.append([0],n_pieces.cumsum()).astype(int)
    #print('work_pieces',n_pieces)

    file_part = [file_list[n_offset[i]:n_offset[i+1]] for i in range(n_cpu)]

    result=[None]*n_cpu
    for rank in range(n_cpu):
        result[rank] = pool.apply_async(dataProcessList, args=(file_part[rank], m_frac, G, r_bin, average_mode))

    # Step 3: Don't forget to close
    pool.close()

    lagri=[]
    corei=[]
    time_profilei=[]
    for i in range(n_cpu):
        lagri.append(result[i].get()[0])
        corei.append(result[i].get()[1])
        time_profilei.append(result[i].get()[2])
    lagr = joinLagrangian(*lagri)
    core = join(*corei)
    time_profile=time_profilei[0]
    for key in time_profile.keys():
        for i in range(1,n_cpu):
            time_profile[key] += time_profilei[i][key]/n_cpu
    return lagr, core, time_profile

