import numpy as np
import multiprocessing as mp
from .base import *
from .data import *
from .lagrangian import *
from .escaper import *
import time
import os


def dataProcessOne(file_path, lagr, core, esc, time_profile, **kwargs): 
    m_frac = lagr.initargs['mass_fraction']
    G=1.0
    r_bin=0.1
    average_mode='sphere'
    read_flag=False

    if ('G' in kwargs.keys()): G=kwargs['G']
    if ('r_max_binary' in kwargs.keys()): r_bin=kwargs['r_max_binary']
    if ('average_mode' in kwargs.keys()): average_mode=kwargs['average_mode']
    if ('read_data' in kwargs.keys()): read_flag=kwargs['read_data']

    fp = open(file_path, 'r')
    header=fp.readline()
    file_id, n_glb, t = header.split()
    fp.close()
    
    if (not read_flag):
        start_time = time.time()
        #print('Loadfile')
        snap=np.loadtxt(file_path, skiprows=1)
        particle=Particle(snap, **kwargs)
        read_time = time.time()

        # find binary
        #print('Find pair')
        kdtree,single,binary=findPair(particle,G,r_bin,True)
        find_pair_time = time.time()
    
        # get cm, density
        #print('Get density')
        cm_pos, cm_vel=core.calcDensityAndCenter(particle,kdtree)
        #print('cm pos:',cm_pos,' vel:',cm_vel)
        get_density_time = time.time()

        #print('Correct center')
        particle.correctCenter(cm_pos, cm_vel)

        # r2
        particle.calcR2()
        # rc

        #print('Core radius')
        rc = core.calcCoreRadius(particle)
        #print('rc: ',rc)

        core.addTime(float(t))
        core.size+=1

        n_frac=m_frac.size+1
        cm_vel=np.array([0,0,0]) # avoid kinetic energy jump 
        single.correctCenter(cm_pos, cm_vel)
        binary.correctCenter(cm_pos, cm_vel)
        center_and_r2_time = time.time()

        single.savetxt(file_path+'.single')
        binary.savetxt(file_path+'.binary')
    else:
        start_time = time.time()
        
        single = Particle(**kwargs)
        p1 = Particle(**kwargs)
        p2 = Particle(**kwargs)
        binary = Binary(p1,p2)

        single.loadtxt(file_path+'.single')
        binary.loadtxt(file_path+'.binary')

        # read from core data
        rc = core.rc[core.time==float(t)]

        read_time = time.time()

        find_pair_time = read_time
        get_density_time = read_time
        center_and_r2_time = read_time
    
    if ('r_escape' in kwargs.keys()): 
        esc.rcut = kwargs['r_escape']
        esc.findEscaper(float(t),single,binary,G)
    
    #print('Lagrangian radius')
    lagr.calcOneSnapshot(float(t), single, binary, rc, average_mode)

    if (not 'r_escape' in kwargs.keys()):
        rhindex=np.where(m_frac==0.5)[0]
        esc.calcRCutIsolate(lagr.all.r[-1,rhindex])
        esc.findEscaper(float(t),single,binary,G)
    lagr_time = time.time()

    time_profile['read'] += read_time-start_time
    time_profile['find_pair'] += find_pair_time-read_time
    time_profile['density'] += get_density_time-find_pair_time
    time_profile['center_core'] += center_and_r2_time-get_density_time
    time_profile['lagr'] += lagr_time-center_and_r2_time

#    return time_profile

def dataProcessList(file_list, **kwargs):
    """ process lagragian calculation for a list of file snapshots
    file_list: file path list
    """
    lagr=LagrangianMultiple(**kwargs)
    time_profile=dict()
    time_profile['read'] = 0.0
    time_profile['find_pair'] = 0.0
    time_profile['density'] = 0.0   
    time_profile['center_core'] = 0.0   
    time_profile['lagr'] = 0.0
    core=Core()
    if ('read_data' in kwargs.keys()):
        if (kwargs['read_data'] == True):
            core_filename=kwargs['filename_prefix']+'.core'
            core.loadtxt(core_filename)
    esc=Escaper()
    for path in file_list:
        #print(' data:',path)
        dataProcessOne(path, lagr, core, esc, time_profile, **kwargs)
    if (len(file_list)>0):
        for key, item in time_profile.items():
            item /= len(file_list)
    return lagr, core, esc, time_profile


def parallelDataProcessList(file_list, n_cpu=int(0), **kwargs):
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
        result[rank] = pool.apply_async(dataProcessList, (file_part[rank],), kwargs)

    # Step 3: Don't forget to close
    pool.close()
    pool.join()

    lagri=[]
    corei=[]
    esci=[]
    time_profilei=[]
    for i in range(n_cpu):
        lagri.append(result[i].get()[0])
        corei.append(result[i].get()[1])
        esci.append(result[i].get()[2])
        time_profilei.append(result[i].get()[3])
    lagr = join(*lagri)
    core = join(*corei)
    esc  = joinEscaper(*esci)
    time_profile=time_profilei[0]
    for key in time_profile.keys():
        for i in range(1,n_cpu):
            time_profile[key] += time_profilei[i][key]/n_cpu
    return lagr, core, esc, time_profile

