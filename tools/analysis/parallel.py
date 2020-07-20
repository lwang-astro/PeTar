import numpy as np
import multiprocessing as mp
from .base import *
from .data import *
from .lagrangian import *
from .escaper import *
from .bse import *
import time
import os


def dataProcessOne(file_path, result, time_profile, read_flag, **kwargs): 
    """Process one snapshot.

    Find binaries of one snapshot, calculate Lagrangian radii, find the system core and find escapers.
    
    Parameters
    ----------
    file_path: list
        The pathes of snapshots
    result: dict
        The results, keys: lagr, core, esc, [bse]
    time_profile: dict
        The CPU (wallclock) time for each parts of calculations
    read_flag: bool
        If true, read single, binary snapshots and core data instead of calculating them
    kwargs: dict ()
        Keywords arguments:
            G: gravitational constant (1.0)
            r_max_binary: maximum separation to detect binaries (0.1)
            average_mode: mode in calculating lagrangian radii (sphere)
            mass_fraction: an 1D numpy.ndarray to indicate the mass fractions to calculate lagrangian radii.
                               Default is np.array([0.1, 0.3, 0.5, 0.7, 0.9])
            interrupt_mode: PeTar interrupt mode: base, bse, none. If not provided, type is none 
    """
    lagr = result['lagr']
    esc  = result['esc']

    m_frac = lagr.initargs['mass_fraction']
    G=1.0
    r_bin=0.1
    average_mode='sphere'

    if ('G' in kwargs.keys()): G=kwargs['G']
    if ('r_max_binary' in kwargs.keys()): r_bin=kwargs['r_max_binary']
    if ('average_mode' in kwargs.keys()): average_mode=kwargs['average_mode']

    fp = open(file_path, 'r')
    header=fp.readline()
    file_id, n_glb, t = header.split()
    fp.close()
    
    if (not read_flag):
        start_time = time.time()

        core = result['core']
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
        
        core = result['core_read']

        single = Particle(**kwargs)
        p1 = Particle(**kwargs)
        p2 = Particle(**kwargs)
        binary = Binary(p1,p2)

        if os.path.getsize(file_path+'.single')>0:
            single.loadtxt(file_path+'.single')
        if os.path.getsize(file_path+'.binary')>0:
            binary.loadtxt(file_path+'.binary')

        # read from core data
        rc = core.rc[core.time==float(t)]

        read_time = time.time()

        find_pair_time = read_time
        get_density_time = read_time
        center_and_r2_time = read_time
    
    
    if ('esc' in result.keys()) & ('r_escape' in kwargs.keys()) : 
        esc.rcut = kwargs['r_escape']
        esc.findEscaper(float(t),single,binary,G)
    
    #print('Lagrangian radius')
    lagr.calcOneSnapshot(float(t), single, binary, rc, average_mode)

    if ('esc' in result.keys()) & (not 'r_escape' in kwargs.keys()):
        rhindex=np.where(m_frac==0.5)[0]
        esc.calcRCutIsolate(lagr.all.r[-1,rhindex])
        esc.findEscaper(float(t),single,binary,G)
    lagr_time = time.time()

    time_profile['read'] += read_time-start_time
    time_profile['find_pair'] += find_pair_time-read_time
    time_profile['density'] += get_density_time-find_pair_time
    time_profile['center_core'] += center_and_r2_time-get_density_time
    time_profile['lagr'] += lagr_time-center_and_r2_time

    if ('bse' in result.keys()):
        bse = result['bse']
        bse.findEvents(float(t),single,binary)
        bse_time = time.time()
        time_profile['bse'] += bse_time - lagr_time

#    return time_profile

def dataProcessList(file_list, read_flag, **kwargs):
    """ process lagragian calculation for a list of file snapshots

    Parameters
    ----------
    file_list: list
        file path list
    read_flag: bool
        indicate whether to read single, binary and core data instead of calculating 
    kwargs: dict
        keyword arguments:
            filename_prefix: filename prefix for output data (data)
            G: gravitational constant (1.0)
            r_max_binary: maximum separation to detect binaries (0.1)
            average_mode: mode in calculating lagrangian radii (sphere)
            mass_fraction: an 1D numpy.ndarray to indicate the mass fractions to calculate lagrangian radii.
                               Default is np.array([0.1, 0.3, 0.5, 0.7, 0.9])
            interrupt_mode: PeTar interrupt mode: base, bse, none. If not provided, type is none 
    """
    result = dict()
    result['lagr']=LagrangianMultiple(**kwargs)
    result['esc']=Escaper()

    time_profile=dict()
    time_profile['read'] = 0.0
    time_profile['find_pair'] = 0.0
    time_profile['density'] = 0.0   
    time_profile['center_core'] = 0.0   
    time_profile['lagr'] = 0.0

    if (read_flag):
        core_filename=kwargs['filename_prefix']+'.core'
        result['core_read']=Core()
        result['core_read'].loadtxt(core_filename)
    else:
        result['core'] = Core()

    if ('interrupt_mode' in kwargs.keys()): 
        interrupt_mode=kwargs['interrupt_mode']
        if (interrupt_mode=='bse'):
            result['bse'] = BSEEvent()
            time_profile['bse'] = 0.0

    for path in file_list:
        #print(' data:',path)
        dataProcessOne(path, result, time_profile, read_flag, **kwargs)

    if (len(file_list)>0):
        for key, item in time_profile.items():
            item /= len(file_list)
    return result, time_profile


def parallelDataProcessList(file_list, n_cpu=int(0), read_flag=False, **kwargs):
    """ parellel process lagragian calculation for a list of file snapshots

    Parameters
    ----------
    file_list: list
        file path list
    n_cpu: int
        number of CPU cores to run parallelly
    read_flag: bool
        indicate whether to read single, binary and core instead of calculating 
    kwargs: dict
        keyword arguments:
            filename_prefix: filename prefix for output data (data)
            G: gravitational constant (1.0)
            r_max_binary: maximum separation to detect binaries (0.1)
            average_mode: mode in calculating lagrangian radii (sphere)
            mass_fraction: an 1D numpy.ndarray to indicate the mass fractions to calculate lagrangian radii.
                               Default is np.array([0.1, 0.3, 0.5, 0.7, 0.9])
            interrupt_mode: PeTar interrupt mode: base, bse, none. If not provided, type is none 
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
        result[rank] = pool.apply_async(dataProcessList, (file_part[rank], read_flag,), kwargs)

    # Step 3: Don't forget to close
    pool.close()
    pool.join()

    time_profile_all=[]
    result_all=dict()
    for i in range(n_cpu):
        resi = result[i].get()[0]
        for key in ['lagr','core','esc','bse']:
            if (key in resi.keys()):
                if (not key in result_all.keys()):
                    result_all[key]=[]
                result_all[key].append(resi[key])
        time_profile_all.append(result[i].get()[1])

    result_gether=dict()
    for key in result_all.keys():
        if (key != 'esc'):
            result_gether[key] = join(*result_all[key])
        else:
            result_gether[key] = joinEscaper(*result_all[key])

    time_profile=dict()

    for key in time_profile_all[0].keys():
        time_profile[key] = 0.0
        for i in range(n_cpu):
            time_profile[key] += time_profile_all[i][key]/n_cpu

    return result_gether, time_profile

