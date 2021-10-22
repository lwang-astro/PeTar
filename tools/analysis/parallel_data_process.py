import numpy as np
import multiprocessing as mp
from .base import *
from .data import *
from .lagrangian import *
from .escaper import *
from .bse import *
from .functions import *
from .external import *
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
        The results, keys: lagr, core|core_read, esc_single, esc_binary, [bse]
        If read_flag = True, core_read is needed, else core is needed
        If interrupt_mode = bse, mobse, BSE based stellar evolution is needed
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
            interrupt_mode: PeTar interrupt mode: base, bse, mobse, none. If not provided, type is none 
            snapshot_format: snapshot format: ascii or binary (ascii)
    """
    lagr = result['lagr']
    esc_single  = result['esc_single']
    esc_binary  = result['esc_binary']

    m_frac = lagr.initargs['mass_fraction']
    G=1.0
    snapshot_format='ascii'
    r_bin=0.1
    average_mode='sphere'
    simple_binary=True
    external_mode='none'
    find_multiple=False

    if ('G' in kwargs.keys()): G=kwargs['G']
    if ('r_max_binary' in kwargs.keys()): r_bin=kwargs['r_max_binary']
    if ('average_mode' in kwargs.keys()): average_mode=kwargs['average_mode']
    if ('simple_binary' in kwargs.keys()): simple_binary=kwargs['simple_binary']
    if ('snapshot_format' in kwargs.keys()): snapshot_format=kwargs['snapshot_format']
    if ('external_mode' in kwargs.keys()): external_mode=kwargs['external_mode']
    if ('find_multiple' in kwargs.keys()): find_multiple=kwargs['find_multiple']

    header = PeTarDataHeader(file_path, **kwargs)
    
    if (not read_flag):
        start_time = time.time()

        core = result['core']
        #print('Loadfile')
        #snap=np.loadtxt(file_path, skiprows=1)
        #particle=Particle(snap, **kwargs)
        particle=Particle(**kwargs)
        if (snapshot_format=='ascii'): particle.loadtxt(file_path, skiprows=1)
        elif (snapshot_format=='binary'): 
            if (external_mode!='none'):
                particle.fromfile(file_path, offset=HEADER_OFFSET_WITH_CM)
            else:
                particle.fromfile(file_path, offset=HEADER_OFFSET)
        else: raise ValueError('Snapshot format unknown, should be binary or ascii, given', snapshot_format)
        time_profile['read'] += time.time() - start_time
        start_time = time.time()

        # find binary
        #print('Find pair')
        kdtree,single,binary=findPair(particle,G,r_bin,True,simple_binary)

        time_profile['find_pair'] += time.time() - start_time
        start_time = time.time()
    
        # get cm, density
        #print('Get density')
        cm_pos, cm_vel=core.calcDensityAndCenter(particle,kdtree)

        # add global offset
        if (external_mode!='none'): 
            core.pos[-1] += header.pos_offset
            core.vel[-1] += header.vel_offset
        #print('cm pos:',cm_pos,' vel:',cm_vel)
        time_profile['density'] += time.time() - start_time
        start_time = time.time()

        #print('Correct center')
        particle.correctCenter(cm_pos, cm_vel)

        # r2
        particle.calcR2()
        # rc

        #print('Core radius')
        rc = core.calcCoreRadius(particle)
        #print('rc: ',rc)

        core.addTime(header.time)
        core.size+=1

        n_frac=m_frac.size+1
        cm_vel=np.array([0,0,0]) # avoid kinetic energy jump 
        single.correctCenter(cm_pos, cm_vel)
        binary.correctCenter(cm_pos, cm_vel)
        time_profile['center_core'] += time.time() - start_time
        start_time = time.time()
        
        if (find_multiple): 
            single_t, binary_t, triple_t, quadruple_t = findMultiple(single,binary,G,r_bin,simple_binary)
            single_t.savetxt(file_path+'.single')
            binary_t.savetxt(file_path+'.binary')
            triple_t.savetxt(file_path+'.triple')
            quadruple_t.savetxt(file_path+'.quadruple')
        else:
            single.savetxt(file_path+'.single')
            binary.savetxt(file_path+'.binary')

        time_profile['save_data'] += time.time() - start_time

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
        rc = core.rc[core.time==header.time]

        time_profile['read'] += time.time() - start_time

    start_time = time.time()
    # calculate central external potential and subtract that from particle pot_ext
    pot_ext = 0
    if (external_mode != 'none'):
        pot_ext = calcCenterPotExt(particle, rc)
        single.pot -= pot_ext
        single.pot_ext -= pot_ext
        binary.p1.pot -= pot_ext
        binary.p2.pot -= pot_ext
        binary.p1.pot_ext -= pot_ext
        binary.p2.pot_ext -= pot_ext
    time_profile['calc_pot'] += time.time() - start_time
    start_time = time.time()

    if ('r_escape' in kwargs.keys()):
        rcut = kwargs['r_escape']
        es_cut = 0

        if (rcut == 'tidal'):
            if (external_mode!='none'): 
                tidal = result['tidal']
                r_gal = np.sqrt(np.sum(core.pos[-1]*core.pos[-1]));
                rcut = tidal.calcTidalSphere(header.time, particle.mass, particle.r2, pot_ext, r_gal, G);
            else:
                raise ValueError('Escape radius is set to tidal radius but the external mode is off')
        if ('e_escape' in kwargs.keys()): 
            es_cut = kwargs['e_escape']

        single = esc_single.findEscaper(header.time, single, rcut, es_cut)
        binary = esc_binary.findEscaper(header.time, binary, rcut, es_cut)
        time_profile['escaper'] += time.time() - start_time
        start_time = time.time()


    #print('Lagrangian radius')
    lagr.calcOneSnapshot(header.time, single, binary, rc, average_mode)

    time_profile['lagr'] += time.time() - start_time
    start_time = time.time()

    if (not 'r_escape' in kwargs.keys()):
        rhindex=np.where(m_frac==0.5)[0]
        rcut = calcRCutIsolate(lagr.all.r[-1,rhindex])
        esc_single.findEscaper(header.time, single, rcut)
        esc_binary.findEscaper(header.time, binary, rcut)

        time_profile['escaper'] += time.time() - start_time
        start_time = time.time()

    if ('bse_status' in result.keys()):
        bse = result['bse_status']
        bse.findEvents(header.time,single,binary)
        time_profile['bse'] += time.time() - start_time

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
    result['esc_single']=SingleEscaper(**kwargs)
    result['esc_binary']=BinaryEscaper(**kwargs)
    result['tidal']=Tidal(**kwargs)

    time_profile=dict()
    for key in ['read','find_pair','density','center_core','save_data','calc_pot','lagr','escaper','bse']:
        time_profile[key] = 0.0

    if (read_flag):
        core_filename=kwargs['filename_prefix']+'.core'
        result['core_read']=Core()
        result['core_read'].loadtxt(core_filename)
    else:
        result['core'] = Core()

    if ('interrupt_mode' in kwargs.keys()): 
        interrupt_mode=kwargs['interrupt_mode']
        if ('bse' in interrupt_mode):
            result['bse_status'] = BSEStatus()

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
            interrupt_mode: PeTar interrupt mode: base, bse, mobse, none. If not provided, type is none 
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
        for key in ['lagr','core','esc_single','esc_binary','bse_status','tidal']:
            if (key in resi.keys()):
                if (not key in result_all.keys()):
                    result_all[key]=[]
                result_all[key].append(resi[key])
        time_profile_all.append(result[i].get()[1])

    result_gether=dict()
    for key in result_all.keys():
        result_gether[key] = join(*result_all[key])

    for key in ['esc_single','esc_binary']:
        result_gether[key].removeDuplicate()

    time_profile=dict()

    for key in time_profile_all[0].keys():
        time_profile[key] = 0.0
        for i in range(n_cpu):
            time_profile[key] += time_profile_all[i][key]/n_cpu

    return result_gether, time_profile

