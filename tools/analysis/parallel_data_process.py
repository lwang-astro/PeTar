import numpy as np
import multiprocessing as mp
from sdar.base import *
from sdar.functions import *
from .data import *
from .lagrangian import *
from .escaper import *
from .bse import *
from .external import *
import time
import os


def _saveRealtimeData(data, filename, output_format):
    """Save updated analysis data immediately after one snapshot is processed."""
    if (data.size == 0):
        return

    if (output_format == 'ascii'):
        with open(filename, 'a') as f:
            data[data.size-1:data.size].savetxt(f)
    elif (output_format == 'binary'):
        with open(filename, 'ab') as f:
            data[data.size-1:data.size].tofile(f)
    elif (output_format == 'npy'):
        with open(filename, 'wb') as f:
            data.save(f)
    else:
        raise ValueError('Output format %s is not supported, should be ascii, binary or npy' % output_format)


def _prepareRealtimeSaveFiles(filename_prefix, save_keys, output_format, write_mode):
    """Prepare output files for real-time saving."""
    if (write_mode == 'a'):
        return

    file_mode = 'w'
    if (output_format in ['binary', 'npy']):
        file_mode = 'wb'

    for key in save_keys:
        with open(filename_prefix + '.' + key, file_mode):
            pass


def _getParallelRealtimePrefix(filename_prefix, rank):
    """Generate per-worker prefix for parallel real-time saving."""
    return filename_prefix + '.parallel.' + str(rank)


def dataProcessOne(file_path, result, time_profile, 
                   read_flag=False, r_max_binary=0.1, 
                   average_mode='sphere', simple_binary=True, 
                   snapshot_format='binary', output_format='binary', 
                   find_multiple=False, r_escape=None, e_escape=None, **kwargs): 
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
    read_flag: bool (False)
        If true, read single, binary snapshots and core data instead of calculating them
    r_max_binary: float (0.1)
        maximum separation to detect binaries (0.1)
    average_mode: str (sphere)
        mode in calculating lagrangian radii (sphere)
    simple_binary: bool (True)
        whether to use simple binary detection (True)
    snapshot_format: str (binary)
        input snapshot format: ascii or binary (binary)
    output_format: str (binary)
        output data format: ascii, binary, npy (binary)
    find_multiple: bool (False)
        whether to find multiple systems (False)
    r_escape: float or string (None)
        escape radius, if set, escaper will be detected. 
        If set to 'tidal', the tidal radius will be calculated and used as escape radius
    e_escape: float or string (None)
        escape energy, if set, escaper will be detected. 
        If set to 'bound_noext', the escapers are defined as those with positive energy 
        when the external potential is not considered, 
        thus the external potential will be subtracted when calculating energy
    kwargs: dict ()
        Keywords arguments:
            G: float (1.0)
               gravitational constant (1.0)
            interrupt_mode: string (none)
               PeTar interrupt mode (set in configure): base, bse, mobse, none
               This option indicates whether columns of stellar evolution exist
            external_mode: string (none)
               PeTar external mode (set in configure): galpy, agama, none 
               This option indicates whether the column of externa potential exist
            use_mpfrc: bool (False)
               If true, add three columns of pos_high indicating the high-precision parts of position
            collect_sp_acc: bool (False)
               If true, the superparticle acceleration is collected and the column acc_sp exists
            float_type: type (np.float64)
                floating point data type
    """
    lagr = result['lagr']
    esc_single  = result['esc_single']
    esc_binary  = result['esc_binary']

    m_frac = lagr.initargs['mass_fraction']
    G=1.0
    external_mode='none'
    m_ext=None

    if ('G' in kwargs.keys()): G=kwargs['G']
    if ('external_mode' in kwargs.keys()): external_mode=kwargs['external_mode']
    if ('m_ext' in result.keys()): m_ext=result['m_ext']

    start_time = time.time()
    header = PeTarDataHeader(file_path, snapshot_format=snapshot_format, **kwargs)
    particle=Particle(**kwargs)
    if (snapshot_format=='ascii'): particle.loadtxt(file_path, skiprows=1)
    elif (snapshot_format=='binary'): 
        if (external_mode!='none'):
            particle.fromfile(file_path, offset=HEADER_OFFSET_WITH_CM)
        else:
            particle.fromfile(file_path, offset=HEADER_OFFSET)
    else: raise ValueError('Snapshot format unknown, should be binary or ascii, given', snapshot_format)
    time_profile['read'] += time.time() - start_time

    detect_single_binary = True # whether to detect single and binary files

    # read from core data
    core = result['core']
    if (read_flag):
        start_time = time.time()
        core_read = result['core_read']    
        tsel = (core_read.time==header.time)
        if (output_format=='npy'):
            single_file_available = (os.path.getsize(file_path+'.single.npy')>0)
            binary_file_available = (os.path.getsize(file_path+'.binary.npy')>0)
        else:
            single_file_available = (os.path.getsize(file_path+'.single')>0)
            binary_file_available = (os.path.getsize(file_path+'.binary')>0)

        # check whether the core data for given time is available
        # if both core data and single/binary files available, set detect_single_binary to False
        # when find_multiple is used, always set detect_single_binary True because the saved single/binary files miss triple/quadruple components
        if (tsel.sum()>0) & (single_file_available | binary_file_available):
            detect_single_binary = False
        time_profile['read'] += time.time() - start_time

    # detect single/binary/multiple and calculate core data
    if (detect_single_binary):

        #print('Loadfile')
        #snap=np.loadtxt(file_path, skiprows=1)
        #particle=Particle(snap, **kwargs)
        start_time = time.time()

        # find binary
        #print('Find pair')
        kdtree,single,binary=findPair(particle, G, r_max_binary, use_kdtree=True, simple_binary=simple_binary)

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
        #cm_vel=np.array([0,0,0]) # avoid kinetic energy jump 
        single.correctCenter(cm_pos, cm_vel)
        binary.correctCenter(cm_pos, cm_vel)
        time_profile['center_core'] += time.time() - start_time
        start_time = time.time()
        
        if (find_multiple): 
            single_t, binary_t, triple_t, quadruple_t = findMultiple(single,binary,G,r_max_binary,simple_binary)
            if (output_format=='ascii'):
                single_t.savetxt(file_path+'.single')
                binary_t.savetxt(file_path+'.binary')
                triple_t.savetxt(file_path+'.triple')
                quadruple_t.savetxt(file_path+'.quadruple')
            elif (output_format=='binary'):
                single_t.tofile(file_path+'.single')
                binary_t.tofile(file_path+'.binary')
                triple_t.tofile(file_path+'.triple')
                quadruple_t.tofile(file_path+'.quadruple')
            elif (output_format=='npy'):
                single_t.save(file_path+'.single')
                binary_t.save(file_path+'.binary')
                triple_t.save(file_path+'.triple')
                quadruple_t.save(file_path+'.quadruple')
            else:
                raise ValueError('Output format %s is not supported, should be ascii, binary or npy' % output_format)
        else:
            if (output_format=='ascii'):
                single.savetxt(file_path+'.single')
                binary.savetxt(file_path+'.binary')
            elif (output_format=='binary'):
                single.tofile(file_path+'.single')
                binary.tofile(file_path+'.binary')
            elif (output_format=='npy'):
                single.save(file_path+'.single')
                binary.save(file_path+'.binary')
            else:
                raise ValueError('Output format %s is not supported, should be ascii, binary or npy' % output_format)
                

        time_profile['save_data'] += time.time() - start_time

    else:
        start_time = time.time()

        core.append(core_read[tsel])
        rc = core_read.rc[tsel]
        cm_pos = core_read.pos[tsel] - header.pos_offset
        cm_vel = core_read.pos[tsel] - header.vel_offset

        #print('Correct center')
        particle.correctCenter(cm_pos, cm_vel)

        # r2
        particle.calcR2()
        
        single = Particle(**kwargs)
        p1 = Particle(**kwargs)
        p2 = Particle(**kwargs)
        binary = Binary(p1,p2,**kwargs)

        if (output_format=='ascii'): single.loadtxt(file_path+'.single')   
        elif (output_format=='binary'): single.fromfile(file_path+'.single')
        elif (output_format=='npy'): single.load(file_path+'.single.npy')
        else: raise ValueError('Output format %s is not supported, should be ascii, binary or npy' % output_format)

        if (output_format=='ascii'): binary.loadtxt(file_path+'.binary')
        elif (output_format=='binary'): binary.fromfile(file_path+'.binary')
        elif (output_format=='npy'): binary.load(file_path+'.binary.npy')
        else: raise ValueError('Output format %s is not supported, should be ascii, binary or npy' % output_format)

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

    if (r_escape is not None):
        if (r_escape == 'tidal'):
            if (external_mode!='none'): 
                tidal = result['tidal']
                tsel = (core.time == header.time)
                pos_c = core.pos[tsel]
                r_gal = np.sqrt(np.sum(pos_c*pos_c))
                M_galaxy = 0
                if (type(m_ext)==np.ndarray):
                    tsel = (m_ext[:,0] == header.time)
                    M_galaxy= m_ext[tsel,1]
                else:
                    M_galaxy = estimateGalaxyMass(pot_ext, r_gal, G)
                r_escape = tidal.calcTidalSphere(header.time, particle.mass, particle.r2, M_galaxy, pot_ext, r_gal, G);
            else:
                raise ValueError('Escape radius is set to tidal radius but the external mode is off')
        else:
            r_escape = float(r_escape)

        if (e_escape is not None): 
            if (e_escape == 'bound_noext'):
                e_escape = 0
                single.pot -= single.pot_ext
                binary.p1.pot -= binary.p1.pot_ext
                binary.p2.pot -= binary.p2.pot_ext
            else:
                e_escape = float(e_escape)
        else:
            e_escape = 0

        single = esc_single.findEscaper(header.time, single, r_escape, e_escape)
        binary = esc_binary.findEscaper(header.time, binary, r_escape, e_escape)
        time_profile['escaper'] += time.time() - start_time
        start_time = time.time()


    #print('Lagrangian radius')
    lagr.calcOneSnapshot(header.time, single, binary, rc, average_mode)

    time_profile['lagr'] += time.time() - start_time
    start_time = time.time()

    if (r_escape is None):
        rhindex=np.where(m_frac==0.5)[0]
        r_escape = calcREscapeIsolate(lagr.all.r[-1,rhindex])
        esc_single.findEscaper(header.time, single, r_escape)
        esc_binary.findEscaper(header.time, binary, r_escape)

        time_profile['escaper'] += time.time() - start_time
        start_time = time.time()

    if ('bse_status' in result.keys()):
        bse = result['bse_status']
        bse.findEvents(header.time,single,binary)
        time_profile['bse'] += time.time() - start_time

#    return time_profile

def dataProcessList(file_list, read_flag=False, **kwargs):
    """ process lagragian calculation for a list of file snapshots

    Parameters
    ----------
    file_list: list
        file path list
    read_flag: bool (False)
        If true, read single, binary snapshots and core data instead of calculating them
        If find_multiple is used, read_flag will be set to False 
        because the saved single/binary files miss triple/quadruple components
    kwargs: dict ()
        keyword arguments, see dataProcessOne
    """
    result = dict()
    result['lagr']=LagrangianMultiple(**kwargs)
    result['esc_single']=SingleEscaper(**kwargs)
    result['esc_binary']=BinaryEscaper(**kwargs)
    result['tidal']=Tidal(**kwargs)

    realtime_save = kwargs.get('realtime_save', True) & ('filename_prefix' in kwargs.keys())
    realtime_save_mode = kwargs.get('realtime_save_mode', 'w')
    output_format = kwargs.get('output_format', 'binary')
    realtime_save_prefix = kwargs.get('realtime_save_prefix', kwargs.get('filename_prefix', None))

    time_profile=dict()
    for key in ['read','find_pair','density','center_core','save_data','calc_pot','lagr','escaper','bse']:
        time_profile[key] = 0.0

    # when find_multiple is used, the saved single/binary files are not complete and miss triple/quadruple components, thus set read_flag to False
    if ('find_multiple' in kwargs.keys()): 
        find_multiple=kwargs['find_multiple']
        if (find_multiple):
            read_flag = False

    result['core']=Core()
    if (read_flag):
        result['core_read'] = Core()
        core_filename=kwargs['filename_prefix']+'.core'
        result['core_read'].loadtxt(core_filename)

    if ('interrupt_mode' in kwargs.keys()): 
        interrupt_mode=kwargs['interrupt_mode']
        if ('bse' in interrupt_mode):
            result['bse_status'] = BSEStatus()

    if ('read_m_ext' in kwargs.keys()):
        read_m_ext=kwargs['read_m_ext'] # filename of m_ext, (time, m)
        result['m_ext'] = np.loadtxt(read_m_ext)

    realtime_save_keys = ['lagr', 'core']
    if ('bse_status' in result.keys()):
        realtime_save_keys.append('bse_status')
    if (kwargs.get('r_escape', None) == 'tidal'):
        realtime_save_keys.append('tidal')

    if (realtime_save):
        _prepareRealtimeSaveFiles(realtime_save_prefix, realtime_save_keys, output_format, realtime_save_mode)

    for path in file_list:
        size_prev = dict()
        if (realtime_save):
            for key in realtime_save_keys:
                size_prev[key] = result[key].size

        dataProcessOne(path, result, time_profile, read_flag=read_flag, **kwargs)

        if (realtime_save):
            start_time = time.time()
            for key in realtime_save_keys:
                if (result[key].size > size_prev[key]):
                    _saveRealtimeData(result[key], realtime_save_prefix+'.'+key, output_format)
            time_profile['save_data'] += time.time() - start_time

    if (len(file_list)>0):
        for key, item in time_profile.items():
            item /= len(file_list)
    return result, time_profile


def parallelDataProcessList(file_list, n_cpu=int(0), **kwargs):
    """ parellel process lagragian calculation for a list of file snapshots

    Parameters
    ----------
    file_list: list
        file path list
    n_cpu: int
        number of CPU cores to run parallelly
    kwargs: dict
        keyword arguments, see dataProcessOne
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
        kwargs_worker = kwargs.copy()
        if ('filename_prefix' in kwargs_worker.keys()):
            kwargs_worker['realtime_save'] = True
            kwargs_worker['realtime_save_mode'] = 'w'
            kwargs_worker['realtime_save_prefix'] = _getParallelRealtimePrefix(kwargs_worker['filename_prefix'], rank)
        else:
            kwargs_worker['realtime_save'] = False
        result[rank] = pool.apply_async(dataProcessList, (file_part[rank],), kwargs_worker)

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

