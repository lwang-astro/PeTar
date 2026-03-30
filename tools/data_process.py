#!/usr/bin/env python3

import numpy as np
import multiprocessing as mp
import sys
import os
import glob
import petar
import getopt


def getParallelRealtimePrefix(filename_prefix, rank):
    return filename_prefix + '.parallel.' + str(rank)


def clearParallelRealtimeFiles(filename_prefix, save_keys):
    for key in save_keys:
        flist = glob.glob(filename_prefix + '.parallel.*.' + key)
        for fname in flist:
            if os.path.exists(fname):
                os.remove(fname)


def readMixDataByKey(key, filename, output_format, kwargs):
    if (key == 'lagr'):
        data = petar.LagrangianMultiple(**kwargs)
    elif (key == 'core'):
        data = petar.Core()
    elif (key == 'bse_status'):
        data = petar.BSEStatus()
    elif (key == 'tidal'):
        data = petar.Tidal(**kwargs)
    else:
        raise ValueError('Unknown key for reading mix data: %s' % key)

    if (output_format == 'ascii'):
        data.loadtxt(filename)
    elif (output_format == 'binary'):
        data.fromfile(filename)
    elif (output_format == 'npy'):
        data.load(filename)
    else:
        raise ValueError('Output format %s is unknown, should be ascii, binary or npy.' % output_format)
    return data


def writeMixData(data, filename, output_format):
    if (output_format == 'ascii'):
        with open(filename, 'w') as f:
            data.savetxt(f)
    elif (output_format == 'binary'):
        with open(filename, 'wb') as f:
            data.tofile(f)
    elif (output_format == 'npy'):
        with open(filename, 'wb') as f:
            data.save(f)
    else:
        raise ValueError('Output format %s is unknown, should be ascii, binary or npy.' % output_format)


def deduplicateByTime(data):
    if (data.size == 0):
        return data
    if (not hasattr(data, 'time')):
        return data

    order = np.argsort(data.time, kind='stable')
    data_sort = data[order]

    # keep the last entry for duplicated time (the most recently recovered one)
    _, ridx = np.unique(data_sort.time[::-1], return_index=True)
    keep_idx = data_sort.size - 1 - ridx
    keep_idx.sort()
    return data_sort[keep_idx]


def recoverParallelRealtimeFiles(filename_prefix, save_keys, output_format, kwargs, append_mode=False):
    recovered_keys = []
    for key in save_keys:
        flist = glob.glob(filename_prefix + '.parallel.*.' + key)
        if (len(flist) == 0):
            continue

        def get_rank(fname):
            fs = fname.split('.')
            if (len(fs) < 3):
                return -1
            try:
                return int(fs[-2])
            except ValueError:
                return -1

        flist.sort(key=get_rank)

        data_temp = []
        for fname in flist:
            if os.path.getsize(fname) > 0:
                data_temp.append(readMixDataByKey(key, fname, output_format, kwargs))

        if (len(data_temp) == 0):
            continue

        data_merge = data_temp[0]
        if (len(data_temp) > 1):
            data_merge = petar.join(*data_temp)

        key_filename = filename_prefix + '.' + key
        if append_mode and os.path.exists(key_filename) and os.path.getsize(key_filename) > 0:
            data_old = readMixDataByKey(key, key_filename, output_format, kwargs)
            data_merge = petar.join(data_old, data_merge)

        data_merge = deduplicateByTime(data_merge)

        writeMixData(data_merge, key_filename, output_format)
        recovered_keys.append(key)

    if (len(recovered_keys) > 0):
        clearParallelRealtimeFiles(filename_prefix, save_keys)

    return recovered_keys


def getProcessedTimeSet(filename_prefix, output_format, kwargs):
    required_keys = ['lagr', 'core']
    if ('interrupt_mode' in kwargs.keys()):
        if ('bse' in kwargs['interrupt_mode']):
            required_keys.append('bse_status')
    if (kwargs.get('r_escape', None) == 'tidal'):
        required_keys.append('tidal')

    processed_time = None
    for key in required_keys:
        key_filename = filename_prefix + '.' + key
        if (not os.path.exists(key_filename)):
            return set()
        if (os.path.getsize(key_filename) == 0):
            return set()

        data = readMixDataByKey(key, key_filename, output_format, kwargs)
        if (data.size == 0):
            return set()

        time_set = set(np.round(data.time, 12))
        if (processed_time is None):
            processed_time = time_set
        else:
            processed_time = processed_time & time_set

        if (len(processed_time) == 0):
            return set()

    return processed_time


def readSnapshotHeader(path, kwargs):
    header_kwargs = kwargs.copy()
    header_kwargs.setdefault('snapshot_format', 'binary')
    return petar.PeTarDataHeader(path, **header_kwargs)


def filterUnprocessedSnapshots(path_list, filename_prefix, output_format, kwargs):
    processed_time = getProcessedTimeSet(filename_prefix, output_format, kwargs)
    if (len(processed_time) == 0):
        return path_list, 0

    path_rest = []
    n_skip = 0
    for path in path_list:
        header = readSnapshotHeader(path, kwargs)
        tkey = np.round(header.time, 12)
        if (tkey in processed_time):
            n_skip += 1
        else:
            path_rest.append(path)
    return path_rest, n_skip


def filterSnapshotsFromTime(path_list, time_min, kwargs):
    if (time_min is None):
        return path_list, 0

    path_rest = []
    n_skip = 0
    for path in path_list:
        header = readSnapshotHeader(path, kwargs)
        if (header.time >= time_min):
            path_rest.append(path)
        else:
            n_skip += 1
    return path_rest, n_skip

if __name__ == '__main__':

    filename_prefix='data'
    average_mode='sphere'
    ftid_file_flag=False
    n_cpu=0
    write_option='w'
    esc_snapshot_format='binary'
    output_format='binary'
    recover_parallel_only=False
    auto_resume=True
    resume_from_time=None

    def usage():
        print("A tool for post-data processing of a list of snapshot files from petar.")
        print("Functionality:")
        print("   1) Generate new snapshots of singles, binaries, and optionally triples/quadruples for each snapshot file.")
        print("   2) Calculate the density center, core radius, and Lagrangian properties based on the density center, including")
        print("      radii and the corresponding properties inside each radius: number of objects, average masses, mean velocities, and velocity dispersions.")
        print("      Binaries are treated as single objects using their center of the mass.")
        print("   3) Identify single and binary escapers and save them into files.")
        print("Usage: petar.data.process [options] [snapshot path list filename]")
        print("   snapshot path list file: A list of snapshot data paths, each line for one snapshot.")
        print("                            This file can be generated by petar.data.gether.")
        print("Options (default arguments shown in parentheses at the end):")
        print("  -h(--help)                Display help information.")
        print("  -p(--filename-prefix) [S] Prefix of output file names as: [prefix].[lagr|esc_[single|binary]|core] (default: data).")
        print("  -m(--mass-fraction)   [S] Lagrangian radii mass fraction, seperated by ',' without empty spaces (default: 0.1,0.3,0.5,0.7,0.9).")
        print("  -G(--gravitational-constant) [F] Gravitational constant (if interrupt-mode=*bse*: "+str(petar.G_MSUN_PC_MYR)+"; else 1.0).")
        print("  -b(--r-max-binary)    [F] Maximum separation for detecting binaries (default: 0.1).")
        print("  -B(--full-binary)         Calculate detailed binary orbital parameters (including orbital angle and phase) with time-consuming computation;")
        print("                            Without this option, only basic orbital parameters (semi-major axis and eccentricity) are calculated.")
        print("  -M(--multiple)            Detect multiple systems (binaries, triples, and quadruples) and save to snapshot files [snapshot_filename].[single|binary|triple|quadruple];")
        print("                            Without this option, only singles and binaries are detected.")
        print("                            When this option is used, the option '-r' cannot be used to restart data processing.")
        print("  -a(--average-mode)    [S] Lagrangian property average mode; choices: (default sphere).")
        print("                                sphere: average from center to Lagrangian radii")
        print("                                shell (average between two neighboring radii")
        print("  -A(--append)              Append new data to existing data files.")
        print("  -r(--read-data)           Read existing single, binary, and core data to avoid expensive KDTree construction; no argument, disabled by default.")
        print("     --r-escape       [S|F] Distance criterion for escaper; if the value is 'tidal', calculate the tidal radius (only works when external-mode is on); otherwise, it is a constant escape distance criterion. If not given, it is 20 times the half-mass radius.")
        print("     --e-escape       [S|F] Energy criterion for escaper when objects are outside distance criterion; only works when --r-escape is used; if the value is 'bound_noext', calculate bound energy without external potential and remove etot > 0; otherwise etot > mass * e-escape (default: 0.0).")
        print("     --m-ext            [S] Read a table of masses of external potential for each time, used for the calculation of tidal radius (default: not used).")
        print("                            The argument is the filename of the table. The file contains two columns: time, mass.")
        print("  -i(--interrupt-mode)  [S] The interruption mode used in petar; choices: no, base, bse, mobse, bseEmp (default: no).")
        print("  -t(--external-mode)   [S] External mode used in petar; choices: galpy, no (default: no).")
        print("  -P(--use_mpfrc)           Include three columns of high-precision parts of particle position x, y, z.")
        print(f"  -s(--snapshot-format) [S] Input snapshot data format: binary, ascii (default: binary).")
        print("                             Refer to the '-i' option of petar.")
        print("     --esc-snapshot-format [S] Set escaper snapshot data format for reading: binary, ascii, npy (follows -s).")
        print("                               These files (*.esc_single, *.esc_binary) are generated by the previous petar.data.process.")
        print("                               This option is used when the option '--append' is switched on.")
        print(f"  -o(--output-format)  [S] Output data format for single, binary, and multiple snapshots: ascii, binary, npy (default: {output_format}).")
        print("     --recover-parallel-only   Only recover and merge unfinished parallel temporary files, then exit.")
        print("     --no-auto-resume          Disable automatic interruption detection/recovery and resume from remaining snapshots (enabled by default).")
        print("     --resume-from-time    [F] Force processing snapshots from this time and later (used as an additional filter on top of auto-resume).")
        print("  -e(--calc-energy)         Enable the calculation of potential energy and virial ratio -(2*ekin/epot) of different Lagrangian radii.")
        print("  -c(--calc-multi-rc)       Enable the calculation of individual core radius for each group chosen for Lagrangian properties (e.g., single, binary, and star type);")
        print("                            The centers are also recalculated for individual groups (time-consuming computation).")
        print("  -n(--n-cpu)           [I] Number of CPU threads for parallel processing (default: all threads).")
        print("     --add-star-type    [S] Calculate addtional Lagrangian properties for specific types of stars.")
        print("          This argument contains a list of type names, separated by ',' without empty spaces.")
        print("          For each given type name, an additional group of data is added in the Lagrangian data file [prefix].lagr.")
        print("          There are four styles of type names:")
        print("            (1) a single SSE type name")
        print("                Calculate Lagrangian properties for the specified single stellar type.")
        print("                For example, if 'BH' is provided, the Lagrangian properties of black holes are calculated.")
        print("            (2) a combination of different SSE types connected by '_'")
        print("                Calculate Lagrangian properties for the specified multiple stellar types.")
        print("                For example, if 'BH_NS_WD' is provided, the Lagrangian properties for black holes, neutron stars, and white dwarfs are calculated.")
        print("            (3) a single SSE type name with the prefix 'no'")
        print("                Calculate Lagrangian properties excluding the specified stellar type.")
        print("                For example, if 'noBH' is provided, the Lagrangian properties excluding black holes are calculated.")
        print("            (4) two types (can be any combination of styles 1-3) are given as '[type 1]__in__[type 2]'")
        print("                Calculate Lagrangian properties excluding the radii of type 1 within the sphere or shell defined by the Lagrangian radii of type 2.")
        print("                For example, if 'BH__in__all' is provided, Lagrangian properties of black holes are calculated within the shell or sphere defined by the Lagrangian radii of all stars.")
        print("                Note that the two type names (excluding 'all') should also be added separately in the list.")
        print("                For example, if 'BH__in__MS' is included, 'BH' and 'MS' should both be added together as: BH,MS,BH__in__MS.")
        print("          - All these styles can be combined and calculated simultaneously, e.g., '--add-star-type NS_BH,MS,NS_BH__in__MS,noBH'.")
        print("          - When this option is used, to read the [prefix].lagr file using petar.LagrangianMultiple, the consistent keyword argument 'add_star_type' should be used.")
        print("            For example, if '--add-star-type BH,MS' is used in petar.data.process, petar.LagrangianMultiple should include the keyword argument 'add_star_type=['BH','MS']'.")
        print("            The corresponding class member names are 'BH' and 'MS'.")
        print("          - The SSE star type names are shown below:")
        print("              LMS:  Deeply or fully convective low mass MS star [0]")
        print("              MS:   Main Sequence star [1]")
        print("              HG:   Hertzsprung Gap [2]")
        print("              GB:   First Giant Branch [3]")
        print("              CHeB: Core Helium Burning [4]")
        print("              FAGB: First Asymptotic Giant Branch [5]")
        print("              SAGB: Second Asymptotic Giant Branch [6]")
        print("              HeMS: Main Sequence Naked Helium star [7]")
        print("              HeHG: Hertzsprung Gap Naked Helium star [8]")
        print("              HeGB: Giant Branch Naked Helium star [9]")
        print("              HeWD: Helium White Dwarf [10]")
        print("              COWD: Carbon/Oxygen White Dwarf [11]")
        print("              ONWD: Oxygen/Neon White Dwarf [12]")
        print("              NS:   Neutron Star [13]")
        print("              BH:   Black Hole [14]")
        print("              SN:   Massless Supernova [15]")
        print("     --add-mass-range   [S] Calculate addtional Lagrangian properties for specific mass ranges of objects.")
        print("          This argument contains a list of mass ranges, separated by ',' without empty spaces.")
        print("          For each given mass range, an additional group of data is added in the Lagrangian data file [prefix].lagr.")
        print("          There are two styles of mass ranges:")
        print("            (1) [minimum mass]_[maximum mass]")
        print("                Calculate Lagrangian properties for the objects in a mass range.")
        print("                For example, if '0.08_1' is given, Lagrangian properties are calculated by selecting objects with masses from 0.08 to 1.0.")
        print("                Be aware that the minimum mass must be > 0.")
        print("            (2) [mass range 1]__in__[mass range 2]")
        print("                Calculate Lagrangian properties of the mass range 1 within the sphere or shell defined by the Lagrangian radii of mass range 2.")
        print("                The [mass range 1] and [mass range 2] have the same syntax of the style (1).")
        print("                For example, if '1_150__in__0.08_1' is given, Lagrangian properties with masses from 1 to 150 are calculated within the sphere or shell defined by the Lagrangian radii with masses from 0.08 to 1.0.")
        print("                To use this style,  both [mass range 1] and [mass range 2] should be added simultaneously.")
        print("          - All these styles can be combined, similar to '--add-star-type'.")
        print("          - When this option is used, to read the [prefix].lagr file using petar.LagrangianMultiple, the consistent keyword argument 'add_mass_range' should be used.")
        print("            For example, if '--add-mass-range 0.08_1,1_150,0.08_1__in__1_150' is used in petar.data.process, petar.LagrangianMultiple should include the keyword argument 'add_mass_range=['0.08_1','1_150','0.08_1__in__1_150'.")
        print("            The corresponding class member names are 'mass_0.08_1', 'mass_1_150', and 'mass_0.08_1__in__1_150'.")
        print("Important notes:")
        print("  1) Ensure correct options are set for '-i', '-t', and '-G' to read snapshots accurately and calculate Kepler orbital parameters correctly.")
        print("     When using the compiled SSE/BSE stellar evolution package, use '-i bse'. Note that even if SSE/BSE is compiled but switched off during petar usage, '-i bse' is still required.")
        print("     Similarly, when the Galpy external potential support is compiled, use '-i galpy' regardless of whether external potential is set in petar options during simulation.")
        print("     Make sure to set the correct value for '-G' based on the units used during petar usage.")
        print("  2) If data is written in BINARY format during petar simulation, use '-s binary'.")
        print("  3) '--add-star-type' functionality is only available when SSE/BSE is used.")
    try:
        shortargs = 'p:m:G:b:MBAea:rt:i:Ps:o:cn:h'
        longargs = ['mass-fraction=','multiple','gravitational-constant=','r-max-binary=','full-binary','average-mode=', 'filename-prefix=','read-data','calc-energy','r-escape=','append','e-escape=','external-mode=','interrupt-mode=','use-mpfrc','snapshot-format=','output-format=','m-ext=','add-star-type=','add-mass-range=','calc-multi-rc','n-cpu=','recover-parallel-only','no-auto-resume','resume-from-time=','help']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        kwargs=dict()
        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-p','--filename-prefix'):
                filename_prefix = arg
            elif opt in ('-m','--mass-fraction'):
                kwargs['mass_fraction'] = np.array([float(x) for x in arg.split(',')])
            elif opt in ('-M','--multiple'):
                kwargs['find_multiple'] = True
            elif opt in ('-G','--gravitational-constant'):
                kwargs['G'] = float(arg)
            elif opt in ('-b','--r-max-binary'):
                kwargs['r_max_binary'] = float(arg)
            elif opt in ('-B','--full-binary'):
                kwargs['simple_mode'] = False
            elif opt in ('-a','--average-mode'):
                kwargs['average_mode'] = arg
            elif opt in ('-A','--append'):
                write_option='a'
            elif opt in ('-e','--calc-energy'):
                kwargs['calc_energy']=True
            elif opt in ('-n','--n-cpu'):
                n_cpu = int(arg)
            elif opt in ('-i','--interrupt-mode'):
                kwargs['interrupt_mode'] = arg
            elif opt in ('-t','--external-mode'):
                kwargs['external_mode'] = arg
            elif opt in ('-P','--use_mpfrc'):
                kwargs['use_mpfrc'] = True
            elif opt in ('-s','--snapshot-format'):
                kwargs['snapshot_format'] = arg
                if (not 'esc_snapshot_format' in kwargs.keys()):
                    kwargs['esc_snapshot_format'] = arg
                    esc_snapshot_format = arg
            elif opt in ('--esc-snapshot-format'):
                kwargs['esc_snapshot_format'] = arg
                esc_snapshot_format = arg
            elif opt in ('-o','--output-format'):
                kwargs['output_format'] = arg
                output_format = arg
            elif opt in ('-r','--read-data'):
                kwargs['read_flag'] = True
            elif opt in ('-c','--calc-multi-rc'):
                kwargs['calc_multi_rc']=True
            elif opt in ('--r-escape'):
                if (arg=='tidal'): 
                    kwargs['r_escape'] = arg
                    ftid_file_flag = True
                else: kwargs['r_escape'] = float(arg)
            elif opt in ('--e-escape'):
                kwargs['e_escape'] = arg
            elif opt in ('--m-ext'):
                kwargs['read_m_ext'] = arg
            elif opt in ('--add-star-type'):
                kwargs['add_star_type'] = [x for x in arg.split(',')]
            elif opt in ('--add-mass-range'):
                kwargs['add_mass_range'] = [x for x in arg.split(',')]
            elif opt == '--recover-parallel-only':
                recover_parallel_only = True
            elif opt == '--no-auto-resume':
                auto_resume = False
            elif opt == '--resume-from-time':
                resume_from_time = float(arg)
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError:
        print('getopt error!')
        usage()
        sys.exit(1)

    filename = remainder[0]

    if (not 'G' in kwargs.keys()):
        if ('interrupt_mode' in kwargs.keys()):
            if ('bse' in kwargs['interrupt_mode']): kwargs['G'] = 0.00449830997959438 # pc^3/(Msun*Myr^2)

    fl = open(filename,'r')
    file_list = fl.read()
    path_list = file_list.splitlines()
    fl.close()

    kwargs['filename_prefix'] = filename_prefix
    kwargs['realtime_save'] = (n_cpu == 1)
    kwargs['realtime_save_mode'] = write_option

    recover_keys=['lagr','core','bse_status']
    if (ftid_file_flag): recover_keys.append('tidal')

    if (auto_resume or recover_parallel_only):
        recovered = recoverParallelRealtimeFiles(filename_prefix, recover_keys, output_format, kwargs, append_mode=True)
        if (len(recovered) > 0):
            print('Recovered parallel temporary files for keys:', ','.join(recovered))
        else:
            print('No parallel temporary files found to recover.')

    if (auto_resume):
        path_list, n_skip = filterUnprocessedSnapshots(path_list, filename_prefix, output_format, kwargs)
        if (n_skip > 0):
            write_option = 'a'
            kwargs['realtime_save_mode'] = write_option
            print('Auto-resume: detected %d processed snapshots, continue with %d remaining snapshots.' % (n_skip, len(path_list)))
            if (len(path_list) > 0):
                header_next = readSnapshotHeader(path_list[0], kwargs)
                print('Auto-resume start from snapshot:', path_list[0], 'time:', header_next.time)

    if (resume_from_time is not None):
        path_list, n_skip_time = filterSnapshotsFromTime(path_list, resume_from_time, kwargs)
        if (n_skip_time > 0):
            write_option = 'a'
            kwargs['realtime_save_mode'] = write_option
        print('Resume-from-time: %.12g, skipped %d snapshots, remaining %d snapshots.' % (resume_from_time, n_skip_time, len(path_list)))
        if (len(path_list) > 0):
            header_next = readSnapshotHeader(path_list[0], kwargs)
            print('Resume-from-time start snapshot:', path_list[0], 'time:', header_next.time)

    if (recover_parallel_only):
        print('Recover-only mode enabled, stop without new snapshot processing.')
        sys.exit(0)

    if (len(path_list) == 0):
        print('No snapshot needs processing, stop.')
        sys.exit(0)

    for key, item in kwargs.items(): print(key,':',item)

    result=dict()
    time_profile=dict()
    if (n_cpu==1):
        result,time_profile = petar.dataProcessList(path_list, **kwargs)
    else:
        result,time_profile = petar.parallelDataProcessList(path_list, n_cpu, **kwargs)

    fout_list=['lagr','core','bse_status']
    if (ftid_file_flag): fout_list.append('tidal')

    realtime_saved_keys = set()
    if (kwargs.get('realtime_save', False)):
        realtime_saved_keys.update(fout_list)

    for key in fout_list:
        if key in result.keys():
            key_filename  = filename_prefix + '.' + key
            if (key in realtime_saved_keys) and (result[key].size > 0):
                print (key,"data is updated during processing in file:",key_filename)
                continue
            if (output_format == 'npy') and (write_option == 'a'):
                result_mix = result[key]
                if os.path.exists(key_filename) and (os.path.getsize(key_filename) > 0):
                    data_old = readMixDataByKey(key, key_filename, output_format, kwargs)
                    result_mix = petar.join(data_old, result[key])
                    result_mix = deduplicateByTime(result_mix)
                writeMixData(result_mix, key_filename, output_format)
                print (key,"data is saved in file:",key_filename)
                continue
            write_mode = write_option
            if (output_format in ['binary','npy']):
                write_mode += 'b'
            with open(key_filename, write_mode) as f:
                if (output_format=='ascii'):
                    result[key].savetxt(f)
                elif (output_format=='binary'):
                    result[key].tofile(f)
                elif (output_format=='npy'):
                    result[key].save(f)
                else:
                    raise ValueError('Output format %s is unknown, should be ascii, binary or npy.' % output_format)
                print (key,"data is saved in file:",key_filename)

    for key in ['esc_single','esc_binary']:
        if key in result.keys():
            key_filename  = filename_prefix + '.' + key
            result_mix = result[key]
            if (write_option=='a'):
                data_read = None
                if (key == 'esc_single'):
                    data_read=petar.SingleEscaper(**kwargs)
                else:
                    data_read=petar.BinaryEscaper(**kwargs)
                if os.path.getsize(key_filename)>0:
                    if (esc_snapshot_format=='ascii'):
                        data_read.loadtxt(key_filename)
                    elif (esc_snapshot_format=='binary'):
                        data_read.fromfile(key_filename)
                    elif (esc_snapshot_format=='npy'):
                        data_read.load(key_filename)
                    else:
                        raise ValueError('Escape snapshot format %s unknown, should be ascii, binary or npy.' % esc_snapshot_format)                                            
                    result_mix = petar.join(data_read,result[key])
                    result_mix.removeDuplicate()
                
            if (output_format=='ascii'):
                with open(key_filename, 'w') as f:
                    result_mix.savetxt(f)
            elif (output_format=='binary'):
                with open(key_filename, 'wb') as f:
                    result_mix.tofile(f)
            elif (output_format=='npy'):
                with open(key_filename, 'wb') as f:
                    result_mix.save(f)
            else:
                raise ValueError('Output format %s is unknown, should be ascii, binary or npy.' % output_format)
            print ("%s data is saved in file: %s" % (key,key_filename))

    if (n_cpu>1):
        clearParallelRealtimeFiles(filename_prefix, fout_list)

    print ('CPU time profile:')
    for key, item in time_profile.items():
        print (key,item,)
