#!/usr/bin/env python3

import numpy as np
import multiprocessing as mp
import sys
import os
import petar
import getopt

if __name__ == '__main__':

    filename_prefix='data'
    average_mode='sphere'
    read_flag=False
    n_cpu=0
    write_option='w'

    def usage():
        print("A tool for processing a list of snapshot data to detect binaries,")
        print("   calculate the density center, the core radius, the Langragian radii (using the density center)")
        print("   and the corresponding properties inside each radius: number of objects, average masses, mean velocities and velocity dispersions")
        print("   binaries are counted as single objects using the c.m. properties")
        print("Usage: petar.data.process [options] data_filename")
        print("data_filename: A list of snapshot data path, each line for one snapshot")
        print("option:")
        print("  -h(--help): help")
        print("  -p(--filename-prefix): prefix of output file names for: [prefix].[lagr|esc.[single|binary]|core] (data)")
        print("  -m(--mass-fraction): Lagrangian radii mass fraction (0.1,0.3,0.5,0.7,0.9)")
        print("  -G(--gravitational-constant): Gravitational constant (if interrupt-mode=(mo)bse: ",petar.G_MSUN_PC_MYR,"; else 1.0)")
        print("  -b(--r-max-binary): maximum sepration for detecting binaries (0.1)")
        print("  -B(--full-binary): calculate full binary orbital parameters (simple_mode=False in Binary class), this option increases computing time")
        print("  -a(--average-mode): Lagrangian properity average mode, choices: sphere: average from center to Lagragian radii; shell: average between two neighbor radii (sphere)")
        print("  -A(--append): append new data to existing data files")
        print("  -r(--read-data): read existing single, binary and core data to avoid expensive KDTree construction, no argument, disabled in default")
        print("     --r-escape: a constant escape distance criterion, in default, it is 20*half-mass radius")
        print("     --e-escape: escape energy criterion, only work together with --r-escape, in default, it is 0.0")
        print("  -i(--interrupt-mode): the interruption mode used in petar, choices: no, base, bse, mobse (no)")
        print("  -t(--external-mode): external mode used in petar, choices: galpy, no (no)")
        print("  -s(--snapshot-format): snapshot data format: binary, ascii (ascii)")
        print("  -n(--n-cpu): number of CPU threads for parallel processing (all threads)")
        print("     --add-star-type: calculate Lagrangian radii and properties for specific types of stars.")
        print("          This argument contain a list of type names, separated by ',' (no space)")
        print("          For each given type name, an additional group of data are added in the Lagrangian data file [prefix].lagr.")
        print("          There are four styles of type names:")
        print("            (1) a single SSE type name (see help(petar.SSEType))")
        print("                For example, if 'BH' is given, one additional class member (type is class Lagrangian), BH, is added.")
        print("                The Lagrangian radii are determined by only counting BHs, ")
        print("                then the properties (average mass, velocity...) of BHs are calculated using these Lagrangian radii.")
        print("            (2) a combination of different SSE types connected by '_'.")
        print("                This will include mutliple SSE types in the calculation.")
        print("                For example, if 'BH_NS_WD' is given, the additional class member, BH_NS_WD,")
        print("                count BHs, NSs and WDs together to calculate the Lagrangian properties.")
        print("            (3) a single SSE type name with the prefix 'no'")
        print("                This will exclude the given SSE type name in the calculation.")
        print("                For exmaple, if 'noBH' is given, the additional class member, noBH, ")
        print("                count all stars except BHs to calculate the Lagrangian properties.")
        print("            (4) two types (can be any case of the style 1-3) are given by '[type 1]__in__[type 2]'")
        print("                Then, the properties of type 1 stars within the spheres or shells of ")
        print("                Lagragian radii of type 2 stars are calculated. The Lagrangian radii of type 1 are not calculated.")
        print("                For example, if 'BH__in__all' is given, Lagrangian properties of BHs are calculated ")
        print("                within the shell or sphere of Lagrangian radii of all stars (instead of Lagrangian radii of BHs).")
        print("                Notice that the two type names (except 'all') should also be added separately in the list.")
        print("                In the case of 'BH__in__all', add_star_type must contain 'BH', i.e. add_star_type=['BH','BH__in__all', ...].")
        print("                In another example, 'BH__in__MS', add_star_type=['BH','MS','BH__in__MS',...].")
        print("          The SSE star type names are shown below:")
        print("              LMS: deeply or fully convective low mass MS star [0]")
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
        print("          All these styles can be combined, e.g. '--add-star-type NS_BH,MS,NS_BH__in__MS,noBH'")
        print("          When this option is used, to read the generated lagrangian data by using petar.LagrangianMultiple in Python3,")
        print("          the consistent keyword argument 'add_star_type' should be used.")
        print("          For example, when '--add-star-type BH,MS' is used, the initialization of petar.LagrangianMultiple should be:")
        print("          lagr = petar.LagrangianMultiple(add_star_type=['BH','MS']) .")
        print("Important note: 1) users should be careful to set the consistent '-i' or -'G' options in order to correctly calculate the Kepler orbital parameters of binaries.")
        print("                2) when data are written in BINARY format, '-s binary' should be used.")
        print("                3) '--add-star-type' only works when the interrupt mode is 'bse' or 'mobse'.")
    try:
        shortargs = 'p:m:G:b:BAa:rt:i:s:n:h'
        longargs = ['mass-fraction=','gravitational-constant=','r-max-binary=','full-binary','average-mode=', 'filename-prefix=','read-data','r-escape=','append','e-escape=','external-mode=','interrupt-mode=','snapshot-format=','add-star-type=','n-cpu=','help']
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
            elif opt in ('-G','--gravitational-constant'):
                kwargs['G'] = float(arg)
            elif opt in ('-b','--r-max-binary'):
                kwargs['r_max_binary'] = float(arg)
            elif opt in ('-B','--full-binary'):
                kwargs['simple_binary'] = False
            elif opt in ('-a','--average-mode'):
                kwargs['average_mode'] = arg
            elif opt in ('-A','--append'):
                write_option='a'
            elif opt in ('-n','--n-cpu'):
                n_cpu = int(arg)
            elif opt in ('-i','--interrupt-mode'):
                kwargs['interrupt_mode'] = arg
            elif opt in ('-t','--external-mode'):
                kwargs['external_mode'] = arg
            elif opt in ('-s','--snapshot-format'):
                kwargs['snapshot_format'] = arg
            elif opt in ('-r','--read-data'):
                read_flag = True
            elif opt in ('--r-escape'):
                kwargs['r_escape'] = float(arg)
            elif opt in ('--e-escape'):
                kwargs['e_escape'] = float(arg)
            elif opt in ('--add-star-type'):
                kwargs['add_star_type'] = [x for x in arg.split(',')]
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

    kwargs['filename_prefix'] = filename_prefix

    for key, item in kwargs.items(): print(key,':',item)

    fl = open(filename,'r')
    file_list = fl.read()
    path_list = file_list.splitlines()

    result=dict()
    time_profile=dict()
    if (n_cpu==1):
        result,time_profile = petar.dataProcessList(path_list, read_flag, **kwargs)
    else:
        result,time_profile = petar.parallelDataProcessList(path_list, n_cpu, read_flag, **kwargs)

    for key in ['lagr','core','bse_status']:
        if key in result.keys():
            key_filename  = filename_prefix + '.' + key
            with open(key_filename, write_option) as f:
                result[key].savetxt(f)
                print (key,"data is saved in file:",key_filename)
     

    if (write_option=='a'):
        if 'esc_single' in result.keys():
            key_filename  = filename_prefix + '.esc_single'
            data_read=petar.SingleEscaper(**kwargs)
            if os.path.getsize(key_filename)>0:
                data_read.loadtxt(key_filename)
                result_mix=petar.join(data_read,result['esc_single'])
                result_mix.removeDuplicate()
                with open(key_filename, 'w') as f:
                    result_mix.savetxt(f)
                    print ("esc_single data is saved in file:",key_filename)
            else:
                with open(key_filename, 'w') as f:
                    result['esc_single'].savetxt(f)
                    print ("esc_single data is saved in file:",key_filename)
        if 'esc_binary' in result.keys():
            key_filename  = filename_prefix + '.esc_binary'
            data_read=petar.BinaryEscaper(**kwargs)
            if os.path.getsize(key_filename)>0:
                data_read.loadtxt(key_filename)
                result_mix=petar.join(data_read,result['esc_binary'])
                result_mix.removeDuplicate()
                with open(key_filename, 'w') as f:
                    result_mix.savetxt(f)
                    print ("esc_binary data is saved in file:",key_filename)
            else:
                with open(key_filename, 'w') as f:
                    result['esc_binary'].savetxt(f)
                    print ("esc_binary data is saved in file:",key_filename)

    else:
        for key in ['esc_single','esc_binary']:
            if key in result.keys():
                key_filename  = filename_prefix + '.' + key
                with open(key_filename, write_option) as f:
                    result[key].savetxt(f)
                    print (key,"data is saved in file:",key_filename)
            

    print ('CPU time profile:')
    for key, item in time_profile.items():
        print (key,item,)
