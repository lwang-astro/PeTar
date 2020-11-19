#!/usr/bin/env python3

import numpy as np
import multiprocessing as mp
import sys
import petar
import getopt

if __name__ == '__main__':

    filename_prefix='data'
    average_mode='sphere'
    read_flag=False
    n_cpu=0

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
        print("  -G(--gravitational-constant): Gravitational constant (if interrupt-mode=bse: ",petar.G_MSUN_PC_MYR,"; else 1.0)")
        print("  -b(--r-max-binary): maximum sepration for detecting binaries (0.1)")
        print("  -B(--full-binary): calculate full binary orbital parameters (simple_mode=False in Binary class), this option increases computing time")
        print("  -a(--average-mode): Lagrangian properity average mode, choices: sphere: average from center to Lagragian radii; shell: average between two neighbor radii (sphere)")
        print("  -r(--read-data): read existing single, binary and core data to avoid expensive KDTree construction, no argument, disabled in default")
        print("     --r-escape: a constant escape distance criterion, in default, it is 20*half-mass radius")
        print("     --e-escape: escape energy criterion, only work together with --r-escape, in default, it is 0.0")
        print("  -i(--interrupt-mode): the interruption mode used in petar, choices: no, base, bse (no)")
        print("  -t(--external-mode): external mode used in petar, choices: galpy, no (no)")
        print("  -s(--snapshot-format): snapshot data format: binary, ascii (ascii)")
        print("  -n(--n-cpu): number of CPU threads for parallel processing (all threads)")
        print("     --add-star-type: calculate Lagrangian radii and properties for specific (SSE) types of stars.")
        print("          A list of SSE star type names should be provided. The choices of names are:")
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
        print("          The values in [] are the type indices defined by SSE/BSE.")
        print("          For each given name, two additional groups of data are added in the Lagrangian data file [prefix].lagr: [star_type_name], [star_type_name]_cross. The definitions are:")
        print("              [star_type_name]: the Lagrangian radii and corresponding properties by only counting this type of stars.")
        print("                                One binary is counted once using its c.m. position and velocity, if two components are both this type of star, the total mass is used.")
        print("                                If one component is this type of star, its mass and c.m. position and velocity are used.")
        print("              [star_type_name]: using the Lagrangian radii of all stars to calculate the Lagrangian properties of this type of stars.")
        print("          Based on the type names, there are three styles of arguments to set the list of [star_type_name]:")
        print("             case 1: a list of names separated by ',': MS,BH,NS ")
        print("                     Calculate Lagrangian radii and properties of MS, BH, NS, respectively.")
        print("                     In the Lagrangian data file [prefix].lagr, 6 additional class members are added: MS, MS_cross, BH, BH_cross, MS, MS_cross.")
        print("             case 2: a combination of names connected by '_': BH_NS_WD]")
        print("                     Calculate Lagrangian radii and propertis by counting BHs, NSs and WDs together; additional class members are BH_NS_WD, BH_NS_WD_cross.")
        print("             case 3: a name with the prefix 'no': noBH")
        print("                     Calculate Lagrangian radii and properties by counting all stars except for BHs; additional class members are noBH, noBH_cross.")
        print("          All three styles can be combined, e.g. '--add-star-type NS_BH,MS,noWD' gives NS_BH, NS_BH_cross, MS, MS_cross, noWD, noWD_cross.")
        print("          To read the [prefix].lagr in Python3, the consistent keyword argument 'add_star_type' should be used in the initialization of petar.LagrangianMultiple.")
        print("          For the above example, the keyword argument is add_star_type=['NS_BH','MS','noWD'].")
        print("          PS: this option only works when the interrupt mode is 'bse'.")
        print("Important note: 1) users should be careful to set the consistent '-i' or -'G' options in order to correctly calculate the Kepler orbital parameters of binaries.")
        print("                2) when data are written in BINARY format, '-s binary' should be used.")

    try:
        shortargs = 'p:m:G:b:Ba:rt:i:s:n:h'
        longargs = ['mass-fraction=','gravitational-constant=','r-max-binary=','full-binary','average-mode=', 'filename-prefix=','read-data','r-escape=','e-escape=','external-mode=','interrupt-mode=','snapshot-format=','add-star-type=','n-cpu=','help']
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
            if (kwargs['interrupt_mode']=='bse'): kwargs['G'] = 0.00449830997959438 # pc^3/(Msun*Myr^2)

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

    for key in ['lagr','core','bse_status', 'esc_single', 'esc_binary']:
        if key in result.keys():
            key_filename  = filename_prefix + '.' + key
            result[key].savetxt(key_filename)
            print (key,"data is saved in file:",key_filename)
     
    print ('CPU time profile:')
    for key, item in time_profile.items():
        print (key,item,)
