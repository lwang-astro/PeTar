#!/usr/bin/env python3
import sys
import getopt
import numpy as np
import warnings

def getPotInstance(pot_name):
    """ Try to get an instance of a give pot_name
    Return
    ----------
    pot_module: module object of pot_name, if it is a combined list, return None
    pot_instance: module instance, if it is not 3D or not available, return None

    """
    pot_module = None
    pot_instance = None
    if (pot_name in dir(galpy.potential)) & ('Potential' in pot_name):
        pot_module = galpy.potential.__getattribute__(pot_name)
        if (type(pot_module) == list): 
            pot_instance = pot_module
            pot_module = None
        elif (type(pot_module) == type):
            # get instance
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    pot_instance = pot_module()
                except (ValueError, TypeError, AttributeError, RuntimeWarning): 
                    pot_module = None
                    pot_instance = None
        else:
            pot_instance = pot_module
            pot_module = type(pot_module)
            
        if (pot_instance != None):
            # remove 2D models
            if (galpy.potential._dim(pot_instance)!=3):
                pot_instance = None
            # remove potential without c support
            if (not _check_c(pot_instance)):
                pot_instance = None

    return pot_module, pot_instance

def savePotTypeArg(config_filename, n_pot, pot_type, pot_arg):
    """ save type argument of a potential

    Parameters
    ----------
    config_filename: string
        file name of output
    n_pot: int
        number of types
    pot_type: np.ndarray.astype(int)
        type list
    pot_arg: np.ndarray
        argument list
    """
    with open(config_filename,'w') as f:
        f.write("0 %d\n" % n_pot)
        for item in pot_type:
            f.write("%d " % item)
        f.write("\n")
        for item in pot_arg:
            f.write("%.14g " % item)
        f.write("\n")

def printPotTypeArg(pot_name, pot_module, pot_instance, print_front_offset=0, print_long_list=False):
    """ print the petar --type-arg options for a given potential
    Format of --type-arg is 
    type1,type2,...:arg1,arg2,... or type1:arg1-1,arg1-2,...|type2:arg2-1,arg2-2,...
    
    Parameters
    ----------
    pot_name: string
        Name of potential
    pot_module: class object
        The corresponding Potential class of a given potential, if the potential is a list of potential instances, input None
    pot_instance: potential instance
        The potential instance or a list of potential instances

    Return
    ----------
    type_arg: petar --type-arg option argument for the given potential
    npot:     number of potential models
    pot_type: a list of potential types used in the given potential
    pot_arg:  a list of potential argments used in the given potential
    """
    name_format='{:30}'
    if (pot_module == None) & (type(pot_instance) == list):
        print(" "*print_front_offset,  name_format.format(pot_name), "a combination of %d potentials:" % len(pot_instance))
        type_arg=""
        npot=0
        pot_type=np.array([]).astype(int)
        pot_arg=np.array([])
        for pot_sub in pot_instance:
            type_arg_sub, n_sub, type_sub, arg_sub = printPotTypeArg(type(pot_sub).__name__, None, pot_sub, print_front_offset+4, print_long_list)
            if (type_arg != ""):
                type_arg += "|"
            type_arg += type_arg_sub
            npot += n_sub
            pot_type = np.append(pot_type,type_sub)
            pot_arg = np.append(pot_arg,arg_sub)
        print(" "*(print_front_offset+4), name_format.format("Combination: "),type_arg)
        return type_arg, npot, pot_type, pot_arg
    elif (pot_instance != None):
        pot_data = _parse_pot(pot_instance)
        npot, pot_type, pot_args = pot_data[:3]

        if (pot_type.size>0):
            pot_args_str=''
            for i in range(pot_args.size-1):
                pot_args_str+=str(pot_args[i])+','
            if (pot_args.size>0): pot_args_str+=str(pot_args[-1])
            type_arg = ''
            if (pot_type.size==1):
                type_arg = str(pot_type[0])+':'+pot_args_str
            else:
                for i in range(pot_type.size-1):
                    type_arg += str(pot_type[i])+','
                type_arg += str(pot_type[-1])+':'+pot_args_str

        if (pot_args.size>12) & (not print_long_list):
            print(" "*print_front_offset, name_format.format(pot_name), "long argument list, require a configure file")
        else:
            print(" "*print_front_offset, name_format.format(pot_name), type_arg)
        return type_arg, npot, pot_type, pot_args
    return "", 0, np.array([]).astype(int), np.array([])


def printPotTitle():
    name_format='{:30}'
    print(name_format.format("Potential name"), " Options for 'petar --type-arg' with default potential parameters")

def printSpliter():
    print("----------------------------------------------------------------------------------------------")

def listPot():
    printPotTitle()
    printSpliter()
    pot_list=[p for p in dir(galpy.potential)] 
    for pot_name in pot_list:
        pot_module, pot_instance = getPotInstance(pot_name)
        if (pot_instance!=None):
            printPotTypeArg(pot_name, pot_module, pot_instance)

if __name__ == '__main__':

    print_help_flag=False
    config_filename=""

    try:
        import galpy
    except ImportError:
        print("galpy is not found, please check whether it is included in the PYTHONPATH")
        sys.exit(2)

    from galpy.orbit.integrateFullOrbit import _parse_pot
    from galpy.potential.Potential import _check_c

    def usage():
        print("A too to show the description of the galpy options for PeTar and supported Galpy potential")
        print("Usage: petar.galpy.help [potential name]")
        print("       Without potential name, a description to setup --galpy-type-arg and --galpy-conf-file is shown,")
        print("       including the list of supported potentials.")
        print("       With a potential name, the description of the potential is shown.")
        print("Options: ")
        print("    -o [string]: output the default type-arg of the given potential to a configure file;")
        print("                 the argument of this option is the filename")
        print("    -h (--help): help")

    try:
        shortargs = 'o:dhl'
        longargs = ['help']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        kwargs=dict()
        for opt,arg in opts:
            if opt in ('-o'):
                config_filename=arg
            elif opt in ('-h','--help'):
                usage()
                sys.exit(1)

    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)


    if (len(remainder)>0):
        pot_name = remainder[0]
        pot_module, pot_instance = getPotInstance(pot_name)
        if (pot_instance!=None):
            printPotTitle()
            type_arg, n_pot, pot_type, pot_arg=printPotTypeArg(pot_name, pot_module, pot_instance, 0, True)
            printSpliter()
            if (config_filename!=""): 
                print("Save type arguments of %s to file %s." % (pot_name,config_filename))
                savePotTypeArg(config_filename, n_pot, pot_type, pot_arg)
            printSpliter()                                                     
            print("Class definition of %s from Galpy:" % pot_name)
            if (type(pot_instance)==list):
                for pot_sub in pot_instance:
                    print(pot_sub.__doc__)
                    print(pot_sub.__init__.__doc__)
            else:
                print(pot_instance.__doc__)
                print(pot_instance.__init__.__doc__)

            # additional comments
            if (pot_name == 'KeplerPotential'):
                print("Notice that the Kepler potential is generated from the PowerSphericalPotential (type index: 7) with alpha = 3")
                print("Thus, the second argument (alpha) must be 3 ")
        else:
            print(pot_name," is not found in supported Galpy potential list. Available potentials are:")
            listPot()
    else:
        print("Here are the detail of description for options of PeTar-Galpy interface to setup units, galpy potential types and arguments.")
        print("Units: In default, no coversion of units between input potential arguments (or PeTar particle data) and Galpy.")
        print("       Thus, if input/PeTar use astronomical units (Msun, pc, Myr), Galpy units are the same.")
        print("       Notice that mass of potential arguments in Galpy is bound with the gravitational constant (G*M).")
        print("       For example, Mass of 1 Msun corresponds to GM of 0.0044983099795944 pc^3/Myr^2 in Galpy potential arguments.")
        print("--galpy-rscale, --galpy-vscale:")
        print("       Users can manually define conversion factor from input/PeTar units to Galpy units for distance and velocity units.")
        print("       To use the Bovy units of Galpy (distance: 8 kpc, velocity: 220 km/s) and astronomical units in PeTar (pc, pc/Myr), the scale factors are:")
        print("           --galpy-rscale 0.000125 --galpy-vscale 0.004444510098549")
        print("       Other scale factors can be derived:")
        print("           time scale          = 0.028124584538757   [Solar orbital period/2 pi] / Myr")
        print("           acceleration scale  = 0.15802936012884    [Bovy acceleration unit] / [pc/Myr^2]")
        print("           potential pscale    = 1.9753670016104e-05 [Bovy potential unit] / [pc^2/Myr^2]")
        print("           G*M scale           = 2.4692087520131e-09 [Bovy gravitational constant * mass unit] / [pc^3/Myr^2]")
        print("Potential types and arguments: Galpy c interface defines the type indices for supported potentials (listed at the end of this manual).")
        print("       Each potential has mutliple number of arguments. To setup types and arguments, there are three ways: --galpy-set (pre-defined), --galpy-type-arg (simple potentials) and --galpy-conf-file (complex potentials)")
        print("--galpy-set: ")
        print("    Pre-defined potentials, right now only two choices:")
        print("    MWPotential2014 : Milkyway potential described in Galpy paper (Bovy 2015)")
        print("    MWPotentialEvolve:[configure filename] : Evolving Milkyway potential based on cosmological evolution (Gomez 2010)")
        print("        This set requires a configure file, the format:")
        print("            Time [value(Myr)] a [value] H0 [value(km/s/Mpc)]")
        print("            Omega [energy radiation matter]")
        print("            dt_scale [value] dt_max [value(Myr)]")
        print("            Halo [m_vir(Msun) c ac]")
        print("            Disk [M(Msun) ra(pc) rb(pc)]")
        print("            Bulge [rho(Msun/pc^(3-alpha)) alpha rcut(pc)]")
        print("        a (scale factor), H0 (Hubble constant) and Omega are used to obtain a at a given time by integrating adot(t)")
        print("        Omega are three normazied densities in FRW cosmological model")
        print("        dt_scale: scale of time to integrate a, if -1, evolve backwards, if 1, evolve forwards")
        print("        dt_max: maximum time step for integrating adot")
        print("        NFW Halo (3): virial mass, concentration, formation epoch of halo (Wechsler 2002, Zhao 2003)")
        print("        Miyamoto Nagai disk (3): mass, a, b")
        print("        Power Spherical bulge with expotential cutoff (3): rho0, power law index, cutoff radius")
        print("--galpy-type-arg: ")
        print("    This is for a quick set-up of potentials with a string combining types and arguments.")
        print("    The format for the string have two styles: (no space in the middle):")
        print("       1) [type1]:[arg1-1],[arg1-2],**|[type2]:[arg2-1],[arg2-2],**")
        print("       2) [type1],[type2],..:[arg1-1],[arg1-2],[arg2-1],**")
        print("       where '|' split different potential types;")
        print("             ':' separates the type indices and the argument list;")
        print("             ',' separates the items of types or arguments in their lists.")
        print("       For example:")
        print("          1) --galpy-type-arg 15:0.0299946,1.8,0.2375|5:0.7574802,0.375,0.035|9:4.85223053,2.0")
        print("          2) --galpy-type-arg 15,5,9:0.0299946,1.8,0.2375,0.7574802,0.375,0.035,4.85223053,2.0")
        print("          both can generate the MWPotential2014 from Galpy (same to --galpy-set MWPotential2014)")
        print("       Users can use Galpy (_parse_pot function, see its online manual) to find and configure the types and arguments of potentials.")
        print("       The defaults values of types and arguments for supported potentials can be found at the end.")
        print("       The potentials defined in --galpy-type-arg and --galpy-set are both added.")
        print("       Thus, don't repeat the same potential sets in the two options.")
        print("--galpy-conf-file: ")
        print("       The configure file containing the parameters for time- and space-dependent potentials.")
        print("       Users can add, modify and remove potentials at a given time.")
        print("       For each potential, the origin (zero) point can be fixed at the galactic frame or follow the motion of particle system.")
        print("       The potential can be also treated as moving object feeling the force from the other potentials and the particle system.")
        print("       The configure file contain a group of blocks with update times.")
        print("       The format of one block like:")
        print("           Time [value (Galpy unit)] Task [task name]")
        print("           [parameters]")
        print("       Each type of parameter follows its label, like 'Time [value]'.")
        print("       Here 'Time' is the update time of potentials during the simulation.")
        print("       The detail of potentials depends on the following 'Task' and parameters.")
        print("       Be careful for the unit conversion of time and parameters between PeTar and Galpy.")
        print("       Don't rename, delete or modify the configure file before the simulation ends.")
        print("       The format of [parameters] depends on the task name (add, remove, update) shown as following")
        print("")
        print("       For the task 'add':")
        print("               Nset [number]")
        print("               [set 0 content]")
        print("               [set 1 content]")
        print("               ...")
        print("           Nset indicates the number of new potential sets, each set shares the same central position and velocity.")
        print("           Each of the following set content contains a few lines of parameters: ")
        print("               Set [index]")
        print("               Ntype [number] Mode [index]")
        print("               GM [mass] Pos [x y z] Vel [vx vy vz]")
        print("               Type [type1 type2 ...]")
        print("               Arg [arg1-1 arg1-2 ... arg2-1 ...]")
        print("               Nchange [number] Index [Arg_index1, Arg_index2, Arg_index3...]")
        print("               *ChangeMode [mode1, mode2, mode3...]")
        print("               *ChangeRate [rate1, rate2, rate3...]")
        print("           The definition of each parameter:")
        print("               Ntype: number of potential types")
        print("                     If there is no argument, the 4th line (Type ...) is still necessary (empty line).")
        print("               Mode: an index to indicate the reference frame for the position and velocity of the potential center:")
        print("                     0: The galactic frame")
        print("                     1: The particle-system frame; also follow the motion of its center")
        print("                     2: The galactic frame, but move based on the forces from other potentials and from the particle system")
        print("               GM: gravitational constant * mass of potential, used to obtain acceleration from particle-system")
        print("               Pos, Vel: the initial position and velocity of the potential center [Galpy unit] in the reference frame depended on the Mode.")
        print("               Type, Arg: types and arguments of each potential, similar to those in --galpy-type-arg.")
        print("                          The number of values in Type line must be Ntype.")
        print("               The three lines following 'Nchange' are options for evolving potential parameters during the simulation.")
        print("                  Nchange: number of evolving arguments, if 0, the following two lines (marked with '*') are not needed.")
        print("                  Arg_index: the indice of arguments to change during simulation;")
        print("                             counting from 0 based on the potential argument list of the local set ('Arg' line).")
        print("                             if the index is -1, change GM instead of potential arguments.")
        print("                  ChangeMode: the changing mode for each Arg_index: ")
        print("                              1: change the argument linearly;")
        print("                              2: change the argument exponentially.")
        print("                  ChangeRate: the changing rate for each Arg_index, based on ChangeMode:")
        print("                              linearly (1): change value per Galpy time unit;")
        print("                              exponentially (2): coefficient 'a' in the expotential delay or increase form: (exp^{a*t})")
        print("           For example, the configure file of MWPotential2014 with a linearly-mass-decrease Plummer potential following the motion of the particle system:")
        print("                 Time 0.0 Task add                    #[Update at time 0; task is 'add']")
        print("                 Nset 2                               #[2 potential sets]")
        print("                 Set 0                                #[Set 0 (MWPotential2014)]")
        print("                 Ntype 3 Mode 0                       #[3 potential types for MilkyWay; Mode 0 (galactic frame)]")
        print("                 GM 0.0 Pos 0.0 0.0 0.0 Vel 0.0 0.0 0.0 #[GM, x, y, z, vx, vy, vz (galactic frame)]")
        print("                 Type 15 5 9                          #[3 type indices (MWPotential2014)]")
        print("                 Arg 251.63858935563147 1.8 1899.9999999999998 306770418.38588977 3000.0 280.0 1965095308.192175 16000.0 #[All arguments for MWPotential2014]")
        print("                 Nchange 0                            #[No changing arguments)]")
        print("                 Set 1                                #[Set 1 (Plummer)")
        print("                 Ntype 1 Mode 1                       #[1 type; Mode 1 (particle-system frame)")
        print("                 GM 4.498309979594381 Pos 0.0 0.0 0.0 Vel 0.0 0.0 0.0  #[GM, x, y, z, vx, vy, vz (particle-system frame)]")
        print("                 Type 17                              #[Type index for Plummer]")
        print("                 Arg 4.498309979594381 1.0            #[Two arguments for Plummer (1000 Msun, 1 pc)")
        print("                 Nchange 2 Index -1 0                 #[2 Changing argument, index -1 is GM and 0 is the first argument (GM of Plummer)]")
        print("                 ChangeMode 1 1                       #[Linearly change]")
        print("                 ChangeRate -0.2 -0.2                 #[Reduce GM and 0.2 pc^3/Myr^2 per Myr]")
        print("           Here the units are astronomical units, where GM: [pc^3/Myr^2], distance: pc and velocity: pc/Myr;")
        print("           The GM of potential set 1 and the GM of the Plummer (first potential argument in Set 1) are the same.")
        print("           Thus, it is need to change both (index -1 and 0) with the same mode and rate.")
        print("           Notice that the comments after the symbol # here are for reference, don't write them in the configure file.")
        print("       For task 'update':")
        print("           The format is the same as task 'add'.")
        print("           Here Nset is the number of potential sets for update.")
        print("           For parameters of each set, only the definition of Mode has a small change:")
        print("           For update, the mode has four choices instead: 0, 1, 2 and -2:")
        print("               0: The potential center is shifted to the new position and velocity (following line) at the galactic frame;")
        print("               1: The potential center is shifted to the new position and velocity at the particle-system frame;")
        print("               2: The (moving) potential center is shifted to the new position and velocity at the galactic frame;")
        print("              -2: The (moving) potential center keeps the original position and velocity, only update potential arguments.")
        print("           For example, stop the mass decrease for the Plummer model at time 1 referring to the above instance:")
        print("               Time 5.0 Task update")
        print("               Nset 1 Index 1      #[only need to change set 1 (Plummer)]")
        print("               Set 1")
        print("               Ntype 1 Mode 1")
        print("               GM 3.498309979594381 Pos 0.0 0.0 0.0 Vel 0.0 0.0 0.0")
        print("               Type 17")
        print("               Arg  3.498309979594381 1.0")
        print("               Nchange 0           #[Stop changing mass]")
        print("       For task 'remove':")
        print("               Nset [number] Index [set_index1 set_index2 ...]")
        print("           Here Nset is the number of potential sets to be removed, followed by the indices of removing sets.")
        print("           For example, remove the Plummer potential at time 2 following the previous instance:")
        print("               Time 6.0 Task remove")
        print("               Nset 1 Index 1")
        print("Users can use --galpy-set, --galpy-type-arg and --galpy-conf-file together.")
        print("Here are the supported list of Potentials and their default type indices and arguments:")
        listPot()
