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
    pot_instance=None
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
        npot, pot_type, pot_args= _parse_pot(pot_instance)

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
        print("A too to show the description of Galpy potential, help to set the parameters of --type-arg in petar commander")
        print("Usage: petar.galpy.help [potential name]")
        print("       This will show the help of the given potential and also the type index for --type-arg option.")
        print("       Without a potential name, all supported potential names in Galpy will be listed.")
        print("       For a potential requiring a long argument list, the type and arguments are not printed, use '-l' to show that.")
        print("       it is suggested to use a configure file to read arguments.")
        print("Options: ")
        print("    -d         : call Python help function to show the document of a give potential")
        print("                 This option only works when a potential name is provided.")
        print("    -o [string]: Output the type-arg to a configure file; the argument of this option is the filename")
        print("    -h (--help): help")

    try:
        shortargs = 'o:dhl'
        longargs = ['help']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        kwargs=dict()
        for opt,arg in opts:
            if opt in ('-d'):
                print_help_flag=True
            elif opt in ('-o'):
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
            type_arg, n_pot, pot_type, pot_arg=printPotTypeArg(pot_name, pot_module, pot_instance, 0)
            printSpliter()
            if (config_filename!=""): 
                print("Save type arguments of %s to file %s." % (pot_name,config_filename))
                savePotTypeArg(config_filename, n_pot, pot_type, pot_arg)
            printSpliter()                                                     
            if (print_help_flag):
                print("Class definition of %s from Galpy:" % pot_name)
                if (type(pot_instance)==list):
                    for pot_sub in pot_instance:
                        print(pot_sub.__doc__)
                        print(pot_sub.__init__.__doc__)
                else:
                    print(pot_instance.__doc__)
                    print(pot_instance.__init__.__doc__)
        else:
            print(pot_name," is not found in supported Galpy potential list. Available potentials are:")
            listPot()
    else:
        listPot()
