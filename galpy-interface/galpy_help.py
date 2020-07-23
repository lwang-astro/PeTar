#!/usr/bin/env python3
import sys
import getopt
import numpy as np

if __name__ == '__main__':

    print_help_flag=False
    print_long_list_flag=False

    try:
        import galpy
    except ImportError:
        print("galpy is not found, please check whether it is included in the PYTHONPATH")
        sys.exit(2)

    from galpy.orbit.integrateFullOrbit import _parse_pot

    def usage():
        print("A too to show the description of Galpy potential, help to set the parameters of --type-arg in petar commander")
        print("Usage: petar.galpy.help [potential name]")
        print("       this will show the help of the given potential and also the type index for --type-arg option")
        print("       without arguments, all available potential names in Galpy will be listed")
        print("Options: ")
        print("    -d         : call Python help function to show the document of a give potential")
        print("    -l         : print the type-arg for the given potential even it has a long argument list")
        print("    -h (--help): help")

    try:
        shortargs = 'dhl'
        longargs = ['help']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        kwargs=dict()
        for opt,arg in opts:
            if opt in ('-d'):
                print_help_flag=True
            if opt in ('-l'):
                print_long_list_flag=True
            elif opt in ('-h','--help'):
                usage()
                sys.exit(1)

    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)


    def getPotInstance(pot_name):
        """ Try to get an instance of a give pot_name
        Return
        ----------
        pot_module: module object of pot_name, if it is a combined list, return None
        pot_instance: module instance, if it is not 3D or not available, return None

        """
        pot_module = None
        pot_instance=None
        if ('Potential' in pot_name):
            pot_module = galpy.potential.__getattribute__(pot_name)
            if (type(pot_module) == list): 
                pot_instance = pot_module
                pot_module = None
            elif (type(pot_module) == type):
                # get instance
                try:
                    pot_instance = pot_module()
                except (ValueError, TypeError, AttributeError, RuntimeWarning): 
                    pot_module = None
                    pot_instance = None
            else:
                pot_instance = pot_module
                pot_module = type(pot_module)
                
            # remove 2D models
            if (pot_instance != None):
                if (galpy.potential._dim(pot_instance)!=3):
                    pot_instance = None

        return pot_module, pot_instance

    name_format='{:30}'
    def printPotItem(pot_name, pot_module, pot_instance, print_front_offset=0, print_long_list=False):
        if (pot_module == None) & (type(pot_instance) == list):
            print(" "*print_front_offset,  name_format.format(pot_name), "a combination of %d potentials:" % len(pot_instance))
            type_arg=[]
            for pot_sub in pot_instance:
                type_arg.append(printPotItem(type(pot_sub).__name__, None, pot_sub, print_front_offset+4, print_long_list))
            type_arg_all=""
            for i in range(len(type_arg)-1):
                type_arg_all+=type_arg[i]+"|"
            type_arg_all+=type_arg[-1]
            print(" "*(print_front_offset+4), name_format.format("Combination: "),type_arg_all)
        elif (pot_instance != None):
            npot, pot_type, pot_args= _parse_pot(pot_instance)

            if (pot_args.size>10):
                if (not print_long_list):
                    print(" "*print_front_offset, name_format.format(pot_name), "long argument list, require a configure file")
                    return ""
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
                print(" "*print_front_offset, name_format.format(pot_name), type_arg)
                return type_arg
        return ""

    def printPotTitle():
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
                printPotItem(pot_name, pot_module, pot_instance)

    if (len(remainder)>0):
        pot_name = remainder[0]
        pot_module, pot_instance = getPotInstance(pot_name)
        if (pot_instance!=None):
            printPotTitle()
            printPotItem(pot_name, pot_module, pot_instance, 0, print_long_list_flag)
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
            print(pot_name," is not found in Galpy potential list. Available potentials are:")
            listPot()
    else:
        listPot()
