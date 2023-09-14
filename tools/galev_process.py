#!/usr/bin/env python3

import numpy as np
import sys
import os
import petar
import getopt

if __name__ == '__main__':

    filename_prefix='data'
    snapshot_format='ascii'
    output_format='ascii'
    data_type='origin'
    mode = 'generate'
    interrupt_mode='none'
    external_mode='none'
    galev_filter = ['Johnson', 'SDSS', 'HST', 'CSST', 'Gaia']
    galev_mode = 'abs_mag'
    FeH = 0

    def usage():
        print("A tool for obtaining snapshots for galev")
        print("Usage: petar.galev.process [options] data_filename")
        print("data_filename: A list of snapshot data path, each line for one snapshot")
        print("option:")
        print("  -h(--help): help")
        print("  -m(--mode): function mode:", mode)
        print("              generate: generate galev input files from snapshots")
        print("              convert: convert format of galev snapshots")
        print("  -d(--data-type): data type of snapshots, multiple types can be combined:",data_type)
        print("       Support type list:")
        print("           origin: original snapshot from petar simulations")
        print("           single: single snapshots from petar.data.process")
        print("           binary: (physical) binary snapshots from petar.data.process")
        print("       For multiple types, combine the names by ',' such as 'origin,single,binary'")
        print("  -s(--snapshot-format): the reading snapshot format for both modes, ascii, binary, npy: ",snapshot_format)
        print("  -i(--interrupt-mode): the interruption mode used in petar, choices: none, base, bse, mobse: ", interrupt_mode)
        print("  -t(--external-mode): external mode used in petar, choices: galpy, none:", external_mode)
        print("  -o(--output-format): output data format for galev snapshot for convert mode:", output_format)
        print("  -f(--filter)       : galev filter list, seperate by comma:",",".join(galev_filter))
        print("  -g(--galev-mode)   : galev mode of values with choices: abs_mag, app_mag, abs_flux, abs_flux:", galev_mode)
        print("  -F(--Fe-H): [Fe/H] value of stars:", FeH)
    try:
        shortargs = 'hm:d:s:i:o:t:F:f:g:'
        longargs = ['help','mode=','data-type=','snapshot-format=','interrupt-mode=','external-mode=','output-format=','Fe-H=','filter=','galev-mode=']
        opts,remainder= getopt.getopt(sys.argv[1:], shortargs, longargs)

        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-m','--mode'):
                mode = arg
            elif opt in ('-d','--data-type'):
                data_type = arg
            elif opt in ('-i','--interrupt-mode'):
                interrupt_mode = arg
            elif opt in ('-t','--external-mode'):
                external_mode = arg
            elif opt in ('-s','--snapshot-format'):
                snapshot_format = arg
            elif opt in ('-o','--output-format'):
                output_format = arg
            elif opt in ('-F','--Fe-H'):
                FeH = float(arg)
            elif opt in ('-f','--filter'):
                galev_filter = [x for x in arg.split(',')]
            elif opt in ('-g','--galev-mode'):
                galev_mode = arg
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError:
        print('getopt error!')
        usage()
        sys.exit(1)

    filename = remainder[0]

    fl = open(filename,'r')
    file_list = fl.read()
    path_list = file_list.splitlines()

    if (mode == 'generate'):
        f = open("snapshot_files.lst","w")

    data_type_list = data_type.split(',')
    for file_path in path_list:
        for dtype in data_type_list:
            if (mode == 'generate'):
                if (dtype == 'origin'):
                    fpath = file_path
                else:
                    fpath = '%s.%s' % (file_path,dtype)
                
                print('generate %s.galev' % fpath)
                if (dtype == 'single') | (dtype == 'origin'):
                    data = petar.Particle(interrupt_mode = interrupt_mode, external_mode = external_mode)
                elif (dtype == 'binary'):
                    data = petar.Binary(member_particle_type=petar.Particle, interrupt_mode = interrupt_mode, external_mode = external_mode)
     
                if (dtype == 'origin'):
                    if (snapshot_format == 'ascii'):
                        data.loadtxt(fpath, skiprows=1)
                    elif (snapshot_format == 'binary') | (snapshot_format == 'npy'):
                        if (external_mode!='none'):
                            data.fromfile(fpath, offset=petar.HEADER_OFFSET_WITH_CM)
                        else:
                            data.fromfile(fpath, offset=petar.HEADER_OFFSET)
                else:
                    if (snapshot_format == 'ascii'):
                        data.loadtxt(fpath)
                    elif (snapshot_format == 'binary'):
                        data.fromfile(fpath)
                    elif (snapshot_format == 'npy'):
                        data.load(fpath+'.npy')
                    else:
                        raise ValueError('Error: snapshot format %s is not supported' % snapshot_format)
     
                read_binary = (dtype == 'binary')
                petar.toGalevSnap(fpath+'.galev', data, FeH, read_binary = read_binary)
     
                f.write(fpath+'.galev\n')
                f.flush()
            else:
                mag = petar.GalevMag(mode = galev_mode, filter = galev_filter)
                if (dtype == 'origin'):
                    fpath = file_path+'.galev.mag'
                else:
                    fpath = file_path+'.'+dtype+'.galev.mag'
                print('convert %s' % fpath)
                if (snapshot_format == 'ascii'):
                    mag.loadtxt(fpath, skiprows=4)
                elif (snapshot_format == 'binary'):
                    mag.fromfile(fpath)
                elif (snapshot_format == 'npy'):
                    mag.load(fpath)
                else:
                    raise ValueError('read format %s is not supported' % read_format)
                if (output_format=='ascii'):
                    mag.savetxt(fpath)
                elif (output_format=='binary'):
                    mag.tofile(fpath)
                elif (output_format=='npy'):
                    mag.save(fpath)
                
    if (mode == 'generate'):
        f.close()
