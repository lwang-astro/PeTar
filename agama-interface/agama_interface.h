#pragma once
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <string>
#include <getopt.h>
#include "../src/io.hpp"
#include "../src/astro_units.hpp"
#include "potential_factory.h"

//! IO parameters for Agama manager
class IOParamsAgama{
public:
    IOParamsContainer input_par_store;
    IOParams<std::string> config_filename; // configure file name
    IOParams<double> rscale; 
    IOParams<double> vscale; 

    bool print_flag;

    IOParamsAgama(): input_par_store(),
                     config_filename(input_par_store, "__NONE__", "agama-conf-file", "A configure file of agama potential"),
                     rscale(input_par_store, 1.0, "agama-rscale", "Radius scale factor from unit of the input particle data (IN) to Agama distance unit"),
                     vscale(input_par_store, 1.0, "agama-vscale", "Velocity scale factor from unit of the input particle data (IN) to Agama velocity unit"),
                     print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] print_format_info: print format information
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const bool print_format_info=true, const int opt_used_pre=0) {
        static int agama_flag=-1;
        const struct option long_options[] = {
            {config_filename.key, required_argument, &agama_flag, 0}, 
            {rscale.key,     required_argument, &agama_flag, 1}, 
            {vscale.key,     required_argument, &agama_flag, 2}, 
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };
        
        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        std::string fname_par;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (agama_flag) {
                case 0:
                    config_filename.value = optarg;
                    if(print_flag) config_filename.print(std::cout);
                    opt_used+=2;
                    break;
                case 1:
                    rscale.value = atof(optarg);
                    if(print_flag) rscale.print(std::cout);
                    opt_used+=2;
                    break;
                case 2:
                    vscale.value = atof(optarg);
                    if(print_flag) vscale.print(std::cout);
                    opt_used+=2;
                    break;
                default:
                    break;
                }
                break;
            case 'p':
                fname_par = optarg;
                if(print_flag) {
                    std::string fagama_par = fname_par+".agama"; 
                    FILE* fpar_in;
                    if( (fpar_in = fopen(fagama_par.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", fagama_par.c_str());
                        abort();
                    }
                    input_par_store.readAscii(fpar_in);
                    fclose(fpar_in);
                }
                opt_used+=2;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
                input_par_store.mpi_broadcast();
                PS::Comm::barrier();
#endif
                break;
            case 'h':
                if(print_flag){
                    std::cout<<"----- Agama options: -----"<<std::endl;
                    input_par_store.printHelp(std::cout, print_format_info);
                    std::cout<<std::endl;
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }
        
        if(print_flag) std::cout<<"----- Finish reading input options of Agama -----\n";

        return opt_used;
    }

};

//! A class to manager the API to Agama
class AgamaManager{
public:
    double rscale;
    double vscale;
    double tscale;
    double fscale;
    double pscale;
    double gmscale;
    std::string conf_fname;
    potential::PtrPotential agama_potential;

    AgamaManager(): rscale(1.0), 
                    vscale(1.0), 
                    tscale(1.0), 
                    fscale(1.0), 
                    pscale(1.0), 
                    gmscale(1.0), 
                    conf_fname("__NONE__"),
                    agama_potential() {}

    //! print reference to cite
    static void printReference(std::ostream & fout, const int offset=4) {
        for (int i=0; i<offset; i++) fout<<" ";
            fout<<"Agama: Vasiliev E., 2019, MNRAS, 482, 1525"<<std::endl;
    }    
    //！

    //! initialization function
    /*!
        @param[in] _input: input parameters
        @param[in] _time: current system time
        @param[in] _print_flag: if true, printing information to std::cout
    */
    void initial(const IOParamsAgama& _input, const double _time, const bool _print_flag=false) {
        // unit scale
        rscale = _input.rscale.value;
        vscale = _input.vscale.value;
        tscale = rscale/vscale;
        fscale = vscale*vscale/rscale;
        pscale = vscale*vscale;
        gmscale = pscale*rscale; 
        conf_fname = _input.config_filename.value;

        agama_potential = potential::readPotential(conf_fname.c_str());
        //    units::ExternalUnits(
        //        units::InternalUnits(units::Kpc, units::Kpc/units::kms),
        //        units::Kpc, units::Kpc/units::kms, units::Msun*input_parameters.gravitational_constant.value));
    }

    //! calculate acceleration and potential at give position
    /*!
      @param[out] acc: [3] acceleration to return
      @param[out] pot: potential to return 
      @param[in] _time: time in input unit
      @param[in] gm: G*mass of particles [input unit]
      @param[in] pos_g: position of particles in the galactic frame [input unit]
      @param[in] pos_l: position of particles in the particle system frame [input unit]
     */
    void calcAccPot(double* acc, double& pot, const double _time, const double gm, const double* pos_g, const double* pos_l) {    
        pot = 0;
        acc[0] = acc[1] = acc[2] = 0.0;

        coord::GradCar grad;
        agama_potential->eval(coord::PosCar(pos_g[0]*rscale, pos_g[1]*rscale, pos_g[2]*rscale), &pot, &grad, NULL, _time*tscale);
        acc[0] = -grad.dx/fscale;
        acc[1] = -grad.dy/fscale;
        acc[2] = -grad.dz/fscale;
        pot /= pscale;
    }

};