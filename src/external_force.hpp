#pragma once

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <getopt.h>
#include "io.hpp"
#include "Common/Float.h"

//! IO parameters manager for external perturbation in hard integration
/*! For initializing the COMMON block variables from the commander option.
  The description of each parameter is also provided.
 */
class IOParamsExternalHard{
public:
    IOParamsContainer input_par_store;
    IOParams<long long int> mode; // option to switch perturbation
    IOParams<double> alpha; // coefficient of drag force

    bool print_flag;
    IOParamsExternalHard(): input_par_store(),
                            mode  (input_par_store, 0, "external-hard-mode", "external hard mode: 0, not used; 1, simple velocity depended friction"),
                            alpha (input_par_store, 1.0, "external-hard-alpha",  "coefficient of velocity depended friction"),
                            print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const int opt_used_pre=0) {
        static int ext_flag=-1;
        const struct option long_options[] = {
            {mode.key,    required_argument, &ext_flag, 0},  
            {alpha.key,   required_argument, &ext_flag, 1},  
            {"help",      no_argument,       0, 'h'},
            {0,0,0,0}
        };

        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        std::string fname_par;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-z:p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (ext_flag) {
                case 0:
                    mode.value = atoi(optarg);
                    if(print_flag) mode.print(std::cout);
                    opt_used+=2;
                    break;            
                case 1:
                    alpha.value = atof(optarg);
                    if(print_flag) alpha.print(std::cout);
                    opt_used+=2;
                    break;            
                default:
                    break;
                }
                break;
            case 'p':
                fname_par = optarg;
                if(print_flag) {
                    std::string fgalpy_par = fname_par+".exthard"; 
                    FILE* fpar_in;
                    if( (fpar_in = fopen(fgalpy_par.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", fgalpy_par.c_str());
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
                    std::cout<<"External perturbation for hard integration, options:"<<std::endl;
                    input_par_store.printHelp(std::cout, 2, 10, 23);
                    std::cout<<std::endl;
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }
        
        if(print_flag) std::cout<<"----- Finish reading input options of external perturbation for hard integration -----\n";

        return opt_used;
    }    
};

class ExternalHardForce{
public:
    bool is_used; // indicator whether external force is used
    Float alpha; // coefficient of drag force

    ExternalHardForce(): is_used(false), alpha(0.0) {}

    //! initial parameters for perturbation
    void initial(const IOParamsExternalHard& _input, const bool _print_flag=false) {
        is_used = bool(_input.mode.value);
        alpha = _input.alpha.value;
    }

    //! External force for one particle in hard part
    /*!
      @param[out] _force: acceleration and jerk
      @param[in] _particle: particle data
      
      Return: the next integration time step
    */
    template<class Tf, class Tp> 
    Float calcAccJerkExternal(Tf& _force, const Tp& _particle){
        
        _force.acc0[0] += -alpha*_particle.vel[0];
        _force.acc0[1] += -alpha*_particle.vel[1];
        _force.acc0[2] += -alpha*_particle.vel[2];
        _force.acc1[0] += _force.acc0[0];
        _force.acc1[1] += _force.acc0[1];
        _force.acc1[2] += _force.acc0[2];
        
        return NUMERIC_FLOAT_MAX;
    }

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        return true;
    }    

    //! print parameters
    void print(std::ostream & _fout) const{
    }    

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(this, sizeof(*this),1,_fp);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }        
};
