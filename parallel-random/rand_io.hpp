#pragma once
#include <iostream>
#include <getopt.h>
#include <chrono>
#include "../src/io.hpp"

//! IO parameters manager for random generator
class IOParamsRand{
public:
    IOParamsContainer input_par_store;
    IOParams<long long int> seed;
    IOParams<std::string> seedfile;
    IOParams<std::string> fname_par;

    bool print_flag;

    IOParamsRand(): input_par_store(),
                    seed  (input_par_store, 0,  "rand-seed",   "Random number seed (positive integer); suppressed when --rand-seedfile is provided; if used, the single integer seed is used to generate multiple seeds for each pair of OpenMP thread and MPI processor", "current cpu time"),
                    seedfile(input_par_store, "__NONE__","rand-seedfile","Name for a file contain random seeds of all threads and MPI processors; For restart the simulation, the randseeds file can be used to restore all seeds", "not used"),
                    fname_par(input_par_store, "input.par", "p", "Input parameter file for random seed (this option should be used first before any other options)",NULL,false),
                    print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] print_format_info: if true, print the format information
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const bool print_format_info=true, const int opt_used_pre=0) {
        static int rand_flag=-1;
        const struct option long_options[] = {
            {seed.key,   required_argument, &rand_flag, 1}, 
            {seedfile.key,   required_argument, &rand_flag, 2}, 
            {"help",     no_argument,       0, 'h'},
            {0,0,0,0}
        };


        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-z:p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (rand_flag) {
                case 1:
                    seed.value = atoi(optarg);
                    if(print_flag) seed.print(std::cout);
                    opt_used+=2;
                    break;
                case 2:
                    seedfile.value = optarg;
                    if(print_flag) seedfile.print(std::cout);
                    opt_used+=2;
                    break;
                default:
                    break;
                }
                break;
            case 'p':
                fname_par.value = optarg;
                if(print_flag) {
                    std::string frand_par = fname_par.value+".rand"; 
                    FILE* fpar_in;
                    if( (fpar_in = fopen(frand_par.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", frand_par.c_str());
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
                    std::cout<<"----- Parallel random generator options: -----"<<std::endl;
                    input_par_store.printHelp(std::cout, print_format_info);
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }

        if (seed.value == 0) {
            auto now = std::chrono::high_resolution_clock::now();
            auto duration = now.time_since_epoch();
            auto millis = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
            seed.value = static_cast<long long int>(millis);
        }
            
        if(print_flag) std::cout<<"----- Finish reading input options of Random generator -----\n";
        return opt_used;
    }    
};
