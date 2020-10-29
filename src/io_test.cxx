#include <iostream>
#include <cstdio>
#include<particle_simulator.hpp>
#include "io.hpp"
#include <getopt.h>

int main(int argc, char** argv) {

    IOParamsContainer input_par_store;
    IOParams<PS::S64> alpha   (input_par_store, 1, "alpha", "Alpha parameter");
    IOParams<PS::S64> beta    (input_par_store, 2, "beta",  "Beta parameter");
    IOParams<PS::F64> gamma_one  (input_par_store, 3.0, "gamma-one", "Gamma 1 parameter");
    IOParams<PS::F64> gamma_two  (input_par_store, 3.0, "gamma-two", "Gamma 2 parameter");
    IOParams<std::string> fname_par  (input_par_store, "io.test.par", "p","IO test parameter file name");    

    static int long_flag=-1;
    static struct option long_options[] = {
        //{"group-data-format",     required_argument, &long_flag, 0},
        {"alpha",      required_argument, &long_flag, 0},
        {"beta",       required_argument, &long_flag, 1},
        {"gamma-one",  required_argument, &long_flag, 2},
        {"gamma-two",  required_argument, &long_flag, 3},
        {"help",             no_argument, 0, 'h'},        
        {0,0,0,0}
    };
    
    int opt_used = 0;
    int copt;
    int option_index;
    optind = 0; // reset getopt
    bool print_flag = true;

    while ((copt = getopt_long(argc, argv, "p:h", long_options, &option_index)) != -1) 
        switch (copt) {
        case 0:
            switch (long_flag) {
            case 0:
                alpha.value = atoi(optarg);
                if(print_flag) alpha.print(std::cout);
                opt_used += 2;
                break;
            case 1:
                beta.value = atoi(optarg);
                if(print_flag) beta.print(std::cout);
                opt_used += 2;
                break;
            case 2:
                gamma_one.value = atof(optarg);
                if(print_flag) gamma_one.print(std::cout);
                opt_used += 2;
                break;
            case 3:
                gamma_two.value = atof(optarg);
                if(print_flag) gamma_two.print(std::cout);
                opt_used += 2;
                break;
            default:
                break;
            }
            break;
        case 'p':
            fname_par.value = optarg;
            if(print_flag) {
                fname_par.print(std::cout);
                FILE* fpar_in;
                if( (fpar_in = fopen(fname_par.value.c_str(),"r")) == NULL) {
                    fprintf(stderr,"Error: Cannot open file %s.\n", fname_par.value.c_str());
                    abort();
                }
                input_par_store.readAscii(fpar_in);
                fclose(fpar_in);
            }
            opt_used += 2;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
            input_par_store.mpi_broadcast();
            PS::Comm::barrier();
#endif
            break;
        case 'h':
            if(print_flag){
                std::cout<<"Usage: petar.io.test [option]\n"
                         <<"Options:"<<std::endl;
                input_par_store.printHelp(std::cout, 2,8,20);
                std::cout<<" -h(--help):   print help"<<std::endl;
            }
            return -1;
        case '?':
            opt_used +=2;
            break;
        default:
            break;
        }
    
    if (print_flag) std::cout<<"Save input parameters to file "<<fname_par.value<<std::endl;
    FILE* fpar_out;
    if( (fpar_out = fopen(fname_par.value.c_str(),"w")) == NULL) {
        fprintf(stderr,"Error: Cannot open file %s.\n", fname_par.value.c_str());
        abort();
    }
    input_par_store.writeAscii(fpar_out);
    fclose(fpar_out);
    
    std::cout<<alpha<<std::endl;
    std::cout<<beta<<std::endl;
    std::cout<<gamma_one<<std::endl;
    std::cout<<gamma_two<<std::endl;
    std::cout<<fname_par<<std::endl;

    std::cout<<"Stored key: values\n";
    input_par_store.print(std::cout);

    return 0;
}
