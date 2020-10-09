#include <iostream>
#include <particle_simulator.hpp>
#include <getopt.h>
#include "soft_ptcl.hpp"
#include "io.hpp"

typedef PS::ParticleSystem<FPSoft> SystemSoft;

int main(int argc, char *argv[]){

    // IO parameters
    IOParamsContainer input_par_store;
    IOParams<PS::S32> direction         (input_par_store, 0, "If 0: BINARY to ASCII; 1 : ASCII to BINARY");
    IOParams<PS::S32> write_mode        (input_par_store, 0, "If 0: replace the snapshot data; 1: create a new file with a suffix of [.B|.A]");
    //IOParams<PS::S32> group_data_format (input_par_store, 0, "When read ASCII mode, if 0: treat group_data as 64bit Integer; 1: treat group_data as artificial particles (old version)");
    IOParams<std::string> fname_list    (input_par_store, "data.snap.lst","The filename of a file containing the list of snapshot data pathes");

    static int long_flag=-1;
    static struct option long_options[] = {
        //{"group-data-format",     required_argument, &long_flag, 0},
        {"help",                  no_argument, 0, 'h'},        
        {0,0,0,0}
    };
        
    int opt_used = 0;
    int copt;
    int option_index;
    optind = 0; // reset getopt
    bool print_flag = true;

    while ((copt = getopt_long(argc, argv, "-i:w:h", long_options, &option_index)) != -1) 
        switch (copt) {
        case 0:
            switch (long_flag) {
            //case 0:
            //    group_data_format.value = atoi(optarg);
            //    if(print_flag) group_data_format.print(std::cout);
            //    opt_used += 2;
            //    assert(group_data_format.value<=1&&group_data_format.value>=0);
            //    break;
            default:
                break;
            }
            break;
        case 'i':
            direction.value = atoi(optarg);
            if(print_flag) direction.print(std::cout);
            opt_used += 2;
            assert(direction.value>=0&&direction.value<=1);
            break;
        case 'w':
            write_mode.value = atoi(optarg);
            if(print_flag) write_mode.print(std::cout);
            opt_used += 2;
            assert(write_mode.value>=0&&write_mode.value<=1);
            break;
        case 'h':
            if(print_flag){
                std::cout<<"Usage: petar.format.transfer [option] filelist"<<std::endl;
                std::cout<<"       filelist: "<<fname_list<<std::endl;
                std::cout<<"Options:  defaulted values are shown after ':'"<<std::endl
                         <<"   -i:  [I] "<<direction<<std::endl
                         <<"   -w:  [I] "<<write_mode<<std::endl
                    //<<"       --group-data-format:  [I] "<<group_data_format<<std::endl
                         <<"   -h(--help):               print help"<<std::endl;
            }
            return -1;
        case '?':
            opt_used +=2;
            break;
        default:
            break;
        }
    
    // count used options
    opt_used ++;
    //std::cout<<"Opt used:"<<opt_used<<std::endl;
    if (opt_used<argc) {
        fname_list.value =argv[argc-1];
        if(print_flag) std::cout<<"Reading file list: "<<fname_list.value<<std::endl;
    }

    std::fstream fin;
    fin.open(fname_list.value,std::fstream::in);
    if(!fin.is_open()) {
        std::cerr<<"Error: data file "<<fname_list.value<<" cannot be open!\n";
        abort();
    }

    if(print_flag) std::cout<<"----- Finish reading input options -----\n";

    SystemSoft data;
    FileHeader file_header;

    while(true) {
        std::string filename;
        fin>>filename;
        if (fin.eof()) break;
        if (print_flag) std::cout<<"Tranfer: "<<filename<<std::endl;
        // Binary to ASCII
        if (direction.value==0) {
            data.readParticleBinary(filename.c_str(), file_header);
            if (write_mode.value==0) 
                data.writeParticleAscii(filename.c_str(), file_header);
            else
                data.writeParticleAscii((filename+".A").c_str(), file_header);
        }
        else {
            data.readParticleAscii(filename.c_str(), file_header);
            if (write_mode.value==0) 
                data.writeParticleBinary(filename.c_str(), file_header);
            else
                data.writeParticleBinary((filename+".B").c_str(), file_header);
        }
    }

    return 0;
}
