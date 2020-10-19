#include <iostream>
#include <particle_simulator.hpp>
#include <getopt.h>
#include "soft_ptcl.hpp"
#include "io.hpp"

typedef PS::ParticleSystem<FPSoft> SystemSoft;

//! redefine readAscii to read group data format, write in I64 format
class FPSoftWriteArtificial: public FPSoft {
public:
    //void writeAscii(FILE* _fout) const{
    //    ParticleBase::writeAscii(_fout);
    //    fprintf(_fout, "%26.17e %lld ", 
    //            this->r_search, this->id);
    //    group_data.artificial.writeAscii(_fout);
    //    changeover.writeAscii(_fout);
    //    fprintf(fp, "%26.17e %26.17e %26.17e %26.17e %26.17e %lld\n", 
    //            this->acc.x, this->acc.y, this->acc.z,  // 9-11
    //            this->pot_tot, this->pot_soft, this->n_ngb);
    //}    

    void readAscii(FILE* _fin) {
        ParticleBase::readAscii(_fin);
        PS::S64 rcount=fscanf(_fin, "%lf %lld ",
                              &this->r_search, &this->id);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
        group_data.artificial.readAscii(_fin);
        changeover.readAscii(_fin);
        rcount=fscanf(_fin, "%lf %lf %lf %lf %lf %lld\n",
                      &this->acc.x, &this->acc.y, &this->acc.z,  // 9-11
                      &this->pot_tot, &this->pot_soft, &this->n_ngb);
        if (rcount<6) {
            std::cerr<<"Error: Data reading fails! requiring data number is 6, only obtain "<<rcount<<".\n";
            abort();
        }
    }
};

typedef PS::ParticleSystem<FPSoftWriteArtificial> SystemSoftWriteArtificial;

int main(int argc, char *argv[]){

    bool b_to_a_flag=false; // If true: BINARY to ASCII; else : ASCII to BINARY
    bool replace_flag=false; // If true: replace the snapshot data; false: create a new file with a suffix of [.B|.A]
    bool group_data_format_f64_flag=false; // When read ASCII mode, if true: treat group_data as 64bit floating (old version); else: treat group_data as 64bit integer
    std::string fname_list("data.snap.lst"); // The filename of a file containing the list of snapshot data pathes

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

    while ((copt = getopt_long(argc, argv, "brgh", long_options, &option_index)) != -1) 
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
        case 'b':
            b_to_a_flag = true;
            if(print_flag) std::cout<<"ASCII to binary\n";
            opt_used ++;
            break;
        case 'r':
            replace_flag = true;
            if(print_flag) std::cout<<"Replace snapshot data\n";
            opt_used ++;
            break;
        case 'g':
            group_data_format_f64_flag = true;
            if(print_flag) std::cout<<"Group data in snapshot uses 64bit floating\n";
            opt_used ++;
            break;
        case 'h':
            if(print_flag){
                std::cout<<"The tool to transfer the format of snapshot data between BINARY and ASCII\n";
                std::cout<<"Usage: petar.format.transfer [option] filelist"<<std::endl;
                std::cout<<"       filelist: "<<fname_list<<std::endl;
                std::cout<<"Options: "<<std::endl
                         <<"   -b  transfer snapshot data format from BINARY to ASCII. If this option is not used, it is ASCII to BINARY"<<std::endl
                         <<"   -r  Replace snapshot data. If this option is not used, new files are created with a suffix of '.B' or '.A'"<<std::endl
                         <<"   -g  Treat group_data as 64bit floating when read ASCII snapshots\n"
                         <<"        This is from the old version before Sep 4, 2020 (GitHub), or -D GROUP_DATA_WRITE_ARTIFICIAL is used in Makefile\n"
                         <<"        After the data transfer, the data will be in new format (64bit Integer)\n"
                         <<"   -h(--help)   print help"<<std::endl;
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
        fname_list =argv[argc-1];
        if(print_flag) std::cout<<"Reading file list: "<<fname_list<<std::endl;
    }

    std::fstream fin;
    fin.open(fname_list,std::fstream::in);
    if(!fin.is_open()) {
        std::cerr<<"Error: data file "<<fname_list<<" cannot be open!\n";
        abort();
    }

    if(print_flag) std::cout<<"----- Finish reading input options -----\n";

    SystemSoft data;
    SystemSoftWriteArtificial data_wa;
    FileHeader file_header;

    while(true) {
        std::string filename;
        fin>>filename;
        if (fin.eof()) break;
        if (print_flag) std::cout<<"Tranfer: "<<filename<<std::endl;
        // Binary to ASCII
        if (b_to_a_flag) {
            data.readParticleBinary(filename.c_str(), file_header);
            if (replace_flag) 
                data.writeParticleAscii(filename.c_str(), file_header);
            else
                data.writeParticleAscii((filename+".A").c_str(), file_header);
        }
        else {
            if (group_data_format_f64_flag) {
                data_wa.readParticleAscii(filename.c_str(), file_header);
                if (replace_flag) 
                    data_wa.writeParticleBinary(filename.c_str(), file_header);
                else
                    data_wa.writeParticleBinary((filename+".B").c_str(), file_header);
            }
            else {
                data.readParticleAscii(filename.c_str(), file_header);
                if (replace_flag) 
                    data.writeParticleBinary(filename.c_str(), file_header);
                else
                    data.writeParticleBinary((filename+".B").c_str(), file_header);
            }
        }
    }

    return 0;
}
