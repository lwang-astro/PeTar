#include <iostream>
#include <particle_simulator.hpp>
#include <getopt.h>
#include "soft_ptcl.hpp"
#include "io.hpp"
#include "status.hpp"

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

//! old version without c.m. pos and vel offset
class FileHeaderNoOffset{
public:
    long long int nfile;  // file id
    long long int n_body;
    double time;
    FileHeaderNoOffset(){
        n_body = 0;
        time = 0.0;
    }
    FileHeaderNoOffset(const long long int ni, const long long int n, const double t) {
        nfile = ni;
        n_body = n;
        time = t;
    }

    int readAscii(FILE * fp){
        int rcount=fscanf(fp, "%lld %lld %lf\n", &nfile, &n_body, &time);
        if (rcount<3) {
          std::cerr<<"Error: cannot read header, please check your data file header!\n";
          abort();
        }
        //std::cout<<"Number of particles ="<<n_body<<";  Time="<<time<<std::endl;
        return n_body;
    }

    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %lld %26.17e\n", nfile, n_body, time);
    }

    int readBinary(FILE* fp){
        size_t rcount=fread(this, sizeof(FileHeader), 1, fp);
        if(rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is "<<1<<" bytes, only obtain "<<rcount<<" bytes.\n";
            abort();
        }
        //std::cout<<"Number of particles ="<<n_body<<";  Time="<<time<<std::endl;
        return n_body;
    }

    void writeBinary(FILE* fp) const{
        fwrite(this, sizeof(FileHeader), 1, fp);
    }
};


typedef PS::ParticleSystem<FPSoftWriteArtificial> SystemSoftWriteArtificial;

int main(int argc, char *argv[]){

    bool b_to_a_flag=false; // If true: BINARY to ASCII; else : ASCII to BINARY
    bool replace_flag=false; // If true: replace the snapshot data; false: create a new file with a suffix of [.B|.A]
    bool group_data_format_f64_flag=false; // When read ASCII mode, if true: treat group_data as 64bit floating (old version); else: treat group_data as 64bit integer
    bool read_one_file_flag=false; // If true: the input file is not a list but the filename of one snapshot
    bool write_ascii_header_flag=false; // If true: only output ascii header
    bool add_record_cm_flag=false; // If true; substract center to header
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

    while ((copt = getopt_long(argc, argv, "brgcfHh", long_options, &option_index)) != -1) 
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
            if(print_flag) std::cout<<"BINARY to ASCII\n";
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
        case 'c':
            add_record_cm_flag = true;
            if(print_flag) std::cout<<"Subtract particle center (no-mass-weighted averaged position and velocity) into header\n";
            opt_used ++;
            break;
        case 'f':
            read_one_file_flag = true;
            if(print_flag) std::cout<<"Read one file instead of a filelist\n";
            opt_used ++;
            break;
        case 'H':
            write_ascii_header_flag = true;
            if(print_flag) std::cout<<"Write only header data of snapshot in ASCII format\n";
            opt_used ++;
            break;
        case 'h':
            if(print_flag){
                std::cout<<"The tool to transfer the format of snapshot data between BINARY and ASCII\n";
                std::cout<<"Usage: petar.format.transfer [option] filelist"<<std::endl;
                std::cout<<"       filelist: "<<fname_list<<std::endl;
                std::cout<<"Stellar evolution method: ";
#ifdef STELLAR_EVOLUTION
#ifdef BSE
                std::cout<<"BSE\n";
#elif MOBSE
                std::cout<<"MOBSE\n";
#else
                std::cout<<"Base\n";
#endif
#else
                std::cout<<"None\n";
#endif
#ifdef EXTERNAL_POT_IN_PTCL
                std::cout<<"External potential column exists\n";
#else
                std::cout<<"External potential column not exists\n";
#endif
                std::cout<<"Important: Ensure that the stellar evolution method and external mode used in the snapshots and this tool are consistent.\n"
                         <<"           If the replace option (-r) is used and the methods are not consistent, the data cannot be recovered!"<<std::endl;
                std::cout<<"Options: "<<std::endl
                         <<"   -b  transfer snapshot data format from BINARY to ASCII. If this option is not used, it is ASCII to BINARY"<<std::endl
                         <<"   -r  Replace snapshot data. If this option is not used, new files are created with a suffix of '.B' or '.A'"<<std::endl
                         <<"   -g  Treat group_data as 64bit floating when read ASCII snapshots\n"
                         <<"        This is from the old version before Sep 4, 2020 (GitHub), or -D GROUP_DATA_WRITE_ARTIFICIAL is used in Makefile\n"
                         <<"        After the data transfer, the data will be in new format (64bit Integer)\n"
                         <<"   -c  Substract particle c.m. position and velocity into header\n"
                         <<"        This is used to transfer data before Dec 18, 2020 (GitHub) to new version when external-mode=galpy\n"
                         <<"   -f  Read one snapshot instead of a list, the filelist should be replaced by the filename of the snapshot\n"
                         <<"   -H  Write only header data of snapshots in ascii format with suffix '.H'. This option suppress '-r' so that data are not replaced"<<std::endl
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

    if(print_flag) std::cout<<"----- Finish reading input options -----\n";

    SystemSoft data;
    SystemSoftWriteArtificial data_wa;
    FileHeader file_header;
    FileHeaderNoOffset file_header_no_offset;
    Status status;

    auto transferOneFile = [&] (const std::string& filename) {
        // Binary to ASCII
        if (b_to_a_flag) {
            if (add_record_cm_flag) {
                data.readParticleBinary(filename.c_str(), file_header_no_offset);
                status.calcAndShiftCenterOfMass(&data[0],data.getNumberOfParticleLocal(), 3, true);
                file_header.nfile  = file_header_no_offset.nfile;
                file_header.n_body = file_header_no_offset.n_body;
                file_header.time   = file_header_no_offset.time;
#ifdef RECORD_CM_IN_HEADER
                file_header.pos_offset = status.pcm.pos;
                file_header.vel_offset = status.pcm.vel;
#endif
            }
            else {
                FILE* fin;
                if( (fin = fopen(filename.c_str(),"r")) == NULL) {
                    std::cerr<<"Error: Cannot open file "<<filename<<"!\n";
                    abort();
                }
                file_header.readBinary(fin);
                if (!write_ascii_header_flag) {
                    data.setNumberOfParticleLocal(file_header.n_body);
                    for (int i=0; i<file_header.n_body; i++) 
                        data[i].readBinary(fin);
                }
                //data.readParticleBinary(filename.c_str(), file_header);
            }
            if (write_ascii_header_flag) {
                FILE* fout;
                std::string fname_header = filename+".H";
                if( (fout = fopen(fname_header.c_str(),"w")) == NULL) {
                    std::cerr<<"Error: Cannot open file "<<fname_header<<"!\n";
                    abort();
                }
                file_header.writeAscii(fout);
            }
            else {
                if (replace_flag) 
                    data.writeParticleAscii(filename.c_str(), file_header);
                else
                    data.writeParticleAscii((filename+".A").c_str(), file_header);
            }
        }
        else {
            if (group_data_format_f64_flag) {
                if (add_record_cm_flag) {
                    data_wa.readParticleAscii(filename.c_str(), file_header_no_offset);
                    status.calcAndShiftCenterOfMass(&data_wa[0], data_wa.getNumberOfParticleLocal(), 3, true);
                    file_header.nfile  = file_header_no_offset.nfile;
                    file_header.n_body = file_header_no_offset.n_body;
                    file_header.time   = file_header_no_offset.time;
#ifdef RECORD_CM_IN_HEADER
                    file_header.pos_offset = status.pcm.pos;
                    file_header.vel_offset = status.pcm.vel;
#endif
                }
                else  data_wa.readParticleAscii(filename.c_str(), file_header);
                if (replace_flag) 
                    data_wa.writeParticleBinary(filename.c_str(), file_header);
                else
                    data_wa.writeParticleBinary((filename+".B").c_str(), file_header);
            }
            else {
                if (add_record_cm_flag) {
                    data.readParticleAscii(filename.c_str(), file_header_no_offset);
                    status.calcAndShiftCenterOfMass(&data[0], data.getNumberOfParticleLocal(), 3, true);
                    file_header.nfile  = file_header_no_offset.nfile;
                    file_header.n_body = file_header_no_offset.n_body;
                    file_header.time   = file_header_no_offset.time;
#ifdef RECORD_CM_IN_HEADER
                    file_header.pos_offset = status.pcm.pos;
                    file_header.vel_offset = status.pcm.vel;
#endif
                }
                else  data.readParticleAscii(filename.c_str(), file_header);
                if (replace_flag) 
                    data.writeParticleBinary(filename.c_str(), file_header);
                else
                    data.writeParticleBinary((filename+".B").c_str(), file_header);
            }
        }
    };

    if (read_one_file_flag) transferOneFile(std::string(fname_list));
    else {
        std::fstream fin;
        fin.open(fname_list,std::fstream::in);
        if(!fin.is_open()) {
            std::cerr<<"Error: data file "<<fname_list<<" cannot be open!\n";
            abort();
        }

        while(true) {
            std::string filename;
            fin>>filename;
            if (fin.eof()) break;
            if (print_flag) std::cout<<"Tranfer: "<<filename<<std::endl;
            transferOneFile(filename);
        }
    }
    return 0;
}
