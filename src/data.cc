#include<iostream>
#include<fstream>
#include<cstdio>
#include<unistd.h>
#include<string>
#include<particle_simulator.hpp>
#include"io.hpp"
#include"soft.hpp"
#include"integrate.hpp"
#include"cluster_list.hpp"

class GlobalParams{
public:
    // Energy
    PS::F64 Time;
    PS::S64 N;
    EnergyAndMomemtum energy, energy_error;

    void  writeAscii(FILE* fout) {
        fprintf(fout, "%26.17e %lld ", this->Time, this->N);
        energy.writeAscii(fout);
        energy_error.writeAscii(fout);
        fprintf(fout,"\n");
    }
    
    void writeBinary(FILE* fout) {
        fwrite(&this->Time, sizeof(PS::F64), 2, fout);
        energy.writeBinary(fout);
        energy_error.writeBinary(fout);
    }
};

int main(int argc, char *argv[]){

    // reading parameters
    IOParams<PS::S32> data_format     (1, "Data read(r)/write(w) format BINARY(B)/ASCII(A): r-B/w-A (3); r-A/w-B (2); rw-A (1); rw-B (0)");
    IOParams<std::string> file_list   ("filelist", "reading filelist from file");
    IOParams<std::string> global_file_name ("global.dat", "global parameter data file name");
    
    int c,optcount=0;
    bool reading_flag=false;
    while((c=getopt(argc,argv,"i:p:h")) != -1) {
        switch(c){
        case 'i':
            data_format.value = atoi(optarg);
            std::cerr<<data_format<<std::endl;
            assert(data_format.value>=0||data_format.value<=3);
            optcount +=2;
            break;
        case 'f':
            reading_flag=true;
            file_list.value=optarg;
            std::cerr<<file_list<<std::endl;
            optcount +=2;
            break;
        case 'g':
            global_file_name.value=optarg;
            std::cerr<<global_file_name<<std::endl;
            optcount +=2;
            break;
        case 'h':
            std::cerr<<"Usage: getdata [option] file_list"<<std::endl;
            std::cerr<<"       Option defaulted values are shown\n"<<std::endl;
            std::cerr<<"  -i: [I] "<<data_format<<std::endl;
            std::cerr<<"  -f: [I] "<<file_list<<std::endl;
            std::cerr<<"  -g: [I] "<<global_file_name<<std::endl;
            return 0;
        }
    }
    
    PS::ReallocatableArray<std::string> loop_list;
    if(reading_flag) {
        std::ifstream fin(file_list.value.c_str());
        if (!fin.is_open()) {
            std::cerr<<"Error: Cannot open file "<<file_list.value<<std::endl;
            abort();
        }
        while(true) {
            std::string fname;
            fin>>fname;
            if(fin.eof()) break;
            loop_list.push_back(fname);
        }
    }
    else {
        for(int i=optcount+1; i<argc; i++) {
            loop_list.push_back(argv[i]);
        }
    }

    GlobalParams g0,gnow;

    FILE* fin;
    FILE* gout;

    if ( (gout = fopen(global_file_name.value.c_str(), "w+")) == NULL) {
       std::cerr<<"Error: Cannot open file "<<global_file_name.value<<std::endl;
       abort();
    }

    for(int i=0; i<loop_list.size(); i++) {
        // Reading particles
        std::cout<<"Reading data from "<<loop_list[i].c_str()<<std::endl;
        if ( (fin = fopen(loop_list[i].c_str(),"r")) == NULL) {
            std::cerr<<"Error: Cannot open file "<<loop_list[i]<<std::endl;
            abort();
        }
        FileHeader file_header;
        PS::ReallocatableArray<FPSoft> particle_data;
        
        if(data_format.value==0||data_format.value==3) { // binary format
            file_header.readBinary(fin);
            PS::S32 n_particle = file_header.n_body;
            particle_data.resizeNoInitialize(n_particle);
            for(int i=0; i<n_particle; i++) 
                particle_data[i].readBinary(fin);
        }
        else {
            file_header.readAscii(fin);
            PS::S32 n_particle = file_header.n_body;
            particle_data.resizeNoInitialize(n_particle);
            for(int i=0; i<n_particle; i++) 
                particle_data[i].readAscii(fin);
        }

        fclose(fin);

        // processing
        std::cout<<"Processing data, FileID= "<<file_header.nfile<<", N="<<particle_data.size()<<", Time="<<file_header.time<<std::endl;
        SearchGroup<FPSoft> groups;
        groups.findGroups(particle_data.getPointer(), particle_data.size(), file_header.n_split);
        
        // Correct CM data
        std::cout<<"Correct c.m. velocity, N_groups= "<<groups.getNumOfGroups()<<std::endl;
#ifdef TIDAL_TENSOR
        for (PS::S32 i=0; i<groups.getNumOfGroups(); i++) 
            subtractFcmAndRecoverCMVec(particle_data.getPointer(), groups.getPtclIndex(i), groups.getGroup(i), groups.getGroupN(i), groups.getGroupPertList(i,file_header.n_split));
#endif
//#ifdef HARD_CM_KICK
//        softKickForCM(particle_data.getPointer(), groups.getPtclList(), groups.getNumOfGroups(), groups.getGroupPertList(0,file_header.n_split), 0.5*file_header.dt_soft, file_header.n_split);
//#endif        
        for (PS::S32 i=0; i<groups.getNumOfGroups(); i++) 
            softKickForOneGroup(particle_data.getPointer(), groups.getPtclIndex(i), groups.getGroup(i), groups.getGroupN(i), groups.getGroupPertList(i,file_header.n_split), 0.5*file_header.dt_soft, file_header.n_split);
        
        // Debug
#ifdef DATA_DEBUG
        for(int i=0; i<particle_data.size(); i++) {
            particle_data[i].writeAscii(stdout);
        }
#endif

        // Check energy
        gnow.Time = file_header.time;
        gnow.N = groups.getGroupListSize() + groups.getPtclN() - groups.getNumOfGroups();
        gnow.energy.clear();
        gnow.energy.calc(particle_data.getPointer(), &groups.getPtclList()[groups.getNumOfGroups()], groups.getPtclN()-groups.getNumOfGroups());
        for (PS::S32 i=0; i<groups.getNumOfGroups(); i++) 
            gnow.energy.calc(particle_data.getPointer(), groups.getGroup(i), groups.getGroupN(i));
        if (i==0)  g0 = gnow;
        else gnow.energy_error = gnow.energy - g0.energy;
        
        if(data_format.value==0||data_format.value==2) {
            gnow.writeBinary(gout);
        }
        else {
            gnow.writeAscii(gout);
        }
        
    }
    fclose(gout);
    
    return 0;
}
