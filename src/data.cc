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
//        fprintf(fout,"\n");
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
    IOParams<PS::S32> snp_style       (1, "Output snapshot style: multiple systems are splitted to individual files (0); All data in one files (1)");
    
    int c,optcount=0;
    bool reading_flag=false;
    while((c=getopt(argc,argv,"i:f:g:h")) != -1) {
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
        case 's':
            snp_style.value = atoi(optarg);
            std::cerr<<snp_style<<std::endl;
            optcount +=2;
            break;
        case 'h':
            std::cerr<<"Usage: getdata [option] file_list"<<std::endl;
            std::cerr<<"       Option defaulted values are shown\n"<<std::endl;
            std::cerr<<"  -i: [I] "<<data_format<<std::endl;
            std::cerr<<"  -f: [I] "<<file_list<<std::endl;
            std::cerr<<"  -g: [I] "<<global_file_name<<std::endl;
            std::cerr<<"  -s: [I] "<<snp_style<<std::endl;
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

    if ( (gout = fopen(global_file_name.value.c_str(), "a+")) == NULL) {
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
//#ifdef DATA_DEBUG
//        for(int i=0; i<particle_data.size(); i++) {
//            particle_data[i].writeAscii(stdout);
//        }
//#endif

        // Check energy
        gnow.Time = file_header.time;
        gnow.N = groups.getGroupListSize() + groups.getPtclN() - groups.getNumOfGroups();
        gnow.energy.clear();
        gnow.energy.calc(particle_data.getPointer(), &groups.getPtclList()[groups.getNumOfGroups()], groups.getPtclN()-groups.getNumOfGroups());
        for (PS::S32 i=0; i<groups.getNumOfGroups(); i++) 
            gnow.energy.calc(particle_data.getPointer(), groups.getGroup(i), groups.getGroupN(i));
        if (i==0)  g0 = gnow;
        else gnow.energy_error = gnow.energy - g0.energy;
        gnow.energy_error.relative(g0.energy);
        
        if(data_format.value==0||data_format.value==2) {
            gnow.writeBinary(gout);
#ifdef DATA_DEBUG
            for(int i=0; i<groups.getNumOfGroups(); i++) {
                const PS::S32* glist=groups.getGroup(i);
                for(int j=0; j<groups.getGroupN(i); j++) {
                    particle_data[glist[j]].ParticleBase::writeBinary(gout);
                }
            }
            for(int i=groups.getNumOfGroups(); i<groups.getPtclN(); i++) {
                particle_data[groups.getPtclIndex(i)].ParticleBase::writeBinary(gout);
            }
#endif
            if (snp_style.value==0) { // individual files
                FILE *dout;
                std::string new_file_name = "snp.M1."+std::to_string(file_header.nfile);
                if ( (dout = fopen(new_file_name.c_str(),"w")) == NULL) {
                    std::cerr<<"Error: Cannot open file "<<new_file_name<<std::endl;
                    abort();
                }
                PS::S32 n_group_max = 0;
                PS::S32 n_groups = groups.getNumOfGroups();
                for(PS::S32 i=0; i<groups.getPtclN(); i++) {
                    FPSoft &pi = particle_data[groups.getPtclIndex(i)];
                    pi.Ptcl::writeBinary(dout);
                    fwrite(&pi.pot_tot,sizeof(PS::F64),1,dout);
                    if(i<n_groups && n_group_max<pi.status) n_group_max = pi.status;
                }
                fclose(dout);
                
                PS::ReallocatableArray<FILE*> mout;
                mout.resizeNoInitialize(n_group_max);
                for(PS::S32 i=2; i<=n_group_max; i++) {
                    std::string new_file_name = "snp.M"+std::to_string(i)+"."+std::to_string(file_header.nfile);
                    if ( (mout[i-2] = fopen(new_file_name.c_str(),"w")) == NULL) {
                        std::cerr<<"Error: Cannot open file "<<new_file_name<<std::endl;
                        abort();
                    }
                }
                for(PS::S32 i=0; i<n_groups; i++) {
                    const PS::S32 ni_group = groups.getGroupN(i);
                    const PS::S32 fi = ni_group - 2;

                    FPSoft &pi = particle_data[groups.getPtclIndex(i)];
                    pi.Ptcl::writeBinary(mout[fi]);
                    fwrite(&pi.pot_tot,sizeof(PS::F64),1,mout[fi]);

                    const PS::S32* glist=groups.getGroup(i);
                    for(PS::S32 j=0; j<ni_group; j++) {
                        FPSoft &pj = particle_data[glist[j]];
                        pj.Ptcl::writeBinary(mout[fi]);
                        fwrite(&pj.pot_tot,sizeof(PS::F64),1,mout[fi]);
                    }
                }
                for (PS::S32 i=0;i<n_group_max-1;i++) fclose(mout[i]);
            }
            else {
                FILE *dout;
                std::string new_file_name = "snp."+std::to_string(file_header.nfile);
                if ( (dout = fopen(new_file_name.c_str(),"w")) == NULL) {
                    std::cerr<<"Error: Cannot open file "<<new_file_name<<std::endl;
                    abort();
                }
                for(int i=0; i<groups.getNumOfGroups(); i++) {
                    const PS::S32* glist=groups.getGroup(i);
                    for(int j=0; j<groups.getGroupN(i); j++) {
                        FPSoft &pj = particle_data[glist[j]];
                        pj.Ptcl::writeBinary(dout);
                        fwrite(&pj.pot_tot,sizeof(PS::F64),1,dout);
                    }
                }
                for(int i=groups.getNumOfGroups(); i<groups.getPtclN(); i++) {
                    FPSoft &pi = particle_data[groups.getPtclIndex(i)];
                    pi.Ptcl::writeBinary(dout);
                    fwrite(&pi.pot_tot,sizeof(PS::F64),1,dout);
                }

                fclose(dout);
            }
            
        }
        else {
            gnow.writeAscii(gout);
#ifdef DATA_DEBUG
            for(int i=0; i<gnow.N; i++) 
                particle_data[i].ParticleBase::writeAscii(gout);
//            for(int i=0; i<groups.getNumOfGroups(); i++) {
//                const PS::S32* glist=groups.getGroup(i);
//                for(int j=0; j<groups.getGroupN(i); j++) {
//                    particle_data[glist[j]].ParticleBase::writeAscii(gout);
//                }
//            }
//            for(int i=groups.getNumOfGroups(); i<groups.getPtclN(); i++) {
//                particle_data[groups.getPtclIndex(i)].ParticleBase::writeAscii(gout);
//            }
#endif
            fprintf(gout,"\n");
            if (snp_style.value==0) { // individual files
                FILE *dout;
                std::string new_file_name = "snp.M1."+std::to_string(file_header.nfile);
                if ( (dout = fopen(new_file_name.c_str(),"w")) == NULL) {
                    std::cerr<<"Error: Cannot open file "<<new_file_name<<std::endl;
                    abort();
                }
                PS::S32 n_group_max = 0;
                PS::S32 n_groups = groups.getNumOfGroups();
                for(PS::S32 i=0; i<groups.getPtclN(); i++) {
                    FPSoft &pi = particle_data[groups.getPtclIndex(i)];
                    pi.Ptcl::writeAscii(dout);
                    fprintf(dout, "%26.17e\n", pi.pot_tot);
                    if(i<n_groups && n_group_max<pi.status) n_group_max = pi.status;
                }
                fclose(dout);
                
                PS::ReallocatableArray<FILE*> mout;
                mout.resizeNoInitialize(n_group_max);
                for(PS::S32 i=2; i<=n_group_max; i++) {
                    std::string new_file_name = "snp.M"+std::to_string(i)+"."+std::to_string(file_header.nfile);
                    if ( (mout[i-2] = fopen(new_file_name.c_str(),"w")) == NULL) {
                        std::cerr<<"Error: Cannot open file "<<new_file_name<<std::endl;
                        abort();
                    }
                }
                for(PS::S32 i=0; i<n_groups; i++) {
                    const PS::S32 ni_group = groups.getGroupN(i);
                    const PS::S32 fi = ni_group - 2;

                    FPSoft &pi = particle_data[groups.getPtclIndex(i)];
                    pi.Ptcl::writeAscii(mout[fi]);
                    fprintf(mout[fi], "%26.17e ", pi.pot_tot);

                    const PS::S32* glist=groups.getGroup(i);
                    for(PS::S32 j=0; j<ni_group; j++) {
                        FPSoft &pj = particle_data[glist[j]];
                        pj.Ptcl::writeAscii(mout[fi]);
                        fprintf(mout[fi], "%26.17e ", pj.pot_tot);
                    }
                    fprintf(mout[i],"\n");
                }
                for (PS::S32 i=0;i<n_group_max-1;i++) fclose(mout[i]);
            }
            else {
                FILE *dout;
                std::string new_file_name = "snp."+std::to_string(file_header.nfile);
                if ( (dout = fopen(new_file_name.c_str(),"w")) == NULL) {
                    std::cerr<<"Error: Cannot open file "<<new_file_name<<std::endl;
                    abort();
                }
                for(int i=0; i<groups.getNumOfGroups(); i++) {
                    const PS::S32* glist=groups.getGroup(i);
                    for(int j=0; j<groups.getGroupN(i); j++) {
                        FPSoft &pj = particle_data[glist[j]];
                        pj.Ptcl::writeAscii(dout);
                        fprintf(dout, "%26.17e\n", pj.pot_tot);
                    }
                }
                for(int i=groups.getNumOfGroups(); i<groups.getPtclN(); i++) {
                    FPSoft &pi = particle_data[groups.getPtclIndex(i)];
                    pi.Ptcl::writeAscii(dout);
                    fprintf(dout, "%26.17e\n", pi.pot_tot);
                }

                fclose(dout);
            }
            
        }
        
    }
    fclose(gout);
    
    return 0;
}
