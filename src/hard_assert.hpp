#pragma once

#include <cassert>
#include <string>
#include "hard_ptcl.hpp"
#include "Hermite/hermite_particle.h"
#include "soft_ptcl.hpp"
#ifdef BSE_BASE
#include "../parallel-random/rand_interface.hpp"
#endif

// Hard debug dump for one cluster
class HardDump{
public:
    typedef H4::ParticleH4<PtclHard> PtclH4;
    PS::F64 time_offset;
    PS::F64 time_end;
    PS::S32 n_ptcl;
    PS::S32 n_arti;
    PS::S32 n_group;
    PS::ReallocatableArray<PS::S32> n_member_in_group;
    PS::ReallocatableArray<FPSoft> ptcl_arti_bk;
    PS::ReallocatableArray<PtclH4> ptcl_bk;
    RandomManager rand_manager;
    bool backup_flag;

    HardDump(): time_offset(0), time_end(0), n_ptcl(0), n_arti(0), n_group(0), n_member_in_group(), ptcl_arti_bk(), ptcl_bk(), rand_manager(), backup_flag(false) {}

    //! backup one hard cluster data 
    /*!
       @param[in] _ptcl_local: hard particle backup
       @param[in] _n_ptcl: cluster member number
       @param[in] _ptcl_artifical: artifical particle backup
       @param[in] _n_group: number of groups
       @param[in] _time_offset: offset of time
       @param[in] _time_end: time ending without offset
       @param[in] _n_artificial: artifical particle number
     */
    void backup(PtclH4 * _ptcl_local,
                const PS::S32 _n_ptcl,
                FPSoft* _ptcl_artifical,
                const PS::S32 _n_group,
                const PS::S32* _n_member_in_group,
                const PS::F64 _time_offset,
                const PS::F64 _time_end,
                const PS::S32 _n_artificial) {
        ptcl_bk.resizeNoInitialize(_n_ptcl);
        for (int i=0; i<_n_ptcl; i++) ptcl_bk[i] = _ptcl_local[i];
        time_offset = _time_offset;
        time_end = _time_end;
        n_ptcl =_n_ptcl;
        n_member_in_group.resizeNoInitialize(_n_group);
        for (int i=0; i<_n_group; i++) n_member_in_group[i] = _n_member_in_group[i];
        if (_ptcl_artifical!=NULL) {
            n_arti = _n_group*_n_artificial;
            ptcl_arti_bk.resizeNoInitialize(n_arti);
            for (int i=0; i<n_arti; i++) ptcl_arti_bk[i] = _ptcl_artifical[i];
        }
        else n_arti = 0;
        n_group = _n_group;
        backup_flag = true;
    }                

    //! Dumping one cluster data for debuging
    /* 
       @param[in] _fname: file name to write
     */
    void dumpOneCluster(const char* _fname) {
        std::FILE* fp = std::fopen(_fname,"w");
        if (fp==NULL) {
            std::cerr<<"Error: filename "<<_fname<<" cannot be open!\n";
            abort();
        }
        fwrite(&time_offset, sizeof(PS::F64),1,fp);
        fwrite(&time_end, sizeof(PS::F64),1,fp);
        // hard particles
        fwrite(&n_ptcl, sizeof(PS::S32), 1, fp);
        for(int i=0; i<n_ptcl; i++) ptcl_bk[i].writeBinary(fp);
        // static member
        PS::F64 ptcl_st_dat[4];
        ptcl_st_dat[0] = Ptcl::search_factor;
        ptcl_st_dat[1] = Ptcl::r_search_min;
        ptcl_st_dat[2] = Ptcl::mean_mass_inv;
        ptcl_st_dat[3] = Ptcl::r_group_crit_ratio;
        fwrite(ptcl_st_dat, sizeof(PS::F64),4, fp);
        // artificial 
        fwrite(&n_arti, sizeof(PS::S32),1,fp);
        fwrite(&n_group, sizeof(PS::S32), 1, fp);
        fwrite(n_member_in_group.getPointer(), sizeof(PS::S32), n_group, fp);
        for (int i=0; i<n_arti; i++) ptcl_arti_bk[i].writeBinary(fp);
        fclose(fp);
#ifdef BSE_BASE
        // rand seed
        rand_manager.writeRandSeedLocalBinary(fp);
#endif
        backup_flag = false;
    }

    //! reading one cluster data for debuging
    /* 
       @param[in]  _fname: file name to read
     */
    void readOneCluster(const char* _fname) {
        std::FILE* fp = std::fopen(_fname,"r");
        if (fp==NULL) {
            std::cerr<<"Error: filename "<<_fname<<" cannot be open!\n";
            abort();
        }
        // read time
        size_t rcount = fread(&time_offset, sizeof(PS::F64),1,fp);
        rcount += fread(&time_end, sizeof(PS::F64),1,fp);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
        // read hard particles
        rcount = fread(&n_ptcl, sizeof(PS::S32),1, fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
        if (n_ptcl<=0) {
            std::cerr<<"Error: particle number "<<n_ptcl<<" <=0 !\n";
            abort();
        }
        ptcl_bk.resizeNoInitialize(n_ptcl);
        for(int i=0; i<n_ptcl; i++) ptcl_bk[i].readBinary(fp);
        // static members
        PS::F64 ptcl_st_dat[4];
        rcount = fread(ptcl_st_dat, sizeof(PS::F64),4, fp);
        if (rcount<4) {
            std::cerr<<"Error: Data reading fails! requiring data number is 3, only obtain "<<rcount<<".\n";
            abort();
        }
        Ptcl::search_factor = ptcl_st_dat[0];
        Ptcl::r_search_min  = ptcl_st_dat[1];
        Ptcl::mean_mass_inv = ptcl_st_dat[2];
        Ptcl::r_group_crit_ratio = ptcl_st_dat[3];
        // artifical particles
        rcount = fread(&n_arti, sizeof(PS::S32),1,fp);
        // number of groups
        rcount += fread(&n_group, sizeof(PS::S32), 1, fp);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
        n_member_in_group.resizeNoInitialize(n_group);
        rcount = fread(n_member_in_group.getPointer(), sizeof(PS::S32), n_group, fp);
        if (rcount<(size_t)n_group) {
            std::cerr<<"Error: Data reading fails! requiring data number is "<<n_group<<", only obtain "<<rcount<<".\n";
            abort();
        }
        if (n_arti<0) {
            std::cerr<<"Error: artificial particle number "<<n_arti<<" <0 !\n";
            abort();
        }
        if (n_group<0) {
            std::cerr<<"Error: group number "<<n_group<<" <0 !\n";
            abort();
        }
        // read artifical particles
        if (n_arti>0) {
            ptcl_arti_bk.resizeNoInitialize(n_arti);
            for (int i=0; i<n_arti; i++) ptcl_arti_bk[i].readBinary(fp);
        }
#ifdef BSE_BASE
        // read rand seed
        rand_manager.readRandSeedLocalBinary(fp);
#endif
        fclose(fp);
    }

};

// a list of hard dump for multi threads
class HardDumpList{
public:
    int size;
    int mpi_rank;
    int dump_number;
    HardDump* hard_dump;

    HardDumpList(): size(0), mpi_rank(0), dump_number(0), hard_dump(NULL) {}

    void initial(const int _nthread, const int _rank=0) {
        size = _nthread;
        mpi_rank = _rank;
        hard_dump = new HardDump[_nthread];
    }

    void clear() {
        size = 0;
        if (hard_dump!=NULL) {
            delete[] hard_dump;
            hard_dump=NULL;
        }
    }

    ~HardDumpList() {
        clear();
    }

    void dumpAll(const char *filename) {
        std::string point(".");
        std::string fname_prefix = filename + point + std::to_string(mpi_rank) + point;
        for (int i=0; i<size; i++) {
            std::time_t tnow = std::time(nullptr);
            std::string fname = fname_prefix + std::to_string(i) + point + std::to_string(dump_number++) + point + std::to_string(tnow);
            if (hard_dump[i].backup_flag) {
                hard_dump[i].dumpOneCluster(fname.c_str());
                std::cerr<<"Dump file: "<<fname.c_str()<<std::endl;
            }
        }
    }

    void dumpThread(const char *filename){
        const PS::S32 ith = PS::Comm::getThreadNum();
        std::string point(".");
        std::time_t tnow = std::time(nullptr);
        //std::tm *local_time = localtime(&tnow);
        std::string fname = filename + point + std::to_string(mpi_rank) + point + std::to_string(ith) + point + std::to_string(dump_number++) + point + std::to_string(tnow);
        if (hard_dump[ith].backup_flag) {
            hard_dump[ith].dumpOneCluster(fname.c_str());
            std::cerr<<"Thread: "<<ith<<" Dump file: "<<fname.c_str()<<std::endl;
        }
    }

    HardDump& operator [] (const int i) {
        assert(i<size);
        return hard_dump[i];
    }
};

static HardDumpList hard_dump;

#ifdef HARD_DUMP
#define DATADUMP(expr) hard_dump.dumpThread(expr)
#else
#define DATADUMP(expr) 
#endif

#ifdef HARD_DEBUG
#define ASSERT(expr)                                                    \
    if(!(expr)) {                                                       \
        std::cerr<<"Assertion! "<<__FILE__<<":"<<__LINE__<<": ("<<#expr<<") fail!"<<std::endl; \
        DATADUMP("hard_dump");                                          \
        abort();                                                        \
    }
#else
#define ASSERT(expr)
#endif
