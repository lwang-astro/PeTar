#pragma once

#include <cassert>
#include <string>
#include "hard_ptcl.hpp"
#include "Hermite/hermite_particle.h"
#include "soft_ptcl.hpp"

// Hard debug dump for one cluster
class HardDump{
public:
    typedef H4::ParticleH4<PtclHard> PtclH4;
    PS::F64 time_end;
    PS::S32 n_ptcl;
    PS::S32 n_arti;
    PS::S32 n_group;
    PS::ReallocatableArray<FPSoft> ptcl_arti_bk;
    PS::ReallocatableArray<PtclH4> ptcl_bk;

    //! backup one hard cluster data 
    /*!
       @param[in] _ptcl_local: hard particle backup
       @param[in] _n_ptcl: cluster member number
       @param[in] _ptcl_artifical: artifical particle backup
       @param[in] _n_group: number of groups
       @param[in] _time_end: time ending
       @param[in] _n_split: artifical particle splitting number
     */
    void backup(PtclH4 * _ptcl_local,
                const PS::S32 _n_ptcl,
                FPSoft* _ptcl_artifical,
                const PS::S32 _n_group,
                const PS::F64 _time_end,
                const PS::S32 _n_split) {
        ptcl_bk.resizeNoInitialize(_n_ptcl);
        for (int i=0; i<_n_ptcl; i++) ptcl_bk[i] = _ptcl_local[i];
        time_end = _time_end;
        n_ptcl =_n_ptcl;
        n_arti = _n_group*(2*_n_split+1);
        ptcl_arti_bk.resizeNoInitialize(n_arti);
        for (int i=0; i<n_arti; i++) ptcl_arti_bk[i] = _ptcl_artifical[i];
        n_group = _n_group;
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
        fwrite(&time_end, sizeof(PS::F64),1,fp);
        // hard particles
        fwrite(&n_ptcl, sizeof(PS::S32), 1, fp);
        for(int i=0; i<n_ptcl; i++) ptcl_bk[i].writeBinary(fp);
        // static member
        PS::F64 ptcl_st_dat[3];
        ptcl_st_dat[0] = Ptcl::search_factor;
        ptcl_st_dat[1] = Ptcl::r_search_min;
        ptcl_st_dat[2] = Ptcl::mean_mass_inv;
        fwrite(ptcl_st_dat, sizeof(PS::F64),3, fp);
        // artificial 
        fwrite(&n_arti, sizeof(PS::S32),1,fp);
        fwrite(&n_group, sizeof(PS::S32), 1, fp);
        for (int i=0; i<n_arti; i++) ptcl_arti_bk[i].writeBinary(fp);
        fclose(fp);
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
        size_t rcount = fread(&time_end, sizeof(PS::F64),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
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
        PS::F64 ptcl_st_dat[3];
        rcount = fread(ptcl_st_dat, sizeof(PS::F64),3, fp);
        if (rcount<3) {
            std::cerr<<"Error: Data reading fails! requiring data number is 3, only obtain "<<rcount<<".\n";
            abort();
        }
        Ptcl::search_factor = ptcl_st_dat[0];
        Ptcl::r_search_min  = ptcl_st_dat[1];
        Ptcl::mean_mass_inv = ptcl_st_dat[2];
        // artifical particles
        rcount = fread(&n_arti, sizeof(PS::S32),1,fp);
        // number of groups
        rcount += fread(&n_group, sizeof(PS::S32), 1, fp);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
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
        fclose(fp);
    }

};

// a list of hard dump for multi threads
class HardDumpList{
public:
    int size;
    HardDump* hard_dump;

    void initial(const int _n) {
        size = _n;
        hard_dump = new HardDump[_n];
    }

    void clear() {
        size = 0;
        delete[] hard_dump;
    }

    ~HardDumpList() {
        clear();
    }

    void dump(const char *filename) {
        for (int i=0; i<size; i++) {
            std::string fname = filename + std::string(".") + std::to_string(i);
            hard_dump[i].dumpOneCluster(fname.c_str());
            std::cerr<<"Dump file: "<<fname.c_str()<<std::endl;
        }
    }

    HardDump& operator [] (const int i) {
        assert(i<size);
        return hard_dump[i];
    }
};

static HardDumpList hard_dump;

#ifdef HARD_DUMP
#define DATADUMP(expr) hard_dump.dump("hard_dump");                                        
#else
#define DATADUMP(expr) 
#endif

#ifdef HARD_DEBUG
#define ASSERT(expr)                                                    \
    if(!(expr)) {                                                       \
        std::cerr<<"Assertion! "<<__FILE__<<":"<<__LINE__<<": ("<<#expr<<") fail!"<<std::endl; \
        DATADUMP();                                                     \
        abort();                                                        \
    }
#else
#define ASSERT(expr)
#endif
