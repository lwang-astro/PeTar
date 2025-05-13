#pragma once

#include <cassert>
#include <string>
#include "hard_ptcl.hpp"
#include "Hermite/hermite_particle.h"
#include "soft_ptcl.hpp"
#ifdef BSE_BASE
#include "../parallel-random/rand_interface.hpp"
#endif
#ifdef EXTERNAL_HARD
#ifdef GALPY
#include "../galpy-interface/galpy_interface.h"
#endif
#endif

// Hard debug dump for one cluster
class HardDump{
public:
    typedef H4::ParticleH4<PtclHard> PtclH4;
    PS::F64 time_offset;
    PS::F64 time_end;
    PS::F64 gcm_mass;
    PS::F64vec gcm_pos;
    PS::F64vec gcm_vel;
    PS::S32 n_ptcl;
    PS::S32 n_arti;
    PS::S32 n_group;
    PS::ReallocatableArray<PS::S32> n_member_in_group;
    PS::ReallocatableArray<FPSoft> ptcl_arti_bk;
    PS::ReallocatableArray<PtclH4> ptcl_bk;
#ifdef BSE_BASE
    uint64_t rand_seed[2];
    RandomManager rand_manager;
#endif
    bool backup_flag;

    HardDump(): time_offset(0), time_end(0), gcm_mass(0), gcm_pos(), gcm_vel(), n_ptcl(0), n_arti(0), n_group(0), n_member_in_group(), ptcl_arti_bk(), ptcl_bk(), 
#ifdef BSE_BASE
                rand_seed{0,0}, rand_manager(), 
#endif
                backup_flag(false) {}

    //! backup one hard cluster data 
    /*!
       @param[in] _ptcl_local: hard particle backup
       @param[in] _n_ptcl: cluster member number
       @param[in] _ptcl_artifical: artifical particle backup
       @param[in] _n_group: number of groups
       @param[in] _time_offset: offset of time
       @param[in] _time_end: time ending without offset
       @param[in] _gcm_mass: mass of the system
       @param[in] _gcm_pos: position of the system
       @param[in] _gcm_vel: velocity of the system
       @param[in] _n_artificial: artifical particle number
     */
    void backup(PtclH4 * _ptcl_local,
                const PS::S32 _n_ptcl,
                FPSoft* _ptcl_artifical,
                const PS::S32 _n_group,
                const PS::S32* _n_member_in_group,
                const PS::F64 _time_offset,
                const PS::F64 _time_end,
                const PS::F64 _gcm_mass,
                const PS::F64vec& _gcm_pos,
                const PS::F64vec& _gcm_vel,
                const PS::S32 _n_artificial) {
        ptcl_bk.resizeNoInitialize(_n_ptcl);
        for (int i=0; i<_n_ptcl; i++) ptcl_bk[i] = _ptcl_local[i];
        time_offset = _time_offset;
        time_end = _time_end;
        gcm_mass = _gcm_mass;
        gcm_pos = _gcm_pos;
        gcm_vel = _gcm_vel;
        n_ptcl =_n_ptcl;
        n_member_in_group.resizeNoInitialize(_n_group);
        for (int i=0; i<_n_group; i++) n_member_in_group[i] = _n_member_in_group[i];
        if (_ptcl_artifical!=NULL) {
            n_arti = _n_group*_n_artificial;
            ptcl_arti_bk.resizeNoInitialize(n_arti);
            for (int i=0; i<n_arti; i++) ptcl_arti_bk[i] = _ptcl_artifical[i];
        }
        else {
            n_arti = 0;
            ptcl_arti_bk.resizeNoInitialize(n_arti);
        }
        n_group = _n_group;
#ifdef BSE_BASE
        rand_manager.getRandSeedLocal(rand_seed);
#endif        
        backup_flag = true;
    }                

    //! write one cluster data with BINARY format
    /*! @param[in] _fout: file IO for write
     */
    void writeOneClusterBinary(FILE* fp) const {
        fwrite(&time_offset, sizeof(PS::F64),1,fp);
        fwrite(&time_end, sizeof(PS::F64),1,fp);
        fwrite(&gcm_mass, sizeof(PS::F64),1,fp);
        fwrite(&gcm_pos, sizeof(PS::F64vec),1,fp);
        fwrite(&gcm_vel, sizeof(PS::F64vec),1,fp);
        // hard particles
        fwrite(&n_ptcl, sizeof(PS::S32), 1, fp);
        for(int i=0; i<n_ptcl; i++) ptcl_bk[i].writeBinary(fp);
        // static member
        PS::F64 ptcl_st_dat[5];
        ptcl_st_dat[0] = Ptcl::search_factor;
        ptcl_st_dat[1] = Ptcl::r_search_min;
        ptcl_st_dat[2] = Ptcl::mean_mass_inv;
        ptcl_st_dat[3] = PtclHard::r_group_over_in;
        ptcl_st_dat[4] = PtclHard::r_search_group_over_in;
        fwrite(ptcl_st_dat, sizeof(PS::F64),5, fp);
        // artificial 
        fwrite(&n_arti, sizeof(PS::S32),1,fp);
        fwrite(&n_group, sizeof(PS::S32), 1, fp);
        fwrite(n_member_in_group.getPointer(), sizeof(PS::S32), n_group, fp);
        for (int i=0; i<n_arti; i++) ptcl_arti_bk[i].writeBinary(fp);
#ifdef BSE_BASE
        fwrite(rand_seed, sizeof(uint64_t), 2, fp);
#endif
    }

    //! Dumping one cluster data and random seed for a given file name
    /* 
       @param[in] _fname: file name to write
       @param[in] _append_flag: if true, append data into existing file and do not seperate random seed to different file
       @param[in] _dump_once_flag: if true, only dump once for the same data
     */
    void dumpOneClusterAndSeed(const char* _fname, const bool _append_flag, const bool _dump_once_flag=true) {
        std::FILE* fp;
        if (_append_flag) 
            fp = std::fopen(_fname,"a");
        else
            fp = std::fopen(_fname,"w");
        if (fp==NULL) {
            std::cerr<<"Error: filename "<<_fname<<" cannot be open!\n";
            abort();
        }
        writeOneClusterBinary(fp);
        fclose(fp);
        if (_dump_once_flag) backup_flag = false;
    }

    //! read one cluster data with BINARY format
    /*! @param[in] _fin: file IO for read
     */
    void readOneClusterBinary(FILE* fp) {
        // read time
        size_t rcount = fread(&time_offset, sizeof(PS::F64),1,fp);
        rcount += fread(&time_end, sizeof(PS::F64),1,fp);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
        // read gcm
        rcount = fread(&gcm_mass, sizeof(PS::F64),1,fp);
        rcount += fread(&gcm_pos, sizeof(PS::F64vec),1,fp);
        rcount += fread(&gcm_vel, sizeof(PS::F64vec),1,fp);
        if (rcount<3) {
            std::cerr<<"Error: Data reading fails! requiring data number is 4, only obtain "<<rcount<<".\n";
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
        PS::F64 ptcl_st_dat[5];
        rcount = fread(ptcl_st_dat, sizeof(PS::F64),5, fp);
        if (rcount<4) {
            std::cerr<<"Error: Data reading fails! requiring data number is 3, only obtain "<<rcount<<".\n";
            abort();
        }
        Ptcl::search_factor = ptcl_st_dat[0];
        Ptcl::r_search_min  = ptcl_st_dat[1];
        Ptcl::mean_mass_inv = ptcl_st_dat[2];
        PtclHard::r_group_over_in = ptcl_st_dat[3];
        PtclHard::r_search_group_over_in = ptcl_st_dat[4];
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
        rand_manager.readRandSeedLocalBinary(fp);
#endif
    }

    //! reading one cluster data and random seed for a given filename
    /* 
       @param[in]  _fname: file name to read
     */
    void loadOneClusterAndSeed(const char* _fname) {
        std::FILE* fp = std::fopen(_fname,"r");
        if (fp==NULL) {
            std::cerr<<"Error: filename "<<_fname<<" cannot be open!\n";
            abort();
        }
        readOneClusterBinary(fp);
        fclose(fp);
    }

};

// a list of hard dump for multi threads
class HardDumpList{
public:
    int size;
    int mpi_rank;
    int omp_level;
    int dump_number;
#ifdef EXTERNAL_HARD
#ifdef GALPY
    GalpyManager* galpy_manager;
#endif
#endif
    HardDump* hard_dump;

    HardDumpList(): size(0), mpi_rank(0), omp_level(0), dump_number(0), 
#ifdef EXTERNAL_HARD
#ifdef GALPY
                    galpy_manager(NULL),
#endif
#endif
                    hard_dump(NULL) {}

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

    void dumpAll(const char *filename, 
                 const bool long_suffix_flag=true, 
                 const bool append_flag=false, 
                 const bool dump_once_flag = true, 
                 const bool print_flag=true) {
        std::string point("_");
        for (int i=0; i<size; i++) {
            if (hard_dump[i].backup_flag) {
                std::time_t tnow = std::time(nullptr);
                std::string fname = filename;
                if (long_suffix_flag) fname = filename + point + std::to_string(hard_dump[i].time_offset) + point + std::to_string(mpi_rank) + point + std::to_string(i) + point + std::to_string(dump_number++) + point + std::to_string(tnow);
                hard_dump[i].dumpOneClusterAndSeed(fname.c_str(), append_flag, dump_once_flag);
#ifdef EXTERNAL_HARD
#ifdef GALPY
                galpy_manager->writePotentialPars((fname+".galpy").c_str(), hard_dump[i].time_offset, false);
#endif
#endif            
                if (print_flag) std::cerr<<"Dump file: "<<fname.c_str()<<std::endl;
            }
        }
    }

    void dumpThread(const char *filename, 
                    const bool long_suffix_flag=true, 
                    const bool append_flag=false, 
                    const bool dump_once_flag = true, 
                    const bool print_flag=true){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        const PS::S32 ith = omp_get_ancestor_thread_num(omp_level);
#else
        const PS::S32 ith = 0;
#endif
        //const PS::S32 ith = PS::Comm::getThreadNum();
        if (hard_dump[ith].backup_flag) {
            std::string point("_");
            std::time_t tnow = std::time(nullptr);
            //std::tm *local_time = localtime(&tnow);
            std::string fname = filename;
            if (long_suffix_flag) fname = filename + point + std::to_string(hard_dump[ith].time_offset) + point + std::to_string(mpi_rank) + point + std::to_string(ith) + point + std::to_string(dump_number++) + point + std::to_string(tnow);
            hard_dump[ith].dumpOneClusterAndSeed(fname.c_str(), append_flag, dump_once_flag);
#ifdef EXTERNAL_HARD
#ifdef GALPY
            galpy_manager->writePotentialPars((fname+".galpy").c_str(), hard_dump[ith].time_offset, false);
#endif
#endif            
            if (print_flag) std::cerr<<"Thread: "<<ith<<" Dump file: "<<fname.c_str()<<std::endl;
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
#define DATADUMPAPP(expr) hard_dump.dumpThread(expr, false, true, false, false)
#else
#define DATADUMP(expr) 
#define DATADUMPAPP(expr) 
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
