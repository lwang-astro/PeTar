#pragma once
#include <iostream>
#include <inttypes.h>
#include "rand.hpp"
#include <cstdio>
#include <string>

//! random generator manager              
class RandomManager{
public:

    //! initial random seeds for all threads and processors
    template <class Tio>
    void initialAll(const Tio& _input) {
        if (_input.seedfile.value!="__NONE__")
            readRandSeeds(_input.seedfile.value.c_str());
        else {
            uint64_t seed_i64 = _input.seed.value;
            srand_parallel(&seed_i64);
        }
    }

    //! initial random seeds for local thread 
    template <class Tio>
    void initialLocal(const Tio& _input) {
        if (_input.seedfile.value!="__NONE__")
            readRandSeedLocal(_input.seedfile.value.c_str());
        else {
            uint64_t seed_i64 = _input.seed.value;
            // use srand_parallel in case with local intialization
            srand_parallel(&seed_i64);
        }
    }

    // initial all seeds directly
    void initialAll(uint64_t seed) {
        srand_parallel(&seed);
    }

    //! print all seeds 
    void printSeeds(std::ostream & fout) const{

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        int rank, n_proc;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        // get number of MPI processors (ranks) in MPI_COMM_WORLD
        MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
#else
        int rank = 0;
        int n_proc = 1;
#endif    

        if (rank==0) fout<<"Random seeds of all threads and processors:\n";
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        for (int i=0; i<n_proc; i++) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            if (rank==i) {
#pragma omp parallel 
                {
#pragma omp critical
                    {
                        int i_omp = omp_get_thread_num();
                        fout<<"Rank["<<rank<<"]-Thread["<<i_omp<<"]: "<<RAND_SEED[0]<<" "<<RAND_SEED[1]<<std::endl;
                    }
                }
            }
        }
    }

    //! write seed of local thread
    void writeRandSeedLocal(FILE* fp) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
        {
            fprintf(fp, "%" PRIu64 " ", RAND_SEED[0]);
            fprintf(fp, "%" PRIu64 " ", RAND_SEED[1]);
        }
#else
        fprintf(fp, "%" PRIu64 " ", RAND_SEED[0]);
        fprintf(fp, "%" PRIu64 " ", RAND_SEED[1]);
#endif
    }

    //! write seed of local thread in binary format
    void writeRandSeedLocalBinary(FILE* fp) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
        {
            fwrite(RAND_SEED, sizeof(uint64_t), 2, fp);
        }
#else
        fwrite(RAND_SEED, sizeof(uint64_t), 2, fp);
#endif
    }

    //! read seed of local thread
    void readRandSeedLocal(FILE* fp) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
        {
            int rcount = 0;
            rcount += fscanf(fp, "%" PRIu64 " ", &RAND_SEED[0]);
            rcount += fscanf(fp, "%" PRIu64 " ", &RAND_SEED[1]);
            if(rcount<2) {
                std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                abort();
            }
        }
#else
        int rcount = 0;
        rcount += fscanf(fp, "%" PRIu64 " ", &RAND_SEED[0]);
        rcount += fscanf(fp, "%" PRIu64 " ", &RAND_SEED[1]);
        if(rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
#endif
    }

    //! read seed of local thread
    void readRandSeedLocalBinary(FILE* fp) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
        {
            int rcount = fread(RAND_SEED, sizeof(uint64_t), 2, fp);
            if(rcount<2) {
                std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                abort();
            }
        }
#else
        int rcount = fread(RAND_SEED, sizeof(uint64_t), 2, fp);
        if(rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
#endif
    }

    //! write random seeds from all threads and MPI processors
    void writeRandSeeds(FILE* fp) {
        
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        int n_proc;
        // get number of MPI processors (ranks) in MPI_COMM_WORLD
        MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
        int rank;
        // get current MPI processor id (rank) in MPI_COMM_WORLD
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
        int rank=0;
        int n_proc = 1;
#endif

        int n_omp = 1;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
        n_omp = omp_get_num_threads();
#endif
        int nseeds = 2*n_proc*n_omp;
        uint64_t rand_seed_all[nseeds];
        uint64_t rand_seed_local[2*n_omp];

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
        {
            int i_omp = omp_get_thread_num();
            rand_seed_local[2*i_omp] = RAND_SEED[0];
            rand_seed_local[2*i_omp+1] = RAND_SEED[1];
        }
#else
        rand_seed_local[0] = RAND_SEED[0];
        rand_seed_local[1] = RAND_SEED[1];
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        MPI_Gather(rand_seed_local, 2*n_omp, MPI_UINT64_T,  rand_seed_all, 2*n_omp, MPI_UINT64_T, 0, MPI_COMM_WORLD);
#else
        for (int i=0; i<2*n_omp; i++) 
            rand_seed_all[i] = rand_seed_local[i];
#endif    
    
        // save all seeds
        if (rank==0) {
            for (int i=0; i<nseeds; i++) {
                fprintf(fp, "%" PRIu64 " ", rand_seed_all[i]);
            }
        }
    }

    //! read random seeds for all threads and MPI processors
    void readRandSeeds(FILE* fp) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        int n_proc;
        // get number of MPI processors (ranks) in MPI_COMM_WORLD
        MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
        int rank;
        // get current MPI processor id (rank) in MPI_COMM_WORLD
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
        int rank=0;
        int n_proc = 1;
#endif

        int n_omp = 1;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
        n_omp = omp_get_num_threads();
#endif

        int nseeds = 2*n_proc*n_omp;
        uint64_t rand_seed_all[nseeds];
        uint64_t rand_seed_local[2*n_omp];

        if (rank==0) {
            int rcount = 0;
            for (int i=0; i<nseeds; i++) {
                rcount += fscanf(fp, "%" PRIu64 " ", &rand_seed_all[i]);
            }

            if(rcount<nseeds) {
                std::cerr<<"Error: Data reading fails! requiring data number is "<<nseeds<<", only obtain "<<rcount<<".\n";
                abort();
            }
        }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
        MPI_Scatter(rand_seed_all, 2*n_omp, MPI_UINT64_T, rand_seed_local, 2*n_omp, MPI_UINT64_T,  0, MPI_COMM_WORLD);
#else
        for (int i=0; i<nseeds; i++)
            rand_seed_local[i] = rand_seed_all[i];
#endif    

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
        {
            int i_omp = omp_get_thread_num();
            RAND_SEED[0] = rand_seed_local[2*i_omp];
            RAND_SEED[1] = rand_seed_local[2*i_omp+1];
        }
#else
        RAND_SEED[0] = rand_seed_local[0];
        RAND_SEED[1] = rand_seed_local[1];
#endif
    }

    //! write rand seeds to file
    void writeRandSeeds(const char* _fname) {
        FILE* fin;
        if( (fin = fopen(_fname,"w")) == NULL) {
            fprintf(stderr,"Error: cannot open file %s.\n", _fname);
            abort();
        }
        writeRandSeeds(fin);
        fclose(fin);
    }

    //! write rand seeds to file for local thread
    void writeRandSeedLocal(const char* _fname) {
        FILE* fin;
        if( (fin = fopen(_fname,"w")) == NULL) {
            fprintf(stderr,"Error: cannot open file %s.\n", _fname);
            abort();
        }
        writeRandSeedLocal(fin);
        fclose(fin);
    }

    //! read rand seeds from file for local thread
    void readRandSeedLocal(const char* _fname) {
        FILE* fin;
        if( (fin = fopen(_fname,"r")) == NULL) {
            fprintf(stderr, "Error: cannot open file %s.\n", _fname);
            abort();
        }
        readRandSeedLocal(fin);
        fclose(fin);
    }

    //! read rand seeds from file
    void readRandSeeds(const char* _fname) {
        FILE* fin;
        if( (fin = fopen(_fname,"r")) == NULL) {
            fprintf(stderr, "Error: cannot open file %s.\n", _fname);
            abort();
        }
        readRandSeeds(fin);
        fclose(fin);
    }
};
