#include <iostream>
#include"rand_manager.hpp"

extern "C" {
    uint64_t rand_uint64_();

    double rand_f64_();

    void srand_parallel_(uint64_t* seed);
}

int main (int argc, char** argv) {

    RandomManager rand_manager;

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
    MPI_Init(&argc, &argv);
    int rank, n_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // get number of MPI processors (ranks) in MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
#else
    int rank = 0;
    int n_proc = 1;
#endif    

    uint64_t seed = 1234;

    // test fortran interface
    if (rank==0) std::cout<<"Test fortran random, generate seeds for each MPI and OMP, input seed: "<<seed<<"\n";
    srand_parallel_(&seed);

    if (rank==0) std::cout<<"Generate uint64 and double random samples from each pair of thread and rank, should be different from each other:\n";
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
                uint64_t si64 = rand_uint64_();
                double sf64 = rand_f64_();
#pragma omp critical
                {
                    int i_omp = omp_get_thread_num();
                    std::cout<<"Rank["<<rank<<"] OMP["<<i_omp<<"] seed:"<<RAND_SEED[0]<<" rand_uint64: "<<si64<<" rand_f64: "<<sf64<<std::endl;
                }
            }
        }
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // test cxx interface
    if (rank==0) std::cout<<"Test c++ random, generate seeds for each MPI and OMP, input seed: "<<seed<<"\n";
    srand_parallel(&seed);

    // test read and write

    if (rank==0) std::cout<<"Test write seeds to file 'rand_seeds'\n";

    FILE *fp;
    if( (fp = fopen("rand_seeds","w")) == NULL) {
        fprintf(stderr,"Error: Cannot open file rand_seeds.\n");
        abort();
    }
    rand_manager.writeRandSeeds(fp);
    fclose(fp);

    // test random int
    if (rank==0) std::cout<<"Generate uint64 and f64 samples from each pair of thread and rank, should be the same as those of fortran random:\n";
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
                uint64_t si64 = rand_uint64();
                double sf64 = rand_f64();
#pragma omp critical
                {
                    int i_omp = omp_get_thread_num();
                    std::cout<<"Rank["<<rank<<"] OMP["<<i_omp<<"] seed:"<<RAND_SEED[0]<<" rand_uint64: "<<si64<<" rand_f64: "<<sf64<<std::endl;
                }
            }
        }
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // test read seeds
    if (rank==0) std::cout<<"Test read seeds from file 'rand_seeds'\n";

    if( (fp = fopen("rand_seeds","r")) == NULL) {
        fprintf(stderr,"Error: Cannot open file rand_seeds.\n");
        abort();
    }
    rand_manager.readRandSeeds(fp);
    fclose(fp);

    // test random int after reading previous seeds, show return the same rand seeds and samples
    if (rank==0) std::cout<<"Generate samples from each thread, should return the same seeds and samples as previous:\n";
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
                uint64_t si64 = rand_uint64();
                double sf64 = rand_f64();
#pragma omp critical
                {
                    int i_omp = omp_get_thread_num();
                    std::cout<<"Rank["<<rank<<"] OMP["<<i_omp<<"] seed:"<<RAND_SEED[0]<<" rand_uint64: "<<si64<<" rand_f64: "<<sf64<<std::endl;
                }
            }
        }
    }
    

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
    MPI_Finalize();
#endif

    return 0;
}
