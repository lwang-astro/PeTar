// random generator for OpenMP and MPI parallelization
#pragma once

#include <stdint.h>
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#include <omp.h>
#endif
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
#include <mpi.h>
#endif

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL 
extern thread_local uint64_t RAND_SEED[2]; // random seeds for one thread
#else
extern uint64_t RAND_SEED[2]; // random seeds
#endif

uint64_t rand_xoroshiro(void);

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void srand_jump(void);

/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^32 starting points,
   from each of which jump() will generate 2^32 non-overlapping
   subsequences for parallel distributed computations. */

void srand_long_jump(void);

// C interface
#if defined(__cplusplus)
extern "C" {
#endif

// return one random integer
    uint64_t rand_uint64(void);

// return one rand double
    double rand_f64(void);

// set random seed in parallel
    void srand_parallel(uint64_t *iseed, const int *rank);

#if defined(__cplusplus)
}
#endif