#include "rand.hpp"
#include <cstddef>

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL 
thread_local uint64_t RAND_SEED[2] = {0, 0}; // random seeds for one thread
#else
uint64_t RAND_SEED[2] = {0, 0}; // random seeds
#endif

/*  Written in 2016-2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* This is xoroshiro128+ 1.0, our best and fastest small-state generator
   for floating-point numbers, but its state space is large enough only
   for mild parallelism. We suggest to use its upper bits for
   floating-point generation, as it is slightly faster than
   xoroshiro128++/xoroshiro128**. It passes all tests we are aware of
   except for the four lower bits, which might fail linearity tests (and
   just those), so if low linear complexity is not considered an issue (as
   it is usually the case) it can be used to generate 64-bit outputs, too;
   moreover, this generator has a very mild Hamming-weight dependency
   making our test (http://prng.di.unimi.it/hwd.php) fail after 5 TB of
   output; we believe this slight bias cannot affect any application. If
   you are concerned, use xoroshiro128++, xoroshiro128** or xoshiro256+.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. 

   NOTE: the parameters (a=24, b=16, b=37) of this version give slightly
   better results in our test than the 2016 version (a=55, b=14, c=36).
*/

static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

uint64_t rand_xoroshiro(void) {
    const uint64_t s0 = RAND_SEED[0];
    uint64_t s1 = RAND_SEED[1];
    const uint64_t result = s0 + s1;

    s1 ^= s0;
    RAND_SEED[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
    RAND_SEED[1] = rotl(s1, 37); // c

    return result;
}

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void srand_jump(void) {
    static const uint64_t JUMP[] = { 0xdf900294d8f554a5, 0x170865df4b3201fc };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    for(std::size_t i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
        for(int b = 0; b < 64; b++) {
            if (JUMP[i] & UINT64_C(1) << b) {
                s0 ^= RAND_SEED[0];
                s1 ^= RAND_SEED[1];
            }
            rand_xoroshiro();
        }

    RAND_SEED[0] = s0;
    RAND_SEED[1] = s1;
}

/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^32 starting points,
   from each of which jump() will generate 2^32 non-overlapping
   subsequences for parallel distributed computations. */

void srand_long_jump(void) {
    static const uint64_t LONG_JUMP[] = { 0xd2a98b26625eee7b, 0xdddf9b1090aa7ac1 };

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    for(std::size_t i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
        for(int b = 0; b < 64; b++) {
            if (LONG_JUMP[i] & UINT64_C(1) << b) {
                s0 ^= RAND_SEED[0];
                s1 ^= RAND_SEED[1];
            }
            rand_xoroshiro();
        }

    RAND_SEED[0] = s0;
    RAND_SEED[1] = s1;
}


// C interface
#if defined(__cplusplus)
extern "C" {
#endif

// return one random integer
    uint64_t rand_uint64(void) {
        return rand_xoroshiro();
    }

// return one rand double
    double rand_f64(void) {
        return 0x1.fffffffffffffP-65 * rand_xoroshiro();
    }

// set random seed in parallel
    void srand_parallel(uint64_t *iseed, const int *rank) {

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
        {
            RAND_SEED[0] = *iseed;
            RAND_SEED[1] = 0;

            int i_omp = omp_get_thread_num();
            int n_omp = omp_get_num_threads();

            int n_loop = n_omp * (*rank) + i_omp + 1;
            for (int i=0; i<n_loop; i++) 
                srand_long_jump();
        }
#else
        RAND_SEED[0] = *iseed;
        RAND_SEED[1] = 0;

        int n_loop = *rank + 1;
        for (int i=0; i<n_loop; i++)
            srand_long_jump();
#endif
    }

#if defined(__cplusplus)
}
#endif
