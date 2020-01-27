# PeTar
Particle-particle \& Particle-tree (_P<sup>3</sup>T_) with slow-down time-transformed symplectic integrator (slow-down algorithmic regularization; _SDAR_) code for simulating gravitational _N_-body systems including close encounters and few-body systems.

## Dependence: 
FDPS: https://github.com/FDPS/FDPS

SDAR: https://github.com/lwang-astro/SDAR

Please download these two codes and put in the same directory where the _PeTar_ directory exist, in order to successfully compile the code.

## Method:
### Algorithm of integration: 
The calculation using the kick-drift-kick mode.
For each cycle, three major steps are done:
1. Search clusters and few-body systems
    1. Construct particle-tree including all real particles for neighbor searching.
    2. Search clusters where each real particle stay in the same cluster as its neighbors.
    3. For each cluster, construct binary tree of all real members and detect few-body systems; for stable few-body systems, construct artificial particles (tidal tensor measure points and orbital-averaged-force sample particles). 
2. Kick (Soft) step
    1. Construct particle-tree for all real and artificial particles and calculate long-distance (soft) forces with linear cutoff (an universal cutoff radius).
    2. Correct changeover function using neighbor list of each particle.
3. Drift (Hard) step
    1. For clusters with more than one memeber, Use _SDAR_ _Hermite_ module (4th order Hermite with slow-down AR method) to integrate the orbits.
    2. For one cluster, just do drift.  

### Parallelization:
1. Tree construction is done by _FPDS_, using _MPI_ and _OpenMP_ .
2. Soft force calculation kernel use _SIMD_ (_AVX_, _AVX2_, _AVX512_).
3. Hard calculation use _OpenMP_ for the loop of clusters

### AMUSE API:
Current support API: gravitational dynamics, gravity field
