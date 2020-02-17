# PeTar
Particle-particle \& Particle-tree (_P<sup>3</sup>T_) with slow-down time-transformed symplectic integrator (slow-down algorithmic regularization; _SDAR_) code for simulating gravitational _N_-body systems including close encounters and few-body systems.

The Doxygen document will be provided in doc directory (not yet done)

## Install:
### Dependence:
FDPS: https://github.com/FDPS/FDPS

SDAR: https://github.com/lwang-astro/SDAR

Please download these two codes and put in the same directory where the _PeTar_ directory exist, in order to successfully compile the code.

### Make:
Onces _FPDS_ and _SDAR_ are available, in the root directoy, use 
```
make
```
to create the excutable file _petar_ and _petar.hard.debug_.
_petar_ is the main routine.
_petar.hard.debug_ is used for debugging if _hard\_dump_ files appears when the code crashes.

### Parallelization options
In defaulted, _MPI_, _OpenMP_ and _AVX2_ parallelization methods are used to compile the code.
If users want to change the options, the Makefile should be modified.
The _PATH_ to _FPDS_, _SDAR_, and _CUDA_ libraries and a few options of parallelization features (_use\_xx_) are shown in the head of the Makefile.
Users can modify these parts.

For example, when the comment symbol '#' is removed for the line
```
use_gpu_cuda=yes
```
the GPU feature is on. 
Users should make sure the _CUDA_ compiler _nvcc_ is installed.

In the future, _GNU AutoConf_ will be implemented.

## Use:
The standard way to use the code is
```
[mpiexec -n X] ./petar [options] [particle data filename]
```
where "[mpiexec -n X]" is used when multiple MPI processors are needed and "X" is the number of processors.

All opitions are listed in the help information, which can be seen by use
```
./petar -h
```
Please ignore the error message (a memory allication issue in _FDPS_) after the help information is printed.

The description of the input particle data file is also shown in the help information. 
All snapshots of particle data outputed in the simulation can be used to restart the simulation. 
To restart the simulation with the same configuration of parameters, use
```
./petar -p input.par [snapshot filename]
```
where _input.par_ is automatically generated from the last run (stored in the same diretory of the simulation).

### Important tips:
To avoid segmetantional fault in simulations in the case of large number of particles, make sure to set the OMP_STACKSIZE large enough.
For example, use
```
OMP_STACKSIZE=128M [mpiexec -n X] ./petar [options] [data filename] 
```

A convenient way is to add
```
export OMP_STACKSIZE=128M
```
in the shell configure/initial file (e.g. .bashrc for _bash_) to avoid type "OMP_STACKSIZE=128M" every time.

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
2. Soft force calculation kernel use _SIMD_ (_AVX_, _AVX2_, _AVX512_) or _GPU_ (_CUDA_).
3. Hard calculation use _OpenMP_ for the loop of clusters

### AMUSE API:
Current support API: gravitational dynamics, gravity field, stopping conditions
