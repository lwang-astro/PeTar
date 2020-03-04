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
./configure
```
to check the local enviroment and automatically detect the compilers and features.
Options for configure can be found by 
```
./configure -h
```
#### A few useful options:
1. Install path
    ```
    ./configure --prefix=[Install path]
    ```
   If the code is already installed before and the executable file (petar.\*\*) exists in the $PATH enviroment, the configure automatically use the same directory for installing.
2. Use MPI
   ```
   ./configure --with-mpi=[auto/yes/no]
   auto: automatical detection (default)
   yes: use MPI c++ compiler
   no: no MPI
3. Use OpenMP
    ```
    ./configure --with-omp=[yes/no]
    yes: (default)
4. Use X86 with SIMD
    ```
    ./configure --with-simd=[auto/avx/avx2/avx512dq]
    auto: automatical detection based on local CPU architecture (default)
    avx: use AVX for tree force calculation and tree neighbor counting
    avx2: use AVX2
    avx512dq: use AVX512F and AVX512DQ
    ```
5. Use GPU (CUDA)
    ```
    ./configure --with-cuda=[yes/no]
    yes: use GPU to calculate tree force
    no: no GPU support (default)
    ```

Multiple options should be combined together, for example:
```
./configure --prefix=/opt --with-cuda=yes
```
will install the executable files in /opt and switch on GPU support.

After configure, use 
```
make
make install
```
to compile and install the code.

The excutable file _petar.\*\*_, _petar.hard.debug_ and _petar.init_ will be installed in [Install path]/bin.
1. _petar.\*\*_ is the main routine. Here \*\* are suffixes which represent the feature of the code based on the configure.
2. _petar.hard.debug_ is used for debugging if _hard\_dump_ files appears when the code crashes.
3. _petar.init_ is used to generate the initial particle data file for start the simulation, the input data file should have the format: mass, position(x,y,z), velocity(x,y,z) each line for one particle

The data analysis tools are written in _PYTHON_.
They are installed in [Install path]/include/petar
Please add [Install path]/include to the _PYTHON_ include path in order to import the code.

## Use:
After installation, if the [Install path]/bin is in system $PATH envirnoment, the standard way to use the code is
(Assume the executable file name is _petar.mpi.omp.avx2_, in other cases, please replace it to the corresponding name)
```
[mpiexec -n X] petar.mpi.omp.avx2 [options] [particle data filename]
```
where "[mpiexec -n X]" is used when multiple MPI processors are needed and "X" is the number of processors.

All opitions are listed in the help information, which can be seen by use
```
petar.mpi.omp.avx2 -h
```
Please ignore the error message (a memory allication issue in _FDPS_) after the help information is printed.

The description of the input particle data file is also shown in the help information. 
All snapshots of particle data outputed in the simulation can be used to restart the simulation. 
To restart the simulation with the same configuration of parameters, use
```
[mpiexec -n X] petar.mpi.omp.avx2 -p input.par [snapshot filename]
```
where _input.par_ is automatically generated from the last run (stored in the same diretory of the simulation).

### Data analysis in _PYTHON_
The data analysis library provide the tools to identify binaries; calculate Lagrangian radii and core radii; obtain system energy error and check performance of each parts of the code.
To use the tool, first 
```
import petar
```
in PYTHON script, IPYTHON or Jupyter. 
petar.parallel_data_process_list can be used to calculate Lagrangian parameters using multiple CPU processors.
For large _N_, the data process is quite slow, thus using multiple CPU processors can speed up the process. 
More useful tools will be implemented in the future.

### Important tips:
To avoid segmetantional fault in simulations in the case of large number of particles, make sure to set the OMP_STACKSIZE large enough.
For example, use
```
OMP_STACKSIZE=128M [mpiexec -n X] petar.mpi.omp.avx2 [options] [data filename] 
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
