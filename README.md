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
   
2. Change MPI parallelization options
   ```
   ./configure --with-mpi=[auto/yes/no]
   ```
   - auto: automatical detection MPI, if exists, use MPI compiler, otherwise use non-MPI compiler (default)
   - yes: use MPI c++ compiler
   - no: non-MPI c++ compiler
   
3. Disable OpenMP parallelization
    ```
    ./configure --disable-openmp
    ```
    By default OpenMP is used.
    
4. Use X86 with SIMD
    ```
    ./configure --with-simd=[auto/avx/avx2/avx512dq]
    ```
    - auto: automatical detection based on local CPU architecture (default)
    - avx: use AVX for tree force calculation and tree neighbor counting
    - avx2: use AVX2
    - avx512dq: use AVX512F and AVX512DQ
    
5. Use GPU (CUDA)
    ```
    ./configure --enable-cuda
    ```
    By default GPU is not used.
    
6. Debug mode
    ```
    ./configure --with-debug=[assert/g/no]
    ```
    - assert: switch on all assertion check
    - g: switch on compiler option '-g -O0 -fbounds-check' in order to support debugger such as gdb
    - no: no debugging, optimized performance (default)    

Multiple options should be combined together, for example:
```
./configure --prefix=/opt --enable-cuda
```
will install the executable files in /opt and switch on GPU support.

After configure, use 
```
make
make install
```
to compile and install the code.

The excutable file _petar_, _petar.hard.debug_ and _petar.init_ will be installed in [Install path]/bin.
1. _petar_ is the main routine. It is actually a soft link to _petar.\*\*_, where the suffix represents the feature of the code based on the configure.
2. _petar.hard.debug_ is used for debugging if _hard\_dump_ files appears when the code crashes.
3. _petar.init_ is used to generate the initial particle data file for start the simulation, the input data file should have the format: mass, position(x,y,z), velocity(x,y,z) each line for one particle

The data analysis tools are written in _PYTHON_.
They are installed in [Install path]/include/petar
Please add [Install path]/include to the _PYTHON_ include path in order to import the code.

## Use:
After installation, if the [Install path]/bin is in system $PATH envirnoment, the standard way to use the code is
```
[mpiexec -n X] petar [options] [particle data filename]
```
where "[mpiexec -n X]" is used when multiple MPI processors are needed and "X" is the number of processors.

All opitions are listed in the help information, which can be seen by use
```
petar -h
```
Please ignore the error message (a memory allication issue in _FDPS_) after the help information is printed.

The description of the input particle data file is also shown in the help information. 
All snapshots of particle data outputed in the simulation can be used to restart the simulation. 
To restart the simulation with the same configuration of parameters, use
```
[mpiexec -n X] petar -p input.par [snapshot filename]
```
where _input.par_ is automatically generated from the last run (stored in the same diretory of the simulation).

### Useful tools
There are a few useful tools helping users to generate initial input data, find a proper tree time step to start the simulations and data analysis.
Each of the tools are stored in the user defined install_path/bin.
The name of tools have a common style: _petar.[tool name]_.
For each tool, users can get help how to use it by
```
petar.[tool name] -h
```

#### Initial input data file
PeTar has an internal Plummer model generator (equal mass system) by using Henon Unit and half-mass radius being 1.0.
If users want to use their own initial particle data, _petar.init_ can help to transfer their particle data to _petar_ intput data file.
To use _petar.init_:
```
petar.init [options] [particle data filename]
```
The particle data file should contain 7 columns: mass, position[3], velocity[3] and one particle each row.
All Binaries come first and the two components should be next to each other. This is important to obtain a correct initial velocity dispersion. 
Options can be found by
```
petar.init -h
```

#### Find tree time step
The performance of _petar_ is very sensitive to the tree time step.
Use _petar.find.dt_ can help to find a proper time step in order to obtain the best performance.
This tool is used as
```
petar.find.dt [particle data filename]
```
The performance depends on the initial particle data file.
The tools will try to perform a short simulations with different tree time steps and listed the performance one by one.
Users can decide which step is the best one.
Notice if a too large time step is tested, the tool may have no response for long time. This suggests that the testing time step is a very bad choice. Users can kill the test.
A few options support to change the numbers of OpenMP threads and MPI processors, the starting time step and the number of primordial binaries.

#### Parallel data process
The _petar.data.process_ can be used to process snapshot data to detect binaries and calculate Langragian, core radii, averaged mass and velocity dispersion.
The single and binary data are stored for each snapshot with additional suffix ".single" and ".binary".
The data of Lagrangian, core and escapers are generated in separate files.
The multiple CPU cores are used to process the data since the KDTree neighbor search for calculating density and detect binaries is very slow.
The basic way to use is
```
petar.data.process [snapshot path list filename]
```
Users should provide a file that contains a list of pathes for the snapshot data files.
Notice it is better to sort the path in the increasing order of evolution time.
The data order in Lagrangian, escapers and core data file follows the order in the snapshot path list.
The Lagrangian and core data can be read by _LagrangianMultiple_ and _Core_ modules by using _loadtxt_(filename).
The escaper data (single and binary) can be read by _Particle_ and _Binary_.

### Data analysis in _PYTHON3_
The data analysis library provide the tools to identify binaries; calculate Lagrangian radii and core radii; obtain system energy error and check performance of each parts of the code.
To use the tool, first 
```
import petar
```
in PYTHON script, IPYTHON or Jupyter. 
The structure of each module in the data analysis tool is based on a mixture of _collection.OrderedDict_ and _numpy.ndarray_.
Each member of the class that store the data is a _numpy.ndarray_. 
For example, the member _mass_ in _petar.Particle_ is a 1D _numpy.ndarray_, the size of array is the total number of particles.
All other 1D members in the same class have the same size. For two dimensional data, the shape of array is column_number x size.

By use _module-name.keys()_, the name of members are listed.
The two special members are _size_ and _ncols_.
- _size_ is the size of one member of 1D _numpy.ndarray_ type.
- _ncols_ is the total number of columns of all data members, if a member has a 2 dimension array, such as _pos_ of 3 x _size_, it is counted as 3. 

Here listed the supported modules.
- Modules for reading outputs of _petar_:
    - _Particle_: the basic particle data (snapshot files).
    - _Status_: the global parameter of the system such as energy and number of particles (status files).
    - _Profile_: the wall-clock time of different parts of the code (profile files).
- Modules for data analysis:
    - _Binary_: the basic binary data, store two member particles and the Kepler orbital information. Can be generated by using _findPair_ function.
    - _LagrangianMultiple_: for calculating the Lagrangian radii and corresponding properties (averaged mass, number of stars and velocity dispersions)
    - _Core_: for calculating the density, density center and core radius

There are also several functions.
- _join_: join two same type instances of modules. For example, _join_(particle1, particle2) will generate a new _Particle_ instance that contain both two data. Each member is numpy.append(particle1.member, particle2.member).
- _joinLagrangian_: for join two Lagrangian type of data
- _findPair_: detect binaries of one particle list by using _scipy.cKDTree_
- _parallelDataProcessList_: use mutliple CPU cores to process a list of snapshot files and generate single and binary snapshots, Lagrangian data, core data and escaper data. For large _N_, the data process is quite slow, thus using multiple CPU processors can speed up the process. 

More useful tools will be implemented in the future.

### Important tips:
To avoid segmetantional fault in simulations in the case of large number of particles, make sure to set the OMP_STACKSIZE large enough.
For example, use
```
OMP_STACKSIZE=128M [mpiexec -n X] petar [options] [data filename] 
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
