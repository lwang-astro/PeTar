# PeTar
A particle-particle \& Particle-tree (_P<sup>3</sup>T_) with slow-down time-transformed symplectic integrator (slow-down algorithmic regularization; _SDAR_) code for simulating gravitational _N_-body systems including close encounters and few-body systems.

The Doxygen document will be provided in doc directory (not yet done).

The reference paper on ArXiv: https://arxiv.org/abs/2006.16560

## Install:
### Dependence:
FDPS: https://github.com/FDPS/FDPS

SDAR: https://github.com/lwang-astro/SDAR

Please download these two codes and put in the same directory where the _PeTar_ directory exist, in order to successfully compile the code.

### Environment:
To successfully compile the code, the C++ compiler (e.g. gcc, icpc) needs the support of the C++11 standard. To use SSE/BSE package, a Fortran (77) compiler (e.g. gfortran) is needed and should be possile to provide API to c++ code (Currently ifort is not supported yet). The MPI compiler (e.g. mpic++) is required to use MPI. NVIDIA GPU and CUDA compiler is required to use GPU acceleration. The SIMD support is designed for the GNU and Intel compilers. It is not tested for other compilers, thus the GNU or Intel compilers are suggested to use. 

### Make:
Once _FPDS_ and _SDAR_ are available, in the root directoy, use 
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
    
    This option switch on the SIMD support for force calculation, the _auto_ case check whether the compiler (GNU or Intel) support the SIMD instructions and choose the newest one. Notice that the supported options of the compiler and the running CPU are different. Please check your CPU instruction whether the compiled option is supported or not. If the CPU can support more than the compiler, it is suggested to change or update the compiler to get better performance.
    
5. Use GPU (CUDA)
    ```
    ./configure --enable-cuda
    ```
    By default GPU is not used. To switch on it, make sure the NVIDIA CUDA is installed and consistent with the c++ compiler.
    
6. Debug mode
    ```
    ./configure --with-debug=[assert/g/no]
    ```
    - assert: switch on all assertion check
    - g: switch on compiler option '-g -O0 -fbounds-check' in order to support debugger such as gdb
    - no: no debugging, optimized performance (default)   
   
7. Use stellar evolution
    ```
    ./configure --with-interrupt=[bse]
    ```
    Currently only SSE/BSE is the available stellar evolution package (Now still in test phase). Notice SSE/BSE is a combined package, the option argument "sse" not work, only "bse" switches on both.
    When this option is switched on, the standard alone tool _petar.bse_ will also be compiled and installed.
    This is a c++ based tool which uses the API of the SSE/BSE from Fortran77 to c++. It can be used to evolve a group of single and binary stars with OpenMP parallelization.

Multiple options should be combined together, for example:
```
./configure --prefix=/opt/petar --enable-cuda
```
will install the executable files in /opt/petar (this directory requires root permission) and switch on GPU support.

After configure, use 
```
make
make install
```
to compile and install the code.

The excutable file _petar_, _petar.hard.debug_, _petar.init_, _petar.find.dt_, _petar.data.process_, _petar.movie_ and _petar.bse_ (if SSE/BSE is used) will be installed in [Install path]/bin.
1. _petar_ is the main routine. It is actually a soft link to _petar.\*\*_, where the suffix represents the feature of the code based on the configure.
2. _petar.hard.debug_ is used for debugging if _hard\_dump_ files appears when the code crashes.
3. _petar.[tool name]_ are a group of tools for initialization, optimize the performance and data analysis. The details can be checked in the section [Useful tools](#useful-tools).

The data analysis module are written in _PYTHON_.
They are installed in [Install path]/include/petar
Please add [Install path]/include to the _PYTHON_ include path in order to import the code.

## Use:
After installation, if the [Install path]/bin is in system $PATH envirnoment, the standard way to use the code is
```
[mpiexec -n X] petar [options] [particle data filename]
```
where "[mpiexec -n X]" is used when multiple MPI processors are needed and "X" is the number of processors.

All snapshots of particle data outputed in the simulation can be used to restart the simulation. 
To restart the simulation with the same configuration of parameters, use
```
[mpiexec -n X] petar -p input.par [options] [snapshot filename]
```
where _input.par_ is automatically generated from the last run (stored in the same diretory of the simulation).
If the user change the name of this file and restart, the next new parameter file also follows the new name.
Notice if the users want to use new options, these options must be put after "-p input.par", otherwises they are rewrited by values stored in input.par.

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

### Help information

All opitions are listed in the help information, which can be seen by use
```
petar -h
```
The description of the input particle data file is also shown in the help information. 

### Output
#### Printed information
When _petar_ is running, there are a few information printed in the front:
- The _FDPS_ logo and _PETAR_ information are printed, where the copyright, and references for citing are shown.
- The input parameters are listed if they are modified when the corresponding options of _petar_ are used.
- If stellar evolution package (SSE/BSE) is used, the common block and global parameters are printed.
- The name of dumped files for the input parameters are shown, in default, they are input.par, input.par.hard, input.par.bse (if BSE is used)
    - input.par: input parameters of petar, can be used to restart the simulation from a snapshot.
    - input.par.hard: input paremeters of hard (short-range interaction part; Hermite + SDAR), can be used for _petar.hard.debug_ to test the dumped hard cluster.
    - input.par.bse: the BSE parameters, if BSE is used, this is the necessary file to restart the simulation and also for _petar.hard.debug_.

Then after the line "Finish parameter initialization",
The status of the simulation is updated every output time interval (the option -o).
The cotent style is like:
- Time, number of real particles, all particles (including artificial particles), removed particles and escaped particles; in local (first MPI process) and in global (all MPI process)
- Energy check: two rows are printed. First shows the physical energy; Second shows the slow-down energy (see reference printed in _petar_ commander).
    - Error/Total: relative error of current step
    - Error: absolute error of current step
    - Error_cum: cumulative absolute error
    - Total: total energy
    - Kinetic: kinetic energy
    - Potential: potential energy
    - Modify: total modified energy (due to e.g. stellar evolution)
    - Modify_group: modified energy from SDAR groups
    - Modify_single: modified energy from singles
    - Error_hard: energy error in short-range interaction (hard part)
    - Error_hard_cum: cumulative energy error in hard part
- Angular momentum: error at current step, cumulative error, components in x, y, z directions and value
- System total mass, center position and velocity
- Other information if conditions are triggered. 

#### Output files
There are a few output files:
- [data filename prefix].[index]: the snapshot files, the format is the same as the input data file. In the help information of _petar_ commander, users can find the definitions of the columns and header
- [data filename prefix].esc.[MPI rank]: the escaped particle information, first lines show column definition
- [data filename prefix].group.[MPI rank]: the information of new and end of groups (binary, triple ...), since the items in each row can be different because of different numbers of members in groups, there is no universal column headers. It is suggested to use petar.gether first to separate groups with different numbers of members into different files. Then use the python tool _petar.Group_ to read the new files.
- [data filename prefix].[s/b]se.[MPI rank]: if BSE is switched on, the files record the SSE and BSE events, such as type changes, Supernovae, binary evolution phase changes. Each line contain the definition of values, thus can be directly read.

Before access these files, it is suggested to run _petar.gether_ tool first to gether the separated files with different MPI ranks to one file for convenience.
This tool also separate groups with different number of members in xx.group to individual files with suffixe ".n[number of members in groups]".
    
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
The usage:
```
petar.init [options] [particle data filename]
```
The particle data file should contain 7 columns: mass, position[3], velocity[3] and one particle per row.
All Binaries come first and the two components should be next to each other. Besides, when binaries exist, the option '-b [binary number]' should be added to the _petar_ commander. This is important to obtain a correct initial velocity dispersion, tree time step and changeover radii. 

Notice when stellar evolution is switched on, the corresponding options "-s bse" should be used together to generate the correct initial files. In this case, the astronomical unit set (Solar mass [Msun], parsec [PC] and Million year [Myr]) are suggested to use for the initial data. Notice the velocity unit should be PC/Myr. Then the corresponding mass scaling factor is 1.0 between units of PeTar and SSE/BSE. Besides, the option "-u 1" should be added to the _petar_ commander in order to use this astronomical unit set.

#### Find tree time step
The performance of _petar_ is very sensitive to the tree time step.
_petar.find.dt_ can help to find a proper time step in order to obtain the best performance.
The usage:
```
petar.find.dt [options] [petar data filename]
```
The performance of _petar_ depends on the initial particle data file in petar input format.
The tools will try to perform a short simulations with different tree time steps and listed the performance one by one.
Users can decide which step is the best one.
Notice if a too large time step is tested, the test will have no response for long time, which indicates a bad choices of time step.
In such case, this tool terminates the test and return the best result from the previous test.

A few options are available in this tool to change the numbers of OpenMP threads and MPI processors, the minimum tree time step to start the test.
If other options are used in _petar_ commander, e.g. -b [binary number], -u [unit set], -G [gravitational constant], these options should also be added with the option -a and the content enclosed by "":
```
petar.find.dt [options] -a "[petar options]" [petar data filename]. 
```
For example,
```
petar.find.dt -m 2 -o 4 -a "-b 100 -u 1" input
```
use 2 MPI processes, 4 OpenMP threads per MPI, 100 primordial binaries and unit set of 1 [Myr, PC, M*] to select the best dt.

Notice that _petar_ only accepts a tree time step of 0.5^[integer number]. 
Thus in the test, if the user specify the minimum step size by '-s [value]', the step size will be regularized if it does not satisfy this requirement.  

#### Parallel data process
The _petar.data.process_ can be used to process snapshot data to detect binaries and calculate Langragian, core radii, averaged mass and velocity dispersion.
The single and binary data are stored for each snapshot with additional suffix ".single" and ".binary".
The data of Lagrangian, core and escapers are generated in separate files.
The multiple CPU cores are used to process the data since the KDTree neighbor search for calculating density and detect binaries is very slow.
The basic way to use:
```
petar.data.process [options] [snapshot path list filename]
```
Users should provide a file that contains a list of pathes for the snapshot data files.
Notice it is better to sort the path in the increasing order of evolution time.
The _sort_ tool is very convenient to get the time-sorted list.
For example, if the data files have the prefix name of 'data.' (the defaulted case), use
```
ls data.[0-9]* |sort -n -k 1.6 >snap.lst
```
will find all data files in the current directory, sort files by using the suffix (values after 'data.') with an increasing order and save the list to the file of _snap.lst_.
Here '-n' indicates that the values to sort are floating-point number.
'-k' defines the starting position of the number in the filename.
In this example, 'data.' has 5 characters and '1.6' represents the 6th chararcters in the first word (here only one word exists).

The data order in Lagrangian, escapers and core data file follows the order in the snapshot path list.
The Lagrangian and core data can be read by _LagrangianMultiple_ and _Core_ modules by using _loadtxt_(filename).
The escaper data (single and binary) can be read by _Particle_ and _Binary_.
By the way, the _petar_ code also generates escaper data by using the energy and distance criterion (see help of _petar_).

#### Movie generator
The _petar.movie_ is a covenient tool to generate a movie from the snapshot files.
It can generate the movies of the positions (x,y) of stars (x, y of positions), the HR diagram if stellar evolution (SSE/BSE) is switched on, the 2D distribution of semi-major axis and eccentricity of binaries.
The snapshot file list is required to generat the movie.
The basic way to use:
```
petar.movie [options] [snapshot path list filename]
```
If users want to plot information of binaries, it is better to use petar.data.process first to generate detect binaries with multiple CPU cores. Then the next time binary detection is not necessary to run again (use option '--generate-binary 2').

This tool use either _imageio_ or _matplotlib.animation_ to generate movies. It is suggested to install _imageio_ in order to generate movie using mutliple CPU cores. This is much faster than the _matplotlib.animation_ which can only use one CPU core. On the other hand, The ffmpeg is also suggested to install in order to support several commonly used movie formats (e.g. mp4, avi).

#### gether output files from different MPI ranks
The _petar.gether_ is used to gether output files from different MPI ranks to one file (e.g. xx.group.[MPI rank]).
For group files, it also generate individual files with suffix ".n[number of members in groups]'.

#### SSE/BSE steller evolution tool
The _petar.bse_ will be generated when --with-interrupt=bse is used during the configuration.
This is the standard alone SSE/BSE tool to evolve stars and binaries to a given time.
All SSE/BSE global parameters can be set in options.
The basic way to evolve a group of stars:
```
petar.bse [options] mass1, mass2 ...
```
If no mass is provided, the tool can evolve a group of stars with equal logarithmic mass interval.

The basic way to evolve a group of binaries:
```
petar.bse [options] -b [binary data file]
```
where the binary data file contain 4 values (mass1, mass2, period, eccentricity) per line.
Notice that the first line has only one value which is the number of binaries.
The units of mass and period depend on the option '--mscale' and '--tscale'. In defaulted case, the units set are Msun and Myr. For example, when the tscale are not 1.0, the period in Myr is calculated by the input period value x tscale. 

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

More useful tools will be implemented in the future. The tools/analysis/parallel.py is a good example to learn how to use this analysis module.

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
