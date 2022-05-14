```
    ██████╗ ███████╗████████╗ █████╗ ██████╗ 
    ██╔══██╗██╔════╝╚══██╔══╝██╔══██╗██╔══██╗
    ██████╔╝█████╗     ██║   ███████║██████╔╝
    ██╔═══╝ ██╔══╝     ██║   ██╔══██║██╔══██╗
    ██║     ███████╗   ██║   ██║  ██║██║  ██║
    ╚═╝     ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝
```

PeTar is a N-body code designed to model collisional stellar systems, where multiplicity (binaries, triples ...) and close encounters are important for dynamical evolution. 
It combines three integration methods:
- The Barnes-Hut tree (Barnes & Hut 1986) is used to calculate long-range forces between particles, which are integrated with a second-order symplectic leap-frog integrator.
- The fourth-order Hermite integrator with block time steps (e.g., Aarseth 2003) is applied to integrate the orbits of stars and the centers-of-mass of multiple systems with short-range forces.
- The slow-down algorithmic regularization method (SDAR; Wang, Nitadori & Makino 2020) is used to integrate the multiple systems, such as hyperbolic encounters, binaries and hierarchical few-body systems.

This readme provide a complete and short documentation to describe how to install and use the code. 
Please carefully read it first before asking questions to developers.
More detail of the algorithms are described in Wang et al. (2020; arXiv: https://arxiv.org/abs/2006.16560).
The Doxygen documentation for developers is under preparation.
In the test folder, there are two sample scripts (sample.sh and sample_galpy.sh) showing how to generate initial condition from mcluster, start the simulation and process the data to generate single, binary snapshots, core information and lagrangian radii. 
sample.sh includes binaries and stellar evolution (bse), sample_galpy.sh add Milky Way potential using galpy based on sample.sh.

The main body of PeTar is written in c++ language. The external modules can have different programme languages.
The data analysis tool is written in Python3. Users need the basic knowledge of Python to access the simulation data.
Especially, it is recommended to learn how to use the Python modules `numpy`, `dict` and `matplotlib`, which cover the most functions to read, process and plot data. 

## About the version
The master branch of PeTar is frequently updated. In order to avoid using the unstable version, users can check the version printed by running `petar -h`.
The version format has three types:
- Develop mode: [PeTar commit count]_[SDAR commit count]. In develop mode, the code can be used and should work properly for most conditions, some features are not fully tested. If users found an unphysical result, please report to the developers.
- Test mode: test_[PeTar commit count]_[SDAR commit count]. In test mode, the code is not confirmed to work properly, please don't use it for production.
- Release mode: r[PeTar release number]_r[SDAR release number]. In release mode, the code is a stable version. It is recommended to use. The release version is updated with a low frequency. Right now the first release is under preparing. 

## Content:
- [Install](#install)
    - [Dependence](#dependence)
    - [Environment](#environment)
        - [For supercomputer](#for-supercomputer)
    - [Make](#make)
        - [A few useful options of configure](#a-few-useful-options-of-configure)
            - [Install path](#install-path)
            - [Change MPI parallelization options](#change-mpi-parallelization-options)
            - [Manually choose compilers](#manually-choose-compilers)
            - [Disable OpenMP parallelization](#disable-openmp-parallelization)
            - [Use X86 with SIMD](#use-x86-with-simd)
            - [Use Fugaku A64FX architecture](#use-fugaku-a64fx-architecture)
            - [Use GPU (CUDA)](#use-gpu-cuda)
            - [Debug mode](#debug-mode)
            - [Use stellar evolution](#use-stellar-evolution)
            - [Use Galpy external potential library](#use-galpy-external-potential-library)
            - [Multiple options](#multiple-options)
- [Use](#use)
    - [petar commander](#petar-commander)
    - [CPU threads](#cpu-threads)
    - [MPI processors](#mpi-processors)
    - [GPU devices](#gpu-devices)
    - [Input data files](#input-data-files)
    - [Restart](#restart)
    - [Options](#options)
    - [Performance optimization](#performance-optimization)
        - [Tree time step](#tree-time-step)
        - [Outer changeover radius](#outer-changeover-radius)
        - [Neighbor searching radius](#neighbor-searching-radius)
        - [Multiple group radius](#multiple-group-radius)
        - [Adjust tree time step and radii](#adjust-tree-time-step-and-radii)
    - [Output](#output)
        - [Printed information](#printed-information)
        - [Output files](#output-files)
        - [Units](#units)
             - [PeTar units](#petar-units)
             - [Stellar evolution units](#stellar-evolution-units)
             - [External potential units](#external-potential-units)
             - [Output file units](#output-file-units)
     - [Warning and errors](#warning-and-errors)
         - [Hard energy significant](#hard-energy-significant)
         - [Large step warning](#large-step-warning)
         - [Hard dump with errors](#hard-dump-with-errors)
         - [Crash with assertion](#crash-with-assertion)
    - [Data format update for old versions](#data-format-update-for-old-versions)
    - [Reference](#reference)
    - [Help information](#help-information)
    - [Useful tools](#useful-tools)
         - [Initial input data file](#initial-input-data-file)
         - [Find tree time step](#find-tree-time-step)
         - [Parallel data process](#parallel-data-process)
         - [Movie generator](#movie-generator)
         - [Gether output files from different MPI ranks](#gether-output-files-from-different-mpi-ranks)
         - [Remove data after a given time](#remove-data-after-a-given-time)
         - [Data format transfer](#data-format-transfer)
         - [Input parameter file format update](#input-parameter-file-format-update)
         - [SSE and BSE based steller evolution tool](#sse-and-bse-based-steller-evolution-tool)
         - [Galpy tool](#galpy-tool)
    - [Data analysis in _Python3_](#data-analysis-in-python3)
         - [Reading particle snapshots](#reading-particle-snapshots)
         - [Check reading consistence](#check-reading-consistence)
         - [Get particle information](#get-particle-information)
         - [Make data selection](#make-data-selection)
         - [Using Matplotlib](#using-matplotlib)
         - [Save data](#save-data)
         - [The reference frame and coordinate system transformation](#the-reference-frame-and-coordinate-system-transformation)
         - [Merge two data](#merge-two-data)
         - [Reading binary snapshots](#reading-binary-snapshots)
         - [Reading triple and quadruple snapshots](#reading-triple-and-quadruple-snapshots)
         - [Reading Lagrangian data](#reading-lagrangian-data)
         - [Use Python help to obtain tool manuals](#use-python-help-to-obtain-tool-manuals)
- [Method](#method)
    - [Algorithm of integration](#algorithm-of-integration)
    - [Parallelization methods](#parallelization-methods)
- [AMUSE API](#amuse-api)

## Install
### Dependence
_FDPS_: https://github.com/FDPS/FDPS (please use v6.0 or v7.0; v7.1 may not work)

_SDAR_: https://github.com/lwang-astro/SDAR

These two libraries are necessary to compile _petar_. 
Notice that the newest version of _FDPS_ (v7.1) has an issue that might cause a crash of _petar_ with an assertion error of NaN check.
After git clone of _FDPS_, please checkout the old release v7.0 by 
```
git checkout v7.0
```
in the _FDPS_ directory.

If users want to use external potential, the _Galpy_ interface is available. Users need to install _Galpy_ either by 
```
pip3 install galpy
```
Or download the source code from https://github.com/jobovy/galpy.

If the source codes of these libraries are put in the same directory where the _PeTar_ directory exist, the configure script (see Section [make](#make)) can detect them automatically. Otherwise users need to provide the pathes of them by adding configure options
```
./configure --with-[code name in lower case]-prefix=[code path]
```

### Environment
To successfully compile the code, the C++ compiler (e.g. GNU gcc/g++, Intel icc/icpc, LLVM clang/clang++) needs the support of the C++11 standard. To use SSE/BSE package, a Fortran (77) compiler, GNU gfortran, is needed and should be possile to provide API to the c++ code, i.e., the libgfortran is required. Currently Intel ifort is not supported yet. The MPI compiler (e.g. mpic++) is required to use MPI. NVIDIA GPU and CUDA compiler is required to use GPU acceleration. The SIMD support is tested for the GNU, Intel and LLVM compilers. It is not tested for others, thus these three kinds of compilers are suggested to use. 
The Fugaku ARM A64FX architecture is also supported. 

All compilers should be searchable in the `$PATH` environment. For example, to use OpenMPI C++ compiler, `mpic++` should be available.
This can be checked by typing `mpic++ --version` in the terminal, it should print the version of current MPI C++ compiler. 
If the result suggests that the commander is not found, then users should properly install OpenMPI first. 

To use _Galpy_ and the analysis tools, the _Python3_ should be available. _Galpy_ also requires the _GSL_ library being installed and can be detected in the load library path.

#### For supercomputer
Genenarlly, the supercomputer provides multiple choices of compilers, including different versions of Intel and GNU compilers.
Before install _PeTar_, users should first check the detail how to correctly setup the compilers by reading the manual or ask the administrators of the supercomputer.

### Make
Once _FPDS_ and _SDAR_ are available, in the root directoy, use 
```
./configure
```
to check the local enviroment and automatically detect the compilers and features.
Once it finished, a summary log will be printed. 
Please read it to correctly set up the enviroment variables (the pathes for executable files and Python libraries).

The options for configure can be found by 
```
./configure -h
```
After configure, use 
```
make
make install
```
to compile and install the code.

The excutable files, _petar_ and _petar.[tool name]_, will be installed in [Install path]/bin.
1. _petar_ is the main routine to perform N-body simulations. It is actually a soft link to _petar.\*\*_, where the suffix represents the feature of the code based on the configure.
2. _petar.[tool name]_ are a group of tools for debugging, datafile initialization, performance optimization and data analysis. The details can be checked in the section [Useful tools](#useful-tools).
For all tools, the commander `[tool name] -h` show all options with descriptions.
Users should always check this first to understand how to properly use the tool.

The generated data files from a simulation can be analyzed by using the _Python3_ based data analysis module.
The _Python3_ module is installed in `[Install path]/include/petar`.
Please add `[Install path]/include` to the _Python_ include path (the environment variable, `$PYTHONPATH`) in order to import the code.

#### A few useful options of configure
##### Install path
```
./configure --prefix=[Install path]
```
If the code is already installed before and the executable file (petar.\*\*) exists in the `$PATH` enviroment, the configure automatically use the same directory for installing.
   
##### Change MPI parallelization options
```
./configure --with-mpi=[auto/yes/no]
```
- auto (default): automatical detection MPI, if exists, use MPI compiler, otherwise use non-MPI compiler
- yes: use MPI c++ compiler
- no: non-MPI c++ compiler
   
##### Manually choose compilers 
Configure will detect the C++, C and Fortran compiler in the default `$PATH` environment. If users want to manually choose these compilers, the CXX, CC and FC can be modified, respectively. For example, when users want to use Intel C++ and C compilers with Intel MPI:
```
CXX=mpiicpc CC=mpiicc ./configure
```
where `mpiicpc` is Intel C++ MPI compiler and `mpiicc` is Intel C MPI compiler. If these compilers are not in the `$PATH` environment, the full path is needed, such as:
```
CXX=[Full path of mpiicpc] CC=[Full path of mpiicc] ./configure
```

##### Disable OpenMP parallelization
```
./configure --disable-openmp
```
By default OpenMP is used.
    
##### Use X86 with SIMD
```
./configure --with-simd=[auto/avx/avx2/avx512dq]
```
- auto (default): automatical detection based on local CPU architecture 
- avx: use AVX for tree force calculation and tree neighbor counting
- avx2: use AVX2
- avx512dq: use AVX512F and AVX512DQ
    
This option switch on the SIMD support for force calculation, the _auto_ case check whether the compiler (GNU or Intel) support the SIMD instructions and choose the newest one. Notice that the supported options of the compiler and the running CPU are different. Please check your CPU instruction whether the compiled option is supported or not. If the CPU can support more than the compiler, it is suggested to change or update the compiler to get better performance.

##### Use Fugaku A64FX architecture
```
./configure --with-arch=fugaku
```

The tree force and neighbor search functions using Fugaku A64FX instruction set are supported now. 
Notice that in Fugaku supercomputer, the configure only work on the running nodes. 
Users should launch an interactive job to configure and compile the code.
    
##### Use GPU (CUDA)
```
./configure --enable-cuda
```
By default GPU is not used. To switch on it, make sure the NVIDIA CUDA is installed and consistent with the c++ compiler.
    
##### Debug mode
```
./configure --with-debug=[assert/g/no]
```
- assert: switch on all assertion check
- g: switch on compiler option '-g -O0 -fbounds-check' in order to support debugger such as gdb
- no: no debugging, optimized performance (default)   
   
##### Use stellar evolution
```
./configure --with-interrupt=[bse/mobse/bseEmp]
```
Currently there are three options of stellar evolution packages based on SSE/BSE (Hurley et al. 2000, MNRAS, 315, 543; 2002, MNRAS, 329, 897):
- bse: the updated SSE/BSE version from Banerjee et al. 2020, A&A, 639, A41.
- mobse: the MOSSE/MOBSE from Giacobbo et al. 2018, MNRAS, 474, 2959.
- bseEmp: the udpated SSE/BSE version from Tanikawa et al. 2020, MNRAS, 495, 4170.

Notice that here all SSE/BSE package names (here after [bse_name]) only contain 'bse', but the SSE package is also included.

When this option is switched on, the standalone tool _petar.[bse_name]_ will also be compiled and installed.
This is a c++ based tool which call the stellar evolution functions to evolve single and binary stars to the given age and metallicity. OpenMP parallelization is used to speed up the calculation if a large group of stars and binaries are provided.

Currently, to use the extreme metal poor evolution track of bseEmp, users should create a soft link in the running directory to 'bse-interface/bseEmp/emptrack/ffbonn'.
Otherwises, the simulation will crash with a file IO error.

When SSE/BSE packages are used, users can control whether to switch on stellar evolution during the simulation by using _petar_ option `--stellar-evolution` and '--detect-interrupt' for single and binary evolution, respectively.
When `--stellar-evolution 2` is used, the dynamical tide for binary stars and hyperbolic gravitational wave energy/angular momentum loss for compact binaries are switched on.
But currently the dynamical tide is still an experimental function, it is not confirmed that the result is always physical. 
In default (`--stellar-evolution 1`), dynamical tide is not switched on.

##### Use _Galpy_ external potential library
```
./configure --with-external=galpy
```
The _Galpy_ library is a _Python_ and _c_ based external potential library, which provides a plenty choices of potentials. 
It is also flexible to combine multiple potentials together (require to use _Galpy_ _Python_ interface to generate the instance, see their document in details).

When this option is switched on, the standalone tool _petar.galpy_ and _petar.galpy.help_ will also be compiled and installed.
- _petar.galpy_ is a simple tool to call _Galpy_ c interface to evaluate the acceleration and potentials for a list of particles with a given potential model.
- _petar.galpy.help_ is a tool (python script) to help users to generate the input options for potential models. When users use _Galpy_ _Python_ interface to design a specific potential, this tool also provides a function to convert a _Galpy_ potential instance to an option or a configure file used by _PeTar_.

##### Multiple options
Multiple options should be combined together, for example:
```
./configure --prefix=/opt/petar --enable-cuda
```
will install the executable files in /opt/petar (this directory requires root permission) and switch on GPU support.


## Use:

### petar commander
After installation, if the `[Install path]/bin` is in the envirnoment variable, `$PATH`, the standard way to use the code is
```
[OMP_STACKSIZE=128M] [OMP_NUM_THREADS=N_threads] [mpiexec -n N_mpi] petar [options] [particle data filename]
```
where `[mpiexec -n N_mpi]` indicates the number of MPI processors (`N_mpi`); `OMP_NUM_THREADS=N_threads` indicates each MPI processor use `N_threads` CPU threads.

### CPU threads

To avoid segmetantional fault in simulations in the case of large number of particles, 
make sure to set the `OMP_STACKSIZE` large enough, for example, 128M is a good choice. 
In addition, make sure that the maximum stack size is unlimited by executing the commander `ulimit -s`.
It should return 'unlimited'. If not, execute `ulimit -s unlimited` before using petar.

To have a better performance, `N_threads` should not be too large (generally <=8).

A convenient way to set `OMP_NUM_THREADS`, `OMP_STACKSIZE` and `ulimit -s`:
Add the lines
```
export OMP_STACKSIZE=128M
export OMP_NUM_THREADS=N_threads
ulimit -s unlimited
```
in the shell configure/initial file (e.g. .bashrc for _bash_).
Then, next time when a shell is open or `source ~/.bashrc`, the two environment variables are automatically set and the maximum stack size is set to unlimited.
Thus, `OMP_STACKSIZE=128M OMP_NUM_THREADS=N_threads` can be removed from the _petar_ commmander line.

### MPI processors
In a multi-nodes cumputing cluster, the way to set `N_threads` and `N_mpi` can be different, please refer to the documentation of the computing cluster.

Notice that `N_mpi x N_threads` should be less than the total available CPU threads in the computing facility.

### GPU devices

When GPU exists, each MPI processor will launch one GPU job. 
The modern NVIDIA GPUs support multiple jobs in one GPU.
Thus, it is no problem that `N_mpi` is larger than the number of GPUs.
But when `N_mpi` is too large, the memory of GPU may not be enough and an error of Cuda Memory allocation may appear. 

When Multiple GPUs exist, each MPI processor will use different GPU based on the processor and GPU IDs.
If users want to only use one specific GPU, the environment variable, CUDA_VISIBLE_DEVICES=[GPU index], can be used. 
For example, 
`CUDA_VISIBLE_DEVICES=1 petar [options] [particle data filename]`
will use the second GPU in the computer (the index counts from 0). 
This variable can also be set in the shell initialization file with `export`.

### Input data files

If users have the initial file that stores masses, positions and velocities of each particles per line, generated by tools like _mcluster_ (https://github.com/lwang-astro/mcluster), the _petar.init_ tool (see [Initial input data file](#initial-input-data-file)) can be used to transfer it to the input file of the _petar_ style.

### Restart

All snapshots of particle data generated in the simulation can be used to restart the simulation. 
To restart the simulation with the same configuration of parameters as before, use
```
[OMP**] [mpiexec**] petar -p input.par [options] [snapshot filename]
```
where _input.par_ is automatically generated from the last run (stored in the same diretory of the simulation).
If users change the name of this file and restart, the next new parameter file also uses the new filename.
Notice if users want to use new options in additional to _input.par_, these options must be put after "-p input.par", otherwises they are rewrited by values stored in _input.par_.

In default, after restart,  the snapshot files with the same name will be replaced; but for other output files, new data will be appended to the existing ones (e.g., filenames with the suffixes: esc, group ...);.
If users want to restart from a snapshot in the middle of a simulation, _petar.data.clear_ can be used to remove the data with time > _t_ in the output files.

### Options
For _petar_, two types of options exist: the single-character options starting with '-' and the long options starting with '--'.
It is better for users to first check all single-character options listed by `petar -h`.
Here a few useful options are listed.
-  -u: the unit set of input data. When `-u 1` is used, the initial snapshot data should use the units of Msun, pc, pc/myr. Otherwise the gravitational constant (G) is assumed to be 1 (Henon unit). When use _petar.init_ tool, users should also be careful to set the consistent unit (see `petar.init -h`). There is no scaling of unit in _petar_, only G can be modified. The option `-u 1` set the value of G in the unit set of [Msun, pc, myr] and when stellar evolution and galactic tidal field are used, it also set the unit scaling factors. 
-  -t: the finishing time. Notice that the tree time steps used in petar is the integer power of 0.5. So if the finishing time does not satisify this criterion. The simulation may not finsh at the exact given time.
-  -o: the time interval for outputing data (snapshot and status). Should be also the integer power of 0.5.
-  -s: the tree time step, should be the integer power of 0.5. When this is set, the changeover radius is automatically determined (see paper for reference), unless `-r` is used manually.
-  -r: set the outer boundary of changeover radii. When this is set, tree time step is automatically determined.
-  -a: data appending option. If `-a 1` is used (default), the new data will be appended to existing files after restart, otherwise files are rewritten. Be careful when you want to restart the simulations with this option.
-  -b: the initial binary fraction number. This is important to correctly calculate the velocity dispersion, which is used to automatically determine tree time step and changeover radii. So it must be correctly set when primordial binaries exist. Besides, all primordial binaries should appear first in the input data file (two neighbor lines per pair).
-  -w: output style. Users can control the whether to output snapshot data. If `-w 2` is used, all particle data will be printed in one line together with a status information of the system per output time. This can be convenient for data analysis when N is small. 
-  -i: the format of snapshot data. This can be used to determine whether read or write snapshot data in BINARY or ASCII format.
-  -G: gravitational constant. 

Notice that `-s` and `-r` significantly affect the accuracy and performance of the simulations. 
Users should be careful for these two options. The _petar.find.dt_ tool can help users to find the optimized values for star clusters.
For a deep understanding and a better configuration, users may need to read the reference paper.

When stellar evolution packages (e.g., BSE) and external potential (e.g., Galpy) are used, the corresponding options are also shown in `petar -h`.

### Performance optimization

The performance of _PeTar_ is sensitive to a few important parameters and also depends on the initial condition of particle (stellar) systems. 
To achieve the best performance for a given input model, users need to adjust the following parameters carefully:

#### Tree time step 
 - _petar_ option `-s`

The tree time step is a fixed time step to calculate the long-range (particle-tree) force. 
The long-range force calculation is one of the most expensive computing with O(N log N). Thus, a smaller time step indicate a more expensive computing per physical time unit. 
However, users cannot increase tree time step too much, see the following discussion about the changeover and neighbor searching radii.

#### Outer changeover radius 
 - _petar_ option `-r`
 
The changeover region is the overlap shell between the long-range and short-range interaction. 
Below the inner region, the short-range interactions are calculated by using 4th order Hermite with individual time steps and SDAR method. 
Above the outer region, the long-range interactions are calculated by the 2nd or 4th-order LeapFrog with particle-tree method. In between, both short- and long-range interactions are calculated. 
The inner and outer radii are mass weighted (`m^(1/3)`) for each particle. The default ratio is 0.1 (_petar_ option `--r-ratio`).

Changeover radii and tree time steps should be consistent to ensure a physical result of simulation (See _PeTar_ paper for detail). A simple way to understand this is to check a circular kepler orbit of two particles with the semi-major axis within the changeover region. Inside the changeover region, the forces between the two particles are splitted to the short-range and the long-range forces. The short-range forces are evaluated every Hermite time step while the long-range forces are calculated every tree time step.
The Hermite time step is smaller than the tree time step. 
Thus, the long-range forces are given as velocity kicks to two particles after a few Hermite time step.
If the tree time step is too large, such as only a few steps per Kepler orbit, the time resolution of long-range forces or velocity kicks is too low that the orbit integration is not accurate (non Kepler orbit). Thus, once the changeover ragion is determined, the tree time step should be small enough to ensure at least a few tens of sampling points to calculate the long-range forces (kicks) along the Kepler orbit. 

Without specifying `-s` and `-r`, _petar_ autodetermines the tree time step and the changeover radii assuming the input model is a sphercial symmetric star cluster with a King or a Plummer like density profile. 
Thus, for a more complex input model without A spherical structure, these parameters may need to be determined manually, following the self-consistent rule (tree time step - changeover radius relation) as described in the above Kepler orbit case. 

#### Neighbor searching radius 
 - _petar_ option `--r-search-min`

The changeover region determines the boundary between the short- and the long-range interactions. In _PeTar_, we define another neighbor searching radius for each particle to detect neighbor candicates that are expected to be inside changeover region during the next tree time step. 
The minimum searching radius is slightly larger than the outer changeover radius, while the used searching radius also depends on the velocity of the particle. If the particle has a high velocity, it can travel a large distance during one tree time step, thus its searching radius should be also large. 

The neighbor searching radius is one important parameter that significantly affects the performance. For the parallal computing of short-range integration, _petar_ first use neighbor searching radius to gather nearby particles into individual clusters. For each cluster, only one CPU core is assigned for the short-range (Hermite+SDAR) integration during the next tree time step. The clustering algorithm ensures that for each member in the cluster, all its neighbors are also inside the cluster. Thus, if the neighbor radii are too large, one huge cluster can form where most of particles in the system are inside. Then, only one CPU core is used to integrate this cluster while all other cores are waiting. This serious load balance issue wastes computing power and can significantly reduce the performance. Thus, for each particle system, the neighbor searching radius cannot be too large.

When the changeover radius (`-r`) is determined, _petar_ automatically calculate the neighbor searching radius. If users want to manually determine the neighbor searching radius, the following options are be set:
 - `--r-search-min`: set the minimum neighbor searching radius reference, the final radius is also mass weighted, similar to the changeover radius.
 - `--search-peri-factor`: set the maximum peri-center criterion assuming two neighbors have a Kepler orbit
 - `--search-vel-factor`: the coefficient to determine the velocity dependent addition to neighbor radius (base neighbor radius + cofficient * velocity * tree time step)

#### Multiple group radius 
 - _petar_ option `--r-bin`

The third important radius that influences the performance is the radius to determine a multiple group where the SDAR method is assigned for.
In a dense stellar system, a multiple group, such as a binary, a triple and a quadruple, can frequently appear. It has much shorter orbital periods of inner members compared to those of single stars orbiting in the host particle system.
The SDAR method is the most important algorithm in _petar_ to ensure the accuracy and the efficiency to integrate their orbits.
The binary stellar evolution is also treated inside the SDAR method. 
The criterion to select group members is the group radius (`--r-bin`). This is autodetermined based on the changeover inner radius.
If the radius is too large, the SDAR method can be very expensive because too many members are selected in a multiple group.
Meanwhile, if the radius is too small, some binaries are not integrated accurately since the Hermite integrator has a systematically long-term drift of energy and angular momentum for the periodic motion.

#### Adjust tree time step and radii 

For starting a new simulation, the autodetermined tree time step and the three radii may not be the best choice.
Users can use the tool `petar.find.dt` (see [Find tree time step](#find-tree-time-step)) to select the best tree time step.
This tool only works with the autodetermination of tree time step and changeover radii of _petar_ (see [Outer changeover radius](#outer-changeover-radius)).


When the structure of the particle system signficiantly evolve after a long time, users may want to redetermine the tree time step and the three radii discussed above to improve the performance. If users want to only change the tree time step and let _petar_ to autodetermine the radii, the following options are required to restart the simulation:
```
[OMP**] [mpiexec**] petar -p input.par -s [new tree time step] -r 0 --r-search-min 0 --r-bin 0 [other options] [snapshot filename for restart]
```
Where `-r 0 --r-search-min 0 --r-bin 0` are used to reset all three radii and switch on the autodetermination based on the new tree time step.
Users can also use `petar.find.dt` to select the best restart tree time step (see [Find tree time step](#find-tree-time-step)).


### Output
#### Printed information

When _petar_ is running, there are a few information printed in the front:
- The _FDPS_ logo and _PETAR_ information are printed, where the copyright, versions and references for citing are shown. The version of PeTar is the commit numbers of PeTar and SDAR on GitHub.
- The switched on features (selected during configure), such as stellar evolution packages, external packages, GPU usage.
- The input parameters are listed if they are modified when the corresponding options of _petar_ are used.
- The Unit scaling for _petar_, stellar evolution packages (e.g. _bse_) and external packages (e.g. _galpy_).
- A short parameter list for tree time step and a few radii that influences the performance.
- If the Galpy is used, the galpy potential setup may be printed.
- If the SSE/BSE based stellar evolution package is used, the common block and global parameters are printed.
- The name of dumped files for the input parameters are shown, in default, they are input.par, input.par.hard, input.par.[bse_name] (if a SSE/BSE based package is used)
    - input.par: input parameters of petar, can be used to restart the simulation from a snapshot.
    - input.par.hard: input paremeters of hard (short-range interaction part; Hermite + SDAR), can be used for _petar.hard.debug_ to test the dumped hard cluster.
    - input.par.[bse_name]: the SSE/BSE based package parameters, if a SSE/BSE based package is used, this is the necessary file to restart the simulation and also for _petar.hard.debug_.
    - input.par.galpy: the Galpy parameter, this is used for restart the simulation.

Then after the line "Finish parameter initialization",
The status of the simulation is updated every output time interval (the option -o).
The content of the status has a style like (One example for the output information at time 1 is shown together):
- Time, number of real particles, all particles (including artificial particles), removed particles and escaped particles; in local (first MPI process) and in global (all MPI process)
    - N_real(loc): number of physical particles (stars) in the local MPI processor (rank 0)
    - N_real(glb): number of physical particles (stars) in all MPI processors
    - N_all(loc): number of all particles including physical and artifical ones in the local MPI processor
    - N_all(glb): number of all particles in all MPI processors
    - N_remove(glb): number of removed particles (e.g. zero-mass particles due to mergers and escpaers) 
    - N_escape(glb): number of escaped particles (outside escape criterion)
```
Time: 1  N_real(loc): 1378  N_real(glb): 1378  N_all(loc): 1378  N_all(glb): 1378  N_remove(glb): 0  N_escape(glb): 0
```
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
    - Error_PP: energy error in short-range particle-particle (PP) interaction 
    - Error_PP_cum: cumulative energy error in PP part
```
Energy:       Error/Total           Error       Error_cum           Total         Kinetic       Potential          Modify    Modify_group   Modify_single        Error_PP    Error_PP_cum
Physic:      1.883442e-05       -644.6969       -644.6969   -3.422972e+07    1.846359e+07   -5.269331e+07        1841.216               0     0.008989855    -8.67599e-06    -8.67599e-06
Slowdown:    1.883442e-05       -644.6969       -644.6969   -3.422972e+07    1.846359e+07   -5.269331e+07        1841.216               0     0.008989855    -8.67599e-06    -8.67599e-06
```
- Angular momentum: error at current step, cumulative error, components in x, y, z directions and value
```
Angular Momemtum:  |L|err: 187484.2  |L|err_cum: 187484.2  L: -1.17383e+07   -1.20975e+07    -1.267862e+09  |L|: 1.267974e+09
```
- System total mass, center position and velocity
```
C.M.: mass: 736.5417 pos: 5127.807   -5729.159    7.270318 vel: -165.7037   -150.611    3.006304
```
- Performance information:
    - Tree step number per output interval
    ```
    Tree step number: 512
    ```
    - Wallclock computing time for each part of the code, two lines are printed indicate the minimum and the maximum times among all MPI processors
        - Total: total computing time per tree time step
        - PP_single: time to integrate single particles without neighbors
        - PP_cluster: time to integrate (short-range) particle clusters using Hermite+SDAR method inside one MPI processor
        - PP_cross: time to integrate (short-range particle clusters using Hermite+SDAR method crossing mutliple MPI processors
        - PP_intrpt*: time used during the interruption of integration (only used in special mode)
        - Tree_NB: time to searching neighbors using particle-tree method
        - Tree_Force: time to calculate long-range force using particle-tree method
        - Force_corr: time to correct long-range force due to the usage of changeover functions
        - Kick: time to kick velocity of particles from the long-range forces
        - FindCluster: time to create particle clusters for short-range interactions
        - CreateGroup: time to find mutliple groups in each clusters
        - Domain_deco: time of domain decomposition for global particle tree
        - Ex_Ptcl: time to exchange particles between different domains
        - Output: time to output data and print information
        - Status: time to calculate global status (e.g., energy) of the system
        - Other: time cost in other parts that are not included in above components.
    ```
    **** Wallclock time per step (local): [Min/Max]
  Total        PP_single    PP_cluster   PP_cross     PP_intrpt*   Tree_NB      Tree_Force   Force_corr   Kick         FindCluster  CreateGroup  Domain_deco  Ex_Ptcl      Output     Status       Other
     0.009742    0.0006008   0.00016328   2.8949e-06            0   0.00069683     0.004435   2.6859e-05   3.3679e-05    8.365e-05   0.00011302   2.0602e-06   4.7604e-05   5.4633e-05   1.8945e-08    0.0034471
    0.0097421   0.00060092   0.00016398   3.0875e-06            0   0.00070318    0.0044427   2.7012e-05   3.3787e-05   8.3776e-05   0.00011312   2.1375e-06   4.7723e-05   5.4636e-05   1.9141e-08    0.0034481
    ```
    - FDPS tree force calculation profile (see FDPS document for detail)
    ```
    **** FDPS tree soft force time profile (local):
  Sample_ptcl  Domain_deco  Ex_ptcl      Set_ptcl_LT  Set_ptcl_GT  Make_LT      Make_GT      SetRootCell  Calc_force   Calc_mom_LT  Calc_mom_GT  Make_LET_1   Make_LET_2   Ex_LET_1     Ex_LET_2     Write_back
            0            0            0   5.2026e-05   1.7637e-06   0.00021123   7.0745e-05   8.7232e-06    0.0022072   4.4457e-05    0.0018011   1.3547e-05            0    4.401e-06            0   2.0096e-05
    ```
    - FDPS tree neighbor searching profile 
    ```
    **** Tree neighbor time profile (local):
  Sample_ptcl  Domain_deco  Ex_ptcl      Set_ptcl_LT  Set_ptcl_GT  Make_LT      Make_GT      SetRootCell  Calc_force   Calc_mom_LT  Calc_mom_GT  Make_LET_1   Make_LET_2   Ex_LET_1     Ex_LET_2     Write_back
            0            0            0   3.6579e-05   1.6174e-06   0.00021966   7.2012e-05   7.7332e-06   0.00025779    3.126e-05    2.935e-05   2.6922e-05            0   1.4563e-06            0   1.2669e-05
    ```
    - Number counts per tree time step
        - PP_single: number of single particles
        - PP_cluster: number of particles in clusters within one MPI processor
        - PP_cross: number of particles in clusters crossing mutliple MPI processors
        - PP_intrpt*: number of particles suffering interruption of integration
        - Cluster: umber of clusters within one MPI processor
        - Cross: number of clusters crossing multiple MPI processors
        - AR_step_sum: number of SDAR integration steps
        - AR_tsyn_sum: number of SDAR integration for time sychronization
        - AR_group_N: number of multiple groups
        - Iso_group_N: number of isolated multiple groups (a cluster containing only one group)
        - H4_step_sum: number of Hermite integration time steps
        - H4_no_NB: number of Hermite integration for particles without neighbors
        - Ep-Ep_sum: number of i and j particle interaction during a particle-tree force calculation
        - Ep-Sp_sum: number of i particle and j superparticle interaction during a particle-tree force calculation
    ```
    **** Number per step (global):
  PP_single    PP_cluster   PP_cross     PP_intrpt*   Cluster      Cross        AR_step_sum  AR_tsyn_sum  AR_group_N   Iso_group_N  H4_step_sum  H4_no_NB     Ep-Ep_sum    Ep-Sp_sum
       1345.6       32.436            0            0       14.744            0       7.0586       4.0996            0            0       329.74            0   1.6972e+06        34134
    ```
    - Histogram of number of members in clusters. The first line show the number of members and the second line show the histogram counts
    ```
    **** Number of members in clusters (global):
            1            2            3            4            5            6
       1345.6        13.24      0.68164      0.31055      0.40234      0.10938
    ```

The performance information is very useful to check whether the simulation has a proper setup of tree time step and radii parameters. 
For a reasonable performance, the Tree_Force wallclock time should dominate the total time. If many multiple groups, such as primordial binaries, exist, the PP_cluster can be also time-consuming. But if it is too large (dominate most of the total time), users should consider to radius the changeover, neighbor searching and group radii (see [Performance optimization](#performance-optimization)).

The histogram of number of members in clusters is also useful to identify whether the neighbor searching radius is too large. 
For a low density system, the maximum number of members should be within 20, like the above example (6).
For a high density system or a system with high-velocity particles, the maximum number may be large.
But if it is within a few hundreds and the PP_cluster wallclock time is not significantly large, it is also fine.

#### Output files
Except the printed information from _petar_ commander, there are a few output files as shown in the following table:

| File name            | Content                                                                                                                    |
| :-------------       | ------------------------------------------------------------------------------------------------------------------------   |
| data.[index]         | Snapshot files for each output time. The format is the same as the input data file.                                        |
|                      | Users can find the definitions of the first line (header) and columns by using `petar -h`.                                 |
|                      | [index] is the output order, counting from 0 (initial snapshot). It is not the time unless the output interval is (`-o`) 1 |
| data.esc.[MPI rank]  | The escaped particle, the columns are the same as snapshot files except for an additional column of escaped time at first  |
| data.group.[MPI rank]| The information of new and end of mutiple systems (binary, triple ...) that are detected in the SDAR integration method.   |
|                      | The definition of the multiple system depends on the distance criterion of _petar_ option `--r-bin`.                       |
| data.status          | The global parameters, e.g., energies, angular momentum, number of particles, system center position and velocity          |
| data.prof.rank.[MPI rank]| The performance measurement for different parts of the code during the simulation                                      |

When the SSE/BSE stellar evolution options (--with-interrupt during configure) are used, there are additional files: 

| File name            | Content                                                                                                                    |
| :-------------       | ------------------------------------------------------------------------------------------------------------------------   |
| data.[sse_name].[MPI rank] | The single stellar evolution records, such as type changes and supernovae. Notice that if a star evolve very fast (less than the dynamical integration time step), the internal type changes may not be recorded|
| data.[bse_name].[MPI rank] | The binary stellar evolution records. All binary type changes are recorded, but if 'Warning: BSE event storage overflow!' appears in the simulation, the binary type change is too frequent so that these changes for the corresponding binary is not recored|

When the Galpy is used (--with-external=galpy), for each snapshot file, the galpy parameter file may also be generated in some conditions:
| File name           | Content                                                                                                                     |
| :-------------      | ------------------------------------------------------------------------------------------------------------------------    |  
| data.[index].galpy  | the galpy parameter file for each snapshot, may need to be used for restart                                                 |


Here 'data' is the default prefix of output files, users can change it by _petar_ option `-f`. For example, by using `-f output`, the output files will be 'output.[index]', 'output.esc.[MPI rank]' ...

[MPI rank] indicates which MPI processor outputs the data. Thus, the number of each data files with this suffix is the same as the MPI processors used in the simulation.
For example, when two MPI processors are used, the escape data files will be 'data.esc.0', 'data.esc.1'.
Before access these files, it is suggested to run _petar.data.gether_ tool first to gether the separated files with different MPI ranks to one file for convenience.

_petar.data.gether_ not only gether the files from different MPI processors, but also generate new files that can be accessed by _petar_ Python data anaylsis tool  (See [Data analysis in _Python3_](#data-analysis-in-python3)) .
- For ".group" files, _petar.data.gether_ separates the few-body groups with different number of members to individual files with suffix ".n[number of members in groups]".
- For ".[sse_name]" and ".[bse_name]" files, this tool separates the type changes, supernova kicks and dynamical mergers into separate files (".type_change, .sn_kick and .dynamic_merge").
To check the details of the generated files from _petar.data.gether_ and cooresponding reading method in Python, see [Gether output files from different MPI ranks](#gether-output-files-from-different-mpi-ranks).

Since each file contains a large number of columns, it is recommended to use the Python data analysis tool to access them. 
The data analysis tool is very convenient to pick up the specific parameter (column) and process the data (selection, calculation and plotting), similar to the `dict` and `numpy` in Python.
It also helps to avoid reading wrong columns. 
Thus, the defnition of columns are not listed in the manual and are also not in the header of files.
Please refer to the help information of the corresponding analysis tool for each file.

The raw snapshots files do not include the information for binaries or multiple systems. To identify them and obtain the Lagrangian and core properties,
_petar.data.process_ tool can be used (see [Parallel data process](#parallel-data-process)).

#### Units

There are 1-3 sets of units in _petar_ depending on the used packages:

##### Petar units

The PeTar without additional stellar evolution and external potential packages follows the unit of the input files. There is no unit conversion inside the PeTar part. Only the gravitational constant can be modified to be consistent with the units from the input data file (See [Options](#options)). Thus, the input data must have consistent unit set: the velocity unit must be consistent with the length unit. For example, the length unit is pc, the velocity unit must be pc/[time unit]. It cannot be km/[time unit] because an additional unit conversion from km to pc is always needed in the integration. This only brings complexities and bugs in the code without any obvious benefit. Thus, petar only allows to change gravitationa constant to keep a self-consisten unit system.

Since no unit conversion exists, the output log, snapshot data (petar part), escaper fiels and group files follow the unit set used in the input data with the corresponding gravitational constant.
For example, if the input data has the unit (Myr, pc, pc/myr), the time and energy in the output log also use the same unit set. The potential energy also includes the gravitational constant.

##### Stellar evolution units

The stellar evolution package (e.g., _bse_) is an additional code. It has a different unit system.
Thus, between the PeTar part and the stellar evolution pacakge, there is an unit conversion.
The conversion scaling factors can be manually defined by the _petar_ options `--bse-rscale` ... (see `petar -h` for details).
It is recommended to use the unit set for input data (Msun, pc, pc/Myr) so that the _petar_ option `-u 1` can be used and the petar will calculate the conversion factor automatically.
The stellar evolution part in the snapshot files and the output files " [data filename prefix].[s/b]se.[MPI rank]" follow the same unit system as that of _[bse_name]_. 

##### external potential units

The external potential package (_galpy_) is also an additional code. There is another unit conversion between the PeTar part and it.
The conversion factors can also be manually defined. For _galpy_, the corresponding options are `--galpy-rscale` ... (see `petar -h` for details).
It is also recommended to use the unit set for input data (Msun, pc, pc/Myr) so that the conversion factor is automatically calculated.
In default unit of Galpy refers to the solar motion, where the length unit is the distance of Sun to the Galactic center, and the velocity unit is the solar velocity in the Galactic frame.

##### Output file units

Here is the table to show the corresponding units for the output files (The data filename prefix is "data" as an example).
| Filename                    | Content                                                     | Unit set                                                           |
| :-------------------------- | :---------------------------------------------------------  | :------------------------------------------------------------      |
| _petar_ output log          | printed information from _petar_                            | PeTar unit                                                         |
| data.[index]                | snapshots                                                   | Particle class: PeTar unit + Stellar evolution unit (see `petar -h`) |
| data.esc.[MPI rank]         | escapers                                                    | Time: PeTar unit; particle: particle class                         |
| data.group.[MPI rank]       | multiple systems                                            | Binary parameters: PeTar unit;  Particle members: particle class   |
| data.[bse_name].[MPI rank]  | binary stellar evolution events                             | Stellar evolution unit                                             |
| data.[sse_name].[MPI rank]  | single stellar evolution events                             | Stellar evolution unit                                             |
| data.status                 | global parameters                                           | PeTar unit                                                         |
| data.prof.rank.[MPI rank]   | performance profiling                                       | Time: second (per tree time step)                                  |

Here 'particle class' indicates data structure of one particle (c++ class) in PeTar. 
Depending on the stellar evolution mode, one particle has the mixed data from PeTar part and stellar evolution part.
These two parts are not in the same unit system.
The units for stellar evolution parameters with the prefix "s_" in the particle class can be found by using `petar -h`. 
The other members follow the units of input data file (PeTar unit).
If users are unclear about the units of the output files. They can also check the help information from the python analysis tool (`help(petar.Particle)`) for reading the corresponding files. 

### Warning and errors

In a simulation, there may be warnings and errors appearing in the printing log. 
Here lists the frequent warnings and errors.

#### Hard energy significant

For one particle cluster with short-range interaction (Hermite+SDAR), the relative energy error after one tree time step exceeds the limitation defined by the option `--energy-err-hard` (the default value is 0.0001). 
Then, a warning message with the title "Hard energy significant" appears and a corresponding file "hard_large_energy.\*" is dumped.

In most of cases, this happens when a triple or a quadruple system exists in the cluster.
There two possibilities that cause the large error: 
- The inner binaries of these systems may be very tight with a large slowdown factor. Since slowdown Hamiltonian is not the same as claasical Hamiltonian. The physical energy is not necessary conserved during the integration because the slowdown methods only ensure that the secular motion is correct. 
- The integration step of SDAR method may not be sufficient small for the multiple system, while reducing step sizes also signficantly increase the computing time.
Either case is difficult to solve if users don't want a significant reduction of computing performance. 
But if users care about the specific particles inside this cluster, they can check whether the orbit of integration is acceptable by using the debuging tool _petar.hard.debug_ to read the dumped hard_large_energy.\* file.
The high energy error is usually generated by a small change of semi-major axis of the most tight binary in the system. 
However, such a small error does not significantly affect the global dynamical evolution of the whole system if it the error only happens once.

#### Large step warning

Sometimes the performance of the code significantly drops with a warning of large steps appearing.
Then, a file "dump_large_step.\*" is dumped.
In this case, there are probably a stable multiple system in the particle cluster.
The AR method is used to integrate the multiple system, but the step count is very large (exceeding the step limit, which can also be set in the input option) so that the performance significantly drops.
This situation cannot be avoided sometimes. 
When the stellar evolution is switched on, it can help in some cases when the inner binaries are tight enough to merge. 
There is no better solution to solve the issue.
If multiple CPU cores are used, users can restart the simulations with less CPU cores. 
Sometimes after restarting, the same stable system may not form so that the problem can be avoided.
If it still forms, the parallel computing is killed so it is not necessary to use multiple CPUs. 
Users can use less CPU resources to pass this phase until the system is disrrupted, and then, restart with the original number of CPU cores.

#### Hard dump with errors
When errors appear in the Hermite-SDAR integration, the error message appears and the simulation is halted.
In this case, a file "hard_dump.\*" is dumped.
This usually indicates that a bug exist in the code.
Users can report this issue by contacting the developer via GitHub or email.
In the content, users need to describe the version of PeTar, the configure options, the initial conditions of the simulations and attach the "hard_dump.\*", the input parameter files (starting with "input.par.\*".

Users can also check the details by using the debug tool _petar.hard.debug_ together with the GDB tool, if users prefer to understand the problems themselves. The knowledge of the source codes of SDAR is required to understand the messages from the debug tool.

#### crash with assertion

Sometimes, the code may crash with an assertion information.
One frequent assertion error is `n_jp<=pg.NJMAX`, e.g.
```
CalcForceEpEpWithLinearCutoffSimd::operator()(const EPISoft*, ParticleSimulator::S32, const EPJSoft*, ParticleSimulator::S32, ForceSoft*): Assertion `n_jp<=pg.NJMAX' failed.
```
This error happens when the particle system has a extremely density contrast, or the density center is not near the coordinate origin (zero) point or a particle is extremely far away from the system.
The most distant particles in the system deterine the outmost boxsize when the particle-tree is constructed.
If the maximum boxsize is too large, the minimum tree cell size is also large that too many particles may exist in one cell that exists the limitation of particle-tree algorithm. Thus, an assertion `n_jp<=pg.NJMAX` appears.
The tree cell near the coordinate origin point has the highest resolution, thus if the density center is far from the origin point, this assertion also can happen. 
Therefore, the solution is to remove very distant particles or choose the coordinate origin point properly.  

The second frequent assertion error is `!std::isnan(vbk.x)`, e.g.
```
SystemHard::driveForOneClusterOMP(ParticleSimulator::F64): Assertion `!std::isnan(vbk.x)' failed.
```
This error happens when FDPS version 7.1 is used. It seems a bug or an inconsistent interface exists in FPDS that can easily cause such an assertion for _petar_.
The solution is to use the FDPS version 7.0.

### Data format update for old versions
The data formats of snapshots, input parameter files and a part of output files have been updated in the past.
If users want to use new version of code to read the old version data, the data transfer is possible.

For snapshot data in ASCII format, after [the version on Aug 8, 2020](https://github.com/lwang-astro/PeTar/commit/0592d70875626071e1bd7aa13dbab30165a98309#diff-25a6634263c1b1f6fc4697a04e2b9904ea4b042a89af59dc93ec1f5d44848a26), the output format of group_data.artificial changes from 64bit floating to 64bit integer in order to keep the full information.
The BINARY format is not affected.
In order to read the old snapshot data, it is needed to transfer the data first by 
```
petar.format.transfer -g [other options] [snapshot path list filename]
```
The new data with BINARY format are created, using the same tool with option `-b`, users can transfer the BINARY format back to the new ASCII format.

The formats of input parameter files generated during simulations (including files from _SSE_/_BSE_ and _Galpy_) update on Oct 18, 2020.
Use [`petar.update.par`](#input-parameter-file-format-update) to update the input files in order to restart the simulations with the newer versions of _PeTar_.
After the update, the reading and modification of input parameter files become much easier.

### Reference
Remember to cite the necessary references when you publish the results using _PeTar_. The references are shown in the help function of _petar_ and the begining of the output after a simulation starts.
When a feature imported from an external library is switched on, e.g. (_SSE_/_BSE_, _Galpy_), the corresponding references are automatically added in the output.

### Help information

All options are listed in the help information. This can be checked by using the commander
```
petar -h
```
The description of the input particle data file is also shown in the help information. 
Before using _PeTar_ and its tools, it is suggested to read the help information first to avoid mistakes.
When new versions of _PeTar_ release, the help information always has the corresponding update.

### Useful tools
There are a few useful tools helping users to generate initial input data, find a proper tree time step to start the simulations and data analysis.
Each of the tools are stored in the user defined install_path/bin.
The name of tools have a common style: _petar.[tool name]_.
For each tool, users can get help how to use it by
```
petar.[tool name] -h
```
Be careful that options with the same name may have a different meaning for each tool.

#### Initial input data file
PeTar has an internal Plummer model generator (equal mass system) by using Henon Unit and half-mass radius being 1.0.
If users want to use their own initial particle data, _petar.init_ can help to transfer their particle data to _petar_ intput data file.
The usage:
```
petar.init [options] [particle data filename]
```
The particle data file should contain 7 columns: mass, position[3], velocity[3] and one particle per row.
All Binaries come first and the two components should be next to each other. Besides, when binaries exist, the option '-b [binary number]' should be added to the _petar_ commander. This is important to obtain a correct initial velocity dispersion, tree time step and changeover radii. 

Notice when stellar evolution is switched on, the corresponding options "-s [bse_name]" should be used together to generate the correct initial files. In this case, the astronomical unit set (Solar mass [Msun], parsec [PC] and Million year [Myr]) are suggested to use for the initial data. Notice the velocity unit should be PC/Myr. Then the corresponding mass scaling factor is 1.0 between units of PeTar and SSE/BSE based code. Besides, the option "-u 1" should be added to the _petar_ commander in order to use this astronomical unit set.

Similarly, when external mode (potential) is switched on, the option '-t' should be used to generate correct number of columns.

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
Thus in the test, if the user specify the minimum step size by `-s [value]` (outside `-a`), the step size will be regularized if it does not satisfy this requirement.  
Be careful that some _petar_ options, such as `-o` and `-s`, cannot be used inside `-a` of _petar.init_. 
See the details by using `petar.find.dt -h`.

For restarting a simulation, users many want to find the new tree time step and determine other parameters (radii) automatically, this can be done by
```
petar.find.dt -m 2 -o 4 -a "-p input.par -r 0 --r-search-min 0 --r-bin 0" [restart snapshot filename]
```

#### Parallel data process
The _petar.data.process_ can be used to process snapshot data to detect binaries, triples and quadruples (binary-binary type), and calculate Langragian, core radii, averaged mass and velocity dispersion.
Notice that the tool calculate the core (density) center and use it to obtain Langrangian radii.
The single and binary data are stored for each snapshot with additional suffix ".single" and ".binary".
The triples and quadruples are optional and they are not used in the calculation of Lagrangian properties. 
The data of Lagrangian, core and escapers are generated in separate files.
The multiple CPU cores are used to process the data since the KDTree neighbor search for calculating density and detect binaries is very slow.
The basic way to use:
```
petar.data.process [options] [snapshot path list filename]
```
Users should provide a file that contains a list of pathes for the snapshot data files.
This file can be generated by using the _petar.data.gether_ tool. 
Manually, users can generate the file by using _ls_ and _egrep_ in Linux.
For example, if the snapshot filename prefix is 'data', the commander is
```
ls |egrep '^data.[0-9]+$'
```
Notice it is better to sort the path in the increasing order of evolution time.
The data order in Lagrangian, escapers and core data file follows the order in the snapshot path list.
The _sort_ tool is very convenient to get the time-sorted list:
```
ls |egrep '^data.[0-9]+$' |sort -n -k 1.6 >snap.lst
```
will find all data files in the current directory, sort files by using the suffix (values after 'data.') with an increasing order and save the list to the file of _snap.lst_.
Here '-n' indicates that the values to sort are floating-point number.
'-k' defines which word (separated by space) and which character in the word is the starting position of the number for sorting.
In this example, '1.6' represents the 6th chararcters in the first word is the starting position to recognize the numbers for sorting.
This is because only one word exists per line in the snap.lst file, and the word 'data.\*\*' has 5 characters before the number (index of the snapshot). 

Users should be careful to set the correct options of gravitational constant (`-G`), interrupt mode (`-i`) and external mode (`-t`) for _petar.data.process_ (see `petar.data.process -h`).
This is important to correctly read the snapshots and calculate the Kepler orbital parameters of binaries.
When users apply the astronomical unit set (`-u 1` in the _petar_ commander) in simulations, `-G 0.00449830997959438` should be used for _petar.data.process_.
Or, if the _BSE_ based package is used, the interrupt mode option `-i [bse_name]` can also set the correct value of G.
When the external mode such as _galpy_ is used, the external mode option `-t galpy` is needed.

Here is the table to show the files generated by _petar.data.process_ and the corresponding python modules (class) to read them.
The the default filename prefix 'data' is assumed.

| Filename            | Content                                                              | The Python analysis classes for reading data files  |
| :------------------ | :------------------------------------------------------------------  | :------------------------|
| data.[\*].single     | snapshots for single stars                                           | petar.Particle(interrupt_mode=[\*], external_mode=[\*])| 
| data.[\*].binary     | snapshots for binaries                                               | petar.Binary(member_particle_type=petar.Particle, simple_mode=[\*], G=[\*], interrupt_mode=[\*], external_mode=[\*]) |
| \*data.[\*].triple    | snapshots for triples (if option `-M` is used)                       | petar.Binary(member_particle_type_one=petar.Particle, member_particle_type_two=[petar.Particle, petar.Particle],  simple_mode=[\*], G=[\*], interrupt_mode=[\*], external_mode=[\*])            |
| \*data.[\*].quadruple | snapshots for binary-binary quadruples (if option `-M` is used)      | petar.Binary(member_particle_type=[petar.Particle, petar.Particle],  simple_mode=[\*], G=[\*], interrupt_mode=[\*], external_mode=[\*])             |
| data.lagr           | Lagrangian and core properties for all objects, single, binaries and user defined stellar types (`add_star_type`)  | petar.LagrangianMultiple(mass_fraction=[\*], calc_energy=[\*], external_mode=[\*], add_star_type=[\*]) |
| data.core           | core position, velocity and radius                                   | petar.Core()               |
| data.esc_single     | Single escapers                                                      | petar.SingleEscaper(interrupt_mode=[\*], external_mode=[\*])   |
| data.esc_binary     | Binary escapers                                                      | petar.BinaryEscaper(member_particle_type=petar.Particle, simple_mode=[\*], G=[\*], interrupt_mode=[\*], external_mode=[\*])      |
| \*data.bse_status   |  The evolution of number counts, maximum and averaged masses of different stellar types | petar.BSEStatus() | 
| \*data.tidal         | Tidal radius data (if option `--r-escape tidal` is used)             | petar.Tidal()              |

PS: the arguments [\*] for the keywords in the python class initialization depend on the configure options for compiling (e.g., `interrupt_mode`, `external_mode`) and the options used in _petar.data.process_. These keywords arguments are optional. They are only needed to provide when the default values are not used.
See help of the Python analysis classes and `petar.data.process -h` for details.

How to use the Python analysis classes to read data are discussed in [Data analysis in _Python3_](#data-analysis-in-python3).

The snapshots (single, binary ...) generated by _petar.data.process_ are shifted to the rest frame where the density center is the coordinate origins.
By adding the core position and velocity from 'data.core' at the corresponding time, the positions and velocities in the initial frame or Galactocentric frame (when _galpy_ is used) can be recovered.

When the snapshot files are in BINARY format, the option `-s binary` can be used for _petar.data.process_ to read the snapshot correctly.
Notice that the generated data from _petar.data.process_ are all in the ASCII format.

Notice that the _petar_ code can also remove escapers and stored the data of escapers by using the energy and distance criterion (in files [data filenam prefix].esc.[MPI rank], see [Output files](#output-files)).
All escapers during the simulations are not stored in the originial snapshot files.
The _petar_ code only has a simple constant escape radial criterion.
Instead, the post-process by _petar.data.process_ can calculate the tidal radius and use it to detect escapers.
These detected escapers will be stored in the post-generated escape files: data.esc_single and data.esc_binary. 
Be careful that these escapers are not removed from the post-generated snapshot files: data.[index].single and data.[index].binary.

For the Lagrangian properties, 'data.lagr' include the radius, average mass, number of objects, different components of velocity and dispersions within different Lagrangian radii. 
The mass fraction of Lagrangian radii is 0.1, 0.3, 0.5, 0.7, 0.9 in default.
The property of core radius is added at the last.
There is an option in _petar.data.process_ to define an arbitrary set of mass function.
When the SSE/BSE based stellar evolution package is used, an additional option `--add-star-type` can be used to calculate the Lagrangian properties for specific type of stars.
See `petar.data.process -h` for details.
When `--add-star-type` is used, the reading function should has the consistent keyword arguments.
The example to read 'data.lagr' is shown in [Reading Lagrangian data](#reading-lagrangian-data).

When `--calc-energy` is used, the potential energy, external potential energy and virial ratio for each Lagrangian radii are calculated.
But be careful when the external potential is used, the virial radio may not be properly estimated when the stellar system has no well defined center (disrrupted phase).

#### Movie generator
The _petar.movie_ is a covenient tool to generate a movie from the snapshot files.
It can generate the movies of the positions (x,y) of stars (x, y of positions), the HR diagram if stellar evolution (SSE/BSE) is switched on, the 2D distribution of semi-major axis and eccentricity of binaries.
The snapshot file list is required to generate the movie.
The basic way to use:
```
petar.movie [options] [snapshot path list filename]
```
If users want to plot information of binaries, it is better to use _petar.data.process_ first to generate detect binaries with multiple CPU cores. Then the movie generator do not need to use expensive KDTree function to detect binaries (use option '--generate-binary 2').

This tool use either _imageio_ or _matplotlib.animation_ (_Python_ modules) to generate movies. It is suggested to install the _imageio_ in order to generate movie using mutliple CPU cores. This is much faster than the _matplotlib.animation_ which can only use one CPU core. On the other hand, The _ffmpeg_ library is also suggested to install in order to support several commonly used movie formats (e.g. mp4, avi). Notice that this is a standalone library, not a _Python_ module. Users should install it in the OS system (e.g. via the apt tool in Ubuntu).

Similar to _petar.data.process_, users should also set correct options of gravitational constant (`-G`), interrupt mode (`-i`) and external mode (`-t`).

#### Gether output files from different MPI ranks

When the MPI is used, each MPI processor generate individual data files.
Thus, some output filenames contains the suffix, [MPI rank]. 
The _petar.data.gether_ is used to gether these output files from different MPI ranks to one file.
In addition, it also helps to split some files into individual file so that the python tools can be used to read the data (see [Output files](#output-files)).
The tool also generates a file, "[output prefix].snap.lst" , that contains the list of all snapshot files sorted by time. This can be used as the input for _petar.data.process_ and _petar.movie_.
The basic usage is
```
petar.data.gether [options] [data filename prefix]
```
where [data filename prefix] is the prefix of data files defined by _petar_ option -f ('data' in default).
If the [output prefix] is not given (option -f of _petar.data.gether_), the [output prefix] is the same as [data filename prefix].

When _SSE_/_BSE_ is used and the code version is before Sep 10, 2020, the data with suffix [.dynamic_merge] has three column less than that of the new version.
The corresponding data are dr, t_peri and sd_factor.
This tool will automatically fill the missing columns with zero (from Column 6 to 8)

Here is the table of files that the tool will generate and the corresponding Python anlysis classes to read (see [Data analysis in _Python3_](#data-analysis-in-python3))
| Original files                | Output files                  | Content                                     | Python classes initialization for reading |
| :--------------               | :------------                 | :---------                                  | :-------------------------------------    |
| data.group.[MPI rank]         | data.group                    | all groups                                  | None                                      |
|                               | data.group.n2                 | binary, hyperbolic                          | petar.GroupInfo(N=2)                      |
|                               | data.group.n3                 | triples                                     | petar.GroupInfo(N=3)                      |
|                               | ...                           | multiple systems...                         | ...                                       |
| data.esc.[MPI_rank]           | data.esc                      | escapers                                    | petar.SingleEscaper(interrupt_mode=[\*], external_mode=[\*]) |
| data.[sse_name].[MPI rank]    | data.[sse_name]               | full single stellar evolution records       | None                                      |
|                               | data.[sse_name].type_change   | single stellar evolution type change events | petar.SSETypeChange()                     |
|                               | data.[sse_name].sn_kick       | supernova natal kick of single star         | petar.SSESNKick()                         |
| data.[bse_name].[MPI rank]    | data.[bse_name]               | full binary stellar evolution records       | None                                      |
|                               | data.[bse_name].type_change   | binary tellar evolution type change events  | petar.BSETypeChange()                     |
|                               | data.[bse_name].sn_kick       | supernova natal kick in binaries            | petar.BSESNKick()                         |
|                               | data.[bse_name].dynamic_merge | dynamical driven (hyperbolic) mergers       | petar.BSEDynamicMerge(less_output=[\*])   |

PS: [\*] is the argument that depends on the configure options that used for compiling. 
'data.group[.n\*]' files are not generated by default because they are huge files. To obtain that, it is necessary to add the option '-g' for _petar.data.gether_.

#### Remove data after a given time
The _petar.data.clear_ is used to remove data after a given time for all output files except the snapshots.
The basic usage is 
```
petar.data.clear -t [time] [data filename prefix]
```
All output files will be firstly backup by rename the file name with an additional suffix '.bk'.
If this tool is mistakely used, the files can be recovered.
If the tool is used again, it will check the backup data files.
Thus a new time criterion larger than the previous one can be to include more data in the new files.

#### Data format transfer
The _petar_ can read and write snapshot data in either BINARY or ASCII format.
The BINARY format is compressed (the file size is less than one half of the ASCII one) and much faster to read and write for both _petar_ and data analysis tool, but cannot be directly read by text editor.
It is suggested to use the BINARY format if the simulation generate a large amount of data and users want to use analysis tools to access the data.
In the class _petar.Particle_ and _petar.PeTarDataHeader_ of the _Python3_ analysis module, reading the BINARY format is supported.

On the other hand, it is also possible to transfer the snapshots data between BINARY and ASCII format by using the tool _petar.format.transfer_:
```
petar.format.transfer [options] [snapshot path list filename]
```
The snapshot file list contain pathes of snapshots that users want to transfer.
In default, new files are generated with a suffix '.B' or '.A'.
To replace the file in order to save space, users can use the option "-r".
This tool can also update the old version of snapshot in ASCII format (before Aug 8, 2020) to the new versions (option -g).
Notice that in the old version, the information stored in the group_data (group c.m. mass and velocities in 64bit floating) is lost.
This is not important since the data can be calculated from data processing.
Also it does not affect restart.

Be careful that the version of _petar.format.transfer_ and _petar_ should be consistent (the same configuration of interrupt mode and external mode).
If snapshots are generated by different version of _petar_, _petar.format.transfer_ can fails to read data or provides wrong tranferred data.

#### Input parameter file format update
The formats of input parameter files generated during simulations (including files from _SSE_/_BSE_ and _Galpy_) update on Oct 18, 2020.
In order to use new version of _PeTar_ to restart the simulations using the old version of input parameter files, the update is needed:
```
petar.update.par [options] [input parameter filename] 
```
Using options `-p`, `-b` and `-t`, users can update input parameters from _PeTar_, _SSE_/_BSE_ and _Galpy_, respectively.
Additional options are used for different features chosen in configure.

After update, the new input parameters files are easier to read.
There are three columns in the file defined as (1) data type of argument, (2) option names and (3) values of arguments.
The reference of first two columns are shown by using the commander `petar -h`.
Users can directly modify the values of arguments in the file.
Besides, it is not necessary to list all options there.
Thus, if new options are implemented in the future version, there is no need to update the file again (unless existing option names change).

#### SSE and BSE based steller evolution tool
The _petar.[bse_name]_ will be generated when the SSE/BSE based stellar evolution package (--with-interrupt=[bse_name]) is used during the configuration.
This is the standalone tool to evolve stars and binaries to a given time.
All global parameters needed can be set in options.
The basic way to evolve a group of stars:
```
petar.[bse_name] [options] mass1, mass2 ...
```
If no mass is provided, the tool can evolve a group of stars with equal logarithmic mass interval.

The basic way to evolve a group of binaries:
```
petar.[bse_name] [options] -b [binary data file]
```
where the binary data file contain 4 values (mass1, mass2, period, eccentricity) per line.
Notice that the first line has only one value which is the number of binaries.
The units of mass and period depend on the option '--mscale' and '--tscale'. In defaulted case, the units set are Msun and Myr. For example, when the tscale are not 1.0, the period in Myr is calculated by the input period value x tscale. 

##### Output files
If there is only one star or one binary, the printed information from the commander line show its evolution history of the type changes and supernove events.
If more than one stars or binaries exist, only the final status are printed.
Additional files that record the evolution history of all stars or binaries are generated if the option `-o` is used.
The argument of this option define the filename prefix.
For example, when `-o output` is used, if single stars exist, 'output.[sse_name].type_change' and 'output.[sse_name].sn_kick' are generated;
if binaries exist, 'output.[bse_name].type_change' and 'output.[bse_name].sn_kick' are generated.
The method to read these files are the same as that for the output files from _petar_ (see [Data analysis in _Python3_](#data-analysis-in-python3) and [Gether output files from different MPI ranks](#gether-output-files-from-different-mpi-ranks)).


#### Galpy tool
The _petar.galpy_ and _petar.galpy.help_ will be generated when --with-external=galpy is used during the configuration.

_petar.galpy_ is the standalone tool to calculate accelerations and potentials for a particle list with a given potential model.
The basic way to use it:
```
petar.galpy [options] [particle data filename]
```

_petar.galpy.help_ is a helper to generate --galpy-type-arg options used in _petar_.
By
```
petar.galpy.help
```
A list of supported potentials are shown with the corresponding --galpy-type-arg options (using default arguments).
Users can check the detail of one potential by
```
petar.galpy.help -d [potential name]
```
This will show the definition of this potential from _Galpy_ document.
The option '-o' can be used to generate a configure file of a given potential. This file can be read by the _petar_ commander option: --galpy-conf-file. 

Notice that in _petar_ commander, there are three options used to read the configuration of potential models: --galpy-type-arg, --galpy-set and --galpy-conf-file.
These three options can be used together, i.e., potential models from all three options will be added together to become a multiple potential. 
Thus be careful not to dupplicate potentials. Whether the potential is dupplicated (printed multiple times) can be checked in the output of _petar_.

### Data analysis in _Python3_
The data analysis library provide the tools to identify binaries; calculate Lagrangian radii and core radii; obtain system energy error and check performance of each parts of the code.
To use the tool, first 
```
import petar
```
in _Python_ script, _IPython_ or _Jupyter_. 

The tools contains classes and functions.
The structure of the class type is based on a mixture of _collection.OrderedDict_ and _numpy.ndarray_.
Each member of the class stores one type of data or a subset of data.
In the former case, the member is a 1D or 2D _numpy.ndarray_.
The 1D array is corresponding to one column in the data file
The array index is corresponding to the row index in data file.
If the first few lines are skiped by using the keyword argument (--skiprows), then it counts after these lines.
For example, the member _mass_ in _petar.Particle_ is a 1D _numpy.ndarray_, the size of array is the total number of particles.
All other one-dimensional members in the same class have the same size.

The 2D array contain multiple columns of the data file.
The shape of array is column number x row number.
For example, the member _pos_ in _petar.Particle_ is the 3D positions of particles, thus it has the shape of (size, 3).

The two special members in all classes are _size_ and _ncols_.
- _size_ is the size of one member of 1D _numpy.ndarray_ type.
- _ncols_ is the total number of columns of all data members, if a member has a 2 dimension array, such as _pos_ of 3 x _size_, it is counted as 3. 

Once users create a class instance, e.g. by `[class instance]=petar.[class name]`, 
they can use `[class instance].keys()` to obtain the list of members.

By using ```help([class instance])```, the names, types and definitions of members (keys) are shown.
The type "1D" or "2D" indicates that the member is _numpy.ndarray_. 
Otherwise, it is the name of the class type and the corresponding member is its class instance which contains a subset of data (multiple columns in data files).
Users can use `help([class instance].[member name])` to obtain the definitions of members in the subset of data.

In some classes, the list of members can change depending on the keyword arguments.
_Particle_ is one typical case. In its help information, members are separated by groups and the final member list are a combination based on the keyword arguments.
After the class instance is created, users can always use `[class instance].keys()` to check the actual list of members.

Here is the list of all classes.
- For reading outputs of _petar_ (need to use _petar.data.gether_ to generate data files first):

| Class name     | Description                                  | Keyword arguments (options shown in [])        |           Corresponding file                                |
| :------------  | :---------------------------------           | :----------------------------------------      | :---------------------------------------------------------- |
| petar.PeTarDataHeader| The header (first line) of snapshot data | external potential option: `external_mode=['galpy', 'none']` |         data.[index]                        |
|                |                                              | Data format: `snapshot_format=['binary', 'ascii']`  |                                                         | 
| petar.Particle | The basic (single) particle data             | stellar evolution option: `interrupt_mode=['bse', 'mobse', 'bseEmp', 'none']` | (single) snapshot files, data.[index], data.[index].single  |
|                |                                              | external potential option: `external_mode=['galpy', 'none']`    |                                        |
| petar.Status   | The global parameter (energy, N ...)         |                                                | data.status                                                 |
| petar.Profile  | The performance of different part of codes   | Whether GPU is used: `use_gpu=[True, False]`   | data.prof.rank.[MPI rank]                                   |
| petar.GroupInfo| Mutliple systems (binary, triple ...)        | Number of members: `N=[2, 3, ...]`             | data.group.n[number of members]                             |

- For outputs from _petar (need to use _petar.data.gether_ first) and _petar.[bse_name]_, when the SSE/BSE based stellar evolution mode is switched on:

| Class name           | Description                                  |   Corresponding file      |
| :------------        | :---------------------------------           | :------------------------ |
| petar.SSETypeChange  | Type change records of single stars          | data.[sse_name].type_change  |
| petar.SSESNKick      | Supernove kick events of single stars        | data.[sse_name].sn_kick      |
| petar.BSETypeChange  | Type change records of binary stars          | data.[bse_name].type_change  |
| petar.BSESNKick      | Supernove kick events of binary stars        | data.[bse_name].sn_kick      |
| petar.BSEDynamicMerge| Dynamically driven mergers                   | data.[bse_name].dynamic_merge|
| petar.BSEMerge       | Mergers of binaries, find mergers (using class function `combine`) from petar.BSETypeChange and petar.BSEDynamicMerge ||
| petar.tide           | Dynamical tide or hyperbolic gravitational wave radiation events | data.[bse_name].tide |

- For data generated by using the _petar.data.process_ (see [Parallel data process](#parallel-data-process) for keyword argument description):

| Class name         | Description            | Keyword arguments                              |           Corresponding file                                         |
| :------------      | :----------------------| :----------------------------------------      | :----------------------------------------------------------          |
| petar.SingleEscaper| Single star escapers   | same as petar.Particle                         | data.esc (from _petar_), data.esc_single (from _petar.data.process_) |
| petar.BinaryEscaper| Binary star escapers   | same as petar.Binary                           | data.esc_binary                                                      |
| petar.LagrangianMultiple| Lagrangian and core properties | mass_fraction, calc_energy, external_mode, add_star_type   | data.lagr                                   | 
| petar.Core         | Data of core radius, center position and velocity of the system  |      | data.core                                                            |
| petar.Binary       | Binary and multiple system data | member_particle_type(\_one|\_two), interrupt_mode, external_mode, simple_mode, G | data.[index].binary, data.[index].triple, data.[index].quadruple       |
| petar.BSEStatus    | The evolution of number counts, maximum and averaged masses of different stellar types || data.bse_status                                      |

There are also several useful functions, see `help(Function name)` in Python to check the manual for the arguments in detail.
| Function name      | Description                                                      |
| :---------------   | :--------------------------------------------------------------- |
| petar.join         | Join two data with the same type. For example, `petar.join(particle1, particle2)` will generate a new `petar.Particle` instance that include `particle1` and `particle2`|
| petar.findPair     | Detect binaries of from a particle snapshot by using _scipy.cKDTree_. This is the function used to detect binaries in _petar.data.process_ |
| petar.findMultiple | Detect triples and quadruples (binary-binary) from single and binary data. It is also used in _petar.data.process_                         |
| petar.parallelDataProcessList| Use mutliple CPU cores to process a list of snapshot files and generate single and binary snapshots, Lagrangian data, core data and escaper data. For large _N_, the data process is quite slow, thus using multiple CPU processors can speed up the process. This is the main function used in _petar.data.process_|
| petar.vecDot       | Dot production of two 2D array |
| petar.vecRot       | Rotate a 3D vector array by using Euler angles (using _scipy.spatial.transform.Rotation_) |
| petar.cantorPairing| generate new ID from two IDs, used to obtain the unique binary ID from IDs of two components |
| petar.calcTrh      | Calculate one-component half-mass relaxation time using Spitzer (1987) formula |
| petar.calcTcr      | Calculate half-mass crossing time |
| petar.convergentPointCheck | Calculate proper motions in the frame of convergent point based on the given velocity and calculate the residuals (van Leeuwen F., 2009; Jerabkova T. et al. 2021) |
| petar.coordinateCorrection | Correct c.m. coordinate based on the difference between snapshot center and observational center in the galactocentric frame |


More useful tools will be implemented in the future. The tools/analysis/parallel_data_process.py is a good example to learn how to use the analysis tool.

#### Reading particle snapshots

Here is one example to use the _Particle_ class to do data analysis for a snapshot.
The snapshot generated by _petar_ contain two parts: hearder and data
In ASCII format, the first line is header, the following are data of each particle per line.

For exmaple, the first snapshot at time zero generated by _petar_ is named as 'data.0' in the default case.
To obtain its header information, _petar.PeTarDataHeader_ class can be used:
```
import petar
header=petar.PeTarDataHeader('data.0')
```
Without the keyword argument `snapshot_format`, it is assumed the reading data is in the ASCII format. 
To read the BINARY format:
```
header=petar.PeTarDataHeader('data.0',snapshot_format='binary')
```
After reading, header contain three members: `fid`,`n` and `time`, which represent the file ID, number of particles and time (in unit as that of the input model) of the snapshot respectively.

If `external_mode=galpy`, the header contain additional members: the offset of position and velocity, `pos_offset` and `vel_offset`, which represent the shift of center of the system in the Galactic tidal field (assuming the coordinate origin is the galactic center).

To read the particle content in ASCII format:
```
particle=petar.Particle(interrupt_mode='bse')
particle.loadtxt('data.0',skiprows=1)
```
Here the keyword argument ```interrupt_mode``` is important to set properly in order to read the snapshot correctly.
The column definitions of snapshots depends on the stellar evolution option (--with-interrupt) and the external potential option (--with-external) used during the configure.
The argument 'bse' indicates that the updated SSE/BSE is used so that external columns exist in the snapshots.
Then, here `particle` contains a member `star` with the class type `petar.SSEStarParameter`.
Similary, if external potential is added, one more column of 'pot_ext' is added.

Since the first line in the snapshot file is the header, `skiprows=1` jumps this line when read data.

If the snapshot data is in the BINARY format, the _fromfile_ function should be used instead of _loadtxt_:
```
particle.fromfile('data.0', offset=petar.HEADER_OFFSET)
```
Here the _fromfile_ is similar to _numpy.fromfile_. The keyword argument `dtype` is defined implicitly, users should not change it. The keyword argument `offset` sets the offset to read the data in bytes. 
`petar.HEADER_OFFSET` is the constant to indicate the offset of snapshot header line.
This `offset` is equivalent to `skiprows=1` in the _loadtxt_ function.

#### Check reading consistence

To use the Python analysis tool to read the data, users must be careful to ensure the keyword arguments in the initialization function are consistent with the options used in simulations (_petar_) and data analysis tool (_petar.data.process_).
This is very important to avoid wrong data reading.
For example, if the bse is used during the configure (`--with-interrupt=bse`), additional columns of stellar evolution parameters are stored in snapshot files.
If the initialization of `petar.Particle` misses the keyword argument `interrupt_mode='bse'` as:
```
particle=petar.Particle()
particle.loadtxt('data.0',skiprows=1)
```
The `particle` will not have the correct data.
In this case, a warning message will appear:
```
UserWarning: The reading data shape[1] or the number of columns (34) mismatches the number of columns defined in the class instance (20)! Make sure whether this is intended and whether you properly choose the correct keyword arguments for the class instance initialziation
```
Here the snapshot file 'data.0' has 34 columns (except the first line) that includes the stellar evolution parameters, but `particle` has only 20 columns (ncols=20).
Such a mismatch can cause a wrong assignment of the data columns and the class members.
Thus, whenever it appears, users need to review whether they have forgotten or incorrectly set the keyword arguments for the initialization. 

#### Get particle information 

To get the number of particles (line numbers in snapshots excluding the first line):
```
print(particle.size)
```

#### Make data selection

To make a filter and create a subset of data with particle mass less than 1.0,
```
particle_set = particle[particle.mass<1.0]
```
This is similar to the _getitem_ function of numpy.ndarray.
This filter can also be used for the members:
```
pos_set = particle.pos[particle.mass<1.0]
```
which generates a 2D numpy.ndarray of particle positions with masses below 1.0.

#### Using Matplotlib 

To plot the mass-distance relation of the subset using _matplotlib_,
```
import matplotlib.pyplot as plt
fig, axes=plt.subplots(1,1)
particle_set.calcR2()
axes.plot(np.sqrt(particle_set.r2), particle_set.mass, '.')
```

#### Save data

To save the subset to a file (ASCII format)
```
particle_set.savetxt([file path of new data])
```
Notice that when additional members are added to the particle instance, the saved data will also include the additional column.
In this example, ```particle_set.calcR2()``` generate a new class member, _r2_, (distance square).
Thus the saved data will have an additional column, _r2_, at the end, which does not exists in the original snapshot.
When users read this saved data, they should add the _r2_ member first in order to read columns correctly:
```
import numpy as np
particle_new = petar.Particle(interrupt_mode='bse')
particle_new.addNewMember('r2',np.array([]))
particle_new.loadtxt([file path of saved data])
```

#### The reference frame and coordinate system transformation

The `petar.Particle` and `petar.Core` have the member function to transform the data type to _astropy.coordinate.SkyCoord_.
The _SkyCoord_ is a powerful Python module that can easily transform the reference frame and the coordinate system of particle data.
For example, when the _galpy_ is used, the simulation of _petar_ use the Galactocentric frame with the Cartesian coordinate system.
When no stellar evolution is used, the input unit is astronomical unit (Msun, pc, pc/Myr) and the output format is ASCII,
the script to transform one snapshot (no stellar evolution and output format is ASCII):
```
import petar
import astropy.units as units

# read snapshot data from petar
particle=petar.Particle(external_mode='galpy')
particle.loadtxt([snapshot path],skiprows=1)

# read c.m. offset in the Galactocentric frame
header=petar.PeTarDataHeader([snapshot path], external_mode='galpy')

# Get SkyCoord data type in the Galactocentric frame
particle_new = particle.toSkyCoord(pos_offset=header.pos_offset, vel_offset=header.vel_offset)
```
This script will generate a new data `particle_new` with the type of `astropy.coordinate.SkyCoord`.
Notice thats the solar position (`galcen_distance=8.0*units.kpc, z_sun=15.*units.pc`) and velocity (`galcen_v_sun=CartesianDifferential([10.0,235.,7.]*u.km/u.s)`) are assumed in _galpy_. 
These are different from the default values of `SkyCoord`. 
In the function `toSkyCoord`, the values are consistent with the assumption of _galpy_.

Then, it is convenient to plot the data (RA, DEC) in the ICRS frame by using Matplotlib:
```
import matplotlib.pyplot as plt
fig, axes=plt.subplots(1,1)
ra = particle_new.icrs.ra
dec = particle_new.icrs.dec
axes.plot(ra,dec,'.')
axes.set_xlabel('RA')
axes.set_ylabel('Dec')
```


#### Merge two data

To join two subset of particle data, particle_set1 and particle_set2, to one,
```
particle_merge=petar.join(particle_set1,particle_set2)
```
For example, if the two sets have sizes of 3 and 5, respectively, the new instance, particle_merge, has a size of 8.

#### Reading binary snapshots

To read a binary star snapshot generated from _findPair_ or _petar.data.process_,
Users should be careful for the member particle type, a safe way to read snapshot when _SSE/BSE_ is used:
```
p1 = petar.Particle(interrupt_mode='bse')
p2 = petar.Particle(interrupt_mode='bse')
binary =petar.Binary(p1,p2)
binary.loadtxt([binary data file path])
```
Or 
```
binary = petar.Binary(member_particle_type=petar.Particle, interrupt_mode='bse', G=petar.G_MSUN_PC_MYR)
binary.loadtxt([binary data file path])
```
Here in the 'bse' mode, the gravitational constant, `G`, should be given the value in the unit set of Msun, pc and myr (0.00449830997959438).
Besides, `petar.G_HENON` is the G in the Henon unit (=1), which is the defaulted value in `petar.Binary`.

#### Reading triple and quadruple snapshots

`petar.Binary` can also be used to read triple and quadruple snapshots.
To read the triple snapshots:
```
binary = petar.Binary(member_particle_type_one=petar.Particle, member_particle_type_two=[petar.Particle, petar.Particle])
binary.loadtxt([triple data file path])
```
where `petar.Binary` is used as a binary tree, where the first component is the outer single and the second is the inner binary

Similarly, the binary-binary quadurple data can be read as:
```
binary = petar.Binary(member_particle_type=[petar.Particle, petar.Particle])
binary.loadtxt([quadruple data file path])
```
where `member_particle_type` represent the type (binary) of both components.

#### Reading Lagrangian data

To read the Lagrangian data generated by _petar.data.process_:
```
lagr = petar.LagrangianMultiple()
lagr.loadtxt([data.lagr path])
```
This data file contains three subsets: all, single and binary, corresponding to Lagrangain data of all objects, singles and binaries, respectively.
Each member in these subsets are 2D numpy ndarray. 
In each row, the properties at different Lagrangian radii and core radii are stored.
In default, the mass fraction is (0.1, 0.3, 0.5, 0.7, 0.9), and thus, 6 members exist in each row with the last member being core property.
For example, to access 0.5 mass fraction Lagrangian radius (half-mass radius) and core radius of all objects:
```
rh = lagr.all.r[:,2] # half-mass radius
rc = lagr.all.r[:,5] # core radius
```

Similarly, to access the averaged mass at half-mass radius:
```
mrh = lagr.all.m[:,2]
```

The _petar.data.process_ has an option to include arbitrary set of mass function.
When the SSE/BSE based stellar evolution package is used, an additional option `--add-star-type` can be used to calculate the Lagrangian properties for specific type of stars.
See `petar.data.process -h` for details.
When `--add-star-type` is used, the reading data should has the consistent keyword arguments
For example,
```
petar.data.process -i bse --add-star-type BH,MS 
```
add additional subsets "BH" and "MS" in the lagragian data file.
It is necessary to add the corresponding keyword argument `add_star_type` when reading the data:
```
lagr = petar.LagrangianMultiple(add_star_type=['BH','MS'])
lagr.loadtxt([data.lagr path])
```
If this argument is missing, the data in lagr can be wrong.
Users should be careful to always check whether the _petar.data.process_ options and _LagrangainMultiple_ keyword arguments are consistent.
A useful way to check the consistence is to count the column numbers in data file and reading as described in [Check reading consistence](#check-reading-consistence)

#### Use Python help to obtain tool manuals

This example show how to read and use some of the petar Python classes, the way to use other classes is very similar, only the class member is different.
As described above, the information of class member can be checked by the help function, e.g.
```
help(petar.Particle)
help(petar.GroupInfo)
```

## Method:
### Algorithm of integration: 
The basic cycle of integration the orbits of particles are described in detail in the reference paper.
Here is the short description:
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

### Parallelization methods:
1. Tree construction is done by _FPDS_, using _MPI_ and _OpenMP_ .
2. Soft force calculation kernel use _SIMD_ (_AVX_, _AVX2_, _AVX512_) or _GPU_ (_CUDA_).
3. Hard calculation (Hermite, SDAR) use _OpenMP_ for the loop of clusters.

### AMUSE API:
Current support API: gravitational dynamics, gravity field, stopping conditions.
Now the official AMUSE version has included the petar as a module. 
When _petar_ is updated, it is suggested to update the version in AMUSE as well.

This can be done by use 
```
make distclean
make
```
in the petar module directory: [AMUSE path]/src/amuse/community/petar/.
Be careful that this commander deletes src directory and all original existing files in the petar module directory. 
If you have modified any file, make a backup first!
