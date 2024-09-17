```
    ██████╗ ███████╗████████╗ █████╗ ██████╗ 
    ██╔══██╗██╔════╝╚══██╔══╝██╔══██╗██╔══██╗
    ██████╔╝█████╗     ██║   ███████║██████╔╝
    ██╔═══╝ ██╔══╝     ██║   ██╔══██║██╔══██╗
    ██║     ███████╗   ██║   ██║  ██║██║  ██║
    ╚═╝     ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝
```

PeTar is an N-body code specifically designed for modeling collisional stellar systems, where factors such as multiplicity (binaries, triples, etc.) and close encounters play a crucial role in dynamical evolution. PeTar offers several key advantages:

- **Precise gravitational force modeling**: PeTar does not employ any softening of gravitational force, enabling accurate tracking of the orbital evolution of binaries, triples, and close encounters.

- **Incorporation of single and binary stellar evolution**: Within the N-body simulation, PeTar dynamically evolves the masses, radii, and stellar types of individual stars. It also tracks significant events like supernovae, mass transfer, common envelope interactions, and binary mergers caused by stellar coalescence/collision and gravitational wave emission.

- **Galactic potential inclusion**: PeTar allows for the simulation of tidal effects on stellar systems by incorporating the galactic potential.

- **Parallel computing capabilities**: PeTar leverages multi-CPU processors/threads and GPU acceleration to accelerate simulations. This enables handling of over $10^7$ particles with a $100\%$ binary fraction.

- **Interoperability with other codes**: PeTar can function as a module within other codes, facilitating the simulation of complex stellar environments. This includes compatibility with frameworks like AMUSE and SPH-based hydrodynamical code Asura-bridge.

The core N-body algorithm of PeTar encompasses three distinct methods, each tailored to the separation distance between individual pairs of particles (stars):

- **Barnes-Hut tree method**: This technique, introduced by Barnes and Hut in 1986, is utilized for computing long-range forces between particles. These forces are then integrated using a second-order symplectic leap-frog integrator.

- **Fourth-order Hermite integrator with block time steps**: PeTar employs this integrator, inspired by Aarseth's work in 2003, to precisely integrate the orbits of stars and the centers-of-mass of multiple systems. It is particularly effective for handling short-range forces.

- **Slow-down algorithmic regularization (SDAR) method**: PeTar leverages the SDAR method, developed by Wang, Nitadori, and Makino in 2020, to integrate the dynamics of close-distance multiple systems. This is especially useful for scenarios involving hyperbolic encounters, binaries, and hierarchical few-body systems.

The core implementation of PeTar is predominantly in the C++ programming language. However, external modules within PeTar can be written in various programming languages.
The data analysis tool accompanying PeTar is developed in Python 3. Users are required to have a fundamental understanding of Python to effectively interact with the simulation data.
In particular, users are encouraged to familiarize themselves with key Python modules such as `numpy`, `dict`, and `matplotlib`. These modules encompass a wide range of functions essential for tasks like data reading, processing, and visualization.

This README document serves as a concise yet comprehensive guide detailing the installation and utilization of the code. It is strongly recommended that users thoroughly review this documentation before reaching out to the developers with any queries.

For a deeper understanding of the algorithms employed, additional details can be found in the work by Wang et al. (2020; available on arXiv: https://arxiv.org/abs/2006.16560).
For developers seeking to understand the code structure, please consult the [Doxygen documentation](https://lwang-astro.github.io/PeTar/doc/html/index.html).

After completing the installation process, users can quickly get started by exploring three sample scripts located in the sample folder: [star\_cluster.sh](https://github.com/lwang-astro/PeTar/blob/master/sample/star_cluster.sh), [star\_cluster\_bse.sh](https://github.com/lwang-astro/PeTar/blob/master/sample/star_cluster_bse.sh), and [star\_cluster\_bse\_galpy.sh](https://github.com/lwang-astro/PeTar/blob/master/test/star_cluster_bse_galpy.sh). These scripts provide practical demonstrations of simulating a star cluster using the PeTar code. They cover tasks such as generating initial conditions using `mcluster`, running simulations, and processing data to produce single and binary snapshots, core information, and Lagrangian radii. Here is a brief description of each script:
- star\_cluster.sh: Simulates a star cluster for up to 100 Myr with 1000 stars initially, following the Kroupa (2001) IMF and including 95% primordial binaries (refer to the `mcluster` manual). This simulation uses only gravitational forces.
- star\_cluster\_bse.sh: Similar to sample.sh but includes single and binary stellar evolution (SSE/BSE) with a metallicity of Z=0.02.
- star\_cluster\_bse\_galpy.sh: Builds upon sample_bse.sh by incorporating the Milky Way potential from Galpy's MWPotential2014 (refer to Bovy 2015 for details).

Furthermore, users can access a Jupyter Notebook titled [data\_analysis.ipynb](https://github.com/lwang-astro/PeTar/blob/master/sample/data_analysis.ipynb), which provides examples of data analysis in Python. By running one of the sample scripts, users can subsequently refer to the demonstrations in this notebook to analyze the simulation results. The data analysis module in PeTar offers greater convenience compared to manually parsing the output files. It is advisable to leverage this module instead of crafting reading code from scratch.

# About the version

PeTar code is maintained across multiple branches, and for users seeking stability, it is advisable to obtain the released version. The master branch undergoes regular updates to introduce new features and address bugs. Users can opt for this version if they find the new features beneficial.

Other branches may not function correctly, and users are advised against interacting with them unless they have consulted the developer and comprehended the implementations and existing issues.

Users can retrieve the code version by running `petar -h` or using `petar` to start a simulation. The code version will be presented in the format `[petar version][suffix]_[SDAR version]`, where:
- `[petar version]` represents the commit count of the PeTar code.
- `[SDAR version]` indicates the commit count of the SDAR code.
- `[suffix]` signifies the development mode. If absent, it denotes the master version. If `[suffix]` is `e` or another term, it indicates an experimental or specialized version.

The major modes are as follows:

- **Release mode (r)**: This mode signifies a stable version recommended for general use. Release versions are infrequently updated, with the initial release currently being prepared.
- **Master mode**: This mode is suitable for most conditions and is expected to function correctly. Users can use it for production runs.
- **Develop mode (e)**: In this mode, the code should operate correctly but may contain bugs that lead to unexpected results. It is advisable to consult developers before use.
- **Test mode (test)**: This mode is unverified for proper functioning and should not be employed for production purposes.

The subsequent sections provide detailed explanations of the installation process and usage instructions. The final sections offer a brief overview of the methods employed in the code and introduce the AMUSE API.

# Content:

- [Installation](#installation)
    - [Dependence](#dependence)
        - [Galpy](#galpy)
        - [Code path](#code-path)
    - [Environment requirements](#environment-requirements)
        - [For supercomputer](#for-supercomputer)
    - [Configuration](#configuration)
        - [Installation Path](#installation-path)
        - [Modifying MPI Parallelization Options](#modifying-mpi-parallelization-options)
        - [Manual Compiler Selection](#manual-compiler-selection)
        - [Disabling OpenMP parallelization](#disabling-openmp-parallelization)
        - [Selecting CPU Architecture](#selecting-cpu-architecture)
        - [Selecting SIMD Instructions](#selecting-simd-instructions)
        - [Enabling GPU Acceleration](#enabling-gpu-acceleration)
        - [Debug mode](#debug-mode)
        - [Utilizing Stellar Evolution](#utilizing-stellar-evolution)
        - [Using External Potential](#using-external-potential)
        - [Combining Multiple options](#combining-multiple-options)
     - [Compilation and Installation](#compilation-and-installation)
- [Usage](#usage)
    - [Preparing the initial condition](#preparing-the-initial-condition)
    - [Starting a Simulation](#starting-a-simulation)
    - [Resuming a Simulation](#resuming-a-simulation)
    - [Using OpenMP](#using-openmp)
    - [Using MPI](#using-mpi)
    - [Using GPU](#using-gpu)
    - [Options](#options)
    - [Performance Optimization](#performance-optimization)
        - [Tree time step](#tree-time-step)
        - [Outer changeover radius](#outer-changeover-radius)
        - [Neighbor searching radius](#neighbor-searching-radius)
        - [Multiple group radius](#multiple-group-radius)
        - [Adjusting tree time step and radii](#adjusting-tree-time-step-and-radii)
    - [Output](#output)
        - [Printed Information](#printed-information)
        - [Output Files](#output-files)
        - [Units](#units)
             - [PeTar Units](#petar-units)
             - [Stellar Evolution Units](#stellar-evolution-units)
             - [External Potential Units](#external-potential-units)
             - [Output File Units](#output-file-units)
    - [Troubleshooting](#troubleshooting)
         - [Significant Hard energy](#significant-hard-energy)
         - [Large Step Warning](#large-step-warning)
         - [Hard Dump with Errors](#hard-dump-with-errors)
         - [Hard Debug Tool](#hard-debug-tool)
         - [Crash with Assertion](#crash-with-assertion)
    - [Data Format Update for Older Versions](#data-format-update-for-older-versions)
    - [Useful Tools](#useful-tools)
         - [Initial Input Data File with `petar.init`](#initial-input-data-file)
         - [Determining the Tree Time Step with `petar.find.it`](#determining-the-tree-time-step)
         - [Gathering Output Files with `petar.data.gether`](#gathering-output-files)
         - [Parallel Data Processing with `petar.data.process`](#parallel-data-processing)
         - [Gathering Specified Objects with `petar.get.object.snap`](#gathering-specified-objects)
         - [Movie Generator with `petar.movie`](#movie-generator)
         - [Data Removal after a specified time with `petar.data.clear`](#data-removal-after-a-specified-time)
         - [Data Format Conversion](#data-format-conversion)
         - [Update of Input Parameter File Format](#update-of-input-parameter-file-format)
         - [Stellar Evolution Tool based on SSE and BSE](#stellar-evolution-tool-based-on-sse-and-bse)
         - [Galpy tool](#galpy-tool)
    - [References](#references)
    - [Help information](#help-information)
    - [Python Data Analysis Module](#python-data-analysis-module)
         - [Reading particle snapshots](#reading-particle-snapshots)
         - [Checking Reading Consistency](#checking-reading-consistency)
         - [Obtaining Particle Information](#obtaining-particle-information)
         - [Data selection](#data-selection)
         - [Plotting Data](#plotting-data)
         - [Saving and Loading Data](#saving-and-loading-data)
         - [Reference Frame and Coordinate System Transformation](#reference-frame-and-coordinate-system-transformation)
         - [Merging Two Datasets](#merging-two-datasets)
         - [Class Functions](#class-functions)
         - [Reading Binary Snapshots](#reading-binary-snapshots)
         - [Reading Triple and Quadruple Snapshots](#reading-triple-and-quadruple-snapshots)
         - [Reading Lagrangian Data](#reading-lagrangian-data)
         - [Reading Stellar Evolution Outputs](#reading-stellar-evolution-outputs)
         - [Reading Group Information](#reading-group-information)
         - [Accessing Tool Manuals Using Python Help](#accessing-tool-manuals-using-python-help)
- [Method](#method)
    - [Integration Algorithm Overview](#integration-algorithm-overview)
    - [Parallelization Methods](#parallelization-methods)
- [AMUSE API Integration](#amuse-api-integration)

# Installation

This section describes how to install PeTar, including the required libraries and codes, the computer environment, and the installation options.

## Dependence

PeTar is built upon the _FDPS_ and the _SDAR_ codes. _FDPS_ offers the particle-tree & particle-particle method and MPI and OpenMP parallelization. _SDAR_ is an N-body code library that incorporates the slow-down algorithmic regularization and Hermite integrator for simulating the motion of particle clusters with short-range and close-distance interactions.

Please download the two codes from the following GitHub links:

- _FDPS_: https://github.com/FDPS/FDPS (please use v6.0 or v7.0; v7.1 may not work)

- _SDAR_: https://github.com/lwang-astro/SDAR

The latest version of FDPS (v7.1) has a known issue that could lead to a crash of PeTar with an assertion error related to NaN check. After cloning FDPS from Git, in the FDPS directory, switch to the previous release v7.0 using the command:
```shell
git checkout v7.0
```

### Galpy

To incorporate external galactic potentials in simulations, users can use the _Galpy_ code through an interface integrated into PeTar. To use Galpy, users should install it either by executing 
```shell
pip3 install --user galpy
```
In this scenario, PeTar can automatically detect _Galpy_.
If `pip3` is unavailable, users can mamually download the source code from https://github.com/jobovy/galpy and specify the code path in the configure command (refer to the following guide).

### Code path

If the source codes of these dependent libraries are located in the same directory as the _PeTar_ directory, the configure script (see Section [Compiling the code](#compiling-the-code)) can automatically detect them. Otherwise, users will need to specify their pathes by adding configure options:
```shell
./configure --with-[code_name_in_lower_case]-prefix=[code path] ...
```

For example, if the _PeTar_, _FDPS_ and _SDAR_ codes are placed in the directory with pathes: `/home/username/code/PeTar`, `/home/username/code/FDPS` and `/home/username/code/SDAR`
the automatic detection will function.

In a different scenario, such as when Galpy is used but not installed via `pip3` and the _Galpy_ source code is located in `/home/username/python/Galpy`, users will have to specify its path using:
```shell
./configure --with-galpy-prefix=/home/username/python/Galpy ...
```

## Environment Requirements

To compile the code successfully, the C++ compiler (e.g., GNU gcc/g++, Intel icc/icpc, LLVM clang/clang++) must support at least the C++11 standard.

Using MPI necessitates the MPI compiler (e.g., mpic++). NVIDIA GPU and CUDA compiler are essential for GPU acceleration. SIMD support has been tested for GNU, Intel, and LLVM compilers. Since it hasn't been tested for others, it's recommended to use these three compiler types. The Fugaku ARM A64FX architecture is also compatible.

To use the SSE/BSE stellar evolution package, a Fortran (77) compiler, GNU gfortran, is required. It should be capable of providing an API to the C++ code, i.e., libgfortran is necessary. Intel ifort is currently not supported.

To use Galpy and the analysis tools, Python3 must be installed. Galpy also mandates the GSL library to be installed and detectable in the load library path.

All compilers should be accessible in the `$PATH` environment. For instance, to employ the OpenMPI C++ compiler, `mpic++` needs to be available. This can be verified by entering `mpic++ --version` in the terminal, which should display the version of the current MPI C++ compiler. If the output indicates that the command is not found, users should install OpenMPI correctly.

### For supercomputer

Generally, the supercomputer offers various compiler options, such as different versions of Intel and GNU compilers. Before installing PeTar, users should ensure they correctly set up the compilers by reviewing the manual or consulting the supercomputer administrators.

## Configuration

Once the required libraries such as FPDS and SDAR are accessible, go to the _PeTar_ directory and run the following command:
```shell
[environment variables] ./configure [options]
```
This command will examine the local environment, automatically identify the compilers and features, with `[environment variables]` and `[options]` representing additional options to manage the features. Upon completion of the configuration process, a summary log will be displayed. It is important to review this log carefully to correctly choose the environment variables and options, including the paths for dependent libraries and compilers.

To view the available environment variables and options for configure, use the following command:
```shell
./configure -h
```

A few useful environment variables and options are presented as follows:

### Installation Path

To specify a custom installation path, use the following command:
```shell
./configure --prefix=[Installation path]
```

The default installation path set by the configure script is `/user/local`, which requires administrator permission to access. It is not recommended to install PeTar there unless all users on the machine need to use the PeTar code. To install the code in a different location, users can add the `--prefix` option. For example, to install the code in `/home/username/tools`, users can include the option in the configure command:
```shell
./configure --prefix=/home/username/tools
```

If PeTar has been previously installed and the executable file (`petar`) is already in the `$PATH` environment, configure will automatically use the same directory for installation.
  
### Modifying MPI Parallelization Options

To enable or disable MPI parallelization, use the following command:
```shell
./configure --with-mpi=[choices] 
```
where `[choices]` can be `auto`, `yes`, or `no`:

- auto (default): Automatically detect MPI. If MPI is available, the MPI compiler will be used; otherwise, a non-MPI compiler will be used.
- yes: use the MPI C++ compiler.
- no: use a non-MPI C++ compiler.

### Manual Compiler Selection

By default, configure will detect the C++, C, and Fortran compilers in the `$PATH` environment. If users prefer to manually specify these compilers, they can modify the environment variables `CXX`, `CC`, and `FC` accordingly. For instance, if users wish to use Intel C++ and C compilers with Intel MPI, they can use the following command:
```shell
CXX=mpiicpc CC=mpiicc ./configure
```
Here, `mpiicpc` represents the Intel C++ MPI compiler, and `mpiicc` denotes the Intel C MPI compiler. If these compilers are not in the `$PATH` environment, users must provide the full path, as illustrated below:
```shell
CXX=/home/username/tool/bin/mpiicpc CC=/home/username/tool/bin/mpiicc ./configure
```
In this example, the Intel MPI compilers are installed in `/home/username/tool/bin/`.

For Mac OS users, clang, clang++, and flang compilers can be used instead of GNU compilers, as illustrated below:
```shell
CXX=clang++ CC=clang FC=flang ./configure
```

### Disabling OpenMP Parallelization

By default, PeTar enables multi-threaded OpenMP parallelization. To disable OpenMP, use the following command:
```shell
./configure --disable-openmp
```

### Selecting CPU Architecture

PeTar can utilize SIMD-like instructions to optimize the performance of tree force calculation and tree neighbor counting. The currently supported CPU architectures for enabling this feature are Intel/AMD x86 and Fugaku ARM A64FX. Users can specify the architecture using the following command:

```shell
./configure --with-arch=[choices]
```

where `[choices]` include `x86` and `fugaku`:

- x86 (default): Intel and AMD CPU architecture, commonly used.
- fugaku: Fugaku supercomputer CPU architecture, supporting A64FX instructions.

Please note that on the Fugaku supercomputer, the configuration only applies to the active nodes. Users are advised to initiate an interactive job to configure and compile the code.

### Selecting SIMD Instructions

PeTar supports multiple SIMD versions to enhance the performance of tree force calculation and tree neighbor counting. Users can choose the SIMD version using the following command:

```shell
./configure --with-simd=[choices]
```

where `[choices]` can be `auto`, `avx`, `avx2`, or `avx512`, with the latter offering the highest speed:

- auto (default): automatically detects based on the local CPU architecture
- avx: uses core-avx-i (theoretical speedup of 4x for single precision, 2x for double precision)
- avx2: uses core-avx2 (8x for single, 4x for double)
- avx512: uses skylake-avx512 (AVX512F and AVX512DQ) (16x for single, 8x for double)

Please note that the supported SIMD options of the compiler and the running CPU may differ. Make sure to use the SIMD version supported by the running CPU.

In the case of a supercomputer, the host and computing nodes might feature distinct CPU architectures. The configure script detects the SIMD version based on the local CPU. It is advisable to verify whether the CPU instructions on the computing node support a superior SIMD choice and opt for that during compilation.

### Enabling GPU Acceleration

PeTar supports the utilization of GPUs based on the CUDA language to accelerate tree force calculations as an alternative speed-up method to SIMD acceleration. To enable this feature, use the following command:

```shell
./configure --enable-cuda
```

By default, the GPU is not utilized. To enable it, ensure that NVIDIA CUDA is installed and compatible with the C++ compiler.

### Debug Mode

If the code crashes or a bug is present, users can enable the debugging mode as follows:

```shell
./configure --with-debug=[choices]
```

where `[choices]` can be `assert`, `g`, or `no`:

- assert: enables all assertion checks to detect unexpected behaviors
- g: activates compiler options '-g -O0 -fbounds-check' to support debuggers like gdb
- no: disables debugging for optimized performance (default)

### Utilizing Stellar Evolution

Users can enable stellar evolution for stars and binaries using the following command:

```shell
./configure --with-interrupt=[choices]
```

where `[choices]` include `bse`, `mobse`, and `bseEmp`.
In this option name, 'interrupt' refers to the N-body integration being interrupted by external effects on particles.

There are currently three options for stellar evolution packages based on SSE/BSE (Hurley et al. 2000, MNRAS, 315, 543; 2002, MNRAS, 329, 897):

- bse: the updated SSE/BSE version from Banerjee et al. 2020, A&A, 639, A41.
- mobse: the MOSSE/MOBSE from Giacobbo et al. 2018, MNRAS, 474, 2959.
- bseEmp: the updated SSE/BSE version from Tanikawa et al. 2020, MNRAS, 495, 4170.

It is important to mention that while all SSE/BSE package names only contain 'bse', the SSE package is also encompassed within them. From now on, the SSE/BSE based package will be denoted in a universal form as '[bse_name]'.

Enabling this option will also compile and install the standalone tool _petar.[bse_name]_. This C++ based tool calls stellar evolution functions to evolve single and binary stars to a specified age and metallicity. OpenMP parallelization is utilized to accelerate calculations when handling a large group of stars and binaries.

To use the extreme metal-poor evolution track of bseEmp, users must create a symbolic link in the working directory to either the _ffbonn_ or _ffgeneva_ directory located in 'PeTar/bse-interface/bseEmp/emptrack/', depending on the selected stellar evolution track mode during the execution of PeTar. Failure to do this will lead to a file I/O error, causing the simulation to crash.

When utilizing SSE/BSE packages, users can control whether to activate stellar evolution during the simulation using the `petar` option `--stellar-evolution` and `--detect-interrupt` for single and binary evolution, respectively. When `--stellar-evolution 2` is specified, dynamical tide for binary stars and hyperbolic gravitational wave energy/angular momentum loss for compact binaries are enabled. It's worth noting that the dynamical tide is still an experimental feature, and its results may not always be physically accurate. By default (`--stellar-evolution 1`), dynamical tide remains inactive.

### Using External Potential

Users can incorporate external potential and force into particles by utilizing the following command:

```shell
./configure --with-external=[choices]
```

Currently, `[choices]` offers only one option: `galpy`. Additional options will be introduced in future versions.

The Galpy library is an external potential library based on Python and C, offering a wide range of potential choices. It allows for the flexible combination of multiple potentials (requiring the use of the Galpy Python interface to instantiate, as detailed in their documentation).

Enabling this option will also compile and install the standalone tools `petar.galpy` and `petar.galpy.help`:
- `petar.galpy` is a straightforward tool that utilizes the Galpy C interface to compute acceleration and potentials for a particle list using a specified potential model.
- `petar.galpy.help` is a Python script tool designed to assist users in generating input options for potential models. When designing a specific potential using the Galpy Python interface, this tool also offers a function to convert a Galpy potential instance into an option or a configuration file used by PeTar.

### Combining Multiple Options

When combining multiple options, they should be used together, as shown in the example below:
```shell
./configure --prefix=/opt/petar --enable-cuda
```
This command will install the executable files in /opt/petar (this directory requires root permission) and activate GPU support.

## Compilation and Installation

After configuring, execute the following commands:
```shell
make
make install
```
to compile and install the code.

The executable files, `petar` and `petar.[tool name]`, will be installed in _[Install path]/bin_.
1. `petar` serves as the primary routine for conducting N-body simulations. It is essentially a symbolic link to `petar.**`, with the suffix reflecting the code's features based on the configuration, such as `petar.avx2.bse`.
2. `petar.[tool name]` comprises a set of tools for debugging, initializing data files, optimizing performance, and analyzing data. Detailed information can be found in the section [Useful tools](#useful-tools). For each tool, running `petar.[tool name] -h` will display all available options along with descriptions. Users are advised to refer to this first to ensure correct tool usage.

If _[Install path]/bin_ is added to the environment variable `$PATH`, the executable files can be directly accessed from any directory in the Linux system using the following command:
```shell
export PATH=$PATH:[Install path]/bin
```

Analysis of simulation-generated data files can be performed using the Python3-based data analysis module located in _[Install path]/include/petar_. To import the code, include `[Install path]/include` in the Python include path (the environment variable `$PYTHONPATH`) with the following command:
```shell
export PYTHONPATH=$PYTHONPATH:[Install path]/include
```

These environment variable configuration commands need to be executed each time a new terminal is opened. To automate the loading of these variables, it is advisable to add these commands to the .bashrc file in the case of a bash Linux system.

# Usage:

PeTar is capable of conducting $N$-body simulations and offers a range of tools for tasks such as initializing input files, post-processing to calculate Lagrangian radii, core radius determination, binary, triple, and quadruple detection, movie creation, and a Python package for data analysis. This section provides a detailed overview of these features.

## Preparing the initial condition

To initiate an $N$-body simulation, users must provide the initial conditions of the particle system. Various tools are available for generating initial conditions, such as [MCluster](https://github.com/lwang-astro/mcluster) and AMUSE.

The initial conditions are stored in a data file with a format that consists of 7 columns representing the masses, positions, and velocities of each particle per line.

Subsequently, users can utilize the `petar.init` tool (refer to [Initial Input Data File](#initial-input-data-file)) to convert this data file into a snapshot file following the `petar` style.

## Starting a Simulation

The process of starting an $N$-body simulation using the `petar` command depends on whether MPI and OpenMP support are compiled.

In the straightforward scenario, the standard procedure for commencing an $N$-body simulation is as follows:
```shell
petar [options] [snapshot filename]
```
Here, `[snapshot filename]` represents the filename of a snapshot of the particle system at a specific time, and `[options]` are utilized to regulate the behaviors of the simulations (refer to [Options](#options)).

For a new simulation, the snapshot file stores the initial conditions of a particle system. The file can be generated from the `petar.init` tool. For a restarted simulation, this file is the outputted snapshot from a previous simulation.

## Using OpenMP

When OpenMP is employed, to prevent segmentation faults in simulations with a large number of particles, users need to set the environment variable `OMP_STACKSIZE` to a sufficiently large value. For instance:
```shell
export OMP_STACKSIZE=128M
```

Furthermore, ensure that the maximum stack size is unlimited by executing the command `ulimit -s`. It should return 'unlimited'. If not, run `ulimit -s unlimited` before using `petar`.

The default number of threads is the maximum number supported by the host machine. Users can specify the number of threads (`N_threads`) by setting the environment variable `OMP_NUM_THREADS`. To enhance performance, `N_threads` should generally be kept at or below 8.

These environment variables can be utilized when executing `petar`, as shown in the following example:
```shell
OMP_STACKSIZE=128M OMP_NUM_THREADS=8 petar [options] [snapshot filename]
```

A convenient approach is to set `OMP_NUM_THREADS`, `OMP_STACKSIZE`, and `ulimit -s` in the initial script file of the terminal, such as the .bashrc file for a Bash system:
```shell
export OMP_STACKSIZE=128M
export OMP_NUM_THREADS=N_threads
ulimit -s unlimited
```
Subsequently, when a terminal is opened or `source ~/.bashrc` is executed in an existing terminal, these lines will be automatically executed. In this scenario, there is no need to specify `OMP_STACKSIZE=128M OMP_NUM_THREADS=N_threads` when using the `petar` command.

When utilizing OpenMP, it is important to note that the simulation may not be reproducible, as dynamic parallelism is employed to integrate the short-range interactions.

## Using MPI

When MPI is utilized, an MPI launcher is required to utilize multiple MPI processors. The standard approach is as follows:
```shell
mpiexec -n [N_mpi] petar [options] [snapshot filename]
```
Here, `[N_mpi]` denotes the number of MPI processors.

For optimal performance, OpenMP and MPI can be used together. In the example below, 4 MPI processors and 8 threads per processor are employed:
```shell
OMP_STACKSIZE=128M OMP_NUM_THREADS=8 mpiexec -n 4 petar [options] [snapshot filename]
```
It is important to ensure that `N_mpi x N_threads` is less than the total available CPU threads in the computing facility.

Please note that on a supercomputer, the MPI launcher may not be named `mpiexec`, and the method for setting the number of OpenMP threads may vary. Refer to the documentation of the job system or consult with the administrator to determine the appropriate approach for using MPI and OpenMP in that specific environment.

## Using GPU

When GPU support is enabled, each MPI processor will initiate one GPU job. Modern NVIDIA GPUs can handle multiple jobs simultaneously.
Therefore, it is acceptable for `N_mpi` to exceed the number of GPUs available. However, if `N_mpi` is too large, the GPU memory may become insufficient, leading to a Cuda Memory allocation error. In such cases, utilizing more OpenMP threads and fewer MPI processors is a preferable approach.

In scenarios where multiple GPUs are present, each MPI processor will utilize a different GPU based on the processor and GPU IDs. If users wish to exclusively utilize a specific GPU, they can employ the environment variable `CUDA_VISIBLE_DEVICES=[GPU index]`. For instance:
```shell
CUDA_VISIBLE_DEVICES=1 petar [options] [particle data filename]
```
This command will utilize the second GPU in the system (indexing starts from 0). The `CUDA_VISIBLE_DEVICES` environment variable can also be configured in the initial script file of the terminal.

## Resuming a Simulation

Any snapshot of particle data generated during a simulation can be utilized to resume the simulation at a specific time. To resume the simulation with the same parameter configuration as before, use the following command:
```shell
petar -p input.par [options] [snapshot filename]
```
Here, _input.par_ stores the previous parameter choices used in a simulation, automatically generated from the prior simulation. 

It is possible to modify the options for resumed simulations in two ways. Users can either directly modify _input.par_ to adjust parameters before resuming or specify new parameters in `[options]` within the `petar` command mentioned earlier. It is crucial to place `[options]` after `-p input.par` to prevent them from being overwritten by the parameters stored in _input.par_. For instance, to update the end time of the simulation to 10 after resuming simulations from the snapshot file `data.5`, use the following command:
```shell
petar -p input.par -t 10 data.5 
```

By default, after resuming, the snapshot files with the same name will be replaced. However, for other output files, new data will be appended to the existing ones (e.g., filenames with suffixes like esc, group, etc.).

In cases where a previous simulation did not reach completion and terminated abnormally, resuming the simulation from the last snapshot file may result in some events being repeated in the output data files. For instance, the same stellar evolution event might be recorded twice in a file with the '.sse' suffix. Similar behavior may occur when restarting the simulation from a non-final snapshot.

To prevent duplicate events, users can utilize `petar.data.clear` to remove events with recorded times greater than a specified time criterion before resuming. However, caution should be exercised to avoid deleting useful data inadvertently.

To overwrite the output files instead of appending to them, the option `-a 0` can be included in the `petar` command.

## Options

Users have the ability to specify various parameters in the options of the `petar` command to control aspects such as the number of binaries, time steps, energy error criteria, and more. There are two types of options available: single-character options starting with '-' and long options starting with '--'. It is advisable for users to initially review all single-character options listed by running `petar -h`. Here are some useful options:

-  `-u`: Sets the unit system for input data. Using `-u 1` requires the initial snapshot data to be in units of Msun, pc, and pc/Myr. Otherwise, the gravitational constant (G) is assumed to be 1 (Henon unit). It's essential to ensure consistency in unit settings when utilizing the `petar.init` tool. `petar` does not scale units, only modifying the gravitational constant. The `-u 1` option sets the value of G in the unit set of [Msun, pc, myr], and adjusts unit scaling factors when using stellar evolution and galactic tidal field.

-  `-t`: Specifies the finishing time. Note that the tree time steps used in `petar` are integer powers of 0.5. If the finishing time does not be a multiple of the tree time step (as specified by the `-s` option), the simulation may not end precisely at the specified time.

-  `-o`: Sets the time interval for outputting data (snapshots and status). This interval should also be a multiple of the tree time step.

-  `-s`: Determines the tree time step, which should be an integer power of 0.5. When set, the changeover radius is automatically calculated unless `-r` is manually specified.

-  `-r`: Defines the outer boundary of changeover radii. When set, the tree time step is automatically determined.

-  `-a`: Controls data appending. Using `-a 1` (default) appends new data to existing files after a restart, while other values rewrite the files. Exercise caution when restarting simulations with this option.

-  `-b`: Specifies the initial number of binaries. This is crucial for accurate velocity dispersion calculation, which in turn affects the automatic determination of tree time step and changeover radii. Primordial binaries should be listed first in the input data file (two neighbor lines per pair).

-  `-w`: Sets the output style. With `-w 2`, all particle data is printed in a single line along with system status information per output time, which can be beneficial for data analysis with small N.

-  `-i`: Determines the format of snapshot data, allowing for BINARY or ASCII format selection.

-  `-G`: Specifies the gravitational constant.

It is important to note that `-s` and `-r` significantly impact simulation accuracy and performance. Users should exercise caution with these options. The `petar.find.dt` tool can assist in finding optimized values for star clusters. For a comprehensive understanding and improved configuration, users may need to refer to the reference paper.

When using stellar evolution packages (e.g., BSE) and external potentials (e.g., Galpy), corresponding options are also displayed in `petar -h`.

## Performance Optimization

The performance of PeTar is influenced by several crucial parameters and is also dependent on the initial conditions of particle (stellar) systems. To achieve optimal performance for a given input model, users must carefully adjust the following parameters:

### Tree Time Step

- `petar` option `-s`

The tree time step represents a fixed time interval for computing the long-range (particle-tree) force. Calculating the long-range force is computationally intensive, with a complexity of $O(N \log N)$. Therefore, a smaller time step results in more computationally expensive calculations per unit of physical time. However, it's essential not to increase the tree time step excessively, as discussed further regarding the changeover and neighbor searching radii.

### Outer Changeover Radius
- `petar` option `-r`

The changeover region denotes the overlapping shell between long-range and short-range interactions. Below the inner region, short-range interactions are computed using 4th-order Hermite integration with individual time steps and the SDAR method. Above the outer region, long-range interactions are calculated using 2nd or 4th-order LeapFrog integration with the particle-tree method. In the intermediate region, both short- and long-range interactions are considered. The inner and outer radii are mass-weighted (`m^(1/3)`) for each particle, with the default ratio set at 0.1 (via the `petar` option `--r-ratio`).

Consistency between changeover radii and tree time steps is crucial to ensure accurate simulation results (refer to the PeTar paper for detailed information). A straightforward way to grasp this concept is by examining a circular Kepler orbit of two particles with a semi-major axis within the changeover region. Inside this region, the forces between the particles are divided into short-range and long-range components. Short-range forces are recalculated at each Hermite time step, while long-range forces are computed at every tree time step. Since the Hermite time step is smaller than the tree time step, the long-range forces impart velocity adjustments to the particles after several Hermite time steps. If the tree time step is too large, with only a few steps per Kepler orbit, the time resolution of long-range forces or velocity adjustments becomes insufficient, leading to inaccurate orbit integration (resulting in non-Keplerian orbits). Therefore, once the changeover region is defined, the tree time step should be sufficiently small to ensure at least several tens of sampling points for calculating the long-range forces (velocity adjustments) along the Kepler orbit.

In the absence of explicitly specified `-s` and `-r` values, `petar` automatically determines the tree time step and changeover radii based on the assumption that the input model represents a spherically symmetric star cluster with a King or Plummer-like density profile. For more intricate input models lacking a spherical structure, these parameters may require manual determination or the utilization of the option `--nstep-dt-soft-kepler`, following the self-consistent rule (tree time step - changeover radius relation) as exemplified in the aforementioned Kepler orbit scenario.

### Neighbor Searching Radius

- `petar` option `--r-search-min`

The changeover region delineates the boundary between short- and long-range interactions. In PeTar, an additional neighbor searching radius is defined for each particle to identify neighboring candidates expected to be within the changeover region during the subsequent tree time step. The minimum searching radius is slightly larger than the outer changeover radius, with the actual searching radius also depends on the particle's velocity.
As a particle's velocity increases, it can cover a greater distance within a single tree time step, necessitating a larger searching radius. However, in cases where the velocity is excessively high, the searching radius may become too large. To address this issue, an additional mechanism is implemented to mitigate the number of neighbors: for every neighbor and the particle, a Kepler orbit assumption is made, and the peri-center distance is computed. If the peri-center lies outside a specified criterion, it is disregarded as a neighbor.

The neighbor searching radius is a crucial parameter that significantly impacts performance. For parallel computing of short-range integration, PeTar initially employs the neighbor searching radius to assemble nearby particles into individual clusters. Each cluster assigns only one CPU core for the short-range (Hermite+SDAR) integration in the following tree time step. The clustering algorithm ensures that all neighbors of each cluster member are also within the same cluster. Excessive neighbor radii can lead to the formation of a single massive cluster encompassing most system particles, resulting in one CPU core handling this cluster while others remain idle. This imbalance in workload undermines computational efficiency and can notably reduce performance. Consequently, the neighbor searching radius for each particle system should not be excessively large.

Upon determining the changeover radius (`-r`), `petar` automatically computes the neighbor searching radius. If users prefer manual determination of the neighbor searching radius, the following options can be configured:
- `--r-search-min`: establishes the minimum neighbor searching radius reference, with the final radius also mass-weighted, akin to the changeover radius.
- `--search-peri-factor`: sets the maximum peri-center criterion assuming two neighbors possess a Kepler orbit.
- `--search-vel-factor`: determines the velocity-dependent addition to the neighbor radius (base neighbor radius + coefficient * velocity * tree time step).

### Multiple Group Radius
- `petar` option `--r-bin`

The third pivotal radius influencing performance is the radius used to identify a multiple group where the SDAR method is applied. In dense stellar systems, multiple groups like binaries, triples, and quadruples frequently occur, exhibiting significantly shorter orbital periods for inner members compared to single stars orbiting within the host particle system. The SDAR method in PeTar plays a crucial role in ensuring accuracy and efficiency in integrating their orbits, with binary stellar evolution also addressed within the SDAR method. The criterion for selecting group members is the group radius (`--r-bin`), which is automatically determined based on the changeover inner radius. If this radius is excessively large, the SDAR method becomes resource-intensive due to an excessive number of selected members in a multiple group. Conversely, if the radius is too small, some binaries may not be integrated accurately, as the Hermite integrator exhibits a systematic long-term drift of energy and angular momentum for periodic motion.

### Adjusting Tree Time Step and Radii

When initiating a new simulation, the automatically determined tree time step and radii may not always be the optimal choice for users. To select the most suitable tree time step, users can utilize the `petar.find.dt` tool (refer to [Find tree time step](#find-tree-time-step)). This tool is compatible only with PeTar's autodetermined tree time step and changeover radii (refer to [Outer Changeover Radius](#outer-changeover-radius)).

In cases where the structure of the particle system undergoes significant evolution over an extended period, users may wish to adjust the tree time step and radii mentioned earlier to enhance performance. If users prefer to modify only the tree time step while allowing `petar` to determine the radii automatically, the options in the following example are necessary to restart the simulation:
```shell
petar -p input.par -s [new tree_time_step] -r 0 --r-search-min 0 --r-bin 0 [other options] [snapshot filename for restart]
```
Here, `-r 0 --r-search-min 0 --r-bin 0` are employed to reset all three radii and activate autodetermination based on the new tree time step. Users can also employ `petar.find.dt` to select the optimal restart tree time step (refer to [Find tree time step](#find-tree-time-step)).

## Output
### Printed Information

When `petar` is running, several pieces of information are displayed at the beginning:
1. The FDPS logo and PeTar details are printed as follows, showcasing copyright information, versions, and references for citation. 
    ```
     //==================================\\
     ||                                  ||
     || ::::::: ::::::. ::::::. .::::::. ||
     || ::      ::    : ::    : ::       ||
     || ::::::  ::    : ::::::'  `:::::. ||
     || ::      ::::::' ::      `......' ||
     ||     Framework for Developing     ||
     ||        Particle Simulator        ||
     ||     Version 7.0 (2021/08)        ||
     \\==================================//
    ...
    ```

2. Enabled features (selected during configuration), such as stellar evolution packages, external packages, and GPU utilization, are listed:
    ```
    Use quadrupole moment in tree force calculation
    Use 3rd order tidal tensor method
    ...
    ```

3. Any modified input parameters are displayed when using corresponding `petar` options.
    ```
    Input data unit, 0: unknown, referring to G; 1: mass:Msun, length:pc, time:Myr, velocity:pc/Myr:   1
    Number of primordial binaries for initialization (assuming the binaries ID=1,2*n_bin):   500
    ...
    ```

4. Unit scaling for PeTar, stellar evolution packages (e.g., SSE/BSE), and external packages (e.g., Galpy) is outlined.
    ```
    ----- Unit set 1: Msun, pc, Myr -----
    gravitational_constant = 0.0044983099795944 pc^3/(Msun*Myr^2)
    ----- Unit conversion for BSE -----
     tscale = 1  Myr / Myr
     mscale = 1  Msun / Msun
     rscale = 44353565.919218  Rsun / pc
     vscale = 0.9778131076864  [km/s] / [pc/Myr]
    ```

5. A brief parameter list for tree time step and key radii influencing performance is provided.
    ```
    ----- Parameter list: -----
      Average mass                      = 0.52232350485643
      Mean inner changeover radius      = 0.0021228717795716
      Mean outer changeover radius      = 0.021228717795716
      Mean SDAR group detection radius  = 0.0016982974236573
      Minimum neighbor searching radius = 0.024787679043148
      Velocity dispersion               = 0.60739605289507
      Tree time step                    = 0.001953125
	  Output time step                  = 1
   ```

   The definitions of these parameters are as follows:
    - `Average mass`: average mass of all objects.
    - `Mean inner changeover radius`: mean inner changeover radius.
    - `Mean outer changeover radius`: mean outer changeover radius.
    - `Mean SDAR group detection radius`: the criterion for selecting SDAR group members.
    - `Minimum neighbor searching radius`: the minimum neighbor searching radius reference.
    - `Velocity dispersion`: velocity dispersion of the system.
    - `Tree time step`: tree time step.
	- `Output time step`: snapshot and status output time step.

    These parameters determine the performance of a simulation, see details in [Performance Optimization](#performance-optimization).

6. If Galpy is utilized, the Galpy potential setup information may be printed.
    ```
    Galpy parameters, time: 0 Next update time: 0
    Potential set 1 Mode: 0 GM: 0 Pos: 0 0 0 Vel: 0 0 0 Acc: 0 0 0
    Potential type indice: 15 5 9
    Potential arguments: 251.63858935563 1.8 1900 306770418.38589 3000 280 1965095308.1922 16000
    ```
7. In case the SSE/BSE-based stellar evolution package is employed, common block and global parameters are showcased.
    ```
     ----- SSE/BSE common block parameter list: -----
     value1: neta:  0.50000000000000000       bwind:   0.0000000000000000       hewind:   1.0000000000000000
    ...
    ```

8. Filenames for dumped input parameters are specified.
    ```
    -----  Dump parameter files -----
    Save input parameters to file input.par
    ...
    ```
    By default, these include:
    - `input.par`: Input parameters of `petar`, useful for restarting the simulation from a snapshot.
    - `input.par.hard`: Input parameters of the hard component (short-range interaction part; Hermite + SDAR), utilized for testing the dumped hard cluster with `_petar.hard.debug_`.
    - `input.par.[bse_name]`: Parameters for the SSE/BSE-based package, necessary for restarting the simulation and for `petar.hard.debug` if an SSE/BSE-based package is used.
    - `input.par.galpy`: Galpy parameters for simulation restart purposes.

    After the "Finish parameter initialization" line, the simulation status is updated at each output time interval (defined by the `-o` option). The status content follows a format similar to the example provided below:

9. Time, number of real particles, all particles (including artificial particles), removed particles, and escaped particles; both locally (within the first MPI process) and globally (across all MPI processes):
    - N_real(loc): Number of physical particles (stars) in the local MPI processor (rank 0).
    - N_real(glb): Number of physical particles (stars) across all MPI processors.
    - N_all(loc): Number of all particles, including physical and artificial ones, in the local MPI processor.
    - N_all(glb): Number of all particles across all MPI processors.
    - N_remove(glb): Number of removed particles (e.g., zero-mass particles due to mergers and escapists).
    - N_escape(glb): Number of escaped particles (beyond the escape criterion).

    Example output information at time 1:
    ```
    Time: 1  N_real(loc): 1378  N_real(glb): 1378  N_all(loc): 1378  N_all(glb): 1378  N_remove(glb): 0  N_escape(glb): 0
    ```

10. Energy check: Two rows are printed. The first row displays physical energy, while the second row shows slow-down energy (referenced in the _petar_ commander).
    - Error/Total: Relative error of the current step.
    - Error: Absolute error of the current step.
    - Error_cum: Cumulative absolute error.
    - Total: Total energy.
    - Kinetic: Kinetic energy.
    - Potential: Potential energy.
    - Modify: Total modified energy (e.g., due to stellar evolution).
    - Modify_group: Modified energy from SDAR groups.
    - Modify_single: Modified energy from singles.
    - Error_PP: Energy error in short-range particle-particle (PP) interaction.
    - Error_PP_cum: Cumulative energy error in the PP part.

    ```
    Energy:       Error/Total           Error       Error_cum           Total         Kinetic       Potential          Modify    Modify_group   Modify_single        Error_PP    Error_PP_cum
    Physic:      1.883442e-05       -644.6969       -644.6969   -3.422972e+07    1.846359e+07   -5.269331e+07        1841.216               0     0.008989855    -8.67599e-06    -8.67599e-06
    Slowdown:    1.883442e-05       -644.6969       -644.6969   -3.422972e+07    1.846359e+07   -5.269331e+07        1841.216               0     0.008989855    -8.67599e-06    -8.67599e-06
    ```

11. Angular momentum: Error at the current step, cumulative error, components in x, y, z directions, and value.

    ```
    Angular Momentum:  |L|err: 187484.2  |L|err_cum: 187484.2  L: -1.17383e+07   -1.20975e+07    -1.267862e+09  |L|: 1.267974e+09
    ```

12. System total mass, center position, and velocity.

    ```
    C.M.: mass: 736.5417 pos: 5127.807   -5729.159    7.270318 vel: -165.7037   -150.611    3.006304
    ```

13. Performance information:
    - Number of tree steps per output interval.
    ```
    Tree step number: 512
    ```
    - Wallclock computing time for each part of the code. Two lines are printed indicating the minimum and maximum times among all MPI processors:
        - Total: Total computing time per tree time step.
        - PP_single: Time to integrate single particles without neighbors.
        - PP_cluster: Time to integrate (short-range) particle clusters using the Hermite+SDAR method within one MPI processor.
        - PP_cross: Time to integrate short-range particle clusters using the Hermite+SDAR method across multiple MPI processors.
        - PP_intrpt*: Time used during the interruption of integration (only used in special mode).
        - Tree_NB: Time to search for neighbors using the particle-tree method.
        - Tree_Force: Time to calculate long-range force using the particle-tree method.
        - Force_corr: Time to correct long-range force due to the usage of changeover functions.
        - Kick: Time to kick the velocity of particles from the long-range forces.
        - FindCluster: Time to create particle clusters for short-range interactions.
        - CreateGroup: Time to find multiple groups within each cluster.
        - Domain_deco: Time for domain decomposition for the global particle tree.
        - Ex_Ptcl: Time to exchange particles between different domains.
        - Output: Time to output data and print information.
        - Status: Time to calculate the global status (e.g., energy) of the system.
        - Other: Time cost in other parts not included in the above components.
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
	- Number counts per tree time step:
	    - PP_single: Number of single particles.
	    - PP_cluster: Number of particles in clusters within one MPI processor.
	    - PP_cross: Number of particles in clusters crossing multiple MPI processors.
	    - PP_intrpt*: Number of particles experiencing interruption of integration.
	    - Cluster: Number of clusters within one MPI processor.
	    - Cross: Number of clusters crossing multiple MPI processors.
	    - AR_step_sum: Number of SDAR integration steps.
	    - AR_tsyn_sum: Number of SDAR integrations for time synchronization.
	    - AR_group_N: Number of multiple groups.
	    - Iso_group_N: Number of isolated multiple groups (a cluster containing only one group).
	    - H4_step_sum: Number of Hermite integration time steps.
	    - H4_no_NB: Number of Hermite integrations for particles without neighbors.
	    - Ep-Ep_sum: Number of interactions between particle i and j during a particle-tree force calculation.
	    - Ep-Sp_sum: Number of interactions between particle i and superparticle j during a particle-tree force calculation.
    ```
    **** Number per step (global):
    PP_single    PP_cluster   PP_cross     PP_intrpt*   Cluster      Cross        AR_step_sum  AR_tsyn_sum  AR_group_N   Iso_group_N  H4_step_sum  H4_no_NB     Ep-Ep_sum    Ep-Sp_sum
       1345.6       32.436            0            0       14.744            0       7.0586       4.0996            0            0       329.74            0   1.6972e+06        34134
    ```
    - Histogram of the number of members in clusters. The first line shows the number of members, and the second line shows the histogram counts.
    ```
    **** Number of members in clusters (global):
            1            2            3            4            5            6
       1345.6        13.24      0.68164      0.31055      0.40234      0.10938
    ```

    The performance information is crucial for verifying whether the simulation has been set up with appropriate tree time steps and radii parameters.

    To ensure reasonable performance, the `Tree_Force` wallclock time should primarily contribute to the total time. If numerous multiple groups, such as primordial binaries, are present, `PP_cluster` may also consume time. However, if `PP_cluster`'s time consumption is excessive (dominating most of the total time), users should consider adjusting the changeover radius, neighbor searching radius, and group radii (refer to [Performance Optimization](#performance-optimization)).

    The histogram depicting the number of members in clusters is valuable for determining whether the neighbor searching radius is too large. In a low-density system, the maximum number of members should typically be around 20, as seen in the example provided (6). In a high-density system or a system with high-velocity particles, the maximum number may be higher. Nonetheless, if it remains within a few hundred members and the PP_cluster wallclock time is not excessively large, the setup is acceptable.

### Output Files

In addition to the printed information provided by the `petar` commander, there are several output files detailed in the table below:

| File name            | Content                                                                                                                    |
| :-------------       | ------------------------------------------------------------------------------------------------------------------------   |
| data.[index]         | Snapshot files for each output time. The format mirrors that of the input data file.                                      |
|                      | Users can reference the definitions of the first line (header) and columns using `petar -h`.                               |
|                      | [index] denotes the output order, starting from 0 (initial snapshot). It does not correspond to time unless the output interval is set to 1. |
| data.[index].randseeds| The random seeds for each OpenMP thread, used for restarting purposes. |
| data.esc.[MPI rank]  | Contains information on escaped particles, with columns matching those in snapshot files, and an additional column for the escaped time at the beginning. |
| data.group.[MPI rank]| Provides details on the start and end of multiple systems (e.g., binary, triple ...) identified during SDAR integration. |
|                      | The definition of a multiple system is based on the distance criterion specified in the `petar` option `--r-bin`.           |
|                      | In cases where a multiple system spans multiple tree time steps, the start event may be recorded multiple times during each tree time step, while only one or no end event is recorded. This behavior is a result of the algorithm's design. |
| data.status          | Includes the evolution of global parameters such as energies, angular momentum, particle count, system center position, and velocity. |
|                      | When the `petar` option `-w 2` is utilized, this file includes all particle information similar to the snapshot files. Instead of segregating particles into distinct lines, all particles are consolidated into a single line. |
| data.prof.rank.[MPI rank]| Offers performance measurements for various parts of the code throughout the simulation.                                  |

When utilizing the SSE/BSE stellar evolution options (--with-interrupt during configure), additional files are generated:

| File name            | Content                                                                                                                    |
| :-------------       | ------------------------------------------------------------------------------------------------------------------------   |
| data.[sse_name].[MPI rank] | Contains records of single stellar evolution events, such as type changes and supernovae. Note that if a star evolves rapidly (less than the dynamical integration time step), internal type changes may not be captured. |
| data.[bse_name].[MPI rank] | Records binary stellar evolution events. All binary type changes are logged; however, if 'Warning: BSE event storage overflow!' appears during the simulation, it indicates that binary type changes are too frequent, resulting in some changes not being recorded for the corresponding binary.

When Galpy is employed (--with-external=galpy), under certain conditions, a galpy parameter file may be generated alongside each snapshot file:
| File name           | Content                                                                                                                     |
| :-------------      | ------------------------------------------------------------------------------------------------------------------------    |  
| data.[index].galpy  | Contains the galpy parameter file for each snapshot, which may be necessary for restart purposes.                                                                                   |

In this context, 'data' serves as the default prefix for output files, although users have the flexibility to modify it using the `petar` option `-f`. 
For instance, by specifying `-f output`, the output files will be named 'output.[index]', 'output.esc.[MPI rank]', and so forth.

The term [MPI rank] denotes the MPI processor responsible for outputting the data. Consequently, the number of data files with this suffix aligns with the count of MPI processors employed in the simulation. For instance, with two MPI processors, the escape data files will be labeled 'data.esc.0' and 'data.esc.1'.

Prior to accessing these files, it is advisable to execute the `petar.data.gether` tool to consolidate the individual files generated by different MPI ranks into a single file for ease of use.

The `petar.data.gether` tool not only consolidates files from various MPI processors but also generates new files accessible by the `petar` Python data analysis tool (refer to [Python Data Analysis Module](#python-data-analysis-module)). 
- For ".group" files, `petar.data.gether` separates few-body groups with varying member counts into individual files labeled with the suffix ".n[number of members in groups]".
- In the case of ".[sse\_name]" and ".[bse\_name]" files, this tool segregates type changes, supernova kicks, and dynamical mergers into distinct files (".type\_change", ".sn\_kick", and ".dynamic\_merge").

For a detailed overview of the files generated by `petar.data.gether` and the corresponding Python reading methods, refer to [Gathering Output Files](#gathering-output-files).

Given the extensive number of columns in each file, it is recommended to utilize the Python data analysis tool for data access. This tool simplifies the process of selecting specific parameters (columns) and conducting data operations (selection, calculation, and plotting), akin to using `dict` and `numpy` in Python. 
Moreover, it helps prevent errors in column reading. Consequently, column definitions are not provided in the manual or within the file headers. 
For column specifics, consult the documentation of the relevant analysis tool associated with each file.

The raw snapshot files do not encompass information on binaries or multiple systems. To identify and acquire details on Lagrangian and core properties, the `petar.data.process` tool can be employed (refer to [Parallel Data Processing](#parallel-data-processing)).

#### Snapshot File Format

By default, `petar` outputs snapshots in ASCII format. The advantage of ASCII format is that users can directly read the file using a standard text reader. However, it is recommended to output snapshots in BINARY format by setting the `-i` option of `petar`. For instance, 
```shell
petar -i 2 ...
```
allows reading the input snapshot in ASCII format and outputting snapshots in BINARY format. BINARY format snapshots preserve precision for floating-point numbers, are significantly faster (at least 10 times faster) to read and write by `petar` and PeTar's data analysis module, and result in much smaller file sizes compared to ASCII format, helping to conserve hard disk space.

Users can also convert between ASCII and BINARY formats after simulations using the `petar.format.transfer` tool (refer to [Data Format Conversion](#data-format-conversion)). Additionally, it is possible to convert post-processed files from `petar.data.process` using the `petar.format.transfer.post` tool.

### Units

There are 1-3 sets of units in `petar` depending on the packages utilized:

#### Petar Units

In the absence of additional stellar evolution and external potential packages, PeTar adheres to the units specified in the input files. Within the PeTar framework, there is no internal unit conversion. The only adjustable parameter is the gravitational constant, which can be modified to align with the units defined in the input data file (Refer to [Options](#options)). Therefore, it is important that the input data maintains a consistent unit set: the velocity unit must correspond to the length unit. For instance, if the length unit is pc, the velocity unit should be pc/[time unit]. It is not permissible to use km/[time unit] as this would need an additional unit conversion from kilometers to parsecs during integration. Such conversions introduce unnecessary complexities and potential bugs without offering any tangible benefits. Consequently, PeTar restricts unit modifications to solely adjusting the gravitational constant to uphold a self-consistent unit system.

Given the absence of unit conversions, the output log, snapshot data (PeTar component), escaper files, and group files adhere to the unit set specified in the input data alongside the corresponding gravitational constant. For instance, if the input data utilizes units such as (Myr, pc, pc/Myr), the time and energy values in the output log will also be expressed in the same unit system. The potential energy calculations will account for the gravitational constant within this context.

#### Stellar Evolution Units

The stellar evolution package (e.g., BSE) operates on a distinct unit system compared to PeTar. Consequently, there exists a unit conversion between the PeTar component and the stellar evolution package. Users can manually define the conversion scaling factors using `petar` options such as `--bse-rscale` (refer to `petar -h` for detailed information). It is advisable to adopt a unit set consistent with input data (e.g., in terms of solar masses, parsecs, and parsecs/Myr) to leverage the `petar` option `-u 1`, enabling automatic calculation of the conversion factor by PeTar. The unit system employed in the snapshot files for the stellar evolution segment and the output files "[data filename prefix].[s/b]se.[MPI rank]" align with the unit system specified in '[bse_name]'.

#### External Potential Units

The external potential package (Galpy) introduces yet another unit system, necessitating a separate unit conversion between the PeTar component and Galpy. Users can manually define the conversion factors for Galpy using options like `--galpy-rscale` (refer to `petar -h` for specifics). To streamline the process and minimize complexities, it is advisable to maintain a unit set consistent with the input data (e.g., in terms of solar masses, parsecs, and parsecs/Myr) for Galpy as well, thereby eliminating the need for additional unit conversions. This approach is the default choice.

It is important to note that this differs from the official unit of Galpy, which is based on the solar motion. In the official Galpy unit system, the length unit corresponds to the distance from the Sun to the Galactic center, and the velocity unit represents the solar velocity in the Galactic frame. Therefore, when interpreting suggested argument values of potentials from Galpy, users should calculate these values in PeTar units instead, ensuring consistency across the calculations.

#### Output File Units

Below is a table illustrating the corresponding units for various output files, with "data" as the example data filename prefix:

| Filename                    | Content                                                     | Unit Set                                                           |
| :-------------------------- | :---------------------------------------------------------  | :---------------------------------------------------------------  |
| `petar` output log          | Printed information from `petar`                            | PeTar unit                                                         |
| data.[index]                | Snapshots                                                   | Particle class: PeTar unit + Stellar evolution unit (refer to `petar -h`) |
| data.esc.[MPI rank]         | Escapers                                                    | Time: PeTar unit; Particle: Particle class                         |
| data.group.[MPI rank]       | Multiple systems                                            | Binary parameters: PeTar unit; Particle members: Particle class   |
| data.[bse_name].[MPI rank]  | Binary stellar evolution events                             | Stellar evolution unit                                             |
| data.[sse_name].[MPI rank]  | Single stellar evolution events                             | Stellar evolution unit                                             |
| data.status                 | Global parameters                                           | PeTar unit                                                         |
| data.prof.rank.[MPI rank]   | Performance profiling                                       | Time: Second (per tree time step)                                  |

In the context of the 'particle class', it denotes the data structure of a single particle (C++ class) within PeTar. Depending on the stellar evolution mode, each particle contains a mix of data from the PeTar component and the stellar evolution component, which may not be in the same unit system. The units for stellar evolution parameters prefixed with "s_" in the particle class can be accessed using `petar -h`. Other members adhere to the units specified in the input data file (PeTar unit).

Should users require clarification on the units of the output files, they can also refer to the help information provided by the Python analysis tool (`help(petar.Particle)`) for guidance on interpreting the respective files.

## Troubleshooting

During a simulation, users may encounter warning and error messages in the printed information of `petar`. When such instances arise, a dump file is generated simultaneously to capture the initial conditions for reproducing the detailed issue. The following section provides an overview of these warnings and errors for user reference.

### Significant Hard Energy

In scenarios involving a particle cluster with short-range interaction (utilizing Hermite+SDAR), if the relative energy error exceeds the threshold specified by the `--energy-err-hard` option of `petar` (default value: 0.0001) during one tree time step, a warning message titled "Hard energy significant" is triggered, and a corresponding dump file named "hard\_large\_energy.*" is created.

This issue commonly arises when a triple or quadruple system exists within the particle system. The large error can be attributed to two main possibilities:
- The inner binaries of these systems may exhibit tight configurations with substantial slowdown factors. As the slowdown Hamiltonian differs from the classical Hamiltonian, physical energy conservation is not guaranteed during integration, focusing instead on ensuring correct secular motion.
- The integration step of the SDAR method might not be sufficiently small for multiple systems. While reducing step sizes can enhance accuracy, it also leads to a notable increase in computational time.

Resolving these cases without compromising computational efficiency can be challenging. However, if users are concerned about specific particles within the cluster, they can assess the integration orbit's acceptability by employing the debugging tool `petar.hard.debug` to analyze the dumped "hard\_large\_energy.*" file. Typically, the high energy error stems from minor changes in the semi-major axis of the tightest binary in the system. Nonetheless, such minor discrepancies typically do not significantly impact the overall dynamical evolution of the system if the error occurs sporadically.

### Large Step Warning

At times, the code's performance may significantly deteriorate, accompanied by a warning indicating large steps. Subsequently, a file named "dump\_large\_step.*" is generated. This situation typically arises when a stable multiple system exists within the particle cluster.

The AR method is employed to integrate the multiple system, but the step count becomes excessively large, surpassing the predefined step limit (which can be adjusted in the input option), leading to a notable decline in performance. Unfortunately, this scenario is unavoidable in some instances.

Enabling stellar evolution can offer some assistance, especially when the inner binaries are sufficiently tight to merge. However, there is no definitive solution to address this issue. If multiple CPU cores are utilized, users have the option to restart the simulations with fewer CPU cores. Sometimes, upon restarting, the same stable system may not form, thereby circumventing the problem.

Should the stable system persist even after restarting, terminating parallel computing can be considered to eliminate the need for multiple CPUs. Users can temporarily reduce CPU resources to navigate through this phase until the system is disrupted. Subsequently, they can restart with the original number of CPU cores.

### Hard Dump with Errors

When errors manifest during the Hermite-SDAR integration, an error message is displayed, and the simulation is halted, triggering the creation of a file named "hard_dump.*". This occurrence typically signifies the presence of a bug within the code.

Users encountering this issue are encouraged to report it by contacting the developer either through GitHub or email. In the report, users should provide essential details such as the version of PeTar, the configuration options, the initial simulation conditions, and include the "hard_dump.\*" file and the input parameter files (prefixed with "input.par.\*").

For those inclined to investigate the issue independently, the debug tool `petar.hard.debug` in conjunction with the GDB tool can be utilized. However, a comprehension of the source codes of SDAR is necessary to interpret the messages generated by the debug tool effectively.

### Hard Debug Tool

The aforementioned warnings and errors generate dump files that can be analyzed using the `petar.hard.debug` tool. This tool facilitates the re-execution of the simulation for one tree time step specifically for the isolated hard sub-cluster associated with the warning, aiding in pinpointing the source of the warning or error.

The basic usage of the tool is as follows:
```shell
petar.hard.debug [dump_file_name] > debug.log 
```
Here, `[dump_file_name]` refers to the name of the dump files discussed in earlier sections, such as "hard\_large\_energy.\*" and "dump\_large\_step.\*". By executing the command above, the `petar.hard.debug` tool will display snapshots of particle data per line in the primary output file (debug.log) along with additional information in the printed messages.

To interpret the debug.log file, users can utilize the Python analysis tool `petar.HardData`. Below is a sample script that reads the debug.log file, converts the first two particles into a binary system, and plots the evolution of the semi-major axis:
```python
# Read petar.hard.debug log, where N_particle represents the total particle number in the sub-cluster (n_ptcl) and N_sd denotes the group number (n_group). Obtain these values from the petar.hard.debug output message.
# Ensure correct option arguments for interrupt_mode and external_mode to accurately interpret the debug.log file.
hard = petar.HardData(member_type=petar.Particle, particle_type='hard', interrupt_mode='bse', external_mode='galpy', N_particle=4, N_sd=2)
hard.loadtxt(path+'debug.log', skiprows=1)

# Retrieve column indices and names
hard.getColumnInfo()

# Convert the first two particles into a binary system
b = petar.Binary(hard.particles.p0, hard.particles.p1, G=petar.G_MSUN_PC_MYR)

# Plot time versus semi-major axis
import matplotlib.pyplot as plt
%matplotlib inline
fig, axes = plt.subplots()
axes.plot(hard.time, b.semi)
```

The initial message in the printed output displays the simulation's input parameters, followed by the particle count (n\_ptcl) and group count (n\_group). Subsequent messages provide insights into the formation, exchange, and disruption of groups (e.g., binaries, triples). If stellar evolution is active, interruption events are also indicated. These details are instrumental in comprehending the factors contributing to energy errors, the presence of large steps, and pinpointing error locations.

There is one drawback of the current version of `petar.hard.debug`. When all particles are within a single group, the integration is performed using the pure SDAR method. In this scenario, the debug.log file may remain empty as snapshot files are only generated when the Hermite integrator is utilized.

`petar.hard.debug` proves more beneficial when used in conjunction with compiler debugging tools like the `gdb` tool. With `gdb`, it becomes feasible to track the simulation's progress step by step within the source code to investigate instances of significant energy errors, a high number of steps, or unexpected crashes. To employ `gdb`, the basic procedure is as follows:
```shell
gdb petar.hard.debug
```
Subsequently, in gdb mode, execute:
```
run [dump_file_name] > debug.log
```
This command will run `petar.hard.debug` with the specified dump file name and store the output in debug.log. Users can establish breakpoints in the source code to halt the execution at precise locations and inspect variable values in the vicinity of the source code. This capability proves invaluable for understanding the exact behavior of the `petar` code and diagnosing issues as they arise. For a more comprehensive understanding, users should familiarize themselves with the usage of `gdb` beforehand.

### Crash with Assertion

At times, the code may crash accompanied by an assertion error message. One common assertion error encountered is `n_jp<=pg.NJMAX`, as illustrated below:
```
CalcForceEpEpWithLinearCutoffSimd::operator()(const EPISoft*, ParticleSimulator::S32, const EPJSoft*, ParticleSimulator::S32, ForceSoft*): Assertion `n_jp<=pg.NJMAX' failed.
```
This error typically occurs when the particle system exhibits an extremely high density contrast, or when the density center is distant from the coordinate origin (zero point), or if a particle is significantly far from the system.

The outermost box size during particle-tree construction is determined by the farthest particles in the system. If the maximum box size is excessively large, it results in a large minimum tree cell size, leading to a situation where too many particles may reside in a single cell, surpassing the limits of the particle-tree algorithm and triggering the assertion `n_jp<=pg.NJMAX`. Additionally, the tree cell closest to the coordinate origin point boasts the highest resolution, hence a distant density center from the origin can also trigger this assertion. Therefore, the resolution lies in removing distant particles or appropriately selecting the coordinate origin point.

Another frequent assertion error is `!std::isnan(vbk.x)`, exemplified by:
```
SystemHard::driveForOneClusterOMP(ParticleSimulator::F64): Assertion `!std::isnan(vbk.x)' failed.
```
This error is observed when FDPS version 7.1 is utilized. It appears that a bug or an inconsistent interface in FDPS can lead to such assertions within `petar`. The recommended solution is to revert to using FDPS version 7.0 to mitigate this issue.

## Data Format Update for Older Versions

Over time, the data formats of snapshots, input parameter files, and certain output files have undergone revisions. Users seeking to utilize a newer version of the code to interpret data from older versions can facilitate data transfer.

For snapshot data in ASCII format, a notable change occurred after [the version released on Aug 8, 2020](https://github.com/lwang-astro/PeTar/commit/0592d70875626071e1bd7aa13dbab30165a98309#diff-25a6634263c1b1f6fc4697a04e2b9904ea4b042a89af59dc93ec1f5d44848a26). Specifically, the output format of `group_data.artificial` transitioned from 64-bit floating point to 64-bit integer to preserve complete information. This alteration solely affects the ASCII format, while the BINARY format remains unchanged. To read older snapshot data, a data conversion step is necessary:
```shell
petar.format.transfer -g [other options] [snapshot path list filename]
```
This process generates new data in BINARY format. By employing the same tool with the `-b` option, users can convert the BINARY format back to the updated ASCII format.

The formats of input parameter files produced during simulations (including files from SSE/BSE and Galpy) were updated on Oct 18, 2020. To adapt the input files for use with newer versions of PeTar, users can employ [`petar.update.par`](#input-parameter-file-format-update). Post-update, the reading and modification of input parameter files are significantly improved, enhancing the ease of restarting simulations with the latest PeTar versions.

## Useful Tools

Several handy tools are available to aid users in generating initial input data, determining an appropriate tree time step to commence simulations, and conducting data analysis. These tools are bundled with `petar` and follow a naming convention of `petar.[tool name]`. To access guidance on utilizing each tool, users can employ the following command:
```shell
petar.[tool name] -h
```
It is important to note that options with identical names may hold distinct meanings across various tools.

The subsequent sections provide detailed descriptions of each tool.

### Initial Input Data File

PeTar features an internal Plummer model generator for an equal-mass system, utilizing the Henon Unit with a half-mass radius of 1.0. Should users prefer to employ their own initial particle data, the `petar.init` tool facilitates the conversion of their particle data into a `petar` input data file. The usage is as follows:
```shell
petar.init [options] [particle data filename]
```
The particle data file should consist of 7 columns: mass, position (3 coordinates), velocity (3 components), with each particle represented in a separate row. Binaries should be listed first, with the two components adjacent to each other. When binaries are present, the `-b [binary number]` option must be included in the `petar` command to ensure accurate initialization of velocity dispersion, tree time step, and changeover radii.

If stellar evolution is activated, the corresponding options `-s [bse_name]` should be used concurrently to generate the correct initial files. In such cases, it is recommended to utilize astronomical units (Solar mass [Msun], parsec [pc], and Million years [Myr]) for the initial data. The velocity unit should be specified as pc/Myr, ensuring a mass scaling factor of 1.0 between PeTar units and SSE/BSE-based code. Additionally, the `-u 1` option should be added to the `petar` command to adopt this astronomical unit set.

Similarly, when external mode (potential) is enabled, the `-t` option should be utilized to ensure the correct number of columns is generated.

### Determining the Tree Time Step

The performance of `petar` is highly dependent on the tree time step chosen. To assist in finding the optimal time step for achieving the best performance, `petar.find.dt` can be utilized. The usage is as follows:
```shell
petar.find.dt [options] [petar snapshot filename]
```
The performance of `petar` relies on the initial particle data file in the petar input format. This tool conducts brief simulations with various tree time steps and presents the performance results sequentially. Users can then determine which time step yields the best performance.

It is important to note that if a time step that is too large is tested, the tool may not respond for an extended period, indicating a suboptimal choice of time step. In such cases, the tool will terminate the test and provide the best result from the previous trials.

Several options are available in this tool to adjust the numbers of OpenMP threads and MPI processors, as well as the minimum tree time step to initiate the test. If other options are used in the `petar` command, such as `-b [binary number]`, `-u [unit set]`, or `-G [gravitational constant]`, these options should be included using the `-a` option with the content enclosed in double quotes:
```shell
petar.find.dt [options] -a "[petar options]" [petar data filename]
```
For example:
```shell
petar.find.dt -m 2 -o 4 -a "-b 100 -u 1" input
```
This command uses 2 MPI processes, 4 OpenMP threads per MPI process, 100 primordial binaries, and a unit set of 1 [Msun, pc, Myr] to determine the best time step.

It is worth noting that `petar` only accepts tree time steps that are integer powers of 0.5. Therefore, during testing, if the user specifies the minimum step size using `-s [value]` (outside `-a`), the step size will be adjusted to meet this requirement if necessary. Users should be cautious as some `petar` options, such as `-o` and `-s`, cannot be used within the `-a` option of `petar.init`. Further details can be obtained by using `petar.find.dt -h`.

For users looking to restart a simulation and automatically determine the new tree time step along with other parameters (radii), the following command can be used:
```shell
petar.find.dt -m 2 -o 4 -a "-p input.par -r 0 --r-search-min 0 --r-bin 0" [restart snapshot filename]
```

### Gathering Output Files

In MPI usage, each MPI processor generates individual data files with filenames containing the suffix `[MPI rank]`. To consolidate these output files from different MPI ranks into a single file, the `petar.data.gether` tool is utilized. Additionally, this tool can split the stellar evolution event files and group files into individual components, enabling the use of Python tools for data analysis (refer to [Output files](#output-files)).

Moreover, `petar.data.gether` generates a file named `"[output prefix].snap.lst"` that includes a sorted list of all snapshot files based on their respective timestamps. This file serves as input for both `petar.data.process` and `petar.movie`.

The basic usage of `petar.data.gether` is as follows:
```shell
petar.data.gether [options] [data filename prefix]
```
Here, `[data filename prefix]` represents the prefix of data files specified by the `petar` option `-f` (default is 'data').

Several options can control the gathered data; use `petar.data.gether -h` to review the details of `[options]`.

- `-f`: specifies the filename prefix of the gathered data, defaulting to the same as `[data filename prefix]`.
- `-i`: prompts whether to overwrite existing gathered data from previous runs of `petar.data.gether`; if not provided, old files are replaced.
- `-l`: generates only the snapshot file list.
- `-g`: gathers group files and splits them into separate files based on the number of members in a group.

Below is a table detailing the files generated by the tool and the corresponding Python analysis classes for reading (refer to [Data analysis in Python3](#data-analysis-in-python3)):

| Original files                | Output files                  | Content                                     | Python classes initialization for reading |
| :--------------               | :------------                 | :---------                                  | :-------------------------------------    |
| data.group.[MPI rank]         | data.group                    | all groups                                  | None                                      |
|                               | data.group.n2                 | binaries, hyperbolic systems                 | petar.GroupInfo(N=2)                      |
|                               | data.group.n3                 | triples                                     | petar.GroupInfo(N=3)                      |
|                               | ...                           | multiple systems...                         | ...                                       |
| data.esc.[MPI_rank]           | data.esc                      | escapers                                    | petar.SingleEscaper(interrupt_mode=[\*], external_mode=[\*]) |
| data.[sse_name].[MPI rank]    | data.[sse_name]               | full single stellar evolution records       | None                                      |
|                               | data.[sse_name].type_change   | single stellar evolution type change events | petar.SSETypeChange()                     |
|                               | data.[sse_name].sn_kick       | supernova natal kick of single star         | petar.SSESNKick()                         |
| data.[bse_name].[MPI rank]    | data.[bse_name]               | full binary stellar evolution records       | None                                      |
|                               | data.[bse_name].type_change   | binary stellar evolution type change events | petar.BSETypeChange()                     |
|                               | data.[bse_name].sn_kick       | supernova natal kick in binaries            | petar.BSESNKick()                         |
|                               | data.[bse_name].dynamic_merge | dynamical-driven (hyperbolic) mergers       | petar.BSEDynamicMerge(less_output=[\*])   |

Here, [\*] represents arguments that depend on the configuration options used for compilation.

Note: 'data.group[.n\*]' files are not generated by default due to their large size. To obtain these files, the `-g` option must be added.

When SSE/BSE is employed and the code version is before Sep 10, 2020, data with the suffix ".dynamic\_merge" contains three fewer columns compared to the new version. This tool automatically fills the missing columns with zeros from Column 6 to 8.

### Parallel Data Processing

The `petar.data.process` tool is utilized for analyzing snapshot data, identifying binaries, triples, and quadruples (specifically binary-binary types), and computing Lagrangian radii, core radii, averaged mass, and velocity dispersion. It is important to highlight that the tool calculates the core (density) center and uses it to determine Lagrangian radii. Single and binary data are saved for each snapshot with the additional suffix ".single" and ".binary", respectively. Triples and quadruples are optional and not used in the computation of Lagrangian properties. Data for Lagrangian, core, and escapers is generated in separate files. Multiple CPU cores are utilized for data processing due to the slow nature of the KDTree neighbor search required for density calculation and binary detection.

The basic usage is as follows:
```shell
petar.data.process [options] [snapshot path list filename]
```
Users need to provide `[snapshot path list filename]` for a file containing a list of paths for the snapshot data files. This file can be generated using the `petar.data.gether` tool. Alternatively, users can manually create the file using commands like `ls` and `egrep` in Linux. For instance, if the snapshot filename prefix is 'data', the command would be:
```shell
ls | egrep '^data.[0-9]+$'
```
It is recommended to sort the paths in increasing order of evolution time. The `sort` tool can be used for this purpose, as shown in the example below:
```shell
ls | egrep '^data.[0-9]+$' | sort -n -k 1.6 > snap.lst
```
This command finds all data files in the current directory, sorts them based on the suffix (values after 'data.') in increasing order, and saves the list to the file 'snap.lst'. The `-n` flag specifies that the values to sort are floating-point numbers, and `-k` defines the starting position of the number for sorting.

Users should ensure to set the correct options for gravitational constant (`-G`), interrupt mode (`-i`), and external mode (`-t`) for `petar.data.process`. This is crucial for reading snapshots and calculating Kepler orbital parameters of binaries. For simulations using astronomical units (`-u 1` in the `petar` command), the gravitational constant `-G 0.00449830997959438` should be used for `petar.data.process`. If using a package like SSE/BSE, the interrupt mode option `-i [bse_name]` can set the correct value of G. When an external mode like Galpy is employed, the external mode option `-t galpy` is necessary.

Below is a table showing the files generated by `petar.data.process` and the corresponding Python modules (classes) for reading them. The default filename prefix 'data' is assumed.

| Filename            | Content                                                              | Python analysis classes for reading data files  |
| :------------------ | :------------------------------------------------------------------  | :------------------------|
| data.[\*].single     | snapshots for single stars                                           | petar.Particle(interrupt_mode=[\*], external_mode=[\*])| 
| data.[\*].binary     | snapshots for binaries                                               | petar.Binary(member_particle_type=petar.Particle, simple_mode=[\*], G=[\*], interrupt_mode=[\*], external_mode=[\*]) |
| \*data.[\*].triple    | snapshots for triples (if option `-M` is used)                       | petar.Binary(member_particle_type_one=petar.Particle, member_particle_type_two=[petar.Particle, petar.Particle],  simple_mode=[\*], G=[\*], interrupt_mode=[\*], external_mode=[\*])            |
| \*data.[\*].quadruple | snapshots for binary-binary quadruples (if option `-M` is used)      | petar.Binary(member_particle_type=[petar.Particle, petar.Particle],  simple_mode=[\*], G=[\*], interrupt_mode=[\*], external_mode=[\*])             |
| data.lagr           | Lagrangian and core properties for all objects, singles, binaries, and user-defined stellar types (`add_star_type`)  | petar.LagrangianMultiple(mass_fraction=[\*], calc_energy=[\*], external_mode=[\*], add_star_type=[\*]) |
| data.core           | core position, velocity, and radius                                   | petar.Core()               |
| data.esc_single     | Single escapers                                                      | petar.SingleEscaper(interrupt_mode=[\*], external_mode=[\*])   |
| data.esc_binary     | Binary escapers                                                      | petar.BinaryEscaper(member_particle_type=petar.Particle, simple_mode=[\*], G=[\*], interrupt_mode=[\*], external_mode=[\*])      |
| \*data.bse_status   |  Evolution of number counts, maximum and averaged masses of different stellar types | petar.BSEStatus() | 
| \*data.tidal         | Tidal radius data (if option `--r-escape tidal` is used)             | petar.Tidal()              |

Note: The arguments [\*] for the keywords in the Python class initialization depend on the configure options for compiling (e.g., `interrupt_mode`, `external_mode`) and the options used in `petar.data.process`. These keyword arguments are optional and only needed when default values are not used. Refer to the help of the Python analysis classes and `petar.data.process -h` for more details.

For details on how to use the Python analysis classes to read data, refer to [Data analysis in Python3](#data-analysis-in-python3).

The snapshots generated by `petar.data.process` are shifted to the rest frame where the density center is the coordinate origin. By adding the core position and velocity from 'data.core' at the corresponding time, positions and velocities in the initial frame or Galactocentric frame (when Galpy is used) can be recovered.

When snapshot files are in BINARY format, the option `-s binary` can be used for `petar.data.process` to read the snapshots correctly. It's important to note that the data generated by `petar.data.process` are all in ASCII format.

Additionally, the `petar` code can remove escapers and store the data of escapers using energy and distance criteria (in files [data filename prefix].esc.[MPI rank], see [Output files](#output-files)). All escapers during the simulations are not stored in the original snapshot files. The `petar` code only applies a simple constant escape radial criterion. The post-processing by `petar.data.process` can calculate the tidal radius and detect escapers, which are then stored in the post-generated escape files: "data.esc\_single" and "data.esc\_binary". These escapers are not removed from the post-generated snapshot files: data.[index].single and data.[index].binary.

For Lagrangian properties, 'data.lagr' includes radius, average mass, number of objects, different components of velocity, and dispersions within different Lagrangian radii. The mass fractions of Lagrangian radii are 0.1, 0.3, 0.5, 0.7, and 0.9 by default. The core radius property is added at the end. There is an option in `petar.data.process` to define an arbitrary set of mass functions. When using the SSE/BSE-based stellar evolution package, an additional option `--add-star-type` can be used to calculate Lagrangian properties for specific types of stars. When `--add-star-type` is used, the reading function should have consistent keyword arguments. An example of reading 'data.lagr' is provided in [Reading Lagrangian data](#reading-lagrangian-data).

When `--calc-energy` is used, potential energy, external potential energy, and virial ratio for each Lagrangian radii are calculated. However, when an external potential is used, the virial ratio may not be accurately estimated in disrupted phases.

### Gathering Specified Objects

The `petar.get.object.snap` tool enables the collection of specified objects from a list of snapshots into a single file with a time series. Users can define IDs, stellar types, mass ranges, and a custom Python script to select objects.

For instance, if users wish to extract the trajectories of objects with IDs 1 and 2, this tool can scan all provided snapshots and consolidate the data of these objects into a single file. The following script demonstrates this process:
```shell
petar.get.object.snap -m id 1_2 [snapshot path list filename]
```
Subsequently, a new file named "object.1_2" is created for objects with IDs 1 and 2. By utilizing `petar.Particle` from the data analysis module, users can access this file and analyze the evolution of these two objects throughout the simulation. For more information, refer to 
 the help information of `petar.get.object.snap -h`.

Similar to configuring `petar.data.process`, users should ensure the correct settings for interrupt mode (`-i`), external mode (`-t`), and snapshot file format (`-s`) are in place while using `petar.get.object.snap`. This guarantees the accurate interpretation of the snapshots.

### Movie Generator

The `petar.movie` tool is a convenient utility for creating movies from snapshot files. It can generate movies showcasing the positions (x, y) of stars, the HR diagram if stellar evolution (SSE/BSE) is enabled, and the 2D distribution of semi-major axis and eccentricity of binaries. To generate a movie, a list of snapshot files is required.

The basic usage of `petar.movie` is as follows:
```shell
petar.movie [options] [snapshot path list filename]
```

For plotting binary information, it is recommended to first utilize `petar.data.process` to detect binaries using multiple CPU cores. This preprocessing step eliminates the need for the movie generator to employ the expensive KDTree function for binary detection (use the option '--generate-binary 2').

This tool utilizes either the `imageio` or `matplotlib.animation` Python modules to create movies. Installing `imageio` is advised for faster movie generation using multiple CPU cores, as `matplotlib.animation` can only utilize a single CPU core. Additionally, it is recommended to install the `ffmpeg` library to support various commonly used movie formats such as mp4 and avi. Please note that `ffmpeg` is a standalone library and not a Python module. Users should install it in the operating system (e.g., via the `apt` tool in Ubuntu).

Similar to configuring `petar.data.process`, users should ensure to set the correct options for the gravitational constant (`-G`), interrupt mode (`-i`), external mode (`-t`) and snapshot file format (`-s`) when using `petar.movie`. This ensures to correctly reading the snapshots.

### Data Removal after a specified time

The `petar.data.clear` tool serves the purpose of removing data recorded after a specified time from all output files except the snapshots. This tool is particularly useful for preventing the repetition of events when restarting an abnormally interrupted simulation.

In scenarios where a simulation does not terminate normally, the last output snapshot lags behind other output files that record events, such as '\*.groups' and '\*.bse' files. Consequently, restarting the simulation from this snapshot may lead to the re-recording of the same event in these files.

The basic syntax for utilizing this tool is as follows:
```shell
petar.data.clear -t [time] [data filename prefix]
```
Before data removal takes place, all output files undergo a backup process by renaming them with an additional '.bk' suffix. This ensures that, in the event of accidental tool usage, the original files can be restored from the backups.

If the tool is executed again, it automatically checks for the presence of backup data files and utilizes them to create new files with removed data. Essentially, this process is akin to restoring the original file before executing the data removal.

Please be aware that the backup files only retain the data prior to the most recent modification. If the tool is utilized multiple times and the backup files are overwritten, it is important to note that the original data may not be completely recoverable.

### Data Format Conversion

`Petar` offers the capability to read and write snapshot files in either BINARY or ASCII format. The BINARY format, being compressed (resulting in file sizes less than half of the ASCII format), facilitates faster read and write operations for both `Petar` and data analysis tools. However, it is important to note that BINARY files cannot be directly interpreted using a text editor. It is recommended to opt for the BINARY format when simulations yield a substantial amount of data and users intend to utilize analysis tools for data interpretation. The `Petar.Particle` and `Petar.PeTarDataHeader` classes in the Python3 analysis module support reading both BINARY and ASCII formats.

Moreover, it is feasible to convert snapshot data between BINARY and ASCII formats using the `petar.format.transfer` and `petar.format.transfer.post` tools. 

`petar.format.transfer` is used for snapshot files outputted from `petar`. The basic syntax for converting a list of snapshot files is as follows:
```shell
petar.format.transfer [options] [snapshot path list filename]
```
The snapshot path list contains paths to the snapshots that users wish to convert. By default, new files are generated with a '.B' or '.A' suffix. To overwrite files for space conservation, users can utilize the `-r` option. This tool can also update older versions of snapshots in ASCII format (prior to Aug 8, 2020) to newer versions using the `-g` option. It is important to note that in the older version, certain information stored in the group_data (group center-of-mass mass and velocities in 64-bit floating point) is lost. However, this loss is inconsequential as the data can be recalculated during data processing and does not impact restart operations.

It is essential to ensure that the versions of `petar.format.transfer` and `Petar` align in terms of interrupt mode and external mode configurations. In cases where snapshots are generated by different versions of `Petar`, `petar.format.transfer` may fail to read data or provide incorrect transferred data.

`petar.format.transfer.post` is designed for snapshots generated by `petar.data.process`. The fundamental syntax is as follows:
```shell
petar.format.transfer.post [options] [snapshot path list filename]
```
The snapshot path list contains paths to the snapshots that users intend to convert. These snapshots include single, binary, triple, and quadruple snapshots generated from `petar.data.process`. Users have the flexibility to convert among three formats: ASCII, BINARY, and npy. The npy format corresponds to the Python Numpy data format. It is necessary for users to maintain consistent interrupt mode and external mode settings to ensure successful format conversion.

### Update of Input Parameter File Format

The formats of input parameter files generated during simulations, including files from SSE/BSE and Galpy, were revised on Oct 18, 2020. To utilize the new version of PeTar for restarting simulations using old versions of input parameter files, an update is required:
```shell
petar.update.par [options] [input parameter filename]
```
By employing options such as `-p`, `-b`, and `-t`, users can update input parameters from PeTar, SSE/BSE, and Galpy, respectively. Additional options cater to different features selected in the configuration.

Post-update, the new input parameter files are more user-friendly. They consist of three columns defined as (1) the data type of the argument, (2) option names, and (3) argument values. The reference for the first two columns can be accessed using the command `petar -h`. Users can directly modify the argument values in the file. Furthermore, it is not necessary to list all options in the file. Consequently, if new options are introduced in future versions, there is no necessity to update the file again unless existing option names undergo changes.

### Stellar Evolution Tool based on SSE and BSE

The `petar.[bse_name]` tool is generated when utilizing the SSE/BSE based stellar evolution package (--with-interrupt=[bse_name]) during configuration, where `[bse_name]` can be 'bse', 'bseEmp' and 'mobse'. This tool functions as a standalone application for evolving stars and binaries up to a specified time. All essential global parameters can be configured through options.

To evolve a group of stars using the tool:
```shell
petar.[bse_name] [options] mass1, mass2 ...
```
If individual masses are not provided, the tool can evolve a group of stars with equal logarithmic mass intervals.

For evolving a group of binaries:
```shell
petar.[bse_name] [options] -b [binary data file]
```
In this scenario, the binary data file contains 4 values (mass1, mass2, period, eccentricity) per line. The first line specifically includes the number of binaries. The units of mass and period are contingent on the '--mscale' and '--tscale' options. By default, the units are set to Msun and Myr. For instance, if the tscale is not 1.0, the period in Myr is calculated as the input period value multiplied by tscale.

#### Evolution History and Output Files

When there is only one star or one binary, the command line displays the evolution history, including type changes and supernova events. However, if there are multiple stars or binaries, only the final status is printed.

Additional files that document the evolution history of all stars or binaries can be generated by using the `-o` option. The argument for this option specifies the filename prefix. For instance, when `-o output` is employed:
- For single stars, 'output.[sse\_name].type\_change' and 'output.[sse\_name].sn\_kick' files are created.
- For binaries, 'output.[bse\_name].type\_change' and 'output.[bse\_name].sn\_kick' files are generated.

Reading these files follows the same method as for the output files from `petar` (refer to [Data analysis in Python3](#data-analysis-in-python3) and [Gathering Output Files](#gathering-output-files-with-petar-data-gether)).

### Galpy Tool

When the external potential package Galpy is compiled with `--with-external=galpy` during configuration, both `petar.galpy` and `petar.galpy.help` are compiled simultaneously.

`petar.galpy` functions as a standalone tool for computing accelerations and potentials for a particle list using a specified potential model. The basic usage is:
```shell
petar.galpy [options] [particle data filename]
```
This tool can create mesh points to measure potential and aid in drawing potential contours for a given potential set.

In the `petar` command line interface for utilizing Galpy, three options are employed to configure potential models: `--galpy-set`, `--galpy-type-arg`, and `--galpy-conf-file`. The latter two options necessitate users to define potential indices and corresponding arguments. `petar.galpy.help` provides essential information to assist users in setting up these parameters.

To begin, execute the following command:
```shell
petar.galpy.help
```
This tool will assist users in generating a configuration file and provide a list of potentials with mathematical formulas, indices, and arguments. These potentials have been tested and verified to function correctly within `petar`.

Moreover, `petar.galpy.help` offers a comprehensive list of potentials sourced from the official Galpy documentation. This list includes indices and default arguments, although they have not been extensively tested within `petar`. 

It is crucial to acknowledge that the official Galpy documentation is primarily designed for the Python interface, and certain descriptions may not align perfectly with the C interface. Users may need to consult the Galpy C interface source codes to ensure the proper setup of potential arguments.

For detailed insights into a specific potential, users can utilize:
```shell
petar.galpy.help [potential name]
```
This command fetches the potential definition from the official Galpy documentation (note: not extensively tested). Additionally, the `-o` option enables the creation of a configuration file template of a specific potential. After adjustment of this file, it can be read by the `petar` command line option: `--galpy-conf-file`.

## References

It is essential to include the relevant references when publishing results obtained using PeTar. These references can be found in the help function of `petar` and are displayed at the start of the output following a simulation run. Additionally, when activating a feature imported from an external library such as SSE/BSE or Galpy, the corresponding references are automatically included in the output.

## Help Information

The help information provides a comprehensive list of all available options and can be accessed using the '-h' option for `petar` and its associated tools. For instance:
```shell
petar -h
petar.data.process -h
```

These commands offer detailed descriptions of the input particle data file format and available options. It is recommended to review the help information prior to utilizing `petar` and its tools to mitigate errors. Updates corresponding to new PeTar versions are consistently integrated into the help information, keeping users informed about the latest features and functionalities.

## Python Data Analysis Module

PeTar features a Python module designed for data analysis purposes. This module facilitates the reading and analysis of output files generated by `petar` and its tools, identification of multiple systems, calculation of Lagrangian radii and core radii, analysis of system energy errors, and performance evaluation of different parts of the code.

To utilize this module, start by importing it in a Python script, IPython, or Jupyter notebook:
```python
import petar
```
Ensure that after installing PeTar, the `include` directory of PeTar has been added to the `PYTHONPATH` environment variable to successfully import this module.

The module is structured using Python classes and functions. Each class is responsible for reading a specific output file, while functions are utilized for data analysis tasks.

After reading a file using a specific class, users can access the class members for analysis, such as mathematical calculations, data selection, and plotting data. Each class member represents a physical parameter, and its data type can be either a NumPy array or another class as a subclass. The array can be 1D or 2D, storing single or multiple columns from the output file, respectively. The first dimension of the array corresponds to the row index in the data file.

For instance, the `petar.Particle` class is used for reading snapshot files outputted by `petar` during a simulation. This class includes members like `mass`, represented as a 1D NumPy array where each element records the mass of a particle. The array size is equal to the total number of particles ($N$).

The class also contains a member `pos`, a $N\times3$ 2D NumPy array where the first dimension represents the number of particles, and the second dimension represents the position vector of each particle.

If the SSE/BSE is activated, the class includes a member `star`, a subclass of type `petar.SSEStarParameter` that stores the stellar evolution parameters of particles.

All classes have special members like `size`, `ncols`, `keys`, and `host`. These members signify:
- `size`: the data size, equivalent to the size of a 1D array member.
- `ncols`: the total number of columns read from the output files, summing column counts of all data members. For a 2D array member like `pos` of $N\times3$, it counts as 3. If a member is a subclass, it counts the `ncols` value of that subclass.
- `keys`: the names and data types of all data members
- `host`: when the class instance is a member of another class instance, such as `star` in `petar.Particle`, `host` references the parent class instance, otherwise it is `None`.

To utilize these classes, users must create a class instance by executing:
```python
[class instance] = petar.[class name]([keyword arguments])
```
Here, `[class instance]` represents the name of a class instance, and `[class name]` represents the class name, e.g., `petar.Particle`. 
After creating the class instance, users can employ `[class instance].keys` or `[class instance].getColumnInfo()` to view the actual list of members. The latter function also provide the corresponding column index in the output file for each member.

Using `help([class name])` or `help([variable name])` reveals the names, types, definitions of class members (keys) and the available keyword arguments. For description of members, the type "1D" or "2D" indicates a NumPy array member. Otherwise, it denotes the subclass type name, and the corresponding member is its class instance containing a subset of data (multiple columns in data files). Users can further explore member definitions in the subset of data by using `help([class instance].[member name])`. The description of the `__init__` method in the help information provide the available keyword arguments.

In certain classes, the list of members may vary depending on the keyword arguments used to create the instance. `petar.Particle` exemplifies this scenario. The help information of this class categorizes members into groups, and the final member list is a combination based on the keyword arguments used. 

The table below lists all class names (omit the prefix "petar."), corresponding keyword arguments, and output files.
For convenience, the output file prefix used in the following filenames is 'data'.

- Reading outputs of `petar` (requires using `petar.data.gether`):

| Class name     | Description                                  | Keyword Arguments (options shown in [])        |           Corresponding file                                |
| :------------  | :---------------------------------           | :----------------------------------------      | :---------------------------------------------------------- |
| PeTarDataHeader| Header (first line) of snapshot data | external potential option: `external_mode=['galpy', 'none']` |         data.[index]                        |
|                |                                              | Data format: `snapshot_format=['binary', 'ascii']`  |                                                         | 
| Particle | Particle data             | stellar evolution option: `interrupt_mode=['bse', 'mobse', 'bseEmp', 'none']` | data.[index], data.[index].single  |
|                |                                              | external potential option: `external_mode=['galpy', 'none']`    |                                        |
| Status   | Global parameters (energy, N ...)         |                                                | data.status                                                 |
| Profile  | Performance metrics of code parts   | GPU usage: `use_gpu=[True, False]`   | data.prof.rank.[MPI rank]                                   |
| GroupInfo| Multiple systems (binary, triple ...)        | Number of members in systems: `N=[2, 3, ...]`             | data.group.n[number of members]                             |

- Outputs from `petar`  or `petar.[bse name]` with SSE/BSE stellar evolution:

| Class name           | Description                                  |   Corresponding file      |
| :------------        | :---------------------------------           | :------------------------ |
| SSETypeChange  | Type change records of single stars          | data.[sse\_name].type\_change  |
| SSESNKick      | Supernove kick events of single stars        | data.[sse\_name].sn\_kick      |
| BSETypeChange  | Type change records of binary stars          | data.[bse\_name].type\_change  |
| BSESNKick      | Supernove kick events of binary stars        | data.[bse\_name].sn\_kick      |
| BSEDynamicMerge| Dynamically driven mergers                   | data.[bse\_name].dynamic\_merge|
| BSEMerge       | Mergers of binaries                          | Generated by using the class function `combine` to gather information from instances of `petar.BSETypeChange` and `petar.BSEDynamicMerge` types.|
| tide           | Dynamical tide or hyperbolic gravitational wave events | data.[bse_name].tide |

- Outputs from `petar.data.process`:

| Class name         | Description            | Keyword arguments                              |           Corresponding file                                         |
| :------------      | :----------------------| :----------------------------------------      | :----------------------------------------------------------          |
| SingleEscaper| Single star escapers   | same as `petar.Particle`                         | data.esc (from `petar`), data.esc\_single|
| BinaryEscaper| Binary star escapers   | same as `petar.Binary`                           | data.esc_binary                                                    |
| LagrangianMultiple| Lagrangian and core properties | `mass_fraction`, `calc_energy`, `external_mode`, `add_star_type`, `add_mass_range`, `calc_multi_rc`   | data.lagr                                   | 
| Core         | Core radius, center position and velocity of the system  |      | data.core                                                            |
| Binary       | Binary and multiple system | `member_particle_type`, `member_particle_type_one`, `member_particle_type_two`,`interrupt_mode`, `external_mode`, `simple_mode`, `G` | data.[index].binary, data.[index].triple, data.[index].quadruple       |
| BSEStatus    | Stellar evolution statistics: the evolution of number counts, maximum and averaged masses of different stellar types || data.bse\_status                                      |

For more information on the keyword arguments, refers to the [Parallel data process](#parallel-data-process-with-petar-data-process) section for a detailed description of the keyword arguments.

The module also includes several useful functions that can be explored using `help(Function name)` in Python. The following table summarize the names and functionalities:

| Function name      | Description                                                      |
| :---------------   | :--------------------------------------------------------------- |
| join         | Join two data instances of the same type. Example: `petar.join(particle1, particle2)` creates a new `petar.Particle` instance with combined data|
| findPair     | Detect binaries of from a particle snapshot data using `scipy.cKDTree`|
| findMultiple | Detect triples and binary-binary quadruples from single and binary data.|
| parallelDataProcessList| Process a list  of snapshot files in parallel to generate single and binary snapshots, Lagrangian data, core data and escaper data. This is the main function used in `petar.data.process`|
| vecDot       | Dot product of two 2D arrays |
| vecRot       | Rotate a 3D vector array using Euler angles via `scipy.spatial.transform.Rotation`) |
| cantorPairing| Generate a new ID from two IDs, useful for obtaining unique binary IDs|
| calcTrh      | Calculate one-component half-mass relaxation time using Spitzer (1987) formula |
| calcTcr      | Calculate half-mass crossing time |
| calcTGW      | Calculate the merging timescale in Myr for gravitational waves using Peters (1964) formula |
| calcTKL      | Calculate the Kozai-Lidov oscillation timescale (Antognini, 2015)|
| convergentPointCheck | Calculate proper motions in the frame of the convergent point and determine residuals (van Leeuwen F., 2009; Jerabkova T. et al. 2021) |
| petar.coordinateCorrection | Correct the c.m. coordinate based on the difference between the snapshot center and the observational center in the galactocentric frame |

Additional useful tools will be implemented in the future. The `tools/analysis/parallel_data_process.py` serves as a good example for learning how to utilize the analysis tool.

The following sections provide exmaples of using the classes to read and analyze output files.

### Reading Particle Snapshots

Here is an example of using the `petar.Particle` class for reading a snapshot of particles. Snapshots generated by `petar` consist of two parts: a header and data. In ASCII format, the first line is the header, followed by the data of each particle per line.

For example, the first snapshot at time zero generated by `petar` is named 'data.0' by default. To obtain its header information, the `petar.PeTarDataHeader` class can be used:

```python
import petar
header = petar.PeTarDataHeader('data.0')
```

If the keyword argument `snapshot_format` is not specified, the data is assumed to be in ASCII format. To read the BINARY format:

```python
header = petar.PeTarDataHeader('data.0', snapshot_format='binary')
```

After reading, the header contains three members: `fid`, `n`, and `time`, representing the file ID, number of particles, and time (in the same unit as that of the input model) of the snapshot, respectively.

If `external_mode=galpy`, the header contains additional members: the position and velocity offsets, `pos_offset` and `vel_offset`, representing the shift of the system's center in the Galactic tidal field (assuming the coordinate origin is the galactic center).

To read the particle content in ASCII format:

```python
particle = petar.Particle(interrupt_mode='bse')
particle.loadtxt('data.0', skiprows=1)
```

Here, the keyword argument `interrupt_mode` is crucial for proper snapshot reading. The column definitions in snapshots depend on the stellar evolution option (`--with-interrupt`) and the external potential option (`--with-external`) used during configuration. The argument `bse` indicates that the updated SSE/BSE is used, so external columns exist in the snapshots. In this case, `particle` contains a member `star` with the class type `petar.SSEStarParameter`. Similarly, if an external potential is added, one more column `pot_ext` is included.

Since the first line in the snapshot file is the header, `skiprows=1` is used to skip this line when reading data.

If the snapshot data is in BINARY format, the `fromfile` function should be used instead of `loadtxt`:

```python
particle.fromfile('data.0', offset=petar.HEADER_OFFSET)
```

Here, `fromfile` is similar to `numpy.fromfile`. The keyword argument `dtype` is implicitly defined, and users should not change it. The `offset` keyword sets the byte offset to read the data. `petar.HEADER_OFFSET` is a constant indicating the snapshot header line offset, equivalent to `skiprows=1` in the `loadtxt` function. If an external potential is included with the corresponding configuration option `--with-interrupt=galpy`, the offset should be set to `petar.HEADER_OFFSET_WITH_CM`.

### Checking Reading Consistency

When using the Python analysis tool to read data, users must ensure that the keyword arguments in the initialization function align with the options set in the PeTar configuration during installation and the options in the data analysis tool `petar.data.process`. This is crucial to prevent incorrect data reading.

For instance, if the configuration includes bse with `--with-interrupt=bse`, additional columns of stellar evolution parameters are saved in snapshot files. If the initialization of `petar.Particle` lacks the keyword argument `interrupt_mode='bse'` like this:

```python
particle = petar.Particle()
particle.loadtxt('data.0', skiprows=1)
```

The `particle` object will not contain the correct data. In such cases, a warning message will be displayed:

```python
UserWarning: The reading data shape[1] or the number of columns (34) mismatches the number of columns defined in the class instance (20)! Make sure whether this is intended and whether you properly choose the correct keyword arguments for the class instance initialization
```

Here, the snapshot file 'data.0' has 34 columns (excluding the first line) that include the stellar evolution parameters, while `particle` is set to have only 20 columns (ncols=20). This mismatch can lead to incorrect assignment of data columns and class members. Therefore, when this warning appears, users should review whether they have forgotten or incorrectly set the keyword arguments during initialization.

It's important to note that by default, Python displays the warning only once. If the same part of the code is executed multiple times, the warning may not reappear even if the data reading is incorrect. To always show the warning message, users can execute the following command first:

```python
import warnings
warnings.simplefilter('always')
```

### Obtaining Particle Information

Users can access particle information such as masses via `particle.mass`. Since `particle.mass` is a NumPy array, it allows for flexible mathematical operations. For example, to calculate the average mass of all particles, users can use:

```python
import numpy as np
mave = np.average(particle.mass)
```

To compute the center of mass (c.m.) position of particle systems, users can perform the following calculation:

```python
import numpy as np
cm_pos = np.sum(particle.pos * particle.mass[:, None], axis=0) / particle.mass.sum()
```

To retrieve the number of particles (line numbers in snapshots excluding the header):

```python
print(particle.size)
```

### Data Selection

Users can utilize NumPy-like indexing, slicing, and conditional selection for PeTar class instances.

For example, to select the data of the first particle:

```python
first_particle = particle[0]
```

To slice the data of the first 10 particles:

```python
particle_set = particle[:10]
```

To apply a filter and create a subset of data with particle masses less than 1.0:

```python
particle_set = particle[particle.mass < 1.0]
```

This operation is akin to the `getitem` function of a NumPy array.

The indexing, slicing, and selection can also be applied to the members:

```python
pos_set = particle.pos[particle.mass < 1.0]
```

This generates a 2D NumPy array of particle positions with masses below 1.0.

### Plotting Data

Utilizing the NumPy-style data, users can leverage plotting modules such as Matplotlib to create visualizations with PeTar data.

For instance, to plot the mass-distance relation of the subset data with mass < 1.0 using Matplotlib:

```python
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 1)
particle_set = particle[particle.mass < 1.0]
particle_set.calcR2()
axes.plot(np.sqrt(particle_set.r2), particle_set.mass, '.')
```

### Saving and Loading Data

Users have the option to save the post-processed data to a new file.

For instance, to save the subset data generated in the previous section in ASCII format:

```python
particle_set.savetxt([file path])
```

It's important to note that when additional members are added to the particle instance, the saved data will also include the additional columns. In the given example, `particle_set.calcR2()` generates a new class member, `r2` (distance squared). Therefore, the saved data will have an additional column at the end, which did not exist in the original snapshot. When users read this saved data, they should first add the `r2` member to ensure correct column reading:

```python
import numpy as np
particle_new = petar.Particle(interrupt_mode='bse')
particle_new.addNewMember('r2', np.array([]))
particle_new.loadtxt([file path])
```

Users can also save and load data in BINARY format by using `tofile` instead of `savetxt` and `fromfile` instead of `loadtxt`, respectively:

```python
# save
particle_set.tofile([file path])

# load
particle_new.fromfile([file path])
```

Additionally, users can save and load data in NumPy format:

```python
# save
particle_set.save([file path])

# load
particle_new.load([file path])
```

### Reference Frame and Coordinate System Transformation

The `petar.Particle` and `petar.Core` classes provide a member function to convert the data type to `astropy.coordinates.SkyCoord`. `SkyCoord` is a powerful Python module that facilitates easy transformation of the reference frame and coordinate system of particle data.

For instance, when using Galpy for simulations in `petar`, the Galactocentric frame with the Cartesian coordinate system is typically employed. In scenarios where no stellar evolution is utilized, the input unit is in astronomical units (Msun, pc, pc/Myr), and the output format is ASCII. Below is a script to transform a single snapshot without stellar evolution and with an ASCII output format:

```python
import petar
import astropy.units as units

# Read snapshot data from petar
particle = petar.Particle(external_mode='galpy')
particle.loadtxt([snapshot path], skiprows=1)

# Read center of mass offset in the Galactocentric frame
header = petar.PeTarDataHeader([snapshot path], external_mode='galpy')

# Obtain SkyCoord data type in the Galactocentric frame
particle_new = particle.toSkyCoord(pos_offset=header.pos_offset, vel_offset=header.vel_offset)
```

This script will generate a new data set `particle_new` with the type `astropy.coordinates.SkyCoord`. It's important to note that the solar position (`galcen_distance=8.0*units.kpc, z_sun=15.*units.pc`) and velocity (`galcen_v_sun=CartesianDifferential([10.0, 235., 7.]*u.km/u.s)`) assumed in `galpy` differ from the default values of `SkyCoord`. The values in the `toSkyCoord` function align with the `galpy` assumptions.

Subsequently, plotting the data (RA, DEC) in the ICRS frame using Matplotlib becomes convenient:

```python
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 1)
ra = particle_new.icrs.ra
dec = particle_new.icrs.dec
axes.plot(ra, dec, '.')
axes.set_xlabel('RA')
axes.set_ylabel('Dec')
```

### Merging Two Datasets

To combine two subsets of particle data, `particle_set1` and `particle_set2`, into one, you can use the following command:

```python
particle_merge = petar.join(particle_set1, particle_set2)
```

For instance, if the sizes of the two sets are 3 and 5, respectively, the new instance `particle_merge` will have a size of 8.

### Class Functions

Each class provides a set of useful functions for data management. Below is a table outlining these functions:

| Function Name       | Description                                                        |
| :-----------------  | :----------------------------------------------------------------- |
| `addNewMember(key, member)` | Add a new data member with the name `key` and content `member`   |
| `getColumnInfo()`    | Retrieve column positions (counting from 0), key member names, and data types from the data file being read |
| `printTable(column_format, print_title)` | Print a table with specified column list and formats |
| `append(data)`      | Append data to the current class instance                        |
| `resize(N)`         | Resize the size of the class instance                             |
| `loadtxt(file_path)` | Load data in ASCII format from a file                             |
| `savetxt(file_path)` | Save data in ASCII format to a file                               |
| `load(file_path)`   | Load data in NumPy format from a file                             |
| `save(file_path)`   | Save data in NumPy format to a file                               |
| `fromfile(file_path)` | Read data in BINARY format from a file                            |
| `tofile(file_path)` | Save data in BINARY format to a file                              |

In addition to these functions, PeTar classes also support basic mathematical operators such as addition and subtraction. For example, when dealing with two sets of particles, executing the following command will create a new class instance where each member is the sum of the corresponding members in `particle_one` and `particle_two`:

```python
particle_add = particle_one + particle_two
``` 

### Reading Binary Snapshots

In the preceding sections, we demonstrated how to read and analyze snapshots generated by `petar`. Here is an additional example illustrating how to read binary snapshot files produced by the `petar.findPair` function or `petar.data.process`.

When reading binary data, users need to pay attention to the member particle type. A reliable approach to reading snapshots when SSE/BSE is employed is as follows:

```python
p1 = petar.Particle(interrupt_mode='bse')
p2 = petar.Particle(interrupt_mode='bse')
binary = petar.Binary(p1, p2)
binary.loadtxt([binary data file path])
```

Alternatively, users can specify the member particle type explicitly along with other parameters:

```python
binary = petar.Binary(member_particle_type=petar.Particle, interrupt_mode='bse', G=petar.G_MSUN_PC_MYR)
binary.loadtxt([binary data file path])
```

In the 'bse' mode, the gravitational constant `G` should be provided in the unit set of Msun, pc, and Myr (0.00449830997959438). Additionally, `petar.G_HENON` represents the gravitational constant in the Henon unit (equal to 1), which is the default value in `petar.Binary`.

### Reading Triple and Quadruple Snapshots

`petar.Binary` can also handle triple and quadruple snapshots for reading purposes. To read triple snapshots, you can use the following approach:

```python
binary = petar.Binary(member_particle_type_one=petar.Particle, member_particle_type_two=[petar.Particle, petar.Particle])
binary.loadtxt([triple data file path])
```

In this setup, `petar.Binary` is structured as a binary tree, with the first component representing the outer single particle and the second component representing the inner binary.

Likewise, for binary-binary quadruple data, the process is as follows:

```python
binary = petar.Binary(member_particle_type=[petar.Particle, petar.Particle])
binary.loadtxt([quadruple data file path])
```

Here, `member_particle_type` denotes the type (binary) of both components in the quadruple system.

### Reading Lagrangian Data

To extract and analyze Lagrangian data generated by `petar.data.process`, users can follow these steps:

```python
lagr = petar.LagrangianMultiple()
lagr.loadtxt([data.lagr path])
```

The Lagrangian data file contains three subsets: `all`, `single`, and `binary`, representing the Lagrangian data of all objects, single objects, and binaries, respectively. Each member in these subsets is stored as a 2D NumPy array, with properties at various Lagrangian radii and core radii recorded in each row. By default, the mass fraction is set to (0.1, 0.3, 0.5, 0.7, 0.9), resulting in six members in each row, with the final member representing core properties.

For example, to access the half-mass radius and core radius at a mass fraction of 0.5 for all objects:

```python
rh = lagr.all.r[:,2]  # half-mass radius
rc = lagr.all.r[:,5]  # core radius
```

Similarly, to retrieve the average mass at the half-mass radius:

```python
mrh = lagr.all.m[:,2]
```

When utilizing `petar.data.process`, the `-m` option allows users to include an arbitrary set of mass fractions. If the SSE/BSE-based stellar evolution package is used, the `--add-star-type` option in `petar.data.process` can be used to compute Lagrangian properties for specific types of stars. When employing `--add-star-type`, additional subsets such as "BH" and "MS" are appended to the Lagrangian data file.

To ensure accurate data retrieval, it is essential to maintain consistency between the options used in `petar.data.process` and the keyword arguments in `LagrangianMultiple`. When reading the data, remember to include relevant keyword arguments, such as `add_star_type`.

For instance, if users want to calculate Lagrangian properties for a subset of black holes and main sequence stars after running a simulation with SSE/BSE and Galpy enabled, the command would be:

```shell
petar.data.process -i bse -t galpy --add-star-type BH,MS data.snap.lst
```

Subsequently, in Python to read the `data.lagr` file generated by `petar.data.process`:

```python
lagr = petar.LagrangianMultiple(external_mode='galpy', add_star_type=['BH', 'MS'])
lagr.loadtxt([data.lagr path])
```

Inconsistencies in setting the `external_mode` and `add_star_type` arguments may lead to warnings or crashes during data loading. One effective method to ensure consistency is by verifying the column numbers in the data file and validating the reading process, as described in [Checking Reading Consistency](#checking-reading-consistency).

Additionally, the `--add-mass-range` option in `petar.data.process` can compute Lagrangian properties for a subset of particles with specified mass ranges. For detailed information on `--add-star-type` and `--add-mass-range`, refer to `petar.data.process -h`.

### Reading Stellar Evolution Outputs

When using SSE/BSE-based stellar evolution packages like updated BSE in a simulation with `petar` or isolated stellar evolution with `petar.bse`, files `data.bse.*` and `data.sse.*` are generated (assuming the default filename prefix 'data' is used). After executing `petar.data.gether data`, several new files are created with the suffixes "type\_change", "sn\_kick", and "dynamic\_merge", corresponding to stellar or binary type change, supernova kick and dynamically driven merger events.

To read these files, users can utilize the following command:
```python
path='./'
prefix='data'
sse_type=petar.SSETypeChange()
sse_type.loadtxt(path+prefix+'.sse.type_change')
sse_kick=petar.SSESNKick()
sse_kick.loadtxt(path+prefix+'.sse.sn_kick')
bse_type=petar.BSETypeChange()
bse_type.loadtxt(path+prefix+'.bse.type_change')
bse_kick=petar.BSESNKick()
bse_kick.loadtxt(path+prefix+'.bse.sn_kick')
bse_dyn=petar.BSEDynamicMerge()
bse_dyn.loadtxt(path+prefix+'.bse.dynamic_merge')
```
Here, `path` denotes the simulation directory path, and `prefix` represents the filename prefix, which defaults to 'data'.

After reading these files, users can display the events in a table. For instance, to print binary stellar evolution events:
```python
bse_type.printTable([('type','%4d'),('init.time','%10g'),('init.type1','%11d'),('init.type2','%11d'),
                     ('init.m1','%9f'),('init.m2','%9f'),('init.semi','%10g'),('init.ecc','%9f'),
                     ('final.type1','%12d'),('final.type2','%12d'),('final.m1','%9f'),('final.m2','%9f'),
                     ('final.semi','%11g'),('final.ecc','%10g')])
```
In this snippet, `printTable` is the function used to display selected columns in a table with specified printing formats.

By employing `petar.BSEMerge`, users can compile potential merger events from `BSETypeChange` and `BSEDynamicMerge`:
```python
merger = petar.BSEMerge()
merger.combine(bse_type, bse_dyn)
merger.printTable()
```
It's important to note that the initial status from `BSEMerge` pertains to the first event record of the binary that undergoes a merger.  Furthermore, if a supernova completely removes the material of stars without leaving any remnants, it may be mistakenly identified as a "merger". Users should thoroughly examine the stellar evolution history to exclude such false merger events.

### Reading Group Information

During a simulation, PeTar records information about the formation and dissolution of various systems, including hyperbolic encounters, binaries, triples, quadruples, and more, during SDAR integration. This information is saved into files named "data.group.[MPI rank]", where "data" serves as the default filename prefix. By utilizing `petar.data.gether -g [filename prefix]`, groups with different member counts are separated into individual files named "data.group.n1", "data.group.n2", and so on, with the number indicating the number of members in the multiple systems. The following example illustrates how to read the groups file and analyze the information:

```python
# Load two-body groups (binary, hyperbolic encounters)
g2 = petar.GroupInfo(N=2)
g2.loadtxt(path+'data.group.n2')

# Load 3-body groups
g3 = petar.GroupInfo(N=3)
g3.loadtxt(path+'data.group.n3')

# Load 4-body groups
g4 = petar.GroupInfo(N=4)
g4.loadtxt(path+'data.group.n4')
```

For $N>2$, such as `g3` and `g4` in the above example, they contain members named `bin0`, `bin1`, and so on, indicating pairs of objects assuming a Kepler orbit in a hierarchical order. For instance, `bin0` comprises two objects with the greatest separation in the multiple system, and these objects can be either individual particles or the center of mass of the inner systems.

The following example presents a table illustrating the 3-body system, where 'bin0' represents the outer pair and 'bin1' denotes the inner pair. The center of mass of the inner pair is one of the members of the outer pair, with the ID of the center of mass particle being the lower of the two inner member IDs.

```python
g3.printTable([('bin0.m1','%10.4f'),('bin0.m2','%10.4f'),
               ('bin0.p1.id','%11d'),('bin0.p2.id','%11d'),
               ('bin0.semi','%10.4g'),('bin0.ecc','%10.4g'),
               ('bin1.m1','%10.4f'),('bin1.m2','%10.4f'),
               ('bin1.p1.id','%11d'),('bin1.p2.id','%11d'),
               ('bin1.semi','%10.4g'),('bin1.ecc','%10.4g')])
```

### Accessing Tool Manuals Using Python Help

The previous sections have demonstrated how to read and utilize various petar Python classes. The process for using other classes is quite similar, with the distinction lying in the class members.

To access detailed information about each class and function, Python's `help` function can be incredibly useful. For instance, to obtain manuals for specific classes like `petar.Particle` and `petar.GroupInfo`, you can execute the following commands:

```python
help(petar.Particle)
help(petar.GroupInfo)
```

By utilizing the `help` function in Python, users can gain insights into the functionalities, methods, and attributes of different classes within the petar package. This approach provides a convenient way to explore the capabilities and usage guidelines of each class, enabling users to make the most of the petar tools effectively.


# Method:


## Integration Algorithm Overview

PeTar offers three methods for conducting collisional N-body simulations. The integration of particle orbits follows a basic cycle, outlined below. For more in-depth information, please refer to the reference paper on PeTar.

1. Search Clusters and Few-Body Systems
   1. Construct a particle tree that includes all real particles for neighbor searching purposes.
   2. Identify clusters where each real particle remains in the same cluster as its neighbors.
   3. For each cluster:
   - Build a binary tree comprising all real members and identify few-body systems.
   - For stable few-body systems, create artificial particles such as tidal tensor measurement points and orbital-averaged-force sample particles.

2. Kick (Soft) Step
   1. Create a particle tree for all real and artificial particles and compute long-distance (soft) forces with a linear cutoff (using a universal cutoff radius).
   2. Adjust the changeover function by utilizing the neighbor list of each particle.

3. Drift (Hard) Step
   1. For clusters with more than one member:
   - Utilize the 4th order Hermite integration with the slow-down AR method to integrate the orbits.
   2. For a single cluster:
   - Perform the drift operation.

## Parallelization Methods

PeTar employs various parallelization techniques to enhance the efficiency of simulations for large values of $N$. The parallel methods utilized in different components include:

- **Tree Construction**: FPDS is utilized for tree construction, leveraging MPI and OpenMP for parallelization.
- **Soft Force Calculation Kernel**: SIMD (such as AVX, AVX2, AVX512, A64FX) or GPU (CUDA) technologies are employed for soft force calculations.
- **Hard Calculation (Hermite, SDAR)**: OpenMP is utilized for the loop of clusters in hard calculations.

When appropriate tree time steps and changeover radii are set, the performance of hard calculations can be comparable to that of soft force calculations. In such cases, the use of GPUs may not significantly enhance performance, unlike in the direct N-body method for soft forces. Therefore, for optimal performance with large $N$, it is advisable to utilize more CPU cores rather than a few CPU cores with GPUs. This recommendation is particularly relevant when a substantial number of binaries are present in the simulation.

By leveraging these parallelization methods effectively, PeTar can efficiently handle the computational demands of large-scale N-body simulations, ensuring robust and accelerated performance across various simulation scenarios.

## AMUSE API Integration

AMUSE stands as a robust software suite that amalgamates various codes, ranging from hydrodynamic and N-body simulations to stellar evolution and radiation transfer models.

Presently, PeTar has been integrated as a module within AMUSE, offering support for gravitational dynamics and gravity field computations. This integration enables users to couple PeTar with hydrodynamic codes, facilitating simulations of collisional stellar systems immersed in gas environments.

It is important to note that the stopping condition feature is currently under development and not yet operational. As a result, the merging of binaries and rapid stellar evolution events like supernova kicks are not supported at this time.

The integration of PeTar into the AMUSE API opens up new avenues for conducting complex simulations that combine gravitational dynamics with other physical processes, paving the way for comprehensive studies of stellar systems within diverse astrophysical contexts.

