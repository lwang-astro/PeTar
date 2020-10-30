```
    ██████╗ ███████╗████████╗ █████╗ ██████╗ 
    ██╔══██╗██╔════╝╚══██╔══╝██╔══██╗██╔══██╗
    ██████╔╝█████╗     ██║   ███████║██████╔╝
    ██╔═══╝ ██╔══╝     ██║   ██╔══██║██╔══██╗
    ██║     ███████╗   ██║   ██║  ██║██║  ██║
    ╚═╝     ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝
```

A particle-particle \& Particle-tree (_P<sup>3</sup>T_) with slow-down time-transformed symplectic integrator (slow-down algorithmic regularization; _SDAR_) code for simulating gravitational _N_-body systems including close encounters and few-body systems.

The Doxygen document will be provided in doc directory (not yet done).

The reference paper on ArXiv: https://arxiv.org/abs/2006.16560

## Install:
### Dependence:
_FDPS_: https://github.com/FDPS/FDPS

_SDAR_: https://github.com/lwang-astro/SDAR

These two libraries are necessary to compile _PeTar_. 

If users want to use external potential, the _Galpy_ interface is available. Users need to install _Galpy_ either by 
```
pip3 install galpy
```
Or download the source code from https://github.com/jobovy/galpy.

If the source codes of these libraries are put in the same directory where the _PeTar_ directory exist, the configure script can detect them automatically. Otherwise users need to provide the pathes of them by adding configure options
```
./configure --with-[code name in lower case]-prefix=[code path]
```

### Environment:
To successfully compile the code, the C++ compiler (e.g. GNU gcc/g++, Intel icc/icpc, LLVM clang/clang++) needs the support of the C++11 standard. To use SSE/BSE package, a Fortran (77) compiler, GNU gfortran, is needed and should be possile to provide API to the c++ code, i.e., the libgfortran is required. Currently Intel ifort is not supported yet. The MPI compiler (e.g. mpic++) is required to use MPI. NVIDIA GPU and CUDA compiler is required to use GPU acceleration. The SIMD support is tested for the GNU, Intel and LLVM compilers. It is not tested for others, thus these three kinds of compilers are suggested to use. 

To use _Galpy_ and the analysis tools, the _Python3_ should be available. _Galpy_ also requires the _GSL_ library being installed and can be detected in the load library path.

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
   - auto (default): automatical detection MPI, if exists, use MPI compiler, otherwise use non-MPI compiler
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
    - auto (default): automatical detection based on local CPU architecture 
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
    Currently only SSE/BSE is the available stellar evolution package. Notice SSE/BSE is a combined package, the option argument "sse" not work, only "bse" switches on both.
    
    When this option is switched on, the standalone tool _petar.bse_ will also be compiled and installed.
    This is a c++ based tool which call the SSE/BSE functions to evolve single and binary stars to the given age and metallicity. OpenMP parallelization is used to speed up the calculation if a large group of stars and binaries are provided.

8. Use _Galpy_ external potential library
    ```
    ./configure --with-external=galpy
    ```
    The _Galpy_ library is a _Python_ and _c_ based external potential library, which provides a plenty choices of potentials. 
    It is also flexible to combine multiple potentials together (require to use _Galpy_ _Python_ interface to generate the instance, see their document in details).
    
    When this option is switched on, the standalone tool _petar.galpy_ and _petar.galpy.help_ will also be compiled and installed.
    - _petar.galpy_ is a simple tool to call _Galpy_ c interface to evaluate the acceleration and potentials for a list of particles with a given potential model.
    - _petar.galpy.help_ is a tool (python script) to help users to generate the input options for potential models. When users use _Galpy_ _Python_ interface to design a specific potential, this tool also provides a function to convert a _Galpy_ potential instance to an option or a configure file used by _PeTar_.

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

The data analysis module are written in _Python_.
They are installed in [Install path]/include/petar
Please add [Install path]/include to the _Python_ include path (the environment variable, $PYTHONPATH) in order to import the code.

## Use:
After installation, if the [Install path]/bin is in the envirnoment variable, $PATH, the standard way to use the code is
```
[mpiexec -n X] petar [options] [particle data filename]
```
where "[mpiexec -n X]" is used when multiple MPI processors are needed and "X" is the number of processors.

If users have the initial file that stores masses, positions and velocities of each particles per line, generated by tools like _mcluster_ (https://github.com/lwang-astro/mcluster), the _petar.init_ tool (see [Initial input data file](#initial-input-data-file)) can be used to transfer it to the input file of the _petar_ style.


All snapshots of particle data generated in the simulation can be used to restart the simulation. 
To restart the simulation with the same configuration of parameters as before, use
```
[mpiexec -n X] petar -p input.par [options] [snapshot filename]
```
where _input.par_ is automatically generated from the last run (stored in the same diretory of the simulation).
If users change the name of this file and restart, the next new parameter file also uses the new filename.
Notice if users want to use new options in additional to _input.par_, these options must be put after "-p input.par", otherwises they are rewrited by values stored in _input.par_.

In default, after restart, the new data will be appended to the existing output files (with suffixes of esc, group ...).
If users want to restart not from the last time but from a middle _t_, _petar.data.clear_ can be used to remove the data with time > _t_ in the output files.

### options
For _petar_, two types of options exist: single-character options starting with '-' and long options starting with '--'.
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
in the shell configure/initial file (e.g. .bashrc for _bash_) to avoid typing "OMP_STACKSIZE=128M" every time.

#### data format update
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
The content of the status has a style like:
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
- [data filename prefix].group.[MPI rank]: the information of new and end of groups (binary, triple ...), since the items in each row can be different because of different numbers of members in groups, there is no universal column headers. It is suggested to use petar.data.gether first to separate groups with different numbers of members into different files. Then use the python tool _petar.Group_ to read the new files.
- [data filename prefix].[s/b]se.[MPI rank]: if BSE is switched on, the files record the SSE and BSE events, such as type changes, Supernovae, binary evolution phase changes. Each line contain the definition of values, thus can be directly read.

Before access these files, it is suggested to run _petar.data.gether_ tool first to gether the separated files with different MPI ranks to one file for convenience.
This tool also separate the few-body groups with different number of members in xx.group files to individual files with suffix ".n[number of members in groups]".
    
##### Debug dump files
During a long-term simulation, a large number of files of "hard_large_energy.xx", "dump_large_step.xx" and "hard_dump.xx" may be generated.
In the main output, the corresponding warning messages are also printed.
This files record the clusters of stars initial conditions for the hard integration (Hermite+SDAR) of one tree step.
- If the integration error exceeds the limit set in the input option (see ```petar -h```), the dump file "hard_large_energy.xx" appears.
- If the AR step count is so large (exceeding the step limit, which can also be set in the input option) that the performance significantly drops, the dump file "dump_large_step.xx" appears.
- If an error appears so that the code may behaviour abnormaly later on, "hard_dump.xx" appears and the code terminates.

In the first two cases, the simulation continues. The users can ignore them if the results of simulations are acceptable. But if something is not correct, this files can help to finger out the issues, by using the debug tool _petar.hard.debug_.

In the third case, the simulation is terminated in order to avoid unpredictable behaviours. The "hard_dump.xx" usually indicates that a bug may exist. It is suggested to report the issues to the developers by attaching "input.par.hard", "input.par.bse" and "hard_dump.xx" files.
Users can also check the details by using the debug tool _petar.hard.debug_ together with the GDB tool, if users prefer to understand the problems themselves. The knowledge of the source codes of SDAR is required to understand the messages from the debug tool.
  
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
This file can be generated by using the _petar.data.gether_ tool. 
Manually, users can generate the file by using _ls_ and _egrep_ in Linux.
For example, if the snapshot filename prefix is 'data', the commander is
```
ls |egrep '^data.[0-9]+$'
```
Notice it is better to sort the path in the increasing order of evolution time.
The _sort_ tool is very convenient to get the time-sorted list:
```
ls |egrep '^data.[0-9]+$' |sort -n -k 1.6 >snap.lst
```
will find all data files in the current directory, sort files by using the suffix (values after 'data.') with an increasing order and save the list to the file of _snap.lst_.
Here '-n' indicates that the values to sort are floating-point number.
'-k' defines which word (separated by space) and which character in the word is the starting position of the number for sorting.
In this example, '1.6' represents the 6th chararcters in the first word is the starting position to recognize the numbers for sorting.
This is because only one word exists per line in the snap.lst file, and the word 'data.\*\*' has 5 characters before the number (index of the snapshot). 

The data order in Lagrangian, escapers and core data file follows the order in the snapshot path list.
The Lagrangian and core data can be read by _LagrangianMultiple_ and _Core_ modules by using _loadtxt_(filename).
The escaper data (single and binary) can be read by _SingleEscaper_ and _BinaryEscaper_.
By the way, the _petar_ code can also remove escapers and stored the data of escapers by using the energy and distance criterion (see help of _petar_).

When the snapshot files are in BINARY format, the option `-s binary` can be used to read the snapshot correctly.
Notice the generated data from _petar.data.process_ are all in ASCII format.

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

#### gether output files from different MPI ranks
The _petar.data.gether_ is used to gether output files from different MPI ranks to one file (e.g. xx.group.[MPI rank]).
For group files, it also generate individual files with suffix ".n[number of members in groups]'.
In addition, the tool generates a file, "[output prefix].snap.lst" , that contains the list of all snapshot files sorted by time. This can be used as the input for _petar.data.process_ and _petar.movie_.
The basic usage is
```
petar.data.gether [options] [data filename prefix]
```
If the [output prefix] is not given (option -f), the [output prefix] is the same as [data filename prefix].

When _SSE_/_BSE_ is used and the code version is before Sep 10, 2020, the data with suffix [.dynamic_merge] has three column less than that of the new version.
The corresponding data are dr, t_peri and sd_factor.
This tool will automatically fill the missing columns with zero (from Column 6 to 8)

#### remove data after a given time
The _petar.data.clear_ is used to remove data after a given time for all output files except the snapshots.
The basic usage is 
```
petar.data.clear -t [time] [data filename prefix]
```
All output files will be firstly backup by rename the file name with an additional suffix '.bk'.
If this tool is mistakely used, the files can be recovered.
If the tool is used again, it will check the backup data files.
Thus a new time criterion larger than the previous one can be to include more data in the new files.

#### data format transfer (BINARY - ASCII)
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

#### SSE/BSE steller evolution tool
The _petar.bse_ will be generated when --with-interrupt=bse is used during the configuration.
This is the standalone SSE/BSE tool to evolve stars and binaries to a given time.
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
    - _Particle_: the basic particle data (snapshot files, [data filename prefix].[index]).
    - _Status_: the global parameter of the system such as energy and number of particles ([data filename prefix].status).
    - _Profile_: the wall-clock time of different parts of the code ([data filename prefix].prof.rank.[rank index]).
    - _GroupInfo_: the formation and disruption of few-body groups log ([data filename prefix].group.n[number of members]).
- For outputs when SSE/BSE is switched on (need to use _petar.data.gether_ to generate data files first):
    - _SSETypeChange_: the log of type change of single stars ([data filename prefix].sse.type_change)
    - _SSESNKick_: the log of SNe kick events of single stars ([data filename prefix].sse.sn_kick)
    - _BSETypeChange_: the log of type change of binary stars ([data filename prefix].bse.type_change)
    - _BSESNKick_: the log of SNe kick events of binary stars ([data filename prefix].bse.sn_kick)
- For data generated by using the _petar.data.process_:
    - _SingleEscaper_: single star escapers, this can read both the escapers generated by _petar.data.process_ and the output from _petar_ ([data filename prefix].esc(\_single)).
    - _BinaryEscaper_: binary star escapers ([data filename prefix].esc\_binary).
    - _LagrangianMultiple_: properties related to Lagrangian and core radii: radii, number, averaged mass, velocity dispersion, rotational velocity in x-y plane ([data filename prefix].lagr).
    - _Core_: data of core radius, center position and velocity of the system ([data filename prefix].core). 
    - _Binary_: the basic binary data, store two member particles and the Kepler orbital information. Can be generated by using _findPair_ function.
    - _BSEStatus_: the evolution of number counts, maximum and averaged masses of different stellar types ([data filename prefix].bse_status).

There are also several useful functions.
- _join_: join two same type instances of modules. For example, _join_(particle1, particle2) will generate a new _Particle_ instance that contain both two data. Each member is numpy.append(particle1.member, particle2.member).
- _findPair_: detect binaries of one particle list by using _scipy.cKDTree_.
- _findMultiple_: detect triples and quadruples (binary-binary) from single and binary data (generated from _findPair_).
- _parallelDataProcessList_: use mutliple CPU cores to process a list of snapshot files and generate single and binary snapshots, Lagrangian data, core data and escaper data. For large _N_, the data process is quite slow, thus using multiple CPU processors can speed up the process. This is the main function used in _petar.data.process_.

More useful tools will be implemented in the future. The tools/analysis/parallel_data_process.py is a good example to learn how to use the analysis tool.

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

To read the particle content in ASCII format:
```
particle=petar.Particle(interrupt_mode='bse')
particle.loadtxt('data.0',skiprows=1)
```
Here the keyword argument ```interrupt_mode``` is important to set properly in order to read the snapshot correctly.
The column definitions of snapshots depends on the stellar evolution option (--with-interrupt) used during the configure.
The argument 'bse' indicates that the SSE/BSE is used so that external columns of SSE/BSE parameters exist in the snapshots.
Then, here `particle` contains a member `star` with the class type `petar.SSEStarParameter`.

Since the first line in the snapshot file is the header, `skiprows=1` jumps this line when read data.

If the snapshot data is in the BINARY format, the _fromfile_ function should be used instead of _loadtxt_:
```
particle.fromfile('data.0', offset=petar.HEADER_OFFSET)
```
Here the _fromfile_ is similar to _numpy.fromfile_. The keyword argument `dtype` is defined implicitly, users should not change it. The keyword argument `offset` sets the offset to read the data in bytes. 
`petar.HEADER_OFFSET` is the constant to indicate the offset of snapshot header line.
This `offset` is equivalent to `skiprows=1` in the _loadtxt_ function.

To get the number of particles (line numbers in snapshots excluding the first line):
```
print(particle.size)
```

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

To plot the mass-distance relation of the subset using _matplotlib_,
```
import matplotlib.pyplot as plt
fig, axes=plt.subplots(1,1)
particle_set.calcR2()
axes.plot(np.sqrt(particle_set.r2), particle_set.mass, '.')
```

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

To join two subset of data, p1 and p2, to one,
```
particle_merge=petar.join(p1,p2)
```
For example, if p1 and p2 have sizes of 3 and 5, respectively, the new instance, particle_merge, has a size of 8.

This example show how to read and use Particle class, the way to use other classes is very similar, only the class member is different.
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
