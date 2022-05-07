#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([PeTar], [1.0], [longwang.astro@live.com])
m4_include([autoconf_m4/ax_check_compile_flag.m4])
m4_include([autoconf_m4/ax_prog_cxx_mpi.m4])
m4_include([autoconf_m4/ax_openmp.m4])
#AC_CONFIG_SRCDIR([src/main.cc])
#AC_CONFIG_HEADERS([config.h])

PROG_NAME="petar"

build_linux=no 
build_windows=no
build_mac=no
host_os=`uname`

# Detect the target system
AC_MSG_CHECKING([system type])
case "${host_os}" in
    Linux*)
        build_linux=yes
        ;;
    Cygwin*|Mingw*)
        build_windows=yes
        ;;
    Darwin*)
        build_mac=yes
        ;;
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac
AC_MSG_RESULT([$host_os])

CXXFLAGS_BK=$CXXFLAGS
CFLAGS_BK=$CFLAGS

AC_LANG([C++])

CXXFLAGS=$CXXFLAGS_BK
CFLAGS=$CFLAGS_BK

AC_PROG_AWK
AC_PROG_GREP

# Checks for programs.
# MPI
  # If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
  # If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
  # If --with-mpi=no is used, use a standard C compiler instead.

AC_ARG_WITH([mpi],
            [AS_HELP_STRING([--with-mpi],
                            [compile with MPI (parallelization) support. If none is found, MPI is not used (auto, yes, no). Default: auto])],
            [],
            [with_mpi=auto])


AX_PROG_CXX_MPI([test x"$with_mpi" != xno],
                [use_mpi=yes
                 PROG_NAME=$PROG_NAME".mpi"],
                [use_mpi=no])
#                 CXXNOMPI=`$CXX -show|awk '{print $1}'`],
#                 CXXNOMPI=$CXX])

AS_IF([test x"$use_mpi" == xno],
      [AS_IF([test x"$with_mpi" = xyes],
             [AC_MSG_FAILURE([MPI compiler requested, but could not use MPI.])],
             [AC_MSG_WARN([No MPI compiler found, disable MPI.])])])

# debug
AC_ARG_WITH([debug],
            [AS_HELP_STRING([--with-debug],
                            [Switch on debugging mode (g: enable debugger support; assert: enable assertion; no: none). Default: no])],
            [],
            [with_debug=no])

DEBFLAGS=' -g'

AX_CHECK_COMPILE_FLAG([-O0],
                      [DEBFLAGS=$DEBFLAGS" -O0"])

AX_CHECK_COMPILE_FLAG([-fbounds-check],
                      [DEBFLAGS=$DEBFLAGS" -fbounds-check"])

AS_IF([test x"$with_debug" == xg],
      [PROG_SUFFIX=".g"
       OPTFLAGS=$DEBFLAGS],
      [AS_IF([test x"$with_debug" != xno],
             [PROG_SUFFIX=".d"])
       AX_CHECK_COMPILE_FLAG([-O3],
                             [OPTFLAGS=" -O3"])])


AX_CHECK_COMPILE_FLAG([-Wall],
                      [OPTFLAGS=$OPTFLAGS" -Wall"
                       DEBFLAGS=$DEBFLAGS" -Wall"])

AC_ARG_WITH([cpp-std],
	    [AS_HELP_STRING([--with-cpp-std],
                            [version of c++ standard (11, 17). Default: 17])],
            [],
	    [with_cpp_std=17])

AX_CHECK_COMPILE_FLAG([-std=c++$with_cpp_std],
                      [OPTFLAGS=$OPTFLAGS" -std=c++"$with_cpp_std
                       DEBFLAGS=$DEBFLAGS" -std=c++"$with_cpp_std
		       CUDAFLAGS=$CUDAFLAGS" -std=c++"$with_cpp_std],
                      [AC_MSG_FAILURE([c++$with_cpp_std support is required for the compiler])])

# OpenMP
AC_ARG_ENABLE([openmp],
              [AS_HELP_STRING([--disable-omp],
                              [disable OpenMP support])],
              [use_omp=no],
              [use_omp=yes])

AS_IF([test x"$use_omp" != xno],
      [AX_OPENMP([use_omp=yes
                  PROG_NAME=$PROG_NAME'.omp'],
                 [use_omp=no])],
      [use_omp=no])

# x86
AC_ARG_WITH([arch],
	    [AS_HELP_STRING([--with-arch],
                            [compile with architecture type (x86, fugaku). Default: x86])],
            [],
	    [with_arch=x86])

AC_ARG_WITH([simd],
            [AS_HELP_STRING([--with-simd],
                [compile with x86 architecture support (avx, avx2, avx512). Default: auto])],
            [],
            [with_simd=auto])

AS_IF([test x"$with_arch" == xx86],
      [use_simd=no
       AS_IF([test x"$with_simd" != xauto],
             [AS_CASE($with_simd,
                      [avx],[arch=core-avx-i],
                      [avx2],[arch=core-avx2],
		      [avx512],[arch=skylake-avx512],
		      [AC_MSG_FAILURE([$with_simd is not support by the compiler])])
	      AX_CHECK_COMPILE_FLAG([-march=$arch],
                                    [SIMDFLAGS=" -march=$arch"
                                     SIMD_TYPE=$with_simd
                                     use_simd=yes],[])],
             [m4_foreach([ax_x86_feature],
                         [[core-avx-i,avx], [core-avx2,avx2], [skylake-avx512,avx512]],
                         [arch=m4_argn(1,ax_x86_feature)
                          simd_name=m4_argn(2,ax_x86_feature)
                          AX_CHECK_COMPILE_FLAG([-march=$arch],
                                                [AC_MSG_CHECKING([whether host CPU support $simd_name])
                                                 AS_IF([test x"$build_mac" != xno],
						       [simd_check=`sysctl -a|grep 'machdep.cpu.*features'|grep -i -c $simd_name`],
						       [simd_check=`grep -c $simd_name /proc/cpuinfo`])
                                                 AS_IF([test $simd_check -gt 0],
                                                       [SIMDFLAGS=" -march=$arch"
                                                        SIMD_TYPE="$simd_name"
                                                        AC_MSG_RESULT([yes])
                                                        use_simd=yes],
                                                       [AC_MSG_RESULT([no])])],
                                                [])])])
       AS_IF([test x"$use_simd" != xno],
             [PROG_NAME=$PROG_NAME"."$SIMD_TYPE])],
      [use_simd=no])
      
AS_IF([test x"$with_arch" == xfugaku],
      [AX_CHECK_COMPILE_FLAG(['-Kfast -Nclang'],
                             [OPTFLAGS=$OPTFLAGS" -Kfast -Nclang"])
       PROG_NAME=$PROG_NAME".fugaku"])

AC_ARG_ENABLE([simd-64],
              [AS_HELP_STRING([--enable-simd-64],
                              [enable x86 SIMD 64 bits floating point support for long-distant tree force])],
              [PROG_NAME=$PROG_NAME".64b"
	       use_simd_64=yes],
              [use_simd_64=no])

# QUAD
AC_ARG_ENABLE([quad],
              [AS_HELP_STRING([--disable-quad],
                              [disable quadrupole-moment calculation for super particles])],
              [use_quad=no],
              [use_quad=yes])

# cuda
AC_ARG_ENABLE([cuda],
              [AS_HELP_STRING([--enable-cuda],
                              [enable CUDA (GPU) acceleration support for long-distant tree force])],
              [use_cuda=yes],
              [use_cuda=no])

AC_ARG_WITH([cuda-prefix], 
            [AS_HELP_STRING([--with-cuda-prefix],
	                    [Prefix of your CUDA installation])],
            [],
	    [cuda_prefix="/usr/local/cuda"])

AC_ARG_WITH([cuda-sdk-prefix], 
            [AS_HELP_STRING([--with-cuda-sdk-prefix],
	                    [Prefix of your CUDA samples (SDK) installation])],
            [cuda_sdk_prefix=$withval],
	    [cuda_sdk_prefix="/usr/local/cuda/samples"])

AS_IF([test "x$use_cuda" != xno],
      [AC_CHECK_PROG([NVCC], [nvcc], [nvcc], [none])
       use_cuda=yes
       AS_IF([test "x$NVCC" == xnone], 
             [AC_PATH_PROG([NVCC], [nvcc], [none], [$cuda_prefix/bin])
              AS_IF([test "x$NVCC" == xnone], 
                    [AC_MSG_FAILURE([can't find CUDA compiler nvcc, please check whether nvcc is in environment PATH or use --with-cuda-prefix to provide the PATH of CUDA installation])])],
             [AC_PROG_SED
              cuda_prefix=`which nvcc|sed 's:/bin/nvcc::'`])

       PROG_NAME=$PROG_NAME".gpu"
       AC_MSG_CHECKING([CUDA version])
       CUDA_VERSION=`$NVCC --version|awk -F ',' 'END{print $(NF-1)}'|awk '{print $2}'`
       AC_MSG_RESULT([$CUDA_VERSION])

       AC_CHECK_FILE(["$cuda_prefix/samples/common/inc/helper_cuda.h"],
                     [CUDAFLAGS=$CUDAFLAGS" -I $cuda_prefix/samples/common/inc"],
                     [AC_CHECK_FILE(["$cuda_sdk_prefix/common/inc/helper_cuda.h"],
                                    [CUDAFLAGS=$CUDAFLAGS" -I $cuda_sdk_prefix/common/inc"],
                                    [AC_CHECK_FILE(["$cuda_sdk_prefix/Common/helper_cuda.h"],
                                                   [CUDAFLAGS=$CUDAFLAGS" -I $cuda_sdk_prefix/Common"],
                                                   [ AC_MSG_FAILURE([can't find CUDA sample inc file helper_cuda.h, please provide correct cuda SDK prefix by using --with-cuda-sdk-prefix])])])])

       AC_CHECK_LIB([cudart],
                    [main], 
                    [CUDALIBS=' -lcudart'],
                    [AC_CHECK_FILE(["$cuda_prefix/lib64/libcudart.so"],
                                   [CUDALIBS=" -L$cuda_prefix/lib64 -lcudart"],
                                   [AC_MSG_FAILURE([can't find CUDA library -lcudart, please provide correct cuda PREFIX by using --with-cuda-prefix])])])
       AC_CHECK_LIB([gomp], 
                    [main],
		    [CUDALIBS=$CUDALIBS' -lgomp'],
		    [AC_CHECK_FILE(["$cuda_prefix/lib64/libgomp.so"],
                                   [CUDALIBS=$CUDALIBS" -lgomp"],
                                   [AC_MSG_FAILURE([can't find CUDA library -lgomp, please provide correct cuda PREFIX by using --with-cuda-prefix])])])],
       [use_cuda=no])


# gperftools
AC_ARG_ENABLE([gperf],
              [AS_HELP_STRING([--enable-gperf],
                              [enable gperftools for profiling])],
              [use_gperf=yes],
              [use_gperf=no])

AC_ARG_WITH([gperf-prefix], 
            [AS_HELP_STRING([--with-gperf-prefix],
	                    [Prefix of your gperftools installation])],
            [gperf_prefix=$withval],
	    [gperf_prefix=/opt/gperftools])

AS_IF([test x"$use_gperf" != xno],
      [use_gperf=yes
       AC_CHECK_LIB([profiler],
                    [main],
		    [GPERFLIBS=' -lprofiler'],
		    [AC_CHECK_FILE(["$gperf_prefix/lib/libprofiler.so"],
                                   [GPERFFLAGS=" -I$gperf_prefix/include"
                                    GPERFLIBS=" -L$gperf_prefix/lib -lprofiler"],
				   [AC_MSG_FAILURE([can't find gperftools library -lprofiler, please provide correct gperftools PREFIX by using --with-gperf-prefix])])])
       AC_CHECK_LIB([tcmalloc],
                    [main],
		    [GPERFLIBS=' -ltcmalloc'],
		    [AC_CHECK_FILE(["$gperf_prefix/lib/libtcmalloc.so"],
                                   [GPERFLIBS=$GPERFLIBS" -ltcmalloc"],
				   [AC_MSG_FAILURE([can't find gperftools library -ltcmalloc, please provide correct gperftools PREFIX by using --with-gperf-prefix])])])],
       [use_gperf=no])

# FDPS
AC_ARG_WITH([fdps-prefix], 
            [AS_HELP_STRING([--with-fdps-prefix],
	                    [Prefix of your FDPS installation. Default: ../FDPS])],
            [fdps_prefix=$withval],
	    [fdps_prefix="../FDPS"])

AC_CHECK_FILE(["$fdps_prefix/src/particle_simulator.hpp"],
	      [INCLUDE=$INCLUDE" -I$fdps_prefix/src"],
	      [AC_MSG_FAILURE([can't find FDPS library, please provide correct PREFIX by using --with-fdps-prefix])])

AC_CHECK_FILE(["$fdps_prefix/src/ps_defs.hpp"],
	      [FDPS_VERSION=`egrep 'Version' $fdps_prefix/src/ps_defs.hpp | sed 's:.*Version\ \(\S\+\).*:\1:'`],
	      [AC_MSG_FAILURE([can't find FDPS ps_defs.hpp, please check the completeness of the FDPS library])])

AS_IF([test ${FDPS_VERSION%%.*} -ge 7 ],
      [use_fdps_comm=yes],
      [use_fdps_comm=no])

# SDAR
AC_ARG_WITH([sdar-prefix], 
            [AS_HELP_STRING([--with-sdar-prefix],
	                    [Prefix of your FDPS installation. Default: ../SDAR])],
            [sdar_prefix=$withval],
	    [sdar_prefix="../SDAR"])

AC_CHECK_FILE(["$sdar_prefix/src/Hermite/hermite_integrator.h"],
	      [INCLUDE=$INCLUDE" -I$sdar_prefix/src"],
	      [AC_MSG_FAILURE([can't find SDAR library, please provide correct PREFIX by using --with-sdar-prefix])])

# GALPY
AC_ARG_WITH([galpy-prefix],
            [AS_HELP_STRING([--with-galpy-prefix],
	                    [Prefix of your Galpy installation. Default: search $HOME/.local, ../galpy])],
            [galpy_prefix=$withval],
	    [galpy_prefix="__search__"])

# Leap frog step method
AC_ARG_WITH([step-mode],
	    [AS_HELP_STRING([--with-step-mode],
	                    [Tree step mode: K[0.5]-D[1.0]-K[0.5] (kdk), K[1/6]-D[1/2]-K[4/6]-D[1/2]-K[1/6] (kdkdk). Default: kdk])],
            [PROG_NAME=$PROG_NAME'.'$with_step_mode],
	    [with_step_mode=kdk])

# interruption options
FCFLAGS_BK=$FCFLAGS
AC_ARG_WITH([interrupt],
	    [AS_HELP_STRING([--with-interrupt],
	                    [Interruption detection mode: Default: off])],
            [PROG_NAME=$PROG_NAME'.'$with_interrupt
	     AS_IF([[test x"$with_interrupt" == xbse] || [test x"$with_interrupt" == xmobse] || [test x"$with_interrupt" == xbseEmp]],
	           [use_fortran=yes
		    AC_LANG([Fortran])
                    #AC_PROG_FC([gfortran ifort flang])
		    AS_IF([test x"$FC" == xgfortran],
		          [FCLIBS=$FCLIBS' -lgfortran'],[])
		    AS_IF([test x"$FC" == xflang],
		          [#AX_CHECK_COMPILE_FLAG([-static-flang-libs],
                          #                       [FCFLAGS=$FCFLAGS" -static-flang-libs"])
		          FCLIBS=$FCLIBS' -lflang'],[])
                    FCFLAGS=$FCFLAGS_BK
		    AS_IF([test x"$with_debug" == xg],
                          [AX_CHECK_COMPILE_FLAG([-g],
                                                 [FCFLAGS=$FCFLAGS" -g"])
                           AX_CHECK_COMPILE_FLAG([-O0],
                                                 [FCFLAGS=$FCFLAGS" -O0"])
                           AX_CHECK_COMPILE_FLAG([-finit-real=nan],
                                                 [FCFLAGS=$FCFLAGS" -finit-real=nan"])
                           AX_CHECK_COMPILE_FLAG([-fbounds-check],
                                                 [FCFLAGS=$FCFLAGS" -fbounds-check"])],
                          [AX_CHECK_COMPILE_FLAG([-O3],
                                                 [FCFLAGS=$FCFLAGS" -O3"])])
                    AX_CHECK_COMPILE_FLAG([-Wall],
                                          [FCFLAGS=$FCFLAGS" -Wall"])])],
	    [with_interrupt=off])

# external potential
AC_ARG_WITH([external],
	    [AS_HELP_STRING([--with-external],
	                    [External potential: Default: off])],
            [PROG_NAME=$PROG_NAME'.'$with_external
	     AS_IF([test x"$with_external" == xgalpy],
	           [AC_LANG([C])
		    #AC_PROG_CC([fcc gcc icc cc])
                    use_c=yes
		    AS_IF([test x"$with_debug" == xg],
                          [AX_CHECK_COMPILE_FLAG([-g],
                                                 [COPTFLAGS=$COPTFLAGS" -g"])
                           AX_CHECK_COMPILE_FLAG([-O0],
                                                 [COPTFLAGS=$COPTFLAGS" -O0"])],
                          [AX_CHECK_COMPILE_FLAG([-O3],
                                                 [COPTFLAGS=" -O3"])])
                    AX_CHECK_COMPILE_FLAG([-Wall],
                                          [COPTFLAGS=$COPTFLAGS" -Wall"])
                    AS_IF([test x"$with_arch" == xfugaku],
                          [AX_CHECK_COMPILE_FLAG(['-Kfast -Nclang'],
                                                 [COPTFLAGS=$COPTFLAGS" -Kfast -Nclang"])])
                    AS_IF([test x"$galpy_prefix" == x__search__],
                          [AC_CHECK_PROG([FIND_CHECK], [find], [find])
			   AS_IF([test x"$FIND_CHECK" == xfind],
                                 [galpy_prefix=`find $HOME/.local -type d -name galpy|tail -1`
                                  AS_IF([test x"$galpy_prefix" == x],
                                        [galpy_prefix=`pwd`'/'`find ../ -type d -name galpy|tail -1`])])])
                    AC_CHECK_FILE(["$galpy_prefix/potential/potential_c_ext/galpy_potentials.h"],
                                  [GALPY_PATH=$galpy_prefix],
                                  [AC_MSG_FAILURE([can't find Galpy library, please provide correct PREFIX by using --with-galpy-prefix])])
                    AC_CHECK_FILE(["$galpy_prefix/__init__.py"],
		                  [GALPY_VERSION=`cut -d'"' -f2 $galpy_prefix/__init__.py`],
				  [AC_MSG_FAILURE([can't find Galpy __init__.py, which contains the version number])])
                    AC_CHECK_PROG([GSL_CFLAGS],[gsl-config],[`gsl-config --cflags`],[])
                    AC_CHECK_PROG([GSL_LIBS],[gsl-config],[`gsl-config --libs`],[none])
		    AS_IF([test x"$GSL_LIBS" == xnone],
                          [AC_CHECK_LIB([gslcblas],[main],
                                        [],
                                        [AC_MSG_FAILURE([can't find GSL library, please install it first, e.g. in Ubuntu 20.04, use sudo apt-get install libgsl-dev])])
                           AC_CHECK_LIB([gsl],[main],
                                        [],
                                        [AC_MSG_FAILURE([can't find GSL library, please install it first, e.g. in Ubuntu 20.04, use sudo apt-get install libgsl-dev])],[-lgslcblas])
                           GSL_LIBS='-lgsl -lgslcblas'])])],
	    [with_external=off])
            

# tidal tensor 3rd order mode
AC_ARG_WITH([tidal-tensor],
	    [AS_HELP_STRING([--with-tidal-tensor],
	                    [Tidal tensor accuracy mode: 2nd/3rd: Default: 3rd])],
            [with_tidal_tensor=$withval
	     PROG_NAME=$PROG_NAME'.'tt$with_tidal_tensor],
	    [with_tidal_tensor=3rd])

# orbit particle methods
AC_ARG_WITH([orbit],
	    [AS_HELP_STRING([--with-orbit],
	                    [Orbit particle method for counter force: os: orbit-sampling; pm: pseudoparticle multipole. Default: pm])],
            [PROG_NAME=$PROG_NAME'.'$with_orbit],
	    [with_orbit=pm])


PROG_NAME=$PROG_NAME$PROG_SUFFIX

# others
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

AC_PREFIX_PROGRAM(petar)

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

# version
AC_CHECK_FILE([$sdar_prefix/VERSION],
	      [sdar_version=`cat $sdar_prefix/VERSION`],
	      [AC_MSG_FAILURE([can't find the file VERSION in SDAR directory, please update SDAR])])
AC_CHECK_FILE([VERSION],
	      [petar_version=`cat VERSION`],
	      [AC_MSG_FAILURE([can't find the file VERSION])])
VERSION=$petar_version'_'$sdar_version

AC_SUBST([with_arch])
AC_SUBST([use_simd])
AC_SUBST([use_simd_64])
AC_SUBST([use_quad])
AC_SUBST([use_mpi])
AC_SUBST([use_cuda])
AC_SUBST([use_omp])
AC_SUBST([use_gperf])
AC_SUBST([with_debug])
AC_SUBST([with_step_mode])
AC_SUBST([with_interrupt])
AC_SUBST([with_external])
AC_SUBST([with_tidal_tensor])
AC_SUBST([with_orbit])
AC_SUBST([OPTFLAGS])
AC_SUBST([COPTFLAGS])
AC_SUBST([DEBFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([FCFLAGS])
AC_SUBST([CFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])
AC_SUBST([NVCC])
AC_SUBST([CUDAFLAGS])
AC_SUBST([CUDALIBS])
AC_SUBST([SIMDFLAGS])
AC_SUBST([GPERFFLAGS])
AC_SUBST([GPERFLIBS])
AC_SUBST([GALPY_PATH])
AC_SUBST([GSL_CFLAGS])
AC_SUBST([GSL_LIBS])
AC_SUBST([INCLUDE])
AC_SUBST([PROG_NAME])
AC_SUBST([CXX])
AC_SUBST([CXXNOMPI])
AC_SUBST([VERSION])
AC_SUBST([FDPS_VERSION])
AC_SUBST([use_fdps_comm])
AC_SUBST([FCLIBS])
AC_SUBST([FLIBS])
AC_SUBST([GALPY_VERSION])

AC_MSG_NOTICE([----------------------Configure Summary--------------------])
AC_MSG_NOTICE([--Features:])
AC_MSG_NOTICE([     Version:           $VERSION])
AC_MSG_NOTICE([     FDPS Version:      $FDPS_VERSION])
AC_MSG_NOTICE([     Using MPI:         $use_mpi])
AC_MSG_NOTICE([     Using GPU:         $use_cuda])
AC_MSG_NOTICE([     Using arch:        $with_arch])
AS_IF([test "x$use_simd" != xno],
      [AS_IF([test "x$use_simd_64" != xno],
             [AC_MSG_NOTICE([     SIMD type:         $SIMD_TYPE 64 bits])],
             [AC_MSG_NOTICE([     SIMD type:         $SIMD_TYPE ])])
       AS_IF([test "x$with_simd" == xauto],
             [AC_MSG_NOTICE([  Notice: this is auto-detected based on the host CPU architecture])
              AC_MSG_NOTICE([          If different CPU is used for running, check whether $SIMD_TYPE is also supported])])])
AC_MSG_NOTICE([     Using OpenMP:      $use_omp])
AC_MSG_NOTICE([     Debug mode:        $with_debug])
AC_MSG_NOTICE([     Step mode:         $with_step_mode])
AC_MSG_NOTICE([     interrupt mode:    $with_interrupt])
AC_MSG_NOTICE([     external mode:     $with_external])
AS_IF([test x"$with_external" == xgalpy],
      [AC_MSG_NOTICE([     Galpy version:     $GALPY_VERSION])])
AC_MSG_NOTICE([     tidal tensor mode: $with_tidal_tensor])
AC_MSG_NOTICE([     orbit mode:        $with_orbit])
AC_MSG_NOTICE([     Using quad:        $use_quad])
AC_MSG_NOTICE([--Compilers:])
AC_MSG_NOTICE([     C++ compiler:      $CXX])
AS_IF([test "x$use_c" == xyes],
      [AC_MSG_NOTICE([     C compiler:        $CC])])
AS_IF([test "x$use_fortran" == xyes],
      [AC_MSG_NOTICE([     Fortran compiler:  $FC])])
AS_IF([test "x$use_cuda" != xno],
      [AC_MSG_NOTICE([     CUDA compiler:     $NVCC])
       AC_MSG_NOTICE([     CUDA version:      $CUDA_VERSION])])
AC_MSG_NOTICE([--Targets:])
AC_MSG_NOTICE([     Program name:      $PROG_NAME])
AS_IF([test x"$prefix" != xNONE],
      [AC_MSG_NOTICE([     Install path:      $prefix])
       AC_MSG_NOTICE([-----------------------------------------------------------])
       AC_MSG_NOTICE([  After install, please])
       AC_MSG_NOTICE([      add $prefix/include into your PTHYON include path])
       AC_MSG_NOTICE([      add $prefix/bin into your PATH])],
      [AC_MSG_WARN([Install path is not defined, please provide it by using --prefix=INSTALL_PATH])])
AS_IF([test "x$use_omp" != xno],
      [AC_MSG_NOTICE([      add 'export OMP_STACKSIZE=128M' in your bash configure file (e.g. .bashrc) to avoid segmentational fault])])
AC_MSG_NOTICE([-----------------------------------------------------------])
AC_MSG_NOTICE([Reference paper: see the beginning of the output when using PeTar])
AC_MSG_NOTICE([Questions and bug reports: https://github.com/lwang-astro/PeTar/issues])
AC_MSG_NOTICE([-----------------------------------------------------------])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([bse-interface/Makefile])
AC_CONFIG_FILES([galpy-interface/Makefile])
AC_CONFIG_FILES([src/get_version.hpp])

AC_OUTPUT
