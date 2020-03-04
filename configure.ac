
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([PeTar], [1.0], [longwang.astro@live.com])
m4_include([autoconf_m4/ax_check_compile_flag.m4])
m4_include([autoconf_m4/ax_prog_cxx_mpi.m4])
m4_include([autoconf_m4/ax_openmp.m4])
m4_include([autoconf_m4/ax_gcc_x86_cpu_supports.m4])
#AC_CONFIG_SRCDIR([src/main.cc])
#AC_CONFIG_HEADERS([config.h])

PROG_NAME="petar"

AC_LANG([C++])

AC_PROG_AWK

# Checks for programs.
# MPI
  # If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
  # If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
  # If --with-mpi=no is used, use a standard C compiler instead.

AC_ARG_WITH([mpi],
            [AS_HELP_STRING([--with-mpi],
                            [compile with MPI (parallelization) support. If none is found, MPI is not used. Default: auto])],
            [],
            [with_mpi=auto])


AX_PROG_CXX_MPI([test x"$with_mpi" != xno],
                [use_mpi=yes
                 PROG_NAME=$PROG_NAME".mpi"
                 CXXNOMPI=`$CXX -show|awk '{print $1}'`],
                [use_mpi=no
                 CXXNOMPI=$CXX])

AS_IF([test x"$use_mpi" == xno],
      [AS_IF([test x"$with_mpi" = xyes],
             [AC_MSG_FAILURE([MPI compiler requested, but could not use MPI.])],
             [AC_MSG_WARN([No MPI compiler found, disable MPI.])])])

# debug
AC_ARG_WITH([debug],
            [AS_HELP_STRING([--with-debug],
                            [Switch on debugging mode. Options: g: enable debugger support; assert: enable assertion; no: none. Default: no])],
            [],
            [with_debug=no])


AS_IF([test x"$with_debug" == xg],
      [OPTFLAGS=" -g -O0"
       PROG_SUFFIX=".g"
       AX_CHECK_COMPILE_FLAG([-fbounds-check],
                             [OPTFLAGS=$OPTFLAGS" -fbounds-check"])],
      [AS_IF([test x"$with_debug" != xno],
             [PROG_SUFFIX=".d"])
       AX_CHECK_COMPILE_FLAG([-O3],
                             [OPTFLAGS=" -O3"])])


AX_CHECK_COMPILE_FLAG([-Wall],[OPTFLAGS=$OPTFLAGS" -Wall"])

AX_CHECK_COMPILE_FLAG([-std=c++11],
                      [OPTFLAGS=$OPTFLAGS" -std=c++11"],
                      [AC_MSG_FAILURE([c++11 support is required for the compiler])])

# OpenMP
AC_ARG_WITH([omp],
	    [AS_HELP_STRING([--with-omp],
                            [compile with OpenMP support. Default: yes])],
            [],
	    [with_omp=yes])

AS_IF([test x"$with_omp" != xno],
      [AX_OPENMP([use_omp=yes
                  PROG_NAME=$PROG_NAME'.omp'],
                 [use_omp=no])],
      [use_omp=no])

# x86
AC_ARG_WITH([x86],
	    [AS_HELP_STRING([--with-x86],
                            [compile with x86 feature support. Default: yes])],
            [],
	    [with_x86=yes])

AC_ARG_WITH([simd],
            [AS_HELP_STRING([--with-simd],
                [compile with x86 SIMD support for long-distant tree force. Options: avx, avx2, avx512dq. Default: auto])],
            [],
            [with_simd=auto])

AS_IF([test x"$with_x86" != xno],
      [use_x86=yes
       use_simd=no
       AS_IF([test x"$with_simd" != xauto],
             [AX_CHECK_COMPILE_FLAG([-m$with_simd],
                                    [SIMDFLAGS=" -m$with_simd"
		                     SIMD_TYPE=$with_simd
				     use_simd=yes],
				    [AC_MSG_FAILURE([$with_simd is not support by the compiler])])],
             [m4_foreach_w([ax_x86_feature],
                           [avx avx2 avx512dq],
                           [AX_GCC_X86_CPU_SUPPORTS(ax_x86_feature,
                                                    [SIMDFLAGS=" -m[]ax_x86_feature"
                                                     SIMD_TYPE="[]ax_x86_feature"
                                                     use_simd=yes],
						     [])])])
       PROG_NAME=$PROG_NAME"."$SIMD_TYPE],
      [use_x86=no])

# cuda
AC_ARG_WITH([cuda],
            [AS_HELP_STRING([--with-cuda],
                [compile with CUDA (GPU) acceleration support for long-distant tree force, Default: no])],
            [],
            [with_cuda=no])

AC_ARG_WITH([cuda-prefix], 
            [AS_HELP_STRING([--with-cuda-prefix],
	                    [Prefix of your CUDA installation])],
            [],
	    [cuda_prefix="/usr/local/cuda"])

AC_ARG_WITH([cuda-sdk-prefix], 
            [AS_HELP_STRING([--with-cuda-sdk-prefix],
	                    [Prefix of your CUDA samples (SDK) installation])],
            [],
	    [cuda_sdk_prefix="/usr/local/cuda/samples"])

AS_IF([test "x$with_cuda" != xno],
      [AC_CHECK_PROG([NVCC], [nvcc], [nvcc], [none])
       use_cuda=yes
       AS_IF([test "x$NVCC" == xnone], 
             [AC_PATH_PROG([NVCC], [nvcc], [none], [$cuda_prefix/bin])
              AS_IF([test "x$NVCC" == xnone], 
                    [AC_MSG_FAILURE([Cannot find CUDA compiler nvcc, please check whether nvcc is in environment PATH or use --with-cuda-prefix to provide the PATH of CUDA installation])])],
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
                                    [AC_MSG_FAILURE([Cannot find CUDA sample inc file helper_cuda.h, please provide correct cuda SDK prefix by using --with-cuda-sdk-prefix])])])

       AC_CHECK_LIB([cudart],
                    [main], 
                    [CUDALIBS=' -lcudart'],
                    [AC_CHECK_FILE(["$cuda_prefix/lib64/libcudart.so"],
                                   [CUDALIBS=" -L$cuda_prefix/lib64 -lcudart"],
                                   [AC_MSG_FAILURE([Cannot find CUDA library -lcudart, please provide correct cuda PREFIX by using --with-cuda-prefix])])])
       AC_CHECK_LIB([gomp], 
                    [main],
		    [CUDALIBS=$CUDALIBS' -lgomp'],
		    [AC_CHECK_FILE(["$cuda_prefix/lib64/libgomp.so"],
                                   [CUDALIBS=$CUDALIBS" -lgomp"],
                                   [AC_MSG_FAILURE([Cannot find CUDA library -lgomp, please provide correct cuda PREFIX by using --with-cuda-prefix])])])],
       [use_cuda=no])


PROG_NAME=$PROG_NAME$PROG_SUFFIX

# gperftools
AC_ARG_WITH([gperf],
            [AS_HELP_STRING([--with-gperf],
                [compile with gperftools for profiling, Default: no])],
            [],
            [with_gperf=no])

AC_ARG_WITH([gperf-prefix], 
            [AS_HELP_STRING([--with-gperf-prefix],
	                    [Prefix of your gperftools installation])],
            [],
	    [gperf_prefix=/opt/gperftools])

AS_IF([test x"$with_gperf" != xno],
      [use_gperf=yes
       AC_CHECK_LIB([profiler],
                    [main],
		    [GPERFLIBS=' -lprofiler']
		    [AC_CHECK_FILE(["$gperf_prefix/lib/libprofiler.so"],
                                   [GPREFFLAGS=" -I$gperf_prefix/include"
                                    GPERFLIBS=" -L$gperf_prefix/lib -lprofiler"],
				   [AC_MSG_FAILURE([Cannot find gperftools library -lprofiler, please provide correct gperftools PREFIX by using --with-gperf-prefix])])])
       AC_CHECK_LIB([tcmalloc],
                    [main],
		    [GPERFLIBS=' -ltcmalloc']
		    [AC_CHECK_FILE(["$gperf_prefix/lib/libtcmalloc.so"],
                                   [GPERFLIBS=GPERFLIBS" -ltcmalloc"],
				   [AC_MSG_FAILURE([Cannot find gperftools library -ltcmalloc, please provide correct gperftools PREFIX by using --with-gperf-prefix])])])],
       [use_gperf=no])

# FDPS
AC_ARG_WITH([fdps-prefix], 
            [AS_HELP_STRING([--with-fdps-prefix],
	                    [Prefix of your FDPS installation. Default: ../FDPS])],
            [],
	    [fdps_prefix="../FDPS"])

AC_CHECK_FILE(["$fdps_prefix/src/particle_simulator.hpp"],
	      [INCLUDE=$INCLUDE" -I$fdps_prefix/src"],
	      [AC_MSG_FAILURE([Cannot find FDPS library, please provide correct PREFIX by using --with-fdps-prefix])])

# SDAR
AC_ARG_WITH([sdar-prefix], 
            [AS_HELP_STRING([--with-sdar-prefix],
	                    [Prefix of your FDPS installation. Default: ../SDAR])],
            [],
	    [sdar_prefix="../SDAR"])

AC_CHECK_FILE(["$sdar_prefix/src/Hermite/hermite_integrator.h"],
	      [INCLUDE=$INCLUDE" -I$sdar_prefix/src"],
	      [AC_MSG_FAILURE([Cannot find SDAR library, please provide correct PREFIX by using --with-sdar-prefix])])

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

AC_SUBST([use_x86])
AC_SUBST([use_simd])
AC_SUBST([use_mpi])
AC_SUBST([use_cuda])
AC_SUBST([use_omp])
AC_SUBST([use_gperf])
AC_SUBST([with_debug])
AC_SUBST([OPTFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])
AC_SUBST([NVCC])
AC_SUBST([CUDAFLAGS])
AC_SUBST([CUDALIBS])
AC_SUBST([SIMDFLAGS])
AC_SUBST([GPERFFLAGS])
AC_SUBST([GPERFLIBS])
AC_SUBST([INCLUDE])
AC_SUBST([PROG_NAME])
AC_SUBST([CXX])
AC_SUBST([CXXNOMPI])

AC_MSG_NOTICE([----------------------Configure Summary--------------------])
AC_MSG_NOTICE([--Features:])
AC_MSG_NOTICE([     Using MPI:      $use_mpi])
AC_MSG_NOTICE([     Using GPU:      $use_cuda])
AC_MSG_NOTICE([     Using x86 :     $use_x86])
AS_IF([test "x$use_simd" != xno],
      [AC_MSG_NOTICE([     SIMD type:      $SIMD_TYPE ])
       AS_IF([test "x$with_simd" == xauto],
             [AC_MSG_NOTICE([  Notice: this is auto-detected based on the local CPU architecture])
              AC_MSG_NOTICE([          If different CPU is used for running, check whether $SIMD_TYPE is also supported])])])
AC_MSG_NOTICE([     Using OpenMP:   $use_omp])
AC_MSG_NOTICE([     Debug mode:     $with_debug])
AC_MSG_NOTICE([--Compilers:])
AC_MSG_NOTICE([     C++ compiler:      $CXX])
AS_IF([test "x$use_cuda" != xno],
      [AC_MSG_NOTICE([     CUDA compiler:     $NVCC])
       AC_MSG_NOTICE([     CUDA version:   $CUDA_VERSION])])
AC_MSG_NOTICE([--Targets:])
AC_MSG_NOTICE([     Program name:   $PROG_NAME])
AS_IF([test x"$prefix" != xNONE],
      [AC_MSG_NOTICE([     Install path:   $prefix])
       AC_MSG_NOTICE([-----------------------------------------------------------])
       AC_MSG_NOTICE([  After install, please])
       AC_MSG_NOTICE([      add $prefix/include into your PTHYON include path])
       AC_MSG_NOTICE([      add $prefix/bin into your PATH])],
      [AC_MSG_WARN([Install path is not defined, please provide it by using --prefix=INSTALL_PATH])])
AC_MSG_NOTICE([-----------------------------------------------------------])
AC_MSG_NOTICE([Reference paper: preparing...])
AC_MSG_NOTICE([Questions and bug reports: https://github.com/lwang-astro/PeTar/issues])
AC_MSG_NOTICE([-----------------------------------------------------------])


AC_CONFIG_FILES([Makefile])

AC_OUTPUT