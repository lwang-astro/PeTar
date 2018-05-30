##FDPS PATH
#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../fdps/src/
#PS_PATH = -I../../../fdps/src/
#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../../../../project/fdps/src
#PS_PATH = -I../FDPS/src
#PS_PATH = -I../fdps.svn/
PS_PATH  = -I/home/lwang/code/fdps/src

##ARC PATH
#ARC_PATH= -I../TSARC/include
#ARC_PATH = -I/home/lwang/code/ARC/include
ARC_PATH = -I/home/lwang/GitHub/ARC/include

##Gperftools PATH
#GPERF_PATH = -L../../soft/gperftools-2.6.90/lib

#ROOT_PATH= ${shell pwd -P}
INCLUDE  = -I./src -I../src

#use_k_computer = yes
#use_xc30_naoj = yes
use_x86 = yes

ifeq ($(use_k_computer),yes)
CXX = time mpiFCCpx
CXXFLAGS = -Kfast
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -x32
CXXFLAGS += -Xg
CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
CXXFLAGS += -DFAST_WALK_K
CXXFLAGS += -std=c++11
#CXXFLAGS += -NRtrap
CXXFLAGS += -Nfjcex
CXXFLAGS += -Krestp=all
CXXFLAGS += -DINTRINSIC_K
# profiling
CXXFLAGS += -Ntl_trt
endif

ifeq ($(use_xc30_naoj),yes)
CXX = time CC
CXXFLAGS = -O3 
CXXFLAGS += -Wall
CXXFLAGS += -march=core-avx2
CXXFLAGS += -ffast-math -funroll-loops
CXXFLAGS += -std=c++11
CXXFLAGS += -DINTRINSIC_X86
#CXXFLAGS += -DUSE_GNU_PARALLEL_SORT
endif

ifeq ($(use_x86),yes)
CXX = time g++
#CXX = time icc

#CXX = kinst-ompp mpicxx

#CXX = tau_cxx.sh  -tau_makefile=/opt/tau-2.26.3/x86_64/lib/Makefile.tau-mpi-openmp -tau_options=-optCompInst

#CXX = time mpicxx
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DMPICH_IGNORE_CXX_SEEKC

CXXFLAGS = -g -O0 -fbounds-check
#CXXFLAGS += -O2

CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_DEBUG_PRINT
CXXFLAGS += -Wall
#CXXFLAGS += -march=skylake-avx512
CXXFLAGS += -march=core-avx2
#CXXFLAGS += -mavx
CXXFLAGS += -ffast-math -funroll-loops
CXXFLAGS += -std=c++11
#CXXFLAGS += ${shell gsl-config --cflags}
CXXFLAGS += -DINTRINSIC_X86
#CXXFLAGS += -DUSE_GNU_PARALLEL_SORT
#CXXLIBS += ${shell gsl-config --libs}
endif

CXXFLAGS += -DDIV_FIX
#CXXFLAGS += -DP3T_64BIT
CXXFLAGS += -DUSE_QUAD
CXXFLAGS += -DUSE_SIMD

CXXFLAGS += -D HARD_CM_KICK
CXXFLAGS += -D TIDAL_TENSOR # Must use HARD_CM_KICK together
CXXFLAGS += -D SOFT_PERT
CXXFLAGS += -D SPLIT_MASS
CXXFLAGS += -D PROFILE
CXXFLAGS += -D ARC_SYM
CXXFLAGS += -D ARC_OPT_SYM2
#CXXFLAGS += -D ARC_SYM_SD_PERIOD

SIMD_DEBFLAGS += -DCALC_EP_64bit
SIMD_DEBFLAGS += -DRSQRT_NR_EPJ_X2
SIMD_DEBFLAGS += -DRSQRT_NR_EPJ_X4
SIMD_DEBFLAGS += -DCALC_SP_64bit
SIMD_DEBFLAGS += -DRSQRT_NR_SPJ_X2
SIMD_DEBFLAGS += -DRSQRT_NR_SPJ_X4
SIMD_DEBFLAGS += -DAVX_PRELOAD

#DEBFLAGS += -D ARC_PROFILE
#DEBFLAGS += -D INTEGRATED_CUTOFF_FUNCTION
#DEBFLAGS += -D ARC_DEBUG
#DEBFLAGS += -D ARC_DEBUG_PRINT
#DEBFLAGS += -D ARC_DEEP_DEBUG
#DEBFLAGS += -D ARC_ERROR
DEBFLAGS += -D ARC_DEBUG_DUMP
#DEBFLAGS += -D ARC_WARN
#DEBFLAGS += -D HARD_DEBUG
#DEBFLAGS += -D HARD_CHECK_ENERGY
DEBFLAGS += -D HARD_DEBUG_DUMP
#DEBFLAGS += -D HARD_DEBUG_PRINT
#DEBFLAGS += -D HARD_DEBUG_PROFILE
#DEBFLAGS += -D NAN_CHECK_DEBUG
#DEBFLAGS += -D DATA_DEBUG
#DEBFLAGS += -D FIX_STEP_DEBUG
#DEBFLAGS += -D DEBUG
#DEBFLAGS += -D DEBUG_TEMP
#DEBFLAGS += -D MAIN_DEBUG

ARC_DEBFLAGS += -D ARC_PROFILE -D ARC_DEBUG -D ARC_DEBUG_PRINT -D ARC_DEEP_DEBUG -D ARC_EEROR -D ARC_WARN
HARD_DEBFLAGS+= -D ARC_PROFILE -D ARC_DEBUG -D ARC_ERROR -D ARC_WARN -D HARD_DEBUG -D HARD_CHECK_ENERGY -D HARD_DEBUG_PRINT

ifneq (x$(GPERF_PATH),x)
CXXLIBS += $(GPERF_PATH) -lprofiler -ltcmalloc 
endif

VPATH=./src ./test ../src

SRC = main.cc hard.hpp soft.hpp hard_force.hpp io.hpp kepler.hpp phantomquad_for_p3t_x86.hpp domain.hpp profile.hpp cluster_list.hpp integrate.hpp init.hpp

all: nbody.out

nbody.out: $(SRC)
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(DEBFLAGS) -o $@ $< $(CXXLIBS)

getdata: data.cc
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(DEBFLAGS) -o $@ $< $(CXXLIBS)

ARC_debug.out: chain_debug.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(ARC_DEBFLAGS) -D DEBUG -o $@ $< $(CXXLIBS)

hard_debug.out: hard_debug.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(HARD_DEBFLAGS) -o $@ $< $(CXXLIBS)

group_debug.out: group_debug.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(HARD_DEBFLAGS) -D STABLE_CHECK_DEBUG -o $@ $< $(CXXLIBS)

rsearchtest: rsearchtest.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

test.out: test.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

searchgrouptest: searchgroup.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(HARD_DEBFLAGS) -o $@ $< $(CXXLIBS)

hermitetest: hermite.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(DEBFLAGS) -o $@ $< $(CXXLIBS)

splinetest: spline.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(DEBFLAGS) -o $@ $< $(CXXLIBS)

keplersolvertest: keplersolver.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(DEBFLAGS) -o $@ $< $(CXXLIBS)

keplertest: keplertest.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(DEBFLAGS) -D STABLE_CHECK_DEBUG -o $@ $< $(CXXLIBS)

hardtest: hard.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(HARD_DEBFLAGS) -o $@ $< $(CXXLIBS)

arctest: arc.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(ARC_DEBFLAGS) -o $@ $< $(CXXLIBS)

simd_test.s: simd_test.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(SIMD_DEBFLAGS) -S $< -o $@  $(CXXLIBS)

simd_test: simd_test.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) $(SIMD_DEBFLAGS)  $< -o $@  $(CXXLIBS)

clean:
	rm *.out *.o
cleanall:
	rm *.out *.hpp~ *.cc~ *.h~

run: nbody.out
	mpiexec -n 2 ./nbody.out -i input.dat.14
