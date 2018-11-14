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
ARC_PATH = -I/home/lwang/code/ARC/include
#ARC_PATH = -I/home/lwang/GitHub/ARC/include

##Gperftools PATH
#GPERF_PATH = ../../soft/gperftools-2.6.90

#ROOT_PATH= ${shell pwd -P}
INCLUDE  = -I./src -I../src

#use_k_computer = yes
#use_xc30_naoj = yes
use_x86 = yes
use_mpi = yes
#debug_mode=yes
#use_intel=yes

ifeq ($(use_k_computer),yes)
CXX = time mpiFCCpx
OPTFLAGS = -Kfast
OPTFLAGS += -Nfjcex
OPTFLAGS += -Krestp=all
OPTFLAGS += -Ntl_trt
#OPTFLAGS += -NRtrap

CXXFLAGS = -Kopenmp
CXXFLAGS += -x32
CXXFLAGS += -Xg
CXXFLAGS += -std=c++11
CXXFLAGS += -DINTRINSIC_K

FDPSFLAGS = -DPARTICLE_SIMULATOR_THREAD_PARALLEL
FDPSFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
FDPSFLAGS += -DFAST_ALL_TO_ALL_FOR_K
FDPSFLAGS += -DFAST_WALK_K
# profiling
endif

ifeq ($(use_xc30_naoj),yes)
CXX = time CC
CXXNOMPI = time CC
OPTFLAGS = -O3 -Wall

CXXFLAGS = -ffast-math -funroll-loops
CXXFLAGS += -march=core-avx2
CXXFLAGS += -std=c++11
CXXFLAGS += -DINTRINSIC_X86

#FDPSFLAGS += -DUSE_GNU_PARALLEL_SORT
endif

ifeq ($(use_x86),yes)
ifeq ($(use_mpi),yes)
CXX = time mpicxx
ifeq ($(use_intel),yes)
CXXNOMPI = time icc
else
CXXNOMPI = time g++
endif
#CXX = kinst-ompp mpicxx
#CXX = tau_cxx.sh  -tau_makefile=/opt/tau-2.26.3/x86_64/lib/Makefile.tau-mpi-openmp -tau_options=-optCompInst
FDPSFLAGS = -DPARTICLE_SIMULATOR_MPI_PARALLEL
FDPSFLAGS += -DMPICH_IGNORE_CXX_SEEKC
else
ifeq ($(use_intel),yes)
CXX = time icc
CXXNOMPI = time icc
else
CXX = time g++
CXXNOMPI = time g++
endif
endif

ifeq ($(debug_mode),yes)
OPTFLAGS = -g -O0 -fbounds-check -Wall -D SANITY_CHECK_REALLOCATABLE_ARRAY -D HARD_DEBUG_PRE_DUMP
else
OPTFLAGS = -O2 -Wall 
#OPTFLAGS += -ffast-math -funroll-loops
endif

CXXFLAGS = -std=c++11
CXXFLAGS += -Wall
CXXFLAGS += -fopenmp
#CXXFLAGS += -march=skylake-avx512
CXXFLAGS += -march=core-avx2
#CXXFLAGS += -mavx
CXXFLAGS += -DINTRINSIC_X86

FDPSFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#FDPSFLAGS += -DPARTICLE_SIMULATOR_DEBUG_PRINT
#FDPSFLAGS += -DUSE_GNU_PARALLEL_SORT
endif

CXXFLAGS += -DDIV_FIX
#CXXFLAGS += -DP3T_64BIT
CXXFLAGS += -DUSE_QUAD
CXXFLAGS += -DUSE_SIMD
CXXFLAGS += -D PROFILE
CXXFLAGS += -D HARD_CHECK_ENERGY
#CXXFLAGS += ${shell gsl-config --cflags}
#CXXLIBS += ${shell gsl-config --libs}

#MT_FLAGS += -D HARD_CM_KICK
MT_FLAGS += -D TIDAL_TENSOR # Must use HARD_CM_KICK together
MT_FLAGS += -D SOFT_PERT
MT_FLAGS += -D SPLIT_MASS
MT_FLAGS += -D ARC_SYM
MT_FLAGS += -D ARC_OPT_SYM2
#MT_FLAGS += -D ARC_SYM_SD_PERIOD

SIMD_DEBFLAGS += -DCALC_EP_64bit
SIMD_DEBFLAGS += -DRSQRT_NR_EPJ_X2
SIMD_DEBFLAGS += -DRSQRT_NR_EPJ_X4
SIMD_DEBFLAGS += -DCALC_SP_64bit
SIMD_DEBFLAGS += -DRSQRT_NR_SPJ_X2
SIMD_DEBFLAGS += -DRSQRT_NR_SPJ_X4
SIMD_DEBFLAGS += -DAVX_PRELOAD

#DEBFLAGS += -D ARC_PROFILE
#DEBFLAGS += -D INTEGRATED_CUTOFF_FUNCTION
DEBFLAGS += -D ARC_DEBUG
#DEBFLAGS += -D ARC_DEBUG_PRINT
#DEBFLAGS += -D ARC_DEEP_DEBUG
DEBFLAGS += -D ARC_ERROR
DEBFLAGS += -D ARC_DEBUG_DUMP
#DEBFLAGS += -D ARC_WARN
DEBFLAGS += -D HARD_DEBUG
DEBFLAGS += -D HARD_DEBUG_PRE_DUMP
DEBFLAGS += -D HARD_DEBUG_DUMP
#DEBFLAGS += -D STABLE_CHECK_DEBUG
#DEBFLAGS += -D CLUSTER_DEBUG
#DEBFLAGS += -D HARD_DEBUG_PRINT
#DEBFLAGS += -D HARD_DEBUG_PROFILE
#DEBFLAGS += -D NAN_CHECK_DEBUG
#DEBFLAGS += -D DATA_DEBUG
#DEBFLAGS += -D FIX_STEP_DEBUG
#DEBFLAGS += -D DEBUG
#DEBFLAGS += -D DEBUG_TEMP
#DEBFLAGS += -D MAIN_DEBUG
#DEBFLAGS += -D CORRECT_FORCE_DEBUG

DEBUG_OPT_FLAGS = -g -O0 -fbounds-check -Wall  -D SANITY_CHECK_REALLOCATABLE_ARRAY
ARC_DEBFLAGS += -D ARC_PROFILE -D ARC_DEBUG -D ARC_DEBUG_PRINT -D ARC_EEROR -D ARC_WARN 
ARC_MT_FLAGS += -D ARC_SYM -D ARC_OPT_SYM2 -D TIDAL_TENSOR -D SPLIT_MASS
HARD_DEBFLAGS+= -D ARC_PROFILE -D ARC_DEBUG -D ARC_ERROR -D ARC_WARN -D HARD_DEBUG -D HARD_CHECK_ENERGY -D HARD_DEBUG_PRINT -D ADJUST_GROUP_DEBUG
HARD_MT_FLAGS += -D ARC_SYM -D ARC_OPT_SYM2 -D TIDAL_TENSOR -D SPLIT_MASS

ifneq (x$(GPERF_PATH),x)
CXXLIBS += -L$(GPERF_PATH)/lib -lprofiler -ltcmalloc 
endif

VPATH=./src ./test ../src

SRC = main.cc hard.hpp soft.hpp hard_force.hpp io.hpp kepler.hpp phantomquad_for_p3t_x86.hpp domain.hpp profile.hpp cluster_list.hpp integrate.hpp init.hpp

all: nbody.out

nbody.out: $(SRC)
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(OPTFLAGS) $(CXXFLAGS) $(FDPSFLAGS) $(MT_FLAGS) $(DEBFLAGS) -o $@ $< $(CXXLIBS)

#getdata: data.cc
#	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(OPTFLAGS) $(CXXFLAGS) $(DEBFLAGS) -o $@ $< $(CXXLIBS)

arc_debug.out: arc_debug.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(MT_FLAGS) $(ARC_DEBFLAGS) -D ARC_DEEP_DEBUG -D DEBUG -o $@ $< $(CXXLIBS)

hard_debug.out: hard_debug.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(MT_FLAGS) $(HARD_DEBFLAGS) -o $@ $< $(CXXLIBS)

stable_debug.out: stable_debug.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(MT_FLAGS) $(HARD_DEBFLAGS) -D STABLE_CHECK_DEBUG -o $@ $< $(CXXLIBS)

adjustgrouptest: adjustgroup_test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(MT_FLAGS) $(HARD_DEBFLAGS) -o $@ $< $(CXXLIBS)

#rsearchtest: rsearch_test.cxx
#	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

test.out: test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

searchgrouptest: searchgroup_test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(MT_FLAGS) $(HARD_DEBFLAGS) -o $@ $< $(CXXLIBS)

hermitetest: hermite_test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(HARD_MT_FLAGS) $(HARD_DEBFLAGS) -o $@ $< $(CXXLIBS)

splinetest: spline_test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(DEBFLAGS) -o $@ $< $(CXXLIBS)

keplersolvertest: keplersolver.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(DEBFLAGS) -o $@ $< $(CXXLIBS)

stabletest: stable_test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(DEBFLAGS) -D STABLE_CHECK_DEBUG -o $@ $< $(CXXLIBS)

hardtest: hard_test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(HARD_MT_FLAGS) $(HARD_DEBFLAGS) -o $@ $< $(CXXLIBS)

arctest: arc_test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(ARC_MT_FLAGS) $(ARC_DEBFLAGS) -o $@ $< $(CXXLIBS)

simdtest.s: simd_test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(OPTFLAGS) $(CXXFLAGS) $(SIMD_DEBFLAGS) -S $< -o $@  $(CXXLIBS)

simdtest: simd_test.cxx
	$(CXXNOMPI) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(OPTFLAGS) $(CXXFLAGS) $(SIMD_DEBFLAGS)  $< -o $@  $(CXXLIBS)

clean:
	rm *.out *.o
cleanall:
	rm *.out *.hpp~ *.cc~ *.h~

run: nbody.out
	mpiexec -n 2 ./nbody.out -i input.dat.14
