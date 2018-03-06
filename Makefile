##FDPS PATH
#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../fdps/src/
#PS_PATH = -I../../../fdps/src/
#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../../../../project/fdps/src
PS_PATH = -I../FDPS/src
#PS_PATH  = -I/home/lwang/code/fdps/src

##ARC PATH
ARC_PATH= -I../TSARC/include
#ARC_PATH = -I/home/lwang/code/ARC/include
#ARC_PATH = -I/home/lwang/GitHub/ARC/include

##Gperftools PATH
#GPERF_PATH = -L/opt/gperftools-2.6.1/lib

#ROOT_PATH= ${shell pwd -P}
INCLUDE  = -I./src -I../src

use_k_computer = yes
#use_xc30_naoj = yes
#use_x86 = yes

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
CXX = time mpicxx
#CXX = kinst-ompp mpicxx
#CXX = tau_cxx.sh  -tau_makefile=/opt/tau-2.26.3/x86_64/lib/Makefile.tau-mpi-openmp -tau_options=-optCompInst 
#CXXFLAGS = -g -O0 -fbounds-check
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -O2
CXXFLAGS += -Wall
CXXFLAGS += -march=core-avx2
CXXFLAGS += -ffast-math -funroll-loops
CXXFLAGS += -std=c++11
#CXXFLAGS += ${shell gsl-config --cflags}
CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
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
#CXXFLAGS += -D ARC_PROFILE

#CXXFLAGS += -D INTEGRATED_CUTOFF_FUNCTION
#CXXFLAGS += -D ARC_DEBUG
#CXXFLAGS += -D ARC_DEBUG_PRINT
#CXXFLAGS += -D ARC_ERROR
#CXXFLAGS += -D ARC_WARN
#CXXFLAGS += -D HARD_DEBUG
#CXXFLAGS += -D HARD_DEBUG_ENERGY
#CXXFLAGS += -D HARD_DEBUG_PRINT
#CXXFLAGS += -D HARD_DEBUG_PROFILE
#CXXFLAGS += -D NAN_CHECK_DEBUG
#CXXFLAGS += -D DATA_DEBUG
#CXXFLAGS += -D FIX_STEP_DEBUG
#CXXFLAGS += -D DEBUG
#CXXFLAGS += -D DEBUG_TEMP
#CXXFLAGS += -D MAIN_DEBUG

ifneq (x$(GPERF_PATH),x)
CXXLIBS += $(GPERF_PATH) -lprofiler -ltcmalloc 
endif

VPATH=./src ./test ../src

SRC = main.cc hard.hpp soft.hpp hard_force.hpp io.hpp kepler.hpp phantomquad_for_p3t_x86.hpp domain.hpp profile.hpp cluster_list.hpp integrate.hpp init.hpp

all: nbody.out

nbody.out: $(SRC)
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

getdata: data.cc
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

ARC_debug.out: chain_debug.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -D DEBUG -o $@ $< $(CXXLIBS)

rsearchtest: rsearchtest.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

test.out: test.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

searchgrouptest: searchgroup.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

hermitetest: hermite.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

splinetest: spline.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

keplersolvertest: keplersolver.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

keplertest: keplertest.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -D STABLE_CHECK_DEBUG -o $@ $< $(CXXLIBS)

hardtest: hard.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

arctest: arc.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) $(INCLUDE) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

clean:
	rm *.out *.o
cleanall:
	rm *.out *.hpp~ *.cc~ *.h~

run: nbody.out
	mpiexec -n 2 ./nbody.out -i input.dat.14
