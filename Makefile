#PS_PATH = -I../../fdps/src/
PS_PATH = -I../fdps/src/
#PS_PATH = -I../../../fdps/src/
#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../../../../project/fdps/src
#PS_PATH = -I../FDPS/src
ARC_PATH = -I../ARC/include

#use_k_computer = yes
#use_xc30_naoj = yes
use_x86 = yes

ifeq ($(use_k_computer),yes)
CXX = time mpiFCCpx
#CXXFLAGS = -Kfast
CXXFLAGS += -x32
CXXFLAGS += -Xg
CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
CXXFLAGS += -DFAST_WALK_K
CXXFLAGS += -std=c++11
#CXXFLAGS += -NRtrap
CXXFLAGS += -Nfjcex
CXXFLAGS += -Krestp=all
CXXFLAGS += -DINTRINSIC_K
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
CXXFLAGS = -O3
CXXFLAGS += -Wall
CXXFLAGS += -march=core-avx2
CXXFLAGS += -ffast-math -funroll-loops
CXXFLAGS += -std=c++11
CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
#CXXFLAGS += -DINTRINSIC_X86
#CXXFLAGS += -DUSE_GNU_PARALLEL_SORT
endif

CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DDIV_FIX -DP3T_64BIT
#CXXFLAGS += -DUSE_QUAD 
CXXFLAGS += -DARC_ERROR
CXXFLAGS += -D ARC_WARN
CXXFLAGS += -D SAFETY_CHECK
#CXXFLAGS += -D DEBUG
#CXXFLAGS += -D DEBUG_TEMP
#CXXFLAGS += -D DEBUG_OUTPUT

VPATH=./src ./test

SRC = main.cc hard.hpp soft.hpp hard_force.hpp io.hpp kepler.hpp phantomquad_for_p3t_x86.hpp domain.hpp profile.hpp cluster_list.hpp integrate.hpp rsearch.hpp

all: nbody.out

nbody.out: $(SRC)
	$(CXX) $(PS_PATH) $(ARC_PATH) $(CXXFLAGS) -o $@ $<

ARC_debug.out: chain_debug.cxx force.hpp
	$(CXX) $(PS_PATH) $(ARC_PATH) -I./src $(CXXFLAGS) -D DEBUG -o $@ $<

rsearchtest: rsearchtest.cxx
	$(CXX) $(PS_PATH) $(ARC_PATH) -I./src $(CXXFLAGS) -o $@ $<

clean:
	rm *.out *.o
cleanall:
	rm *.out *.hpp~ *.cc~ *.h~

run: nbody.out
	mpiexec -n 2 ./nbody.out -i input.dat.14
