#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../fdps/
#PS_PATH = -I../../../fdps/src/
#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../../../../project/fdps/src
PS_PATH = -I../FDPS/src
ARC_PATH = -I../ARC/include

#use_k_computer = yes
#use_xc30_naoj = yes
use_x86 = yes

ifeq ($(use_k_computer),yes)
CXX = time mpiFCCpx
#CXXFLAGS = -Kfast
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
endif

ifeq ($(use_xc30_naoj),yes)
CXX = time CC
CXXFLAGS = -O3 
CXXFLAGS += -Wall
CXXFLAGS += -march=core-avx2
CXXFLAGS += -ffast-math -funroll-loops
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
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
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -std=c++11
CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
#CXXFLAGS += -D ARC_WARN
CXXFLAGS += -DREAD_FILE
#CXXFLAGS += -D DEBUG
#CXXFLAGS += -DINTRINSIC_X86
#CXXFLAGS += -DUSE_GNU_PARALLEL_SORT
endif

VPATH=./src

SRC = main.cc hard.hpp class.hpp force.hpp io.hpp kepler.hpp phantomquad_for_p3t_x86.hpp domain.hpp profile.hpp cluster_list.hpp

all:nbody.out

nbody.out:$(SRC)
	$(CXX) -DUSE_QUAD -DDIV_FIX -DP3T_64BIT $(PS_PATH) $(ARC_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm *.out *.o
cleanall:
	rm *.out *.hpp~ *.cc~ *.h~
