#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../fdps/src/
#PS_PATH = -I../../../fdps/src/
#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../../../../project/fdps/src
#PS_PATH = -I../FDPS/src
PS_PATH  = -I/home/lwang/code/fdps/src
ARC_PATH = -I/home/lwang/GitHub/ARC/include
INCLUDE  = -I./src -I../src

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
CXXFLAGS = -g
CXXFLAGS += -Wall
CXXFLAGS += -march=core-avx2
CXXFLAGS += -ffast-math -funroll-loops
CXXFLAGS += -std=c++11
CXXFLAGS += ${shell gsl-config --cflags}
CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
#CXXFLAGS += -DINTRINSIC_X86
#CXXFLAGS += -DUSE_GNU_PARALLEL_SORT
CXXLIBS += ${shell gsl-config --libs}
endif

CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DDIV_FIX -DP3T_64BIT
#CXXFLAGS += -DUSE_QUAD 
CXXFLAGS += -D ARC_DEBUG
CXXFLAGS += -D ARC_ERROR
CXXFLAGS += -D ARC_WARN
CXXFLAGS += -D SAFETY_CHECK
CXXFLAGS += -D HARD_DEBUG
CXXFLAGS += -D HARD_DEBUG_PRINT
CXXFLAGS += -D FIX_STEP_DEBUG
#CXXFLAGS += -D DEBUG
#CXXFLAGS += -D DEBUG_TEMP
#CXXFLAGS += -D DEBUG_OUTPUT

VPATH=./src ./test

SRC = main.cc hard.hpp soft.hpp hard_force.hpp io.hpp kepler.hpp phantomquad_for_p3t_x86.hpp domain.hpp profile.hpp cluster_list.hpp integrate.hpp rsearch.hpp

all: nbody.out

nbody.out: $(SRC)
	$(CXX) $(PS_PATH) $(ARC_PATH) $(CXXFLAGS) -o $@ $< $(CXXLIBS)

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
