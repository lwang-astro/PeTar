# standard amuse configuration include
# config.mk will be made after ./configure has run
# don't run make directly in the PeTar/amuse-interface directory, this is used inside amuse.
AMUSE_DIR?=$(shell amusifier --get-amuse-dir)
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

INCLUDE  += -I src/PeTar/src -I src/SDAR/src -I src/FDPS/src 

MPIFLAGS +=  -DPARTICLE_SIMULATOR_MPI_PARALLEL -DMPICH_IGNORE_CXX_SEEKC

CXXFLAGS = ${INCLUDE} -fPIC -Wall -std=c++11  -fopenmp
CXXFLAGS += -O2 
#CXXFLAGS += -g -O0

# arch
CXXFLAGS += -march=core-avx2
CXXFLAGS += -D INTRINSIC_X86
CXXFLAGS += -D USE_SIMD
CXXFLAGS += -D DIV_FIX

CXXFLAGS += -D PARTICLE_SIMULATOR_THREAD_PARALLEL 
# MPI
CXXFLAGS += -D PARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -D MPICH_IGNORE_CXX_SEEKC

CXXFLAGS += -D SOFT_PERT -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D CLUSTER_VELOCITY
CXXFLAGS += -D USE_QUAD
CXXFLAGS += -D STELLAR_EVOLUTION

CXXFLAGS += -D PROFILE
CXXFLAGS += -D HARD_CHECK_ENERGY
CXXFLAGS += -D TIDAL_TENSOR_3RD

# debug
#CXXFLAGS +=  -D PETAR_DEBUG -D AR_DEBUG -D AR_DEBUG_DUMP -D AR_WARN -D HARD_DEBUG -D HARD_DUMP -D CLUSTER_DEBUG -D ARTIFICIAL_PARTICLE_DEBUG 
CXXFLAGS += -D INTERFACE_DEBUG_PRINT
CXXFLAGS += -D INTERFACE_DEBUG

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o 

FILELIST=interface.cc interface.py interface.h test_interface.cc test_interface.py run_one_cluster_movie.py

CODELIB = 

all: petar_worker __init__.py test_interface.py run_one_cluster_movie.py

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ petar_worker worker_code.cc
#	make -C src clean

distclean: clean
	rm -rf src $(FILELIST)

#$(CODELIB):
#	make -C src all

src:
	mkdir src

src/PeTar: |src
	npath=$(shell pwd); \
	cd src; \
	wget "https://github.com/AICS-SC-GROUP/PeTar/archive/master.zip" ; \
	unzip master.zip; \
	mv PeTar-master PeTar; \
	rm -f master.zip; \
	cd $(npath); \
	cp src/PeTar/amuse-interface/Makefile ./

src/SDAR: |src
	npath=$(shell pwd); \
	cd src; \
	wget "https://github.com/lwang-astro/SDAR/archive/master.zip" ; \
	unzip master.zip; \
	mv SDAR-master SDAR; \
	rm -f master.zip; \
	cd $(npath);

src/FDPS: |src
	npath=$(shell pwd); \
	cd src; \
	wget "https://github.com/FDPS/FDPS/archive/v5.0g.zip" ;\
	unzip v5.0g.zip; \
	mv FDPS-5.0g FDPS; \
	rm -f v5.0g.zip; \
	cd $(npath);



__init__.py:
	touch $@

$(FILELIST): |src/PeTar
	ln -sf src/PeTar/amuse-interface/$@

src/PeTar/src/get_version.hpp: src/PeTar/src/get_version.hpp.in |src/PeTar
	sed 's/@VERSION@/'`cat src/PeTar/VERSION`'_'`cat src/SDAR/VERSION`'/g' src/PeTar/src/get_version.hpp.in >src/PeTar/src/get_version.hpp

test_interface: $(OBJS) test_interface.cc
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $^ -o $@ $(SC_MPI_CLIBS) $(LIBS) $(AM_LIBS)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py petarInterface -o $@

worker_code.h: interface.py 
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py petarInterface -o $@

petar_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS) interface.h |src/PeTar src/SDAR src/FDPS
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(LDFLAGS) $< $(OBJS) $(CODELIB) -o $@ $(SC_MPI_CLIBS) $(LIBS) $(AM_LIBS)

interface.o: interface.cc src/PeTar/src/get_version.hpp 
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -c -o $@ $< 
