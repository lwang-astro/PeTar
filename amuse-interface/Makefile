# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=$(shell amusifier --get-amuse-dir)
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

INCLUDE  += -I src/PeTar/src -I src/ARModule/src -I src/FDPS/src 

MPIFLAGS =  -DPARTICLE_SIMULATOR_MPI_PARALLEL -DMPICH_IGNORE_CXX_SEEKC

CXXFLAGS = ${INCLUDE} -fPIC -O2 -Wall -std=c++11  -fopenmp -march=core-avx2 -DINTRINSIC_X86 -DDIV_FIX -DUSE_QUAD  -DPARTICLE_SIMULATOR_THREAD_PARALLEL -D PROFILE  -D HARD_CHECK_ENERGY -D SOFT_PERT -D AR_TTL -D AR_TTL_SLOWDOWN_INNER -D SLOWDOWN_TIMESCALE -D CLUSTER_VELOCITY -D PETAR_DEBUG -D AR_DEBUG -D AR_DEBUG_DUMP -D AR_WARN -D HARD_DEBUG -D HARD_DUMP -D CLUSTER_DEBUG -D ARTIFICIAL_PARTICLE_DEBUG 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o 

CODELIB = 

all: petar_worker 

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ petar_worker worker_code.cc
	make -C src clean

distclean: clean
	make -C src distclean

$(CODELIB):
	make -C src all

interface.cc:
	bash down.sh 

test_interface:   $(OBJS) test_interface.o
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py petarInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py petarInterface -o $@

petar_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 
