FDPS_PATH = -I../FDPS/src
SDAR_PATH = -I../SDAR/src
GPERF_PATH= /opt/gperftools/2.7
CUDA_PATH = /usr/local/cuda
ROOT_PATH = .

#use_k_computer = yes
#use_xc30_naoj = yes
use_x86 = yes
use_mpi = yes
use_gpu_cuda=yes
use_simd= yes
use_omp = yes
use_gperf=yes
#use_intel=yes
#debug_mode=yes

#MT_FLAGS += -D HARD_CM_KICK
MT_FLAGS += -D USE_QUAD
MT_FLAGS += -D SOFT_PERT
MT_FLAGS += -D AR_TTL
MT_FLAGS += -D AR_SLOWDOWN_INNER
#MT_FLAGS += -D AR_TTL_GT_BINARY_INNER
#MT_FLAGS += -D FIX_CHANGEOVER
MT_FLAGS += -D AR_SLOWDOWN_TIMESCALE
#MT_FLAGS += -D AR_SLOWDOWN_MASSRATIO
MT_FLAGS += -D CLUSTER_VELOCITY
MT_FLAGS += -D HARD_CHECK_ENERGY
#MT_FLAGS += -D KDKDK_2ND
#MT_FLAGS += -D KDKDK_4TH
#MT_FLAGS += -D ONLY_SOFT
#MT_FLAGS += -D STELLAR_EVOLUTION

#DEBFLAGS += -DPARTICLE_SIMULATOR_DEBUG_PRINT
#DEBFLAGS += -D INTEGRATED_CUTOFF_FUNCTION
DEBFLAGS += -D PETAR_DEBUG
DEBFLAGS += -D AR_DEBUG
DEBFLAGS += -D AR_DEBUG_DUMP
#DEBFLAGS += -D AR_DEBUG_PRINT
#DEBFLAGS += -D AR_DEEP_DEBUG
DEBFLAGS += -D AR_WARN
#DEBFLAGS += -D AR_COLLECT_DS_MODIFY_INFO
DEBFLAGS += -D HARD_DEBUG
DEBFLAGS += -D HARD_DUMP
#DEBFLAGS += -D HARD_INTERRUPT_PRINT
#DEBFLAGS += -D STABLE_CHECK_DEBUG
DEBFLAGS += -D CLUSTER_DEBUG
DEBFLAGS += -D ARTIFICIAL_PARTICLE_DEBUG
#DEBFLAGS += -D CLUSTER_DEBUG_PRINT
#DEBFLAGS += -D HARD_CLUSTER_PRINT
#DEBFLAGS += -D HARD_DEBUG_PRINT
#DEBFLAGS += -D HARD_DEBUG_PROFILE
#DEBFLAGS += -D NAN_CHECK_DEBUG
#DEBFLAGS += -D DATA_DEBUG
#DEBFLAGS += -D FIX_STEP_DEBUG
#DEBFLAGS += -D DEBUG
#DEBFLAGS += -D MAIN_DEBUG
#DEBFLAGS += -D CORRECT_FORCE_DEBUG

SIMD_DEBFLAGS += -DRSQRT_NR_EPJ_X2
#SIMD_DEBFLAGS += -DRSQRT_NR_EPJ_X4
#SIMD_DEBFLAGS += -DCALC_SP_64bit
#SIMD_DEBFLAGS += -DRSQRT_NR_SPJ_X2
#SIMD_DEBFLAGS += -DRSQRT_NR_SPJ_X4
SIMD_DEBFLAGS += -DAVX_PRELOAD

DEBUG_OPT_FLAGS = -g -O0 -fbounds-check -Wall  -D SANITY_CHECK_REALLOCATABLE_ARRAY
HARD_DEBFLAGS+= -D AR_DEBUG -D AR_DEBUG_DUMP -D AR_WARN -D HARD_DEBUG -D HARD_CHECK_ENERGY -D HARD_DEBUG_PRINT -D ADJUST_GROUP_DEBUG -D HERMITE_DEBUG -D AR_COLLECT_DS_MODIFY_INFO
HARD_MT_FLAGS += -D AR_TTL -D AR_SLOWDOWN_INNER -D AR_SLOWDOWN_TIMESCALE -D HARD_CHECK_ENERGY -D STABLE_CHECK_DEBUG_PRINT -D ARTIFICIAL_PARTICLE_DEBUG


#------------------------------------
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

#------------------------------------
ifeq ($(use_xc30_naoj),yes)
CXX = time CC
CXXNOMPI = time CC
OPTFLAGS = -O3 -Wall

CXXFLAGS = -ffast-math -funroll-loops
CXXFLAGS += -march=core-avx2
CXXFLAGS += -std=c++11
CXXFLAGS += -DINTRINSIC_X86

#FDPSFLAGS += -DUSE_GNU_PARALLEL_SORT
endif # xc30

#------------------------------------
ifeq ($(use_x86),yes)

ifeq ($(use_mpi),yes) 
CXX = time mpicxx

FDPSFLAGS = -DPARTICLE_SIMULATOR_MPI_PARALLEL
FDPSFLAGS += -DMPICH_IGNORE_CXX_SEEKC
else # no mpi

ifeq ($(use_intel),yes)
CXX = time icc
else # no intel
CXX = time g++
endif # intel

endif # mpi

ifeq ($(use_intel),yes)
CXXNOMPI = time icc
else # no intel
CXXNOMPI = time g++
endif # intel

ifeq ($(debug_mode),yes)
OPTFLAGS = -g -O0 -fbounds-check -Wall -D SANITY_CHECK_REALLOCATABLE_ARRAY 
else # no debug
OPTFLAGS = -O2 -Wall 
#OPTFLAGS += -ffast-math -funroll-loops
endif # debug

CXXFLAGS = -std=c++11
CXXFLAGS += -Wall
#CXXFLAGS += -march=skylake-avx512
CXXFLAGS += -march=core-avx2
#CXXFLAGS += -mavx
CXXFLAGS += -DINTRINSIC_X86

ifeq ($(use_omp),yes)
CXXFLAGS += -fopenmp
CXXFLAGS += -D PARTICLE_SIMULATOR_THREAD_PARALLEL
endif # omp

endif # end x86
#---------------------------------------------

#---------------------------------------------
ifeq ($(use_gpu_cuda),yes)
CUDAFLAGS = -D PARTICLE_SIMULATOR_GPU_MULIT_WALK_INDEX
FDPSFLAGS += $(CUDAFLAGS)
NVCC = nvcc -std=c++11 -Xcompiler="$(OPTFLAGS) $(CXXFLAGS) $(CUDAFLAGS) $(MT_FLAGS)"
CXXLIBS += -L$(CUDA_PATH)/lib64 -lcudart -lgomp
OBJS += force_gpu_cuda.o
CXXFLAGS += -D USE_GPU
CXXFLAGS += -D GPU_PROFILE
endif # gpu cuda

#-------------------------------------
ifeq ($(use_simd),yes)
CXXFLAGS += -D DIV_FIX
#CXXFLAGS += -D P3T_64BIT
CXXFLAGS += -D USE_SIMD
BIN_NAME =$(BIN_NAME).simd
endif  # simd
#-------------------------------------

ifneq (x$(use_gperf),x)
CXXFLAGS += -I$(GPERF_PATH)/include -D GPERF_PROFILE
CXXLIBS += -L$(GPERF_PATH)/lib -lprofiler -ltcmalloc 
endif # gperf
#-------------------------------------

CXXFLAGS += -D PROFILE

VPATH=${ROOT_PATH}/src ${ROOT_PATH}/test

SRC = main.cc ${shell ls ${ROOT_PATH}/src/*.hpp} 
INCLUDE  = -I${ROOT_PATH}/src $(FDPS_PATH) $(SDAR_PATH)

all: petar petar.hard.debug

petar:  $(SRC) $(OBJS)
	$(CXX) $(INCLUDE) $(OPTFLAGS) $(CXXFLAGS) $(FDPSFLAGS) $(MT_FLAGS) $(DEBFLAGS) -o $@ $< $(OBJS) $(CXXLIBS)

petar.hard.debug: hard_debug.cxx
	$(CXXNOMPI) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(MT_FLAGS) $(HARD_DEBFLAGS) -D HARD_DEBUG_PRINT_TITLE -D STABLE_CHECK_DEBUG -o $@ $< $(CXXLIBS)

petar.hard.test: hard_test.cxx
	$(CXXNOMPI) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS) $(HARD_MT_FLAGS) $(HARD_DEBFLAGS) -o $@ $< $(CXXLIBS)

petar.simd.test: simd_test.cxx $(OBJS)
	$(CXXNOMPI) $(INCLUDE) $(OPTFLAGS) $(CXXFLAGS) $(SIMD_DEBFLAGS) -DPARTICLE_SIMULATOR_THREAD_PARALLEL $^ -o $@  $(CXXLIBS)

petar.tt.test: tidal_tensor_test.cxx
	$(CXXNOMPI) $(INCLUDE) $(DEBUG_OPT_FLAGS) $(CXXFLAGS)  $< -o $@  $(CXXLIBS)

force_gpu_cuda.o: force_gpu_cuda.cu
	$(NVCC) $(INCLUDE) -I$(CUDA_PATH)/samples/common/inc -c $< -o $@ 

clean:
	rm *.out *.o 
