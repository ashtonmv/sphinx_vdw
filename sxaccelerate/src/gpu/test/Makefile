ifeq ($(SYS),amd64_sles11)
   # RZG
   GCCPATH=/afs/ipp/.cs/gcc/4.8.3/amd64_sles11/bin/
   CUCFLAGS+=-arch=sm_35 -DSM_35
else
   # ThinkPad
   GCCPATH=/opt/apps/gcc/4.8.3/bin/
endif

# --- number of threads used on the host CPU socket
GPUHOSTTHREADS=4

# --- CUDA compiler
CUC=nvcc
CUCFLAGS+=-O3
CUCFLAGS+=--use_fast_math
CUCFLAGS+=--compiler-bindir=$(GCCPATH)
CUCFLAGS+=-Xcompiler=-Wall
CUCFLAGS+=-Xcompiler=-fopenmp

# C++ host compiler
CXX=$(GCCPATH)g++
CXXFLAGS+=-O3
CXXFLAGS+=-march=native
CXXFLAGS+=-Wall
# the following flag is necessary to enable vectorization of the inner loop
CXXFLAGS+=-funsafe-math-optimizations
#CXXFLAGS+=-fopt-info-vec
#CXXFLAGS+=-fopt-info-vec-missed
CXXFLAGS+=-fopenmp

# --- makefile targets
GPUEXE=sx_gpu

all: $(GPUEXE)

$(GPUEXE): ../SxGpu.cu
	$(CUC) $(CUCFLAGS) -DSTANDALONE -o $@ -lcublas $<

.PHONY: clean
clean:
	rm -f *.o *.so $(GPUEXE)

.PHONY: rungpu
rungpu: $(GPUEXE)
	numactl --cpunodebind=0  env OMP_NUM_THREADS=$(GPUHOSTTHREADS) ./$<
