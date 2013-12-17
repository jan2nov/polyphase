###############################################################
# Makefile for Artemis project -- Polyphase filter
# created by JaN
# contact: jan.novotny@fpf.slu.cz, karel.adamek@fpf.slu.cz
###############################################################

###############################################################
# CUDA_HOME are supposed to be on default position
# and set it in your PATH .bashrc
###############################################################
SDK := /home/novotny/NVIDIA_GPU_Computing_SDK/CUDALibraries/common/inc/
INC := -I$(CUDA_HOME)/include #-I$(SDK)
LIB := -L$(CUDA_HOME)/lib64 -lcudart -lcufft -lfftw3f -lcuda

# use this compilers
GCC = g++
NVCC = nvcc

# name of executables file
ANALYZE = analyze.x

###############################################################
# Basic flags for compilers, one for debug options
# fmad flags used for reason of floating point operation
###############################################################
NVCCFLAGS = -O3 -fmad=false -arch=sm_20 --ptxas-options=-v --use_fast_math -Xcompiler -Wextra
#NVCCFLAGS= -g -G -arch=sm_20 --ptxas-options=-v --use_fast_math -Wall -Wextra -Xcompiler
GCC_OPTS =-O3 -Wall -Wextra 

all:	analyze

main.o:	main.c
	$(GCC) -c main.c $(GCC_OPTS) $(INC)

analyze: polyphase-analyze.o utils_cuda.h timer.h polyphase.o reference.o Makefile
	$(NVCC) -o $(ANALYZE) polyphase-analyze.o polyphase.o reference.o $(LIB) $(NVCCFLAGS)

polyphase-analyze.o: polyphase-analyze.c utils_file.h timer.h data.h
	$(GCC) -c polyphase-analyze.c $(GCC_OPTS) $(INC)

polyphase.o: polyphase.cu utils_cuda.h data.h timer.h
	$(NVCC) -c polyphase.cu $(NVCCFLAGS)

reference.o: reference.c
	$(GCC) -c reference.c $(GCC_OPTS) $(INC)

clean:	
	rm -f *.o *.~ main $(ANALYZE)

run:
	./$(ANALYZE) 1 120 >>gpu_1.dat
	./$(ANALYZE) 1 1200 >>gpu_1.dat
	./$(ANALYZE) 1 6000 >>gpu_1.dat
	./$(ANALYZE) 1 12000 >>gpu_1.dat
	./$(ANALYZE) 1 60000 >>gpu_1.dat
	./$(ANALYZE) 2 120 >>gpu_2.dat
	./$(ANALYZE) 2 1200 >>gpu_2.dat
	./$(ANALYZE) 2 6000 >>gpu_2.dat
	./$(ANALYZE) 2 12000 >>gpu_2.dat
	./$(ANALYZE) 2 60000 >>gpu_2.dat
	./$(ANALYZE) 4 120 >>gpu_4.dat
	./$(ANALYZE) 4 1200 >>gpu_4.dat
	./$(ANALYZE) 4 6000 >>gpu_4.dat
	./$(ANALYZE) 4 12000 >>gpu_4.dat
	./$(ANALYZE) 4 60000 >>gpu_4.dat
	./$(ANALYZE) 1 60000 >>gpu_th.dat
	./$(ANALYZE) 2 60000 >>gpu_th.dat
	./$(ANALYZE) 4 60000 >>gpu_th.dat
	./$(ANALYZE) 8 60000 >>gpu_th.dat
	./$(ANALYZE) 16 60000 >>gpu_th.dat
