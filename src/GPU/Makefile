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
# g++ just because the file write
GCC = g++
NVCC = nvcc

# name of executables file
ANALYZE = analyze.x
STREAM = stream.x

###############################################################
# Basic flags for compilers, one for debug options
# fmad flags used for reason of floating point operation
###############################################################
NVCCFLAGS = -O3 -fmad=false -arch=sm_20 --ptxas-options=-v --use_fast_math -Xcompiler -Wextra
#NVCCFLAGS= -g -G -arch=sm_20 --ptxas-options=-v --use_fast_math -Xcompiler -Wextra
GCC_OPTS =-O3 -Wall -Wextra 



all:	analyze stream

main.o:	main.c
	$(GCC) -c main.c $(GCC_OPTS) $(INC)

stream: polyphase-stream.o utils_cuda.h timer.h stream-test.o reference.o Makefile
	$(NVCC) -o $(STREAM) polyphase-stream.o stream-test.o reference.o $(LIB) $(NVCCFLAGS)

analyze: polyphase-analyze.o utils_cuda.h timer.h polyphase.o reference.o Makefile
	$(NVCC) -o $(ANALYZE) polyphase-analyze.o polyphase.o reference.o $(LIB) $(NVCCFLAGS)

polyphase-analyze.o: polyphase-analyze.c  timer.h data.h
	$(GCC) -c polyphase-analyze.c $(GCC_OPTS) $(INC)

polyphase-stream.o: polyphase-stream.c utils_file.h timer.h data.h
	$(GCC) -c polyphase-stream.c $(GCC_OPTS) $(INC)

stream-test.o: stream-test.cu utils_cuda.h data.h timer.h
	$(NVCC) -c stream-test.cu $(NVCCFLAGS)

polyphase.o: polyphase.cu utils_cuda.h utils_file.h data.h timer.h
	$(NVCC) -c polyphase.cu $(NVCCFLAGS)

reference.o: reference.c
	$(GCC) -c reference.c $(GCC_OPTS) $(INC)

clean:	
	rm -f *.o *.~ *.dat main $(ANALYZE)
	
run:
	rm -f GPU-polyphase.dat
	@printf '#num \t fir \t\t fir+fft \t bandwidth \t flops \t\t HtoD \t\t DtoH\n' > GPU-polyphase.dat
	./$(ANALYZE) 2 10007 > GPU-10.dat
	./$(ANALYZE) 2 15007 > GPU-15.dat
	./$(ANALYZE) 2 20007 > GPU-20.dat
	./$(ANALYZE) 2 30007 > GPU-30.dat
	./$(ANALYZE) 2 40007 > GPU-40.dat
	./$(ANALYZE) 2 50007 > GPU-50.dat
	./$(ANALYZE) 2 60007 > GPU-60.dat

help:
	@printf '\n\n-----------------------------------------------------\n'
	@printf '\t\tWelcome, to the help\n\n'
	@printf 'For compiling all type: \t\t make\n'
	@printf 'For compiling analyze: \t\t\t make analyze\n'
	@printf 'For compiling stream: \t\t\t make stream\n'
	@printf 'To create data file with times: \t make run\n'
	@printf '-----------------------------------------------------\n\n'
