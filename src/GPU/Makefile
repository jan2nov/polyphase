SDK := /home/novotny/NVIDIA_GPU_Computing_SDK/CUDALibraries/common/inc/
INC := -I$(CUDA_HOME)/include -I$(SDK) -I
LIB := -L$(CUDA_HOME)/lib64 -lcudart -lcufft -lfftw3f -lcuda

NVCCFLAGS := -g -G -arch=sm_20 -DCUDA_DEVICE=1 --ptxas-options=-v --use_fast_math

all:	main

main:	main.cu
	nvcc main.cu -o main $(INC) $(NVCCFLAGS) $(LIB)

clean:	
	rm -f main