default:
	@cat README	
HEADERS = log.h su.h uthash.h vector.h
triSYCL_DIR=$(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)
#SYCL FLAGS
SYCLFLAGS = -O3 -lm -std=c++1y -lOpenCL -I$(triSYCL_DIR)/include \
	-I$(triSYCL_DIR)/tests/common -fopenmp -openmp
FLAGS= -O3 -fopenmp
c=g++
sycl: stencil-sycl.cpp
	$c stencil-sycl.cpp $(FLAGS) -osycl.x
omp: stencil27.cpp
	$c stencil27.cpp $(FLAGS) $(SYCLFLAGS) -oomp.x
cl: stencil-ocl.cpp
	$c stencil-ocl.cpp -g -lm -lOpenCL -ocl.x
ompn: stencil27-normal.cpp
	$c stencil27-normal.cpp $(FLAGS) -o ompn.x
clean:
	rm -f *.x
