 #/***************************************************************************
 #   Copyright (C) 2015,2016 by Edson Borin, Hercules Cardoso da Silva, Ian Liu Rodrigues
 #   Authors: Edson Borin - edson@ic.unicamp.br
 #	      Hercules Cardoso da Silva - hercules.cardoso.silva1@gmail.com
 #	      Ian Liu Rodrigues - ian.liu@ggaunicamp.com 
 #   
 #    This file is part of SYCL-OpenCL-OpenMP-Benchmark.
 #
 #   SYCL-OpenCL-OpenMP-Benchmark is free software: you can redistribute it and/or modify
 #   it under the terms of the GNU General Public License as published by
 #   the Free Software Foundation, either version 3 of the License, or
 #   (at your option) any later version.
 #
 #   Foobar is distributed in the hope that it will be useful,
 #   but WITHOUT ANY WARRANTY; without even the implied warranty of
 #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #   GNU General Public License for more details.
 #
 #   You should have received a copy of the GNU General Public License
 #   along with Foobar.  If not, see <http://www.gnu.org/licenses/>. 

default:
	@cat README	
nvFlags= -I/usr/local/cuda-7.0/include -L/usr/lib/nvidia-current
HEADERS = log.h su.h uthash.h vector.h
triSYCL_DIR=$(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)

#AMD FLAGS
#FLAGS=-DINTEL -lm -O3 -DGPU -I/opt/AMDAPPSDK-3.0/include/CL/ -L/opt/AMDAPPSDK-3.0/lib/x86_64 -lOpenCL
#SYCL FLAGS
FLAGS = -O3 -lm  -L/opt/AMDAPPSDK-3.0/lib/x86_64 -lOpenCL -std=c++1y -I$(triSYCL_DIR)/include \
	-I$(triSYCL_DIR)/tests/common -fopenmp

#FLAGS FLAGS = -g -lm -lOpenCL -std=c++1y -I$(triSYCL_DIR)/include -I$(triSYCL_DIR)/tests/common -fopenmp
#FLAGS= -DINTEL -DCPU -lOpenCL -O3 -lm
#debug flags
#FLAGS= -DINTEL -DGPU   -lOpenCL -O3
INPUTDATA=simple-synthetic.su
c=g++-4.9
cmp.sycl.x: $(HEADERS) semblance.cpp cmp.cpp
	$c -I ./ $(EXTRA) -o cmp.opencl.x semblance.cpp su.cpp interpol.cpp cmp.cpp log.c  $(FLAGS)	
run: cmp.x
	./cmp.x 1.98e-7 1.77e-6 101 600 0.002 $(INPUTDATA)
clean:
	rm -f *.x c.su cmp.coher.su cmp.stack.su
