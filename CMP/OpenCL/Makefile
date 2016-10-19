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
FLAGS= -DLINUX -DCPU	
INPUTDATA=simple-synthetic.su
c= gcc
cmp.gcc.x: $(HEADERS) semblance.c log.c
	$c -lOpenCL  $(FLAGS)   -I ./ -O3 semblance.c log.c -o cmp.gcc.x

cmp.clang.x: $(HEADERS) semblance.c log.c
	$c -lOpenCL $(FLAGS) -I ./ -O3 semblance.c log.c -o cmp.clang.x

cmp.opencl.x: $(HEADERS) semblance-ocl.c log.c
	$c   -O3  $(FLAGS) -lOpenCL -I ./ $(EXTRA) -o cmp.opencl.x semblance-ocl.c log.c

run: cmp.x
	./cmp.x 1.98e-7 1.77e-6 101 600 0.002 $(INPUTDATA)

clean:
	rm -f *.x c.su cmp.coher.su cmp.stack.su
