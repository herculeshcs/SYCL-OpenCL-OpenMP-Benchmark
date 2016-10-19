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

# C Compiler
CC=gcc
# Compiler flags
CFLAGS= -fopenmp -O3 -std=gnu99 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -fopenmp
LINK_FLAGS=-fopenmp -lm 

INCLUDES=-I ./
OBJECTS=cmp.o gather.o interpol.o log.o semblance.o su.o

cmp.x: $(OBJECTS)
	$(CC) $(OBJECTS) $(LINK_FLAGS) -o cmp.x

%.o : %.c
	$(CC) $(INCLUDES) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) cmp.x c.su cmp.coher.su cmp.stack.su
