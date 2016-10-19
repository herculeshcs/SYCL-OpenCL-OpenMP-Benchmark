 /***************************************************************************
 *   Copyright (C) 2015,2016 by Edson Borin, Hercules Cardoso da Silva, Ian Liu Rodrigues
 *   Authors: Edson Borin - edson@ic.unicamp.br
 *	      Hercules Cardoso da Silva - hercules.cardoso.silva1@gmail.com
 *	      Ian Liu Rodrigues - ian.liu@ggaunicamp.com 
 *   
 *    This file is part of SYCL-OpenCL-OpenMP-Benchmark.
 *
 *   SYCL-OpenCL-OpenMP-Benchmark is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Foobar.  If not, see <http://www.gnu.org/licenses/>. 
 ****************************************************************************/	
#ifndef VECTOR_H__
#define VECTOR_H__

#include <stdlib.h>

#define vector_t(type) struct {int len, cap; type t; type *a;}
#define vector_init(v) ((v).len = (v).cap = 0, (v).a = NULL)
#define vector_push(v, x) do { \
	if ((v).len == (v).cap) { \
		(v).cap *= 2, (v).cap++; \
		(v).a = realloc((v).a, sizeof(__typeof((v).t))*(v).cap); \
	} \
	(v).a[(v).len++] = x; \
} while (0)
#define vector_get(v, i) ((v).a[i])
#define vector_set(v, i, x) ((v).a[i]=x)

#endif /* VECTOR_H__ */
