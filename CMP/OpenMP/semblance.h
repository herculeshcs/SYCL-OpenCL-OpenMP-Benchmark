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

#ifndef SEMBLANCE_H__
#define SEMBLANCE_H__

#include <vector.h>
#include <su.h>

typedef struct fast_aperture fast_aperture_t;
struct fast_aperture{
	float ap_m, ap_h, ap_t;
	vector_t(fast_su_trace_t*) traces;
};
typedef struct aperture aperture_t;

struct aperture {
	float ap_m, ap_h, ap_t;
	vector_t(su_trace_t*) traces;
};

float semblance_2d(fast_aperture_t *ap,float*  vC,int t0s, float m0x, float m0y,float *stack,float idt, float t0, int tau, int w, int nc, float * Copt);
 void inline_su_get_midpoint(su_trace_t *tr, float *mx, float *my);
void inline_su_get_halfoffset(su_trace_t *tr, float *hx, float *hy);
static float inline_get_scalco(su_trace_t *tr);

#endif /* SEMBLANCE_H__ */
