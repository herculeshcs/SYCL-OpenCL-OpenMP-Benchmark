/*
 *  Copyright 2014-2015 Ian Liu Rodrigues <ian.liu88@gmail.com>
 *
 *  This file is part of CMP.
 *
 *  CMP is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CMP is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CMP.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SEMBLANCE_H__
#define SEMBLANCE_H__

#include <vector.h>
#include <su.h>

typedef struct aperture aperture_t;

struct aperture {
	float ap_m, ap_h, ap_t;
	vector_t(su_trace_t*) traces;
};

float semblance_2d(aperture_t *ap,
		float A, float B, float C,
		int t0s, float m0x, float m0y,
		float *stack);

#endif /* SEMBLANCE_H__ */
