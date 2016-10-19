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

#include "semblance.h"

#include <math.h>
#include <string.h>
#include <vector.h>
#include <utils.h>
#include <interpol.h>

static float time2d(float A, float B, float C, float t0,
		float m0x, float m0y,
		float mx, float my,
		float hx, float hy)
{
	float _mx = mx - m0x;
	float _my = my - m0y;
	float _m2 = _mx*_mx + _my*_my;
	float t2 = t0 + A * sqrt(_m2);
	t2 *= t2;
	t2 += B * _m2;
	t2 += C * (hx*hx + hy*hy);
	if (t2 < 0)
		return -1;
	else
		return sqrt(t2);
}

static float get_halfoffset(su_trace_t *tr)
{
	float hx = (float)(tr->gx - tr->sx) / 2;
	float hy = (float)(tr->gy - tr->sy) / 2;
	return sqrt(hx * hx + hy * hy);
}

static float get_midpoint(su_trace_t *tr)
{
	float mx = (float)(tr->gx + tr->sx) / 2;
	float my = (float)(tr->gy + tr->sy) / 2;
	return sqrt(mx * mx + my * my);
}

float semblance_2d(aperture_t *ap,
		float A, float B, float C,
		int t0s, float m0x, float m0y,
		float *stack)
{
	su_trace_t *tr = vector_get(ap->traces, 0);
	float dt = (float) tr->dt / 1000000;
	float idt = 1 / dt;
	float t0 = t0s * dt;
	int tau = MAX((int)(ap->ap_t * idt), 0);
	int w = 2 * tau + 1;
	float num[w], den[w];
	memset(&num[0], 0, sizeof(num));
	memset(&den[0], 0, sizeof(den));
	int M = 0, skip = 0;
	float _stack = 0;
	int count =0;
	for (int i = 0; i < ap->traces.len; i++) {
		tr = vector_get(ap->traces, i);
		float mx, my, hx, hy;
		su_get_midpoint(tr, &mx, &my);
		su_get_halfoffset(tr, &hx, &hy);
		float t = time2d(A, B, C, t0, m0x, m0y, mx, my, hx, hy);
		int it = (int)(t * idt);
		if (it - tau >= 0 && it + tau+1 < tr->ns) {
			for (int j = 0; j < w; j++) {
				int k = it + j - tau;
				float v = interpol_linear(k, k+1,
						tr->data[k], tr->data[k+1],
						t*idt + j - tau);
				num[j] += v;
				den[j] += v*v;
				_stack += v;
			}
			M++;
		} else if (++skip == 2) {
			goto error;
		}

	}

	float sem = 0;
	float aux = 0;
	for (int j = 0; j < w; j++) {
		sem += num[j] * num[j];
		aux += den[j];
	}

	if (stack) {
		_stack /= M*w;
		*stack = _stack;
	}

	return sem / aux / M;

error:
	return 0;
}
