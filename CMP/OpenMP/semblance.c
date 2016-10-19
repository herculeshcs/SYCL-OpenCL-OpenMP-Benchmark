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
#include "semblance.h"

#include <math.h>
#include <string.h>
#include <vector.h>
#include <utils.h>
#include <interpol.h>
static float time2d(float C, float t0,float m0x, float m0y,float mx, float my,float hx, float hy)
{
	float t2 = t0*t0;
	t2 += C * (hx*hx + hy*hy);
	if (t2 < 0)
		return -1;
	else
		return sqrt(t2);
}

static float get_halfoffset(su_trace_t *tr)
{
	float hx = (float)(tr->gx - tr->sx) /2;
	float hy = (float)(tr->gy - tr->sy) /2;
	return sqrt(hx * hx + hy * hy);
}

static float get_midpoint(su_trace_t *tr)
{
	float mx = (float)(tr->gx + tr->sx) /2;
	float my = (float)(tr->gy + tr->sy) /2;
	return sqrt(mx * mx + my * my);
}
float semblance_2d(fast_aperture_t *ap,float* vC,int t0s, float m0x, float m0y,float *stack,float idt, float t0, int tau, int w, int nc, float * Copt)
{	
	fast_su_trace_t *tr = vector_get(ap->traces, 0);
	float num[nc][w], den[nc][w];
	int M[nc], vskip[nc];
	char error[nc];
	memset(&num[0], 0, sizeof(num));
	memset(&den[0], 0, sizeof(den));
	memset(&M[0],0,sizeof(M));
	memset(&vskip[0],0,sizeof(vskip));
	memset(&error[0],0,sizeof(error));
	float _stack = 0;
	int len = ap->traces.len;
	float t[len][nc];

	#pragma vector
	for(int i =0;i<len;i++)
	{
		for(int k=0;k<nc;k++)
		{
			float C = vC[k];
			tr = vector_get(ap->traces, i);
			float mx, my, hx, hy;
			fast_su_get_midpoint(tr, &mx, &my);
			fast_su_get_halfoffset(tr, &hx, &hy);
			t[i][k] = time2d(C, t0, m0x, m0y, mx, my, hx, hy);
		}
	}
	for (int i = 0; i < len; i++) {
		for(int k =0;k<nc;k++)
		{
			int it = (int)(t[i][k] * idt);
			int ittau = it -tau;
			if (ittau >= 0 && it + tau < tr->ns) {
				tr = vector_get(ap->traces, i);
				float * data = tr->data;
				float v;
				int k1;
				for (int j = 0; j < w; j++) {
					k1 = ittau + j;
					v = interpol_linear(k1, k1+1,data[k1], data[k1+1],t[i][k]*idt + j - tau);
					 num[k][j] += v;
					 den[k][j] += v*v;
					}
					M[k]++;
			} else if (++vskip[k] == 2) {
				 error[k]=1; //return 0;
			}
		}
	}	
	float sem = 0;
	float aux = 0;
	float smax = 0;
	float auxC=0;
	float auxStack=0;
	for(int k =0;k<nc;k++)
	{
		sem =0;
		aux = 0;	
		_stack=0;
		#pragma simd reduction(+:sem) reduction(+:aux) reduction(+:_stack)
		for (int j = 0; j < w; j++) {
			sem += num[k][j] * num[k][j];
			aux += den[k][j];
			_stack+= num[k][j];
		}
		float s = sem /(aux*M[k]);
		if(s>smax && error[k]==0)
		{
			smax=s;
			auxStack=_stack;
			auxC=vC[k];
			if (stack) {
				auxStack /= M[k]*w;
					
			}
		}
	}
	*Copt= auxC;
	if(stack)
	*stack = auxStack;
	return smax;

}
