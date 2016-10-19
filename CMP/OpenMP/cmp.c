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
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <vector.h>
#include <uthash.h>
#include <interpol.h>
#include <semblance.h>
#include <log.h>
#include <su.h>
#include <sys/time.h>

double mysecond()
{
        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
struct cdp_traces {
	int cdp;
	vector_t(su_trace_t) traces;
	UT_hash_handle hh;
};

float getmax_C(fast_aperture_t *ap, int t0s, float c0, float c1, int nc, float *sem, float *stack, float idt, int tau,int w, float dt)
{
	float Copt;
	float smax = 0;
	float _stack = 0;
	float aux = (float) 1.0f/nc;
	float aux1 = c1-c0;
	float aux3  = aux1*aux;
	float t0 = t0s * dt;
	float vC[nc];
	float m0x, m0y;
	fast_su_get_midpoint(ap->traces.a[0], &m0x, &m0y);
	for (int i = 0; i < nc; i++) {
		vC[i]= c0+(aux3)*i;
	}
		smax = semblance_2d(ap,vC, t0s, m0x, m0y, stack,idt,t0,tau,w,nc,&Copt);
	
	if (sem)
		*sem = smax;

	return Copt;
}

int main(int argc, char *argv[])
{
	if (argc != 8) {
		fprintf(stderr, "Usage: %s C0 C1 NC APH TAU INPUT NUMBEROFTHREADS\n", argv[0]);
		exit(1);
	}
	char namef[100];
	sprintf(namef,"time_1_%s.txt\n",argv[7]);
	FILE * timeFile = fopen(namef,"a+");
	float c0 = atof(argv[1]);
	float c1 = atof(argv[2]);
	int nc = atoi(argv[3]);
	int aph = atoi(argv[4]);
	float tau = strtof(argv[5], NULL);
	char *path = argv[6];
	FILE *fp = fopen(path, "r");

	su_trace_t tr;
	struct cdp_traces *cdp_traces = NULL;
	while (su_fgettr(fp, &tr)) {
		float hx, hy;
		su_get_halfoffset(&tr, &hx, &hy);
		if (hx*hx + hy*hy > aph*aph)
			continue;
		int cdp = tr.cdp;
		struct cdp_traces *val;
		HASH_FIND_INT(cdp_traces, &cdp, val);
		if (!val) {
			val = malloc(sizeof(struct cdp_traces));
			val->cdp = cdp;
			vector_init(val->traces);
			HASH_ADD_INT(cdp_traces, cdp, val);
		}
		vector_push(val->traces, tr);
	}

	FILE *C_out = fopen("c.su", "w");
	FILE *S_out = fopen("cmp.coher.su", "w");
	FILE *stack = fopen("cmp.stack.su", "w");
	float progress_max = 1.0f / (HASH_COUNT(cdp_traces) - 1);
	int k =0;
	
	int cap = 4096;	
	struct cdp_traces *it = (struct cdp_traces*) malloc(sizeof(struct cdp_traces)*cap);	
	struct cdp_traces *iter;
	int i;
 	double   t= mysecond();
	for (i = 0, iter = cdp_traces; iter; iter = iter->hh.next,i++) {
		if(i == cap-1) {
			cap = cap*2;
			it = realloc(it, sizeof(struct cdp_traces)*cap);
		}
		it[i] = *iter;
	}
	int size = i;
	#pragma omp parallel
	{
		#pragma omp  for schedule(dynamic) nowait
		for (int i = 0; i < size; i++) {
			struct cdp_traces *iter = &it[i];
			su_trace_t *trs = iter->traces.a;
			su_trace_t ctr, str, stacktr;
			memcpy(&ctr, &trs[0], SU_HEADER_SIZE);
			ctr.offset = 0;
			ctr.sx = ctr.gx = (trs[0].sx + trs[0].gx) >> 1;
			ctr.sy = ctr.gy = (trs[0].sy + trs[0].gy) >> 1;

			memcpy(&str, &ctr, SU_HEADER_SIZE);
			memcpy(&stacktr, &ctr, SU_HEADER_SIZE);
			su_init(&ctr);
			su_init(&str);
			su_init(&stacktr);

			fast_aperture_t ap;
			ap.ap_m = 0;
			ap.ap_h = aph;
			ap.ap_t = tau;
			su_trace_t *tr = &vector_get(iter->traces, 0);
			float dt = (float) tr->dt / 1000000;
			float idt = (float) 1.0f/ dt;
	
			int tau; 
			int temptau= (int)(ap.ap_t * idt);	
			 if(temptau > 0)
			{
				tau = temptau;
			}
			else
			{
				tau=0;
			}
			int w = 2 * tau + 1;
			vector_init(ap.traces);
			vector_alloc(ap.traces,iter->traces.len);
			for (int i = 0; i < iter->traces.len; i++)
			{
				vector_push_otim(ap.traces,convert_to_fast(&vector_get(iter->traces,i)));
			}
			for (int t0 = 0; t0 < trs[0].ns; t0++) {
				float sem, stk;
				float C = getmax_C(&ap, t0, c0, c1, nc, &sem, &stk,idt,tau,w,dt);
				ctr.data[t0] = C;
				str.data[t0] = sem;
				stacktr.data[t0] = stk;
			}
+			#pragma omp critical 
			{
				su_fputtr(C_out, &ctr);
				su_fputtr(S_out, &str);
				su_fputtr(stack, &stacktr);
			}
			log_progress(++k * progress_max, "Processing CDP %d", ctr.cdp);
		}
	}
	double end= mysecond();
	fprintf(timeFile,"%f\n",end-t);
	printf("\n");
	return 0;
}
