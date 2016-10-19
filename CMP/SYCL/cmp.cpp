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
#define SQR(x) ((x)*(x))
#ifndef MAX
# define MAX(a, b) ((a)>(b)?(a):(b))
#endif
#ifndef MIN
# define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <vector.h>
#include <uthash.h>
#include <semblance.h>
#include <log.h>
#include <su.h>
#include <sys/time.h>
#include <CL/sycl.hpp>
 using namespace cl::sycl;

static float interpol_linear(float x0, float x1, float y0, float y1, float x)
{
	return (y1 - y0) * (x - x0) / (x1 - x0) + y0;
}
 struct cdp_traces {
 	int cdp;
 	vector_t(su_trace_t) traces;
 	UT_hash_handle hh;
 };

 float getmax_C(aperture_t *ap, int t0s, float c0, float c1, int nc, float *sem, float *stack)
 {
 	float Copt;
 	float smax = 0;
 	float _stack = 0;
 	for (int i = 0; i < nc; i++) {
 		float C = c0 + (c1 - c0) * i / nc;
 		float m0x, m0y;
 		su_get_midpoint(ap->traces.a[0], &m0x, &m0y);
 		float s = semblance_2d(ap, 0, 0, C, t0s, m0x, m0y, &_stack);
 		if (s > smax) {
 			smax = s;
 			Copt = C;
 			*stack = _stack;
 		}
 	}
 	if (sem)
 		*sem = smax;
 	return Copt;
 }

double mysecond()
{
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

 static float get_trace_scalco(float trace_scalco)
 {
 	if (trace_scalco == 0) return 1;
 	else if (trace_scalco > 0)  return trace_scalco;
 	else return 1.0f / trace_scalco;
 }
 static float time2d(float C, float t0,
 	float m0x, float m0y,
 	float mx, float my,
 	float hx, float hy)
 {
  float t2 = t0;
  t2 *= t2;
  t2 += C * (hx*hx + hy*hy);
  if (t2 < 0)
  	return -1;
  else
  	return sqrt(t2);
}
int main(int argc, char *argv[])
{


	if (argc != 8) {
		fprintf(stderr, "Usage: %s C0 C1 NC APH TAU INPUT host\n", argv[0]);
		exit(1);
	}

	float c0 = atof(argv[1]);
	float c1 = atof(argv[2]);
	int nc = atoi(argv[3]);
	int aph = atoi(argv[4]);
	float tau = strtof(argv[5], NULL);
	char *path = argv[6];
	FILE *fp = fopen(path, "r");
	queue myQueue;
	su_trace_t tr;
	struct cdp_traces *cdp_traces = NULL;
	int min_ns = -1;
	int max_ns = -1;
	while (su_fgettr(fp, &tr)) {
		float hx, hy;
		su_get_halfoffset(&tr, &hx, &hy);
		if (hx*hx + hy*hy > aph*aph)
			continue;
		 if (min_ns < 1) min_ns = tr.ns; // First trace.
   		 if (max_ns < 1) max_ns = tr.ns; // First trace.
   		 if (min_ns > tr.ns) min_ns = tr.ns;
   		 if (max_ns < tr.ns) max_ns = tr.ns;
   		 int cdp = tr.cdp;
   		 struct cdp_traces *val;
   		 HASH_FIND_INT(cdp_traces, &cdp, val);
   		 if (!val) {
   		 	val = (struct cdp_traces*) malloc(sizeof(struct cdp_traces));
   		 	val->cdp = cdp;
   		 	vector_init(val->traces);
   		 	HASH_ADD_INT(cdp_traces, cdp, val);
   		 }
   		 vector_push(val->traces, tr);
   		}

	  long unsigned int max_ntrs = 1; // maximum traces per gather
	  void * it;
	  int numberOfCDPs=0;
	  struct cdp_traces** gathers = (struct cdp_traces**) malloc(sizeof(struct cdp_traces*) * HASH_COUNT(cdp_traces));
	  for (it = cdp_traces; it; it = ((struct cdp_traces *)it)->hh.next) {
	  	gathers[numberOfCDPs++] =(struct cdp_traces *)it;
	  	if (max_ntrs < ((struct cdp_traces *)it)->traces.len) 
	  		max_ntrs = ((struct cdp_traces *)it)->traces.len;

	  }
  	long unsigned int ng = 512;
  	long unsigned int mns = (long unsigned int) max_ns;
  	buffer <float,3> cube({mns,max_ntrs,ng});
  	buffer <int,2> gx({max_ntrs,ng});
  	buffer <int,2> gy({max_ntrs,ng});
  	buffer <int,2> sx({max_ntrs,ng});
  	buffer <int,2> sy({max_ntrs,ng});
  	buffer <short,2> scalco({max_ntrs,ng});
  	buffer <int,1> Dntrs(ng);
  	buffer <int,1> tr0_dt(ng);
  	buffer <float,2> ctr({mns,ng});
  	buffer <float,2> stk({mns,ng});
  	buffer <float,2> sem({mns,ng});
  	auto DNtrs = Dntrs.get_access<access::write, access::host_buffer>();
  	auto Dgx = gx.get_access<access::write, access::host_buffer>();
  	auto Dgy = gy.get_access<access::write, access::host_buffer>();
  	auto Dsx = sx.get_access<access::write, access::host_buffer>();
  	auto Dsy = sy.get_access<access::write, access::host_buffer>();
  	auto Dscalco = scalco.get_access<access::write, access::host_buffer>();
  	auto Dcube = cube.get_access<access::write, access::host_buffer>();
  	auto Dctr = ctr.get_access<access::write, access::host_buffer>();
  	auto Dtr0_dt = tr0_dt.get_access<access::write,access::host_buffer> ();	
  	FILE *C_out = fopen("cmp.c.su", "w");
  	FILE *S_out = fopen("cmp.coher.su", "w");
  	FILE *stack = fopen("cmp.stack.su", "w");
  	float progress_max = 1.0f / (HASH_COUNT(cdp_traces) - 1);
  	int k = 0;
  	int total_gathers = numberOfCDPs;
  	int ngathers;
	su_trace_t ctr_t[ng];
  	for(int global_gather = 0; global_gather < total_gathers;){
  		ngathers = MIN(ng, total_gathers - global_gather);
  		int local_gather_idx;
  		for(int local_gather_idx=0;local_gather_idx<ngathers;local_gather_idx++,global_gather++)
  		{
  			su_trace_t *trs  = gathers[global_gather]->traces.a;
  			int         ntrs = gathers[global_gather]->traces.len;
  			 memcpy(&ctr_t[local_gather_idx], &trs[0], SU_HEADER_SIZE);
	 		ctr_t[local_gather_idx].offset = 0;
	 		ctr_t[local_gather_idx].sx = ctr_t[local_gather_idx].gx = (trs[0].sx + trs[0].gx) / 2;
			ctr_t[local_gather_idx].sy = ctr_t[local_gather_idx].gy = (trs[0].sy + trs[0].gy) / 2;
			su_init(&ctr_t[local_gather_idx]);
  			DNtrs[local_gather_idx]= ntrs;
  			Dtr0_dt[local_gather_idx] = trs[0].dt;
  			for (int i = 0; i < ntrs; i++)
  			{
  				su_trace_t *tr = &(trs[i]);
  				Dgx[i][local_gather_idx] = tr->gx;
  				Dgy[i][local_gather_idx] = tr->gy;
  				Dsx[i][local_gather_idx] = tr->sx;
  				Dsy[i][local_gather_idx] = tr->sy;
  				Dscalco[i][local_gather_idx] = tr->scalco;
  			
  				for(int t0  = 0;t0 < tr->ns;t0++)
  				{
  					Dcube[t0][i][local_gather_idx] = tr->data[t0];
  				} 			
  			}
  		}
  	}
  	const int size =mns;
  	const int groupsize = ngathers;
  	
  	std::cout<<"ngathers  = "<<ngathers<<"\n"; 
  	std::cout<<"ns = "<<mns<<"\n";
  	std::cout<<"c0 = "<<c0<<"c1 = "<<c1<<"\n";
	char namef[100];
	namef[0]=0;
	sprintf(namef,"time_%s.txt",argv[7]);
  	FILE * timeFile = fopen(namef,"a+");
	double t = mysecond();
  	myQueue.submit([&](handler &cgh){
  		auto deviceNtrs = Dntrs.get_access<access::read>(cgh);
  		auto deviceGx = gx.get_access <access::read>(cgh);
  		auto deviceGy = gy.get_access  <access::read> (cgh);
  		auto deviceSx = sx.get_access <access::read> (cgh);
  		auto deviceSy = sy.get_access <access::read> (cgh);
  		auto deviceScalco = scalco.get_access <access::read> (cgh);
  		auto deviceCube = cube.get_access <access::read> (cgh);
  		auto deviceCtr  = ctr.get_access <access::write> (cgh);
  		auto deviceTr0_dt=  tr0_dt.get_access<access::read> (cgh);
      auto deviceStk = stk.get_access<access:write::write> (cgh);
      auto deviceSem = sem.get_access<access:write::write>(cgh);
 cgh.parallel_for<class cmp>(mns, [=] (int index) {
     
     			int wgroup = index;
     			for(int j=0;j<ngathers;j++)
     			{
		        float caux = (c1-c0);
  				int gather_idx =j;
  				float best_sem = 0;
  				float best_stk = 0;
  				float best_ctr = 0;
  				int TAU = tau;
  				float dt = (float) deviceTr0_dt[make_id(gather_idx)] / 1000000;
  				float idt = 1 / dt;
  				int tau1 = MAX((int)((float)TAU * idt), 0);
  				int w = 2 * tau1 + 1;
  				int ntraces = deviceNtrs[make_id(gather_idx)];
  				for( int jc = 0;jc < nc;jc++)
  				{
  					float C = c0 + caux*jc/nc;
  					id<2> pos1(0,gather_idx);
  					float _s = get_trace_scalco(deviceScalco[pos1]);
  					float m0x = _s * (deviceGx[pos1]+deviceSx[pos1])*0.5;
  					float m0y = _s * (deviceGy[pos1]+deviceSy[pos1])*0.5;
  					float t0 = (float) wgroup * (float) dt;
  					int M = 0, skip = 0;
  					float stack_value = 0;
  					float den = 0;
  					int valid;
  					valid =0;

  					float num[16];
  					for (int ii=0; ii<w; ii++)
  					{
  						num[ii]=0;
  					}
  					for(int tr_idx =0;tr_idx<ntraces;tr_idx++)
  					{
  						id <2> index1(tr_idx,gather_idx);
  						float __s = get_trace_scalco(deviceScalco[index1]);
					  	 float t_gx = deviceGx[index1];
					  	 float t_gy = deviceGy[index1];
					  	 float t_sx = deviceSx[index1];
					   float t_sy = deviceSy[index1];
					   float mx = __s * (t_gx + t_sx) * 0.5;
					   float my = __s * (t_gy + t_sy) * 0.5;
					   float hx = __s * (t_gx - t_sx) * 0.5;
					   float hy = __s * (t_gy - t_sy) * 0.5;
					   float t = time2d(C, t0, m0x, m0y, mx, my, hx, hy);
					   int it = (int)(t * idt);
					   if (it-tau >= 0 && it+tau+1<mns) {
					   	for (int j = 0; j < w; j++) {
					   		int k = it + j - tau;         
					         float d1 = deviceCube[make_id(k,tr_idx,gather_idx)]; 
					         float d2 = deviceCube[make_id(k+1,tr_idx,gather_idx)];
					         float v = interpol_linear(k, k+1, d1, d2, t*idt + j - tau);
					         num[j]      += v;
					         den         += v*v;
					         stack_value += v;
					     }
					     M++; 
					 } else if (++skip == 2) {
					 	valid = 0;
					 }
					}
					float numerator = 0;
					for (int j = 0; j < w; j++) {
						numerator += num[j] * num[j];
			      den += denv[j];
			      stack_value += num[j];
					}
			      // compute semblance and stacked value
					stack_value = stack_value / (M*w);
					float s = (numerator / den) / M;
			      // -- semblance_2d for C -- End -- //
					if (valid)
						s= 0;
					if (s > best_sem) {
						best_sem = s;
						best_ctr = C;
						best_stk = stack_value;
					}

				}
			  
	     	   deviceCtr[make_id(wgroup,gather_idx)]=best_ctr;
           deviceStk[make_id(wgroup,gather_idx)]=best_stk;
           deviceSem[make_id(wgroup,gather_idx)]=best_sem;
 	     	}
	     			
			});
});
	
auto Dc_out = ctr.get_access<access::read,access::host_buffer> ();
double end = mysecond();
fprintf(timeFile,"%f\n",end-t);
		for(int i =0;i<ngathers;i++)
		 {
			for (int t0 = 0; t0 < mns; t0++) {
				float C = Dc_out[t0][i];
				ctr_t[i].data[t0] = C;
		}
			       su_fputtr(C_out, &ctr_t[i]);
		}
		printf("\n");

		return 0;
	}
