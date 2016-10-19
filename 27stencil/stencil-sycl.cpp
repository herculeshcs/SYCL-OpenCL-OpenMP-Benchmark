 /***************************************************************************
 *   Copyright (C) 2016 by Edson Borin, Hercules Cardoso da Silva
 *   Authors: Edson Borin - edson@ic.unicamp.br
 *	      Hercules Cardoso da Silva - hercules.cardoso.silva1@gmail.com
 *
 *    Parts of this code are based on the code available in the master's thesis "Investigation of the OpenCL SYCL Programming Model" by Angelos Trigkas.
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
#include <sys/time.h>
#include <CL/sycl.hpp>
#define ITERATIONS 1000
#define SZ 97
int init(double *v1, double * v2, int sz)
{
	int i;
	for(i =0;i<sz;i++)
	{
		v1[i]= 2;
		v2[i]= 3;
	}
}
double mysecond()
{
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
 using namespace cl::sycl;
int main (int argc, char ** argv)
{
	 double *  a3;
	 double *  a4;
	 double fac = 1.1;
	 int n = 10e6;
	 a3 = new double[n+1];
	 a4 = new double[n+1];

	 init(a3,a4,n);
	 int k,i,j;
	 int iter;
	 n=100;
	 char namef[100];
	 char namefcopy[100];
	 sprintf(namef,"timeSYCL_%s_%s.txt\n",argv[1],argv[2]);
	 sprintf(namefcopy,"timeCopy_%s_%s.txt\n",argv[1],argv[2]); 
	 FILE * timeFile = fopen(namef,"a+");
 	 FILE * timeFileCopy = fopen(namefcopy,"a+"); 
	 queue myQueue;
	 int sz = 100;
	 int SIZE= atoi(argv[1]);
	 buffer<double, 1> d_a0(SIZE +1);
	 buffer<double, 1> d_a1(SIZE +1);
	 auto h_a0 = d_a0.get_access<access::write,access::host_buffer>();
	 auto h_a1 = d_a1.get_access<access::write,access::host_buffer>();
	 for(int i =0;i< sz*sz*sz;i++)
	 {
	 	h_a0[i]=1;
	 	h_a1[i]=1;
	 }
	 double computeTime=0;
	 double t = mysecond();
	 myQueue.submit([&](handler &cgh){
	 auto a0 = d_a0.get_access<access::read_write>(cgh);
	 auto a3 = d_a1.get_access<access::read_write>(cgh);
	 sz =  cbrt(SIZE/2);
	 n = sz-2;
	 size_t max = n;
	 for (iter = 0; iter < ITERATIONS; iter++) {
	 double tc = mysecond();
	 cgh.parallel_for < class compute>({max,max,max},[=] (id<3> index) {
		  int i = index.get(0)+1;
		  int j = index.get(1)+1;
		  int k = index.get(2)+1;
				 a3[make_id(i*sz*sz+j*sz+k)] = (
				 a0[make_id(i*sz*sz+(j-1)*sz+k)] + a0[make_id(i*sz*sz+(j+1)*sz+k)] +
				 a0[make_id((i-1)*sz*sz+j*sz+k)] + a0[make_id((i+1)*sz*sz+j*sz+k)] +
				 a0[make_id((i-1)*sz*sz+(j-1)*sz+k)] + a0[make_id((i-1)*sz*sz+(j+1)*sz+k)] +
				 a0[make_id((i+1)*sz*sz+(j-1)*sz+k)] + a0[make_id((i+1)*sz*sz+(j+1)*sz+k)] +
				 a0[make_id(i*sz*sz+(j-1)*sz+(k-1))] + a0[make_id(i*sz*sz+(j+1)*sz+( k-1))] +
				 a0[make_id((i-1)*sz*sz+j*sz+(k-1))] + a0[make_id((i+1)*sz*sz+j*sz+( k-1))] +
				 a0[make_id((i-1)*sz*sz+(j-1)*sz+(k-1))] + a0[make_id((i-1)*sz*sz+(j+1)*sz+(k-1))] +
				 a0[make_id((i+1)*sz*sz+(j-1)*sz+(k-1))] + a0[make_id((i+1)*sz*sz+(j+1)*sz+(k-1))] +
				 a0[make_id(i*sz*sz+(j-1)*sz+(k+1))] + a0[make_id(i*sz*sz+(j+1)*sz+( k+1))] +
				 a0[make_id((i-1)*sz*sz+j*sz+(k+1))] + a0[make_id((i+1)*sz*sz+j*sz+( k+1))] +
				 a0[make_id((i-1)*sz*sz+(j-1)*sz+(k+1))] + a0[make_id((i-1)*sz*sz+(j+1)*sz+(k+1))] +
				 a0[make_id((i+1)*sz*sz+(j-1)*sz+(k+1))] + a0[make_id((i+1)*sz*sz+(j+1)*sz+(k+1))] +
				 a0[make_id(i*sz*sz+j*sz+(k-1))] + a0[make_id(i*sz*sz+j*sz+(k+1))]
			 ) * fac;
	 });
 	 computeTime+= mysecond() - tc; 
	 cgh.parallel_for <class copy> ({max*max*max},([=](int item)
	 {
		a0[make_id(item)] = a3[make_id(item)];
		}));
	  }
	 }); // Command group ends here
	double end = mysecond();
	FILE * f=fopen("OutputSYCL.txt","w");
	auto a0 = d_a0.get_access<access::read_write,access::host_buffer>();
	auto a1 = d_a1.get_access<access::read_write,access::host_buffer>();
	n=97;
	for (i = 1; i < n+1; i++) {
		for (j = 1; j < n+1; j++) {
			for (k = 1; k < n+1; k++) {
	 			 fprintf(f,"%f\n",a1[i*sz*sz+j*sz+k]);
 			}
		}
	 }

 fprintf(timeFile,"%f\n",computeTime);
 fprintf(timeFileCopy,"%f\n",(end -t)-computeTime);
}
