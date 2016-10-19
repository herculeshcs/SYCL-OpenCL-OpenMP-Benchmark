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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#define ITERATIONS 1000
#define SIZE 10e6
#define SZ 97
int init(double *v1, double * v2, int sz)
{
	int i;
	for(i =0;i<sz;i++)
	{
		v1[i]=  1e-7;
		v2[i]=  1e-7;
	}
}
double mysecond()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp,&tzp);
	return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

int main(void )
{

	double *  a0;
	double *  a1;
	double fac = 2.0;
	int n = SIZE;
	int sz  =  SZ;
	a0 = (double*) malloc(sizeof(double)*n+1);
	a1 = (double*) malloc(sizeof(double)*n+1);
	init(a0,a1,n);
	int k,i,j;
	int iter;
	FILE * timeFile = fopen("TimeOpenMP.txt","a+");
	sz =  cbrt(1000000/2);
	printf("sZ = %d\n",sz);
	n = sz-2;
	double t = mysecond();
	for (iter = 0; iter < ITERATIONS; iter++) {
		#pragma omp parallel for  default(none) private(i,j,k) shared(a0, a1, fac, n, sz) collapse(2)
		for (i = 1; i < n+1; i++) {
		 	for (j = 1; j < n+1; j++) {
			 	for (k = 1; k < n+1; k++) {
					 a1[i*sz*sz+j*sz+k] = (
					 a0[i*sz*sz+(j-1)*sz+k] + a0[i*sz*sz+(j+1)*sz+k] +
					 a0[(i-1)*sz*sz+j*sz+k] + a0[(i+1)*sz*sz+j*sz+k] +
					 a0[(i-1)*sz*sz+(j-1)*sz+k] + a0[(i-1)*sz*sz+(j+1)*sz+k] +
					 a0[(i+1)*sz*sz+(j-1)*sz+k] + a0[(i+1)*sz*sz+(j+1)*sz+k] +
					 a0[i*sz*sz+(j-1)*sz+(k-1)] + a0[i*sz*sz+(j+1)*sz+(k-1)] +
					 a0[(i-1)*sz*sz+j*sz+(k-1)] + a0[(i+1)*sz*sz+j*sz+(k-1)] +
					 a0[(i-1)*sz*sz+(j-1)*sz+(k-1)] + a0[(i-1)*sz*sz+(j+1)*sz+(k-1)] +
					 a0[(i+1)*sz*sz+(j-1)*sz+(k-1)] + a0[(i+1)*sz*sz+(j+1)*sz+(k-1)] +
					 a0[i*sz*sz+(j-1)*sz+(k+1)] + a0[i*sz*sz+(j+1)*sz+(k+1)] +
					 a0[(i-1)*sz*sz+j*sz+(k+1)] + a0[(i+1)*sz*sz+j*sz+(k+1)] +
					 a0[(i-1)*sz*sz+(j-1)*sz+(k+1)] + a0[(i-1)*sz*sz+(j+1)*sz+(k+1)] +
					 a0[(i+1)*sz*sz+(j-1)*sz+(k+1)] + a0[(i+1)*sz*sz+(j+1)*sz+(k+1)] +
					 a0[i*sz*sz+j*sz+(k-1)] + a0[i*sz*sz+j*sz+(k+1)]
					 )*fac;
				}
			}
		 }
	 #pragma omp parallel for default(none) private(i,j,k) shared(a0, a1, n, sz) collapse(3)
	for (i = 1; i < n+1; i++) {
		for (j = 1; j < n+1; j++) {
			for (k = 1; k < n+1; k++) {
		 	a0[i*sz*sz+j*sz+k] = a1[i*sz*sz+j*sz+k];
	 			}
			}
		 }
	 } /* send iteration loop */
	 double end = mysecond();
	 fprintf(timeFile,"%f\n",end-t);
		FILE * f=fopen("OutputOpenmp.txt","w");
		for (i = 1; i < n+1; i++) {
			for (j = 1; j < n+1; j++) {
				for (k = 1; k < n+1; k++) {
		 			 fprintf(f,"%.10lf\n",a1[i*sz*sz+j*sz+k]);
	 			}
			}
		 }
		return 0;
}
