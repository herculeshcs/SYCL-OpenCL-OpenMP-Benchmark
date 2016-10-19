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
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#define ITERATIONS 1000
#define MITERATIONS 100
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

int main(int argc, char * argv[] )
{
	if(argc != 3)
	{
		std::cout<<"Usage: "<<argv[0]<<" Size Numeber_of_threads"<<std::endl;
		std::cout<<"Example: "<<agv[0]<<" 1000 16"<<std::endl;
	}

	int Size = atoi(argv[1]);
	int th = atoi(argv[2]);
	printf("size = %d ,threads =  %d \n", Size,th);
	double *  a0;
	double *  a1;
	double fac = 2.0;
	int sz =  cbrt(Size/2);
	printf("sZ = %d\n",sz);
	int n = sz-2;
	a0 = (double*) malloc(sizeof(double)*Size+1);
	a1 = (double*) malloc(sizeof(double)*Size+1);
	init(a0,a1,Size);
	int k,i,j;
	int iter;

	char filename[100];
	filename[0]=0;
	strcat(filename,"TimeOpenMPNormal_");
	strcat(filename,argv[1]);
	strcat(filename,"_");
	strcat(filename,argv[2]);

	char filenameCopy[100];
	filenameCopy[0]=0;
	strcat(filenameCopy,"timeOpenMPNormalCopy_");
	strcat(filenameCopy,argv[1]);
	strcat(filenameCopy,"_");
	strcat(filenameCopy,argv[2]);

	FILE * timeFileCopy = fopen(filenameCopy,"a+");
	FILE * timeFile = fopen(filename,"a+");
	double t = mysecond();
	double t1 = mysecond();
	double parc = 0;
	omp_set_dynamic(0);
	for (iter = 0; iter < ITERATIONS; iter++) {
		double t1 = mysecond();
	 #pragma omp parallel for  default(none) private(i,j,k) shared(a0, a1, fac, n, sz) collapse(3)  num_threads(th) 
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
	parc += mysecond() - t1;
	#pragma omp parallel for default(none) private(i,j,k) shared(a0, a1, n, sz) collapse(3) num_threads(th)
	for (i = 1; i < n+1; i++) {
		for (j = 1; j < n+1; j++) {
			for (k = 1; k < n+1; k++) {
		 	a0[i*sz*sz+j*sz+k] = a1[i*sz*sz+j*sz+k];
	 			}
			}
		 }
	 } /* end iteration loop */
	 double end = mysecond();
	 fprintf(timeFile,"%f\n",parc);
	 fprintf(timeFileCopy,"%f\n",(end -t) -parc);
	 fclose(timeFile);
	 fclode(timeFile);
		return 0;
}
