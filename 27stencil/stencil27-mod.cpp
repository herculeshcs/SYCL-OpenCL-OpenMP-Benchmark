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
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#define ITERATIONS 1000
// function for intilize vectors
int init(double *v1, double * v2, int sz)
{
	for(int i =0;i<sz;i++)
	{
		v1[i]=  1e-7;
		v2[i]=  1e-7;
	}
}
// funcition to measure time
double mySecond()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp,&tzp);
	return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
//Number of interation for interation loop
#define ITERATIONS 1000
int main(int argc, char * argv[] )
{
	if(argc != 3)
	{
		std::cout<<"Usage: "<<argv[0]<<" Size Numeber_of_threads"<<std::endl;
		std::cout<<"Example: "<<agv[0]<<" 1000 16"<<std::endl;
	}

	int  size = atoi(argv[1]);
	int th = atoi(argv[2]);
	double *  a0;
	double *  a1;
	double fac = 2.0;
	int sz =  cbrt(size/2);
	int n = sz-2;
	//alocate memory for a0 and a1
	a0 = (double*) malloc(sizeof(double)*size+1);
	a1 = (double*) malloc(sizeof(double)*size+1);
	// Initialize a0 and a1
	init(a0,a1,size);
	int k,i,j;
	int iter;
	char fileNameCopy[100];
	char fileName[100];
	sprintf(fileName,"timeOpenMP_%s_%s",argv[1],argv[2]);
	sprintf(fileNameCopy,"timeOpenMPCopy_%s_%s",argv[1],argv[2]).
	// Open files for appending time
	FILE * timeFileCopy = fopen(fileNameCopy,"a+");
	FILE * timeFile = fopen(fileName,"a+");
	double t = mySecond();
		double computeTime = 0;
	//start interation loop
	for (iter = 0; iter < ITERATIONS; iter++) {
	 double t1 = mySecond();
	 #pragma omp parallel for  private(i,j,k) schedule(dynamic) shared(a1) num_threads(th)
	 //Compute loop 
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
	 computeTime += mySecond() - t1;
	 #pragma omp parallel for  private(i,j,k) schedule(dynamic) shared(a1) num_threads(th)
	//Copy loop
	for (i = 1; i < n+1; i++) {
		for (j = 1; j < n+1; j++) {
			for (k = 1; k < n+1; k++) {
		 	a0[i*sz*sz+j*sz+k] = a1[i*sz*sz+j*sz+k];
	 			}
			}
		 }
	 } //end iteration loop
	 double end = mySecond();
	 fprintf(timeFile,"%f\n",computeTime);
	 fprintf(timeFileCopy,"%f\n",(end -t) -computeTime);
	 std::cout<<"Compute Time = "<<computeTime<<std::endl;
	 std::cout<<"Copy Time = "<<(end - t) -computeTime<<std::endl;
	 fclose(timeFile);
		return 0;
}
