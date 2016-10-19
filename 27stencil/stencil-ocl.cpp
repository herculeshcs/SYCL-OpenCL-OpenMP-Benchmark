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
#include <CL/cl.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#define ITERATIONS 1000
#define MAX_SOURCE_SIZE 4000
void print_opencl_error(FILE* fh, cl_int err)
{
#define PRINT_ERR(code) case code : fprintf(fh, #code); break
	switch(err) {
		PRINT_ERR(CL_INVALID_PROGRAM);
		PRINT_ERR(CL_INVALID_VALUE);
		PRINT_ERR(CL_INVALID_DEVICE);
		PRINT_ERR(CL_INVALID_BINARY);
		PRINT_ERR(CL_INVALID_BUILD_OPTIONS);
		PRINT_ERR(CL_INVALID_OPERATION);
		PRINT_ERR(CL_COMPILER_NOT_AVAILABLE);
		PRINT_ERR(CL_BUILD_PROGRAM_FAILURE);
		PRINT_ERR(CL_OUT_OF_RESOURCES);
		PRINT_ERR(CL_OUT_OF_HOST_MEMORY);
	default:
		fprintf(fh, "unknown code");
	break;
  };
  return;
}
double mysecond()
{
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
int main (int argc, char **argv)
{
	cl_platform_id platform;
	cl_device_id device;
	cl_context context;
	cl_command_queue command_queue;
	cl_program program;
	cl_kernel kernel_stencil, kernel_copy;
	cl_kernel kernel_compute;
	cl_mem d_a0;
	cl_mem d_a1;
	cl_uint num_platforms;
	cl_uint num_devices;
	cl_int err;  

	clGetPlatformIDs(1, &platform, &num_platforms);
	clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, &num_devices);
	context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
	if(err != CL_SUCCESS)
	{
		fprintf(stderr,"ERROR: clCreaterContext(...) return error code: %d\n",err);
	}
	command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);
	if(err != CL_SUCCESS)
	{
		fprintf(stderr,"ERROR: clCreatCommandQueue(...) return error code: %d\n",err);
		exit(0);
	} 
	char namef[100];
	char nameCopyf[100]; 
	sprintf(namef, "timeCL_%s_%s.txt\n",argv[1],argv[2]); 
	sprintf(nameCopyf,"timeCLCopy_%s_%s.txt\n",argv[1],argv[2]);
	FILE * timeFile = fopen(namef,"a+");
	FILE * timeFileCopy = fopen(nameCopyf,"a+");
	FILE * fp = fopen("kernel.cl", "r");
	if(fp==NULL)
	{
		fprintf(stderr,"ERROR: Not found kernel.cl\n");
		exit(0);
	}
	char * source_str = new char[MAX_SOURCE_SIZE];
	size_t source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
	program = clCreateProgramWithSource(context, 1, (const char **) &source_str, (const size_t *)&source_size, &err);
	if(err != CL_SUCCESS)
	{
		fprintf(stderr,"ERROR: clCreateProgramWithSource return error code: %d\n",err);
		exit(0);
	}
	err= clBuildProgram(program, 1, &device, NULL, NULL, NULL);
	if (err != CL_SUCCESS) {
		size_t len;
		fprintf(stderr,"Error: Failed to build program executable (");
		print_opencl_error(stderr, err);
		fprintf(stderr,")\n");
		err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
		fprintf(stderr,"Error size = %d\n", (int) len);
		char buffer[4000];
		err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, (len+10), buffer, &len);
		if (err != CL_SUCCESS) {
			fprintf(stderr,"Error: Could not read build log. clGetProgramBuildInfo(...) returned ");
			print_opencl_error(stderr, err);
			fprintf(stderr,"\n");
	    }
	    else {
	      fprintf(stderr,"-- BUILD LOG ----------\n%s\n-----------------------\n", buffer);
	    }
	    exit(1);
	  }  
	kernel_stencil = clCreateKernel(program, "stencil", &err);
	if(err != CL_SUCESS)
	{
		fprintf(stderr,"ERROR: clCreateKernel return error code: %d\n",err);
	}  
	  kernel_copy = clCreateKernel(program, "copy", &err);
	if(err != CL_SUCESS)
	{
		fprintf(stderr,"ERROR: clCreateKernel return error code: %d\n",err);
	}
	int Size = atoi(argv[1]);
	d_a0 = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(double)*Size+1, NULL, NULL);
	d_a1 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(double)*Size+1, NULL, NULL);
	int ma = Size;
	double * a0_init = new double[Size+1];
	double * a1 = new double[Size+1];
	for(int i =0;i< ma;i++)
	{
		a0_init[i]=1e-7;
	}
	clEnqueueWriteBuffer(command_queue, d_a0, CL_TRUE, 0, sizeof(double)*Size, (void*)a0_init, 0, NULL, NULL);
	int iter;
	kernel_compute = kernel_stencil;
	double fac = 2;
	cl_ulong  total_time = 0;
	cl_ulong  compute_time = 0;
	for (iter = 0; iter < ITERATIONS; iter++) {
		clSetKernelArg(kernel_compute, 0, sizeof(cl_mem), &d_a0);
		clSetKernelArg(kernel_compute, 1, sizeof(cl_mem), &d_a1);
		clSetKernelArg(kernel_compute, 2, sizeof(cl_mem), &fac);
		int  sz =  cbrt(Size/2);
		int n = sz-2;
		cl_event event;
		size_t global[3] = { n, n, n };
		err=clEnqueueNDRangeKernel(command_queue, kernel_compute, 3, NULL, global, NULL, 0, NULL,&event);
		clWaitForEvents(1, &event);
		cl_ulong time_start, time_end;
		clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
		clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
		total_time+= time_end-time_start;
		compute_time += time_end-time_start;
		if (err) {
			printf("Error: clEnqueueNDRangeKernel(..) returned error code %d (", err);
			print_opencl_error(stderr, err);
			fprintf(stderr, ")\n");  	       
	     }
	 
	clFinish(command_queue);
	clSetKernelArg(kernel_copy, 0, sizeof(cl_mem), &d_a0);
	clSetKernelArg(kernel_copy, 1, sizeof(cl_mem), &d_a1);
	clEnqueueNDRangeKernel(command_queue, kernel_copy, 3, NULL, global, NULL, 0, NULL, &event);
	clWaitForEvents(1, &event);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	total_time += time_end-time_start;

	clFinish(command_queue);
		if (err) {
			printf(" COPY Error: clEnqueueNDRangeKernel(..) returned error code %d (", err);
			print_opencl_error(stderr, err);
			fprintf(stderr, ")\n");  
			exit (1);
	     }
	 }
	clEnqueueReadBuffer(command_queue, d_a0, CL_TRUE, 0, sizeof(double)*Size, a1, 0, NULL, NULL);
	fprintf(timeFile,"%f\n",compute_time/1000000000.0);
	fprintf(timeFileCopy,"%f\n",(total_time-compute_time)/1000000000.0);
}
