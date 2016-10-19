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

#include "su.h"
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

/*------------------------------------------------*/
/* Code to read the wall clock time.              */
#include <sys/time.h>
#include <time.h>
double mysecond()
{
 	struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

#define TIMER(name)				\
  double total_ ## name = 0.0;			\
  double lap_ ## name

#define TIMEP(name,str)				\
  printf("%20s: %0.12f\n", str, total_ ## name)

#define TIMEPP(totaln,name,str)			\
  printf("%20s: %0.12f : %.3f%%\n", str, total_ ## name, 100*(total_ ## name / total_ ## totaln))
#define TIMERSUM(name,name1)\
total_ ## name += total_ ## name1
#define LAP_START(name) lap_ ## name = mysecond()
#define LAP_END(name) total_ ## name += mysecond() - lap_ ## name

/*------------------------------------------------*/

int verbosity_level = 0;
#define VERBOSE(level,msg) do {if(verbosity_level >= level) { fprintf(stderr,msg); } }while(0)

#define SQR(x) ((x)*(x))
#ifndef MAX
# define MAX(a, b) ((a)>(b)?(a):(b))
#endif
#ifndef MIN
# define MIN(a, b) ((a)<(b)?(a):(b))
#endif

void su_init(su_trace_t *tr)
{
	if (tr->ns <= 0)
		return;
	tr->data = malloc(sizeof(float) * tr->ns);
}

void su_free(su_trace_t *tr)
{
	free(tr->data);
}

int su_fgettr(FILE *file, su_trace_t *tr)
{
	size_t size = fread(tr, SU_HEADER_SIZE, 1, file);
	if (size < 1)
		return 0;
	tr->data = malloc(tr->ns * sizeof(float));
	size = fread(tr->data, sizeof(float), tr->ns, file);
	if (size < tr->ns) {
		su_free(tr);
		return 0;
	}
	return 1;
}

int su_gettr(su_trace_t *tr)
{
	return su_fgettr(stdin, tr);
}

int su_fputtr(FILE *file, su_trace_t *tr)
{
	size_t size = SU_HEADER_SIZE;
	if (fwrite(tr, size, 1, file) < 1)
		return 0;
	if (fwrite(tr->data, sizeof(float), tr->ns, file) < tr->ns)
		return 0;
	return 1;
}

int su_puttr(su_trace_t *tr)
{
	return su_fputtr(stdout, tr);
}

int su_get_cdp(su_trace_t *tr)
{
	return tr->cdp;
}

void su_get_source(su_trace_t *tr, float *sx, float *sy)
{
  float s;
  if (tr->scalco == 0)
    s = 1;
  else if (tr->scalco > 0)
    s = tr->scalco;
  else 
    s = 1.0f / tr->scalco;
  *sx = s * tr->sx;
  *sy = s * tr->sy;
}

void su_get_receiver(su_trace_t *tr, float *gx, float *gy)
{
  float s;
  if (tr->scalco == 0)
    s = 1;
  else if (tr->scalco > 0)
    s = tr->scalco;
  else 
    s = 1.0f / tr->scalco;
  *gx = s * tr->gx;
  *gy = s * tr->gy;
}

#include <vector.h>
#include <su.h>

typedef struct aperture aperture_t;

struct aperture {
	float ap_m, ap_h, ap_t;
	vector_t(su_trace_t*) traces;
};

#include <math.h>
#include <string.h>
#include <vector.h>

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <vector.h>
#include <uthash.h>
#include <log.h>
#include <su.h>

struct cdp_traces {
	int cdp;
	vector_t(su_trace_t) traces;
	UT_hash_handle hh;
};

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

int read_txt_file(char* program_src, char* filename, int max)
{
  FILE* fh = fopen(filename, "r");
  if (!fh) {
    fprintf(stderr,"fopen returned error!\n");
    return 1;
  }

  size_t sz = fread((void *) program_src, 1, max, fh);
  if (sz <= 0) {
    fprintf(stderr,"fread returned error!\n");
    return 2;
  }

  program_src[sz] = '\0';

  fclose(fh);

  return 0;
}
#ifdef MAC
#include <OpenCL/opencl.h>
#endif
#ifdef LINUX
#include <CL/cl.h>
#endif

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

static float get_scalco(su_trace_t *tr)
{
	if (tr->scalco == 0)
		return 1;
	if (tr->scalco > 0)
		return tr->scalco;
	return 1.0f / tr->scalco;
}

void su_get_halfoffset(su_trace_t *tr, float *hx, float *hy)
{
	float s = get_scalco(tr);
	*hx = s * (tr->gx - tr->sx) * 0.5;
	*hy = s * (tr->gy - tr->sy) * 0.5;
}

int compare_gathers_by_ntraces(const void* a, const void* b) {
  struct cdp_traces **A = (struct cdp_traces**) a;
  struct cdp_traces **B = (struct cdp_traces**) b;
  return (*B)->traces.len - (*A)->traces.len;
}

int main(int argc, char *argv[])
{
	
 FILE *  cTxt = fopen("c_output","w");
  int totalCoisas = 0;
  TIMER(total);
  TIMER(read_traces);
  TIMER(sort_traces);
  TIMER(ocl_setup);
  TIMER(cmp_processing);
  TIMER(write_bf);
  TIMER(read_bf);
  TIMER(kernel_exec);
  TIMER(memcp);
  TIMER(memcpout);
  TIMER(kernel_time);
  LAP_START(total);
  
  if (argc != 10) {
    fprintf(stderr, "Usage: %s C0 C1 NC APH TAU INPUT VERBOSITY GPU\n", argv[0]);
    exit(1);
  }

  double wpt = mysecond();
  
  char* c0_str = argv[1];
  float c0 = atof(argv[1]);
  char* c1_str = argv[2];
  float c1 = atof(argv[2]);
  char* nc_str = argv[3];
  int nc = atoi(argv[3]);
  int aph = atoi(argv[4]);
  char* tau_str = argv[5];
  float tau = strtof(argv[5], NULL);
  char *path = argv[6];
  char * suf= argv[8];
  char timeFile[100];
  sprintf(timeFile,"time_%s_%s",suf,argv[9]);
  if (suf == NULL)
	{
			suf = (char*) malloc(sizeof(char)*1);
		   suf[0]=0;
	}  
  FILE *  time = fopen(timeFile,"a+");
  FILE *  mem = fopen("mem","a+");
  FILE * sem_s = fopen("sem","a+");
  verbosity_level = atoi(argv[7]);
#ifndef NG
  #define NG 64
  int ng = 64;
#else
  #warning "NG =  " NG
  int ng = NG;
#endif

  VERBOSE(2,"fopen()\n");

  FILE *fp = fopen(path, "r");
  
  /* Read traces and group them by CDP. ( */
  VERBOSE(2,"Reading traces\n");
  LAP_START(read_traces);
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
      val = malloc(sizeof(struct cdp_traces));
      val->cdp = cdp;
      vector_init(val->traces);
      HASH_ADD_INT(cdp_traces, cdp, val);
    }
    vector_push(val->traces, tr);
  }

  fclose (fp);

  /* Verify and set the number of samples per trace: ns */
  if (min_ns != max_ns) {
    fprintf(stderr,"ERROR: There are traces with different number of samples. Min = %d, Max = %d\n", min_ns, max_ns);
    exit(1);
  }
	
  int ns = min_ns;
  fprintf(stderr, "ns = %d\n", ns);
  fprintf(stderr, "ng = %d\n", ng);
  LAP_END(read_traces);
 int v[100];
 int vi=0;
 for(vi=0;vi< 100;vi++)
 {
  v[vi]=0;
 }

  /* Sort gathers by number of traces. This should improve load balance on GPU. */
  int total_gathers = 0;

  struct cdp_traces** gathers = malloc(sizeof(struct cdp_traces*) * HASH_COUNT(cdp_traces));
  fprintf(stderr,"HASH_COUNT = %d\n", HASH_COUNT(cdp_traces));
  int max_ntrs = 1; // maximum traces per gather
  struct cdp_traces *it;
  int numberOfCDPs=0;
  for (it = cdp_traces; it; it = it->hh.next) {
    gathers[total_gathers++] = it;
    numberOfCDPs++;
    v[it->traces.len]++;
    if (max_ntrs < it->traces.len) 
      max_ntrs = it->traces.len;
          
  }
  int aux = atoi(nc_str);

  int numberOfSemblances=0;
  int ttraces=0;
  for(vi=0;vi < 100;vi++)
  {
    numberOfSemblances+=  vi*v[vi]*ns*aux;
    ttraces += vi*v[vi];
  }
  /* ------- OpenCL basic stuff -------- */
  LAP_START(ocl_setup);
  cl_int err;     // error code returned from api calls
  
  //size_t global;  // global domain size for our calculation
  //size_t local;   // local domain size for our calculation
  
  cl_platform_id cpPlatform; // OpenCL platform
  cl_device_id device_id;    // compute device id
  cl_context context;        // compute context
  cl_command_queue commands; // compute command queue
  cl_event prof_event;

  // Connect to a compute device
  VERBOSE(2,"clGetPlatformIDs\n");
  	cl_uint numberOfPlataforms;
  err = clGetPlatformIDs(1, &cpPlatform, &numberOfPlataforms);
  if (err != CL_SUCCESS) {
    fprintf(stderr, "Error: Failed to find a platform!\n");
    return EXIT_FAILURE;
  }
  
  // Get a device of the appropriate type
  VERBOSE(2,"clGetDeviceIDs\n");
   	#ifdef XEONPHI
		 VERBOSE(2,"Use Acceleretor\n");
  		err = clGetDeviceIDs(cpPlatform,CL_DEVICE_TYPE_ACCELERATOR, 1, &device_id, NULL);
  	 #endif
  	 #ifdef GPU 
        VERBOSE(2,"Use GPU\n");
  	   err = clGetDeviceIDs(cpPlatform,CL_DEVICE_TYPE_GPU, 1, &device_id, NULL);
  	 #endif
  	 #ifdef CPU
	VERBOSE(2,"Use CPU\n");
  	 err = clGetDeviceIDs(cpPlatform,CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
  	 #endif
  	   
  if (err != CL_SUCCESS) {
    fprintf(stderr, "Error: Failed to create a device group!, %d\n", err);
    return EXIT_FAILURE;
  }
  
  // Create a compute context
  VERBOSE(2,"clCreateContext\n");
  context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
  if (!context) {
    fprintf(stderr, "Error: Failed to create a compute context!\n");
    return EXIT_FAILURE;
  }
  
  // Create a command commands
  VERBOSE(2,"clCreateCommandQueue\n");
  commands = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);
  if (!commands) {
    fprintf(stderr, "Error: Failed to create a command queue!\n");
    return EXIT_FAILURE;
  }

  /* ----------------------------- */

  /* ----- OpenCL Generate OpenCL Kernel!!! ------ */

  /* -- Creating Program Objects -- */
  char program_src[64*1024]; // TODO: 64k shoulbe be enough -- Revise it in case compute_gather.cl is changed
	// Define kernel constants using opencl preprocessor definitions
  char options[1024];
  #ifndef DEBUG
  sprintf(options,"-D C0=(float)%s -D C1=(float)%s -D NC=%s -D TAU=(float)%s -D MAX_NTRS=%d -D MAXW=%d -D NG=%d", c0_str, c1_str, nc_str, tau_str, max_ntrs, 16 /*maxw*/,NG);
  // Read source into program_src.
  #endif
  #ifdef DEBUG
  sprintf(options,"-D C0=(float)%s -D C1=(float)%s -D NC=%s -D TAU=(float)%s -D MAX_NTRS=%d -D MAXW=%d -D NG=%d -g -s \"/home/hercules_silva/Mestrado/cmp-tool-clean/compute_gather.cl\"		", c0_str, c1_str, nc_str, tau_str, max_ntrs, 16 /*maxw*/,NG);
  #endif	
  VERBOSE(2,"read_txt_file\n");
  if (read_txt_file(program_src, "compute_gather.cl", 64*1024)) {
    fprintf(stderr, "Error: could not rea	d \"compute_gather.cl\" file\n");
    return EXIT_FAILURE;
  }

  // Create the compute program from the source buffer
    strcat(program_src,"//");
    strcat(program_src,options);
  char* program_ptr = program_src;

  VERBOSE(2,"clCreateProgramWithSource\n");
  cl_program program = clCreateProgramWithSource(context, 1, (const char **) &program_ptr, NULL, &err);
  if (!program) {
    fprintf(stderr, "Error: Failed to create opencl program (err = %d)!\n", err);
    return EXIT_FAILURE;
  }
  printf("OpenCL compiler options: %s\n", options);

  // Build the program executable
  VERBOSE(2,"clBuildProgram\n");
  printf("%s\n",options);
  err = clBuildProgram(program, 0, NULL, options, NULL, NULL);
  if (err != CL_SUCCESS) {
    size_t len;
    fprintf(stderr,"Error: Failed to build program executable (");
    print_opencl_error(stderr, err);
    fprintf(stderr,")\n");
    err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
    fprintf(stderr,"Error size = %d\n", (int) len);
    char* buffer = malloc(sizeof(char)*(len+10));
    err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, (len+10), buffer, &len);
    if (err != CL_SUCCESS) {
      fprintf(stderr,"Error: Could not read build log. clGetProgramBuildInfo(...) returned ");
      print_opencl_error(stderr, err);
      fprintf(stderr,"\n");
    }
    else {
      fprintf(stderr,"-- BUILD LOG ----------\n%s\n-----------------------\n", buffer);
    }
    free(buffer);	   
    exit(1);
  }  

  // Create the compute kernel in the program
  cl_kernel kernel;          // compute kernel
  VERBOSE(2,"clCreateKernel\n");
  LAP_START(kernel_time);
  kernel = clCreateKernel(program, "compute_gather", &err);
  LAP_END(kernel_time);
  if (!kernel || err != CL_SUCCESS) {
    fprintf(stderr,"Error: Failed to create compute kernel. clCreateKernel(...) returned (");
    print_opencl_error(stderr, err);
    fprintf(stderr,")\n");
    exit(1);
  }


  VERBOSE(2,"Checking maximum work-group size.\n");
  // Get the maximum work group size for executing the kernel on the device
  size_t max_WG_sz;
  err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(max_WG_sz), &max_WG_sz, NULL);
  if (err != CL_SUCCESS) {
    fprintf(stderr,"Error: Failed to retrieve kernel work group info (CL_KERNEL_WORK_GROUP_SIZE)!\n");
    exit(1);
  }
  if (ng <= max_WG_sz) {
	    fprintf(stderr,"Maximum workgroup size (CL_KERNEL_WORK_GROUP_SIZE) = %zu -- number of gathers per batch = %d\n", max_WG_sz, ng);
  } 
  else {
    fprintf(stderr,"Error: maximum workgroup size (CL_KERNEL_WORK_GROUP_SIZE) = %zu < number of gathers per batch = %d\n", max_WG_sz, ng);
    exit(1);
  }

  /* --- Device memory vectors --- */
  VERBOSE(2,"Allocating memory on device\n");

#define GET_2D_DATA(data,i,j,szj) data[i*szj + j]
#define GET_3D_DATA(data,i,j,k,szj,szk) data[i*(szj*szk) + j*szk + k]

int pos3d(int i, int j, int k, int szj,int szk)
{
	return i*(szj*szk) + j*szk + k;
	
}

#define CREATE_INPUT_BUFFER(name,type,sz)				\
    cl_mem name = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(type) * sz, NULL, NULL); \
    if (!name) {							\
      fprintf(stderr, "Error when creating input buffer named %s\n", #name); \
    }									\
    type* _ ## name = malloc(sizeof(type) * sz);				\
    if (! _ ## name) {							\
      fprintf(stderr, "Error when creating host buffer named _%s\n", #name); \
    }									\
    do {} while(0)
 
#define CREATE_OUTPUT_BUFFER(name,type,sz)				\
    cl_mem name = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(type) * sz, NULL, NULL); \
    if (!name) {							\
      fprintf(stderr, "Error when creating output buffer named %s\n", #name); \
    }									\
    type* _ ## name = malloc(sizeof(type) * sz);				\
    if (!_ ## name) {							\
      fprintf(stderr, "Error when creating host buffer named _%s\n", #name); \
    }									\
    do {} while(0)
  
  CREATE_INPUT_BUFFER(cube, float, ns*max_ntrs*ng);
  CREATE_INPUT_BUFFER(gx, int, max_ntrs*ng);
  CREATE_INPUT_BUFFER(gy, int, max_ntrs*ng);
  CREATE_INPUT_BUFFER(sx, int, max_ntrs*ng);
  CREATE_INPUT_BUFFER(sy, int, max_ntrs*ng);
  CREATE_INPUT_BUFFER(scalco, short, max_ntrs*ng);
  CREATE_INPUT_BUFFER(ntrs, int, ng);
  CREATE_INPUT_BUFFER(tr0_dt, short, ng);
  CREATE_OUTPUT_BUFFER(ctr, float, ns*ng);
  CREATE_OUTPUT_BUFFER(stk, float, ns*ng);
  CREATE_OUTPUT_BUFFER(sem, float, ns*ng);
  /* ----------------------------- */

  /* --- Create output files --- */
  VERBOSE(2,"Create output files\n");
  	char cfile[100];
  	char sfile[100];
  	char stackfile[100];
 	sprintf(cfile,"C%s.su",suf);
 	sprintf(sfile,"cmp.coher%s.su",suf);
  	sprintf(stackfile,"cmp.stack%s.su",suf);
   FILE *C_out = fopen(cfile, "w");
   FILE *S_out = fopen(sfile, "w");
   FILE *stack = fopen(stackfile, "w");
  LAP_END(ocl_setup);
  LAP_START(cmp_processing);
   su_trace_t ctr_t[ng];
   su_trace_t str_t[ng];
   su_trace_t stk_t[ng];
int global_gather_idx;
	int conta=0;
   for ( global_gather_idx = 0; global_gather_idx < total_gathers;) {
     /* Assemble a set of NG gathers (CDP panels) to send to the GPU. */
     VERBOSE(2,"Grouping gathers to send to GPU\n");
     int ngathers = MIN(ng, total_gathers - global_gather_idx);
     int local_gather_idx;
     for (local_gather_idx = 0; local_gather_idx < ngathers; local_gather_idx++, global_gather_idx++) {


	 su_trace_t *trs  = gathers[global_gather_idx]->traces.a;
	 int         ntrs = gathers[global_gather_idx]->traces.len;

	 /* Setup headers of output data. */ 
	 memcpy(&ctr_t[local_gather_idx], &trs[0], SU_HEADER_SIZE);
	 ctr_t[local_gather_idx].offset = 0;
	 ctr_t[local_gather_idx].sx = ctr_t[local_gather_idx].gx = (trs[0].sx + trs[0].gx) / 2;
	 ctr_t[local_gather_idx].sy = ctr_t[local_gather_idx].gy = (trs[0].sy + trs[0].gy) / 2;
   	 memcpy(&str_t[local_gather_idx], &ctr_t[local_gather_idx], SU_HEADER_SIZE);
	 memcpy(&stk_t[local_gather_idx], &ctr_t[local_gather_idx], SU_HEADER_SIZE);
	 su_init(&ctr_t[local_gather_idx]);
	 su_init(&str_t[local_gather_idx]);
	 su_init(&stk_t[local_gather_idx]);

	 /* Init kernel input data. */
	 _ntrs[local_gather_idx] = ntrs;
	 _tr0_dt[local_gather_idx] = trs[0].dt;
	 int i;
	 for (i = 0; i < ntrs; i++) {
	   su_trace_t *tr = &(trs[i]);
	   GET_2D_DATA(_gx,i,local_gather_idx,ngathers) = tr->gx;
	   GET_2D_DATA(_gy,i,local_gather_idx,ngathers) = tr->gy;
	   GET_2D_DATA(_sx,i,local_gather_idx,ngathers) = tr->sx;
	   GET_2D_DATA(_sy,i,local_gather_idx,ngathers) = tr->sy;
	   GET_2D_DATA(_scalco,i,local_gather_idx,ngathers) = tr->scalco;
	   int t0;
	   for (t0=0; t0<tr->ns; t0++)
	   {
	     if(pos3d(t0,i,local_gather_idx,max_ntrs,ngathers)>=ns*max_ntrs*ng)
	     {
	     		exit(0);
	     }
	     GET_3D_DATA(_cube,t0,i,local_gather_idx,max_ntrs,ngathers) = tr->data[t0];
		}	     
	     	
	 }
     }
     /* Send gathers to GPU and calngl kernel. */

#define ENQUEUE_WRITE_BUFFER(name,type,count)				\
     err = clEnqueueWriteBuffer(commands, name, CL_TRUE, 0, sizeof(type) * count, _ ## name, 0, NULL, NULL); \
     if (err != CL_SUCCESS) {						\
       fprintf(stderr, "Error: Failed to enqueue command to write buffer %s!\n", #name); \
       exit(1);								\
     } do {} while(0)

     VERBOSE(3,"Enqueuing commands to write input data to GPU\n");
     
     LAP_START(memcp);
     ENQUEUE_WRITE_BUFFER(cube, float, ns*max_ntrs*ngathers);
     ENQUEUE_WRITE_BUFFER(gx, int, max_ntrs*ngathers);
     ENQUEUE_WRITE_BUFFER(gy, int, max_ntrs*ngathers);
     ENQUEUE_WRITE_BUFFER(sx, int, max_ntrs*ngathers);
     ENQUEUE_WRITE_BUFFER(sy, int, max_ntrs*ngathers);
     ENQUEUE_WRITE_BUFFER(scalco, short, max_ntrs*ngathers);
     ENQUEUE_WRITE_BUFFER(ntrs, int, ngathers);
     ENQUEUE_WRITE_BUFFER(tr0_dt, short, ngathers);
     LAP_END(memcp);
#define SET_KERNEL_PARM(idx,name,type)					\
     err = clSetKernelArg(kernel, idx, sizeof(type), name);		\
     if (err != CL_SUCCESS) {						\
       fprintf(stderr, "Error: Failed to set kernel argument: %s\n", #name); \
       exit(1);								\
     } do {} while(0)
     VERBOSE(3,"Set the arguments to the compute kernel\n");
     SET_KERNEL_PARM( 0, &cube, cl_mem);
     SET_KERNEL_PARM( 1, &gx, cl_mem);
     SET_KERNEL_PARM( 2, &gy, cl_mem);
     SET_KERNEL_PARM( 3, &sx, cl_mem);
     SET_KERNEL_PARM( 4, &sy, cl_mem);
     SET_KERNEL_PARM( 5, &scalco, cl_mem);
     SET_KERNEL_PARM( 6, &ntrs, cl_mem);
     SET_KERNEL_PARM( 7, &tr0_dt, cl_mem);
     SET_KERNEL_PARM( 8, &ngathers, int);
     SET_KERNEL_PARM( 9, &ns, int);
     SET_KERNEL_PARM(10, &ctr, cl_mem);
     SET_KERNEL_PARM(11, &sem, cl_mem);
     SET_KERNEL_PARM(12, &stk, cl_mem);


     /* Wait for GPU to finish, read results back and write to output files. 
      * - Output files need a header. Inherit it from original traces 
      */
     
     VERBOSE(3,"Execute kernel and wait for it to finish\n");
     const size_t global = ngathers * ns;
     const size_t work_group_sz = ngathers;
     err = clEnqueueNDRangeKernel(commands,   // Command queue
				  kernel,     // Kernel to be executed
				  1,          // work_dim: number of dimensions used to specify the global work-items and work-items in the work-group
				  NULL,       // global_work_offset
				  &global,    // global_work_size
				  &work_group_sz, // local_work_size -- size of work group
				  0,          // num_events_in_wait_list
				  NULL,       // event_wait_list
				  &prof_event);      // event
     if (err) {
       fprintf(stderr, "Error: clEnqueueNDRangeKernel(..) returned error code %d (", err);
       print_opencl_error(stderr, err);
       fprintf(stderr, ")\n");	
       exit (1);
     }
     
     // Wait for all commands to complete
     clFinish(commands);
	
     err = clWaitForEvents(1, &prof_event);
     cl_ulong start_time, end_time;
     size_t return_bytes;
     err |= clGetEventProfilingInfo(prof_event,
				   CL_PROFILING_COMMAND_QUEUED,
				   sizeof(cl_ulong),
				   &start_time,
				   &return_bytes);
     err |= clGetEventProfilingInfo(prof_event,
				    CL_PROFILING_COMMAND_END,
				    sizeof(cl_ulong),
				    &end_time,
				    &return_bytes);
     if (err != CL_SUCCESS) {
       fprintf(stderr, "Error: Failed to wait for prof_event!\n");
     }
     total_kernel_exec += ((double)(end_time - start_time))/1000000000;


     VERBOSE(3,"Copy results back from GPU and write to file\n");
#define READ_OUTPUT_BUFFER(name,type,sz)				\
     err = clEnqueueReadBuffer(commands, name, CL_TRUE /* Blocking */,	\
			       0, sizeof(type) * sz, _ ## name, 0, NULL, NULL ); \
     if (err != CL_SUCCESS) {						\
       fprintf(stderr, "Error: Failed to read output buffer: %s\n", #name); \
       exit(1);								\
     } do {} while(0)
	  LAP_START(memcpout);
     READ_OUTPUT_BUFFER(ctr, float, ns*ngathers);
     READ_OUTPUT_BUFFER(stk, float, ns*ngathers);
     READ_OUTPUT_BUFFER(sem, float, ns*ngathers);
     LAP_END(memcpout);
#ifdef CHECK_RESULTS
     /* - Check results ---- */
     
     if (global_gather_idx == ng) {
       /* Printing semblance for local gathers 0, 1 and 2 on second batch of gathers (64). */
       fprintf(stderr,"g0     g1\n");
       int t0;
       for (t0=0; t0<ns; t0++) {
	 float _g0 = GET_2D_DATA(_sem,t0,0,ngathers);
	 float _g1 = GET_2D_DATA(_sem,t0,1,ngathers);
	 float _g2 = GET_2D_DATA(_sem,t0,2,ngathers);
	 float _g3 = GET_2D_DATA(_sem,t0,3,ngathers);
	 float _g4 = GET_2D_DATA(_sem,t0,4,ngathers);
	 fprintf(stderr, "> %.16f %.16f %.16f %.16f %.16f\n", _g0, _g1, _g2, _g3, _g4);
       }
     }
#endif // CHECK_RESULTS


     /* Re-layout data to appropriate arrays. 
     */
     int g;
     for (g=0; g<ngathers; g++) {
		int t0;
    fprintf(cTxt,"@ %d\n",totalCoisas++);
       for (t0=0; t0 < ns; t0++) {
	 ctr_t[g].data[t0] = GET_2D_DATA(_ctr,t0,g,ngathers);
	 str_t[g].data[t0] = GET_2D_DATA(_sem,t0,g,ngathers);
	 stk_t[g].data[t0] = GET_2D_DATA(_stk,t0,g,ngathers);
	 
	 fprintf(cTxt,"D %.16f \n",ctr_t[g].data[t0]);

       }
       /* Write to output files. */
       su_fputtr(C_out, &ctr_t[g]);
       su_fputtr(S_out, &str_t[g]);
       su_fputtr(stack, &stk_t[g]);		
     }
   }
   LAP_END(cmp_processing);

   fclose(C_out);
   fclose(S_out);
   fclose(stack);
   LAP_END(total);
   printf("%20s: %12s : %s\n", "Step", "Time in sec.", "Percentage");
   TIMEP(total, "Total time");
   fprintf(time,"%.2f\n",total_kernel_exec);
   fprintf(sem_s,"%e\n",(float)numberOfSemblances/total_kernel_exec);
   fprintf(mem,"%.2f\n",total_memcp+total_memcpout);
   TIMERSUM(memcp,memcpout);
   TIMEPP(total,read_traces, "Reading traces");
   TIMEPP(total,ocl_setup, "OpenCL Setup");
   TIMEPP(total,cmp_processing, "CMP Processing");
   TIMEPP(total,kernel_exec, "Kernel exec");

  return 0;
}
