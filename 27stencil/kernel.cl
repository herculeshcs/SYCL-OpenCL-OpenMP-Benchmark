 /***************************************************************************
 *   Copyright (C) 2015,2016 by Edson Borin, Hercules Cardoso da Silva, Ian Liu Rodrigues
 *   Authors: Edson Borin - edson@ic.unicamp.br
 *	      Hercules Cardoso da Silva - hercules.cardoso.silva1@gmail.com
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

 __kernel void stencil(__global double *in, __global double *out, double fac) {
	 int i = get_global_id(0)+1;
	 int j = get_global_id(1)+1;
	 int k = get_global_id(2)+1;
	 int sz = get_global_size(0)+2;
	 out[i * sz * sz+j*sz+k] = (
	 in[i*sz*sz+( j-1)*sz+k] + in[i*sz*sz+(j+1)*sz+k] +
	 in[(i-1) * sz * sz+j*sz+k] + in[(i+1) * sz * sz+j*sz+k] +
	 in[(i-1)* sz * sz+(j-1) * sz+k] + in[(i-1) * sz * sz+( j+1)*sz+k] +
	 in[(i+1) * sz * sz+(j-1) * sz+k] + in[(i+1) * sz * sz+(j+1)*sz+k] +
	 in[i * sz * sz+(j-1) * sz+(k-1)] + in[i * sz * sz+(j+1) * sz+(k-1)] +
	 in[(i-1) * sz * sz+j * sz+(k-1)] + in[(i+1) * sz * sz+j* sz+(k-1)] +
	 in[(i-1) * sz * sz+(j-1) * sz+(k-1)] + in[(i-1) * sz * sz+(j+1) * sz+(k-1)] +
	 in[(i+1)*sz*sz+(j-1)*sz+(k-1)] + in[(i+1)*sz*sz+(j+1)*sz+(k-1)] +
	 in[i*sz*sz+( j-1)*sz+(k+1)] + in[i*sz*sz+(j+1)*sz+(k+1)] +
	 in[(i-1)*sz*sz+j*sz+(k+1)] + in[(i+1)*sz*sz+j*sz+(k+1)] +
	 in[(i-1)*sz*sz+(j-1)*sz+(k+1)] + in[(i-1)*sz*sz+(j+1)*sz+(k+1)] +
	 in[(i+1)*sz*sz+(j-1)*sz+(k+1)] + in[(i+1)*sz*sz+(j+1)*sz+(k+1)] +
	 in[i*sz*sz+j*sz+(k-1)] + in[i*sz*sz+j*sz+(k+1)]
	 ) * fac;
 }
__kernel void copy(__global double *a0, __global double *a1) {
   int i = get_global_id(0)+1;
   int j = get_global_id(1)+1;
   int k = get_global_id(2)+1;
   int sz = get_global_size(0)+2;
   a0[i*sz*sz+j*sz+k] = a1[i*sz*sz+j*sz+k];
}
