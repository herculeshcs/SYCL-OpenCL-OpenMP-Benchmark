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
// Constants C0, C1, NC, TAU, MAX_NTRS, and MAXW must be defined during compilation
#define SQR(x) ((x)*(x))
#ifndef MAX
# define MAX(a, b) ((a)>(b)?(a):(b))
#endif
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

static float interpol_linear(float x0, float x1, float y0, float y1, float x)
{
  return (y1 - y0) * (x - x0) / (x1 - x0) + y0;
}
static int pos2d(int i, int j, int szj)
{
  return i*szj + j;
}
static int pos3d(int i, int j, int k, int szj,int szk)
{
 return i*(szj*szk) + j*szk + k;
 }
#define GET_2D_DATA(data,i,j,szj) data[i*szj + j]
#define GET_3D_DATA(data,i,j,k,szj,szk) data[i*(szj*szk) + j*szk + k]
__kernel void compute_gather(
   __global float* cube,   // cube[ns][MAX_NTRS][ng] --  --- use constant memory?
   __global int* gx,       // gx[MAX_NTRS][ng]      -- traces' gx
   __global int* gy,       // gy[MAX_NTRS][ng]      -- traces' gy
   __global int* sx,       // sx[MAX_NTRS][ng]      -- traces' sx
   __global int* sy,       // sy[MAX_NTRS][ng]      -- traces' sy
   __global short* scalco, // scalco[MAX_NTRS][ng]  -- traces' scalco
   __global int* ntrs,     // ntrs[ng]              -- number of traces per gather
   __global unsigned short* tr0_dt,   // tr0_dt[ng]            -- trace0->dt
   int                ng,    // number of gathers
   int                ns,    // number of samples per trace
   __global   float* ctr,    // ctr[ns][ng]     -- best C for every (t0,gather) pair
   __global   float* sem,    // sem[ns][ng]     -- best semblance for every (t0, gather) pair
   __global   float* stk   // stk[ns][ng]     -- stacked value for best semb. for every (t0, gather) pair
)  
{
    // There are ngathers * nsamples work items. They are grouped into nsamples workgroups.
    // Each work group will work on a different t0 (local_t0).
    // Each work item will work on a different gather (local_gather).

    size_t t0_idx = get_group_id(0); // From 0 to ns-1 -- 1 work group per t0 (ns work groups)
    size_t gather_idx = get_local_id(0); // From 0 to ng-1 -- 1 work item per gather : X gathers => X work items per group
    int ntraces = ntrs[gather_idx];
    __private float best_sem = 0;
    __private float best_stk = 0;
    __private float best_ctr = 0;

    // Loop invariant semblance computation
    float dt = (float) tr0_dt[gather_idx] / 1000000;
    float idt = 1 / dt;
    int tau = MAX((int)((float)TAU * idt), 0);
    int w = 2 * tau + 1;
    // Local arrays used by semblance.
    __private float num[MAXW];
    int valid =0;
    float aux=(C1-C0);
    // For each C, compute semblance and pick the best
    for (int i = 0; i < NC; i++) {
      valid=0;
      float C = C0 + aux * i / NC;
      //GET_2D_DATA(scalco,0,gather_idx,ng)
      float _s = get_trace_scalco(scalco[pos2d(0,gather_idx,ng)]); // (scalco[0][gather_idx]);
      //GET_2D_DATA(gx,0,gather_idx,ng)
      float m0x = _s * (gx[pos2d(0,gather_idx,ng)] + sx[pos2d(0,gather_idx,ng)]) * 0.5; // (gx[0][gather_idx] + sx[0][gather_idx]) * 0.5;
      //GET_2D_DATA(gy,0,gather_idx,ng)
      float m0y = _s * (gy[pos2d(0,gather_idx,ng)] + sy[pos2d(0,gather_idx,ng)]) * 0.5; // (gy[0][gather_idx] + sy[0][gather_idx]) * 0.5;
      float t0 = (float) t0_idx * (float) dt;
      int M = 0, skip = 0;
      float stack_value = 0;
      float den = 0;
      // Clear num[w]
      for (int ii=0; ii<w; ii++)
      {
         num[ii]=0;
      }

      // For each trace 
   for (int tr_idx = 0; tr_idx < ntraces; tr_idx++) {
  float __s = get_trace_scalco(scalco[pos2d(tr_idx,gather_idx,ng)]);
  float t_gx = gx[pos2d(tr_idx,gather_idx,ng)]; 
  float t_gy = gy[pos2d(tr_idx,gather_idx,ng)];
  float t_sx = sx[pos2d(tr_idx,gather_idx,ng)];
  float t_sy = sy[pos2d(tr_idx,gather_idx,ng)];
  float mx = __s * (t_gx + t_sx) * 0.5;
  float my = __s * (t_gy + t_sy) * 0.5;
  float hx = __s * (t_gx - t_sx) * 0.5;
  float hy = __s * (t_gy - t_sy) * 0.5;

  float t = time2d(C, t0, m0x, m0y, mx, my, hx, hy);
  int it = (int)(t * idt);
  if (it-tau >= 0 && it+tau+1<ns) {
     for (int j = 0; j < w; j++) {
         int k = it + j - tau;         
         float d1 = cube[pos3d(k,tr_idx,gather_idx,MAX_NTRS,ng)]; 
       float d2 = cube[pos3d(k+1,tr_idx,gather_idx,MAX_NTRS,ng)];
         float v = interpol_linear(k, k+1, d1, d2, t*idt + j - tau);
         num[j]      += v;
         den         += v*v;
         stack_value += v;
     }
           M++; 
        } else if (++skip == 2) {
      // Error
        valid = 1;
  }
      }  // for each trace end 

      float numerator = 0;
      for (int j = 0; j < w; j++) {
      numerator += num[j] * num[j];
      //  den += denv[j];
      //stack_value += num[j];
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
      ctr[pos2d(t0_idx,gather_idx,ng)]=best_ctr;
      sem[pos2d(t0_idx,gather_idx,ng)]=best_sem;
      stk[pos2d(t0_idx,gather_idx,ng)]=best_stk;
}


