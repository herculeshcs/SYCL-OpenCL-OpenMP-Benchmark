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
#ifndef SU_H__
#define SU_H__

#include <stdio.h>

#define SU_HEADER_SIZE ((unsigned long)&(((su_trace_t*)0)->data))
typedef struct fast_su_trace fast_su_trace_t;
struct fast_su_trace{
	unsigned short dt;
	float *data;
	short scalco;
	int sx;
	int sy;
	int gx;
	int gy;
	unsigned short ns;
};

typedef struct su_trace su_trace_t;
struct su_trace {
	int tracl;
	int tracr;
	int fldr;
	int tracf;
	int ep;
	int cdp;
	int cdpt;
	short trid;
	short nvs;
	short nhs;
	short duse;
	int offset;
	int gelev;
	int selev;
	int sdepth;
	int gdel;
	int sdel;
	int swdep;
	int gwdep;
	short scalel;
	short scalco;
	int sx;
	int sy;
	int gx;
	int gy;
	short counit;
	short wevel;
	short swevel;
	short sut;
	short gut;
	short sstat;
	short gstat;
	short tstat;
	short laga;
	short lagb;
	short delrt;
	short muts;
	short mute;
	unsigned short ns;
	unsigned short dt;
	short gain;
	short igc;
	short igi;
	short corr;
	short sfs;
	short sfe;
	short slen;
	short styp;
	short stas;
	short stae;
	short tatyp;
	short afilf;
	short afils;
	short nofilf;
	short nofils;
	short lcf;
	short hcf;
	short lcs;
	short hcs;
	short year;
	short day;
	short hour;
	short minute;
	short sec;
	short timbas;
	short trwf;
	short grnors;
	short grnofr;
	short grnlof;
	short gaps;
	short otrav;
	float d1;
	float f1;
	float d2;
	float f2;
	float ungpow;
	float unscale;
	int ntr;
	short mark;
        short shortpad;
	short unass[14];
	float *data;
};

/*
 * This function allocates `tr->data' array based on `tr.ns', so must be
 * already set!
 */
void su_init(su_trace_t *tr);

/*
 * Frees internal structure.
 */
void su_free(su_trace_t *tr);

/*
 * Reads a trace header and its data from `file'.
 */
int su_fgettr(FILE *file, su_trace_t *tr);

/*
 * Reads a trace header and its data from `stdin'.
 */
int su_gettr(su_trace_t *tr);

/*
 * Writes a trace header and its data to `file'.
 */
int su_fputtr(FILE *file, su_trace_t *tr);

/*
 * Writes a trace header and its data to `stdout'.
 */
int su_puttr(su_trace_t *tr);

/*
 * Returns the CDP key of `tr'.
 */
int su_get_cdp(su_trace_t *tr);

/*
 * Returns the source position of `tr'.
 */
void su_get_source(su_trace_t *tr, float *sx, float *sy);

/*
 * Returns the receiver position of `tr'.
 */
void su_get_receiver(su_trace_t *tr, float *gx, float *gy);

/*
 * Returns the midpoint of `tr'.
 */
void su_get_midpoint(su_trace_t *tr, float *mx, float *my);

/*
 * Returns the half offset of `tr'.
 */
void su_get_halfoffset(su_trace_t *tr, float *hx, float *hy);

/*
 * Returns the midpoint of `tr'.
 */
void fast_su_get_midpoint(fast_su_trace_t *tr, float *mx, float *my);

/*
 * Returns the half offset of `tr'.
 */
void fast_su_get_halfoffset(fast_su_trace_t *tr, float *hx, float *hy);

/*
 * Writes a trace header and its data to `file', in mic.
 */
fast_su_trace_t * convert_to_fast(su_trace_t * a);




#endif
