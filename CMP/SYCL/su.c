/*
 *  Copyright 2014-2015 Ian Liu Rodrigues <ian.liu88@gmail.com>
 *
 *  This file is part of CMP.
 *
 *  CMP is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CMP is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CMP.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "su.h"

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
static float get_scalco(su_trace_t *tr)
{
	if (tr->scalco == 0)
		return 1;	
	if (tr->scalco > 0)
		return tr->scalco;
	return 1.0f / tr->scalco;
}

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
	float s = get_scalco(tr);
	*sx = s * tr->sx;
	*sy = s * tr->sy;
}

void su_get_receiver(su_trace_t *tr, float *gx, float *gy)
{
	float s = get_scalco(tr);
	*gx = s * tr->gx;
	*gy = s * tr->gy;
}

void su_get_midpoint(su_trace_t *tr, float *mx, float *my)
{
	float s = get_scalco(tr);
	*mx = s * (tr->gx + tr->sx) * 0.5;
	*my = s * (tr->gy + tr->sy) * 0.5;
}

void su_get_halfoffset(su_trace_t *tr, float *hx, float *hy)
{
	float s = get_scalco(tr);
	*hx = s * (tr->gx - tr->sx) * 0.5;
	*hy = s * (tr->gy - tr->sy) * 0.5;

}
