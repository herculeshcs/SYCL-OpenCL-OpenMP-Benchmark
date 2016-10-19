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

#include "gather.h"

#include <utils.h>

void gather_get_traces_for_aperture(int *m0s,
		su_trace_t *traces, int ntr,
		aperture_t *ap)
{
	vector_init(ap->traces);
	float x0, y0;
	float x, y;
	su_get_receiver(&traces[*m0s], &x0, &y0);
	su_get_receiver(traces, &x, &y);
	float dist = SQR(x0 - x) + SQR(y0 - y);
	while (ntr-- && dist > SQR(ap->ap_m)) {
		traces++;
		su_get_receiver(traces, &x, &y);
		dist = SQR(x0 - x) + SQR(y0 - y);
		(*m0s)--;
	}
	ntr++;
	if (ntr == 0)
		return;
	while (ntr-- && dist <= SQR(ap->ap_m)) {
		vector_push(ap->traces, traces);
		traces++;
		su_get_receiver(traces, &x, &y);
		dist = SQR(x0 - x) + SQR(y0 - y);
	}
}

