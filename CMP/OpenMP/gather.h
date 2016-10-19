/***************************************************************************
 *
 * Copyright 2015 Ian Liu Rodrigues and Edson Borin
 *
 ***************************************************************************/

#ifndef GATHER_H__
#define GATHER_H__

#include <su.h>
#include <semblance.h>

void gather_get_traces_for_aperture(int *m0s,
		su_trace_t *traces, int ntr,
		aperture_t *ap);

#endif /* GATHER_H__ */
