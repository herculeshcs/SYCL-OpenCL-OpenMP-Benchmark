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

#ifndef VECTOR_H__
#define VECTOR_H__

#include <stdlib.h>

#define vector_t(type) struct {int len, cap; type t; type *a;}
#define vector_init(v) ((v).len = (v).cap = 0, (v).a = NULL)
#define vector_push(v, x) do { \
	if ((v).len == (v).cap) { \
		(v).cap *= 2, (v).cap++; \
		(v).a = (su_trace_t*)realloc((v).a, sizeof(__typeof((v).t))*(v).cap); \
	} \
	(v).a[(v).len++] = x; \
} while (0)
#define vector_get(v, i) ((v).a[i])
#define vector_set(v, i, x) ((v).a[i]=x)
#define vector_push2(v, x) do { \
	if ((v).len == (v).cap) { \
		(v).cap *= 2, (v).cap++; \
		(v).a = (su_trace_t**)realloc((v).a, sizeof(__typeof((v).t))*(v).cap); \
	} \
	(v).a[(v).len++] = x; \
} while (0)
#endif /* VECTOR_H__ */
