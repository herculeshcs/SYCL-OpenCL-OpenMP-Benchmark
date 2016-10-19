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

#include "log.h"

#include <stdio.h>
#include <stdarg.h>

void log_progress(float percentage, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	printf("[%3d%%] ", (int)(percentage * 100));
	vprintf(fmt, ap);
	printf("\033[0K\r");
	fflush(stdout);
	va_end(ap);
}
