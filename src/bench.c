/* routines to evaluate CPU time
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bench.c,v 1.3 2006/09/27 00:19:17 ichiki Exp $
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <stdlib.h> /* exit() */
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "bench.h"

/* return process time difference since the last call in mili-seconds */
long
ptime_ms (void)
{
  static long t0 = 0;

  struct rusage info;
  long tmp = t0;

  if (getrusage (RUSAGE_SELF, & info)) exit (1);
  t0 = (info.ru_utime.tv_sec  + info.ru_stime.tv_sec)  * 1000
    + (info.ru_utime.tv_usec + info.ru_stime.tv_usec) / 1000;
  return (t0 - tmp);
}

/* return process time difference since the last call in mili-seconds */
double
ptime_ms_d (void)
{
  static double t0 = 0;

  struct rusage info;
  double tmp = t0;

  if (getrusage (RUSAGE_SELF, & info)) exit (1);
  t0 = ((double) info.ru_utime.tv_sec
	+ (double) info.ru_stime.tv_sec) * 1000.0
    + ((double) info.ru_utime.tv_usec
       + (double)info.ru_stime.tv_usec) / 1000.0;
  return (t0 - tmp);
}

/* return process time difference since the last call in micro-seconds */
long
ptime_micros (void)
{
  static long t0 = 0;

  struct rusage info;
  long tmp = t0;

  if (getrusage (RUSAGE_SELF, & info)) exit (1);
  t0 = (info.ru_utime.tv_usec + info.ru_stime.tv_usec);
  return (t0 - tmp);
}
