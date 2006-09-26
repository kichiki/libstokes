/* routines to evaluate CPU time
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bench.c,v 1.2 2006/09/26 18:21:16 ichiki Exp $
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
