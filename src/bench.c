/* routines to evaluate CPU time
 * Copyright (C) 1999,2002 Kengo Ichiki <ichiki@pegasus.me.jhu.edu>
 * $Id: bench.c,v 1.1 2002/05/27 22:31:14 ichiki Exp $
 */

#include <stdlib.h> /* exit() */
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "bench.h"

/* benchmark routine */
/*  現プロセスが消費した時間の<差分>をミリ秒単位で返す  */
long
ptime_ms (void)
{
    struct rusage info;
    static long t0 = 0;
    long tmp = t0;

    if (getrusage (RUSAGE_SELF, & info)) exit (1);
    t0 = (info.ru_utime.tv_sec  + info.ru_stime.tv_sec)  * 1000
       + (info.ru_utime.tv_usec + info.ru_stime.tv_usec) / 1000;
    return (t0 - tmp);
}

double
ptime_ms_d (void)
{
  struct rusage info;
  static double t0 = 0;
  double tmp = t0;

  if (getrusage (RUSAGE_SELF, & info)) exit (1);
  t0 = ((double) info.ru_utime.tv_sec
	+ (double) info.ru_stime.tv_sec) * 1000.0
    + ((double) info.ru_utime.tv_usec
       + (double)info.ru_stime.tv_usec) / 1000.0;
  return (t0 - tmp);
}

long
ptime_micros (void)
{
    struct rusage info;
    static long t0 = 0;
    long tmp = t0;

    if (getrusage (RUSAGE_SELF, & info)) exit (1);
    t0 = (info.ru_utime.tv_usec + info.ru_stime.tv_usec);
    return (t0 - tmp);
}
