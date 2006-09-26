/* header file for bench.c --
 * routines to evaluate CPU time
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bench.h,v 1.1 2006/09/26 18:21:31 ichiki Exp $
 */
#ifndef	_BENCH_H_
#define	_BENCH_H_


/* return process time difference since the last call in mili-seconds */
long
ptime_ms (void);

/* return process time difference since the last call in mili-seconds */
double
ptime_ms_d (void);

/* return process time difference since the last call in micro-seconds */
long
ptime_micros (void);

#endif /* !_BENCH_H_ */
