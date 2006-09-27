/* header file for bench.c --
 * routines to evaluate CPU time
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bench.h,v 1.2 2006/09/27 00:19:49 ichiki Exp $
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
