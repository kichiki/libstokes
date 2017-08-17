/* header file for periodicity.c --
 * periodicity handling routines
 * Copyright (C) 1995-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: periodicity.h,v 1.1 2006/10/22 04:05:33 kichiki Exp $
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
#ifndef	_PERIODICITY_H_
#define	_PERIODICITY_H_


/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  x [np * 3] : position to check; the range is [0, l[xyz])
 */
void
check_periodic (struct stokes * sys,
		double *x);

/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  angle [np * 3] : angle to check; the range is [0, 2\pi)
 */
void
check_angle (struct stokes * sys, double *angle);

/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  x [np * 2] : position to check; the range is [0, l[xy])
 */
void
check_periodic_2d (struct stokes * sys,
		   double *x);
/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  angle [np] : angle to check; the range is [0, 2\pi)
 */
void
check_angle_2d (struct stokes * sys, double *angle);


#endif /* !_PERIODICITY_H_ */
