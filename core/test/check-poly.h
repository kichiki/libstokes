/* header file for check-poly.c --
 * test code for polydisperse handling for non-periodic systems
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-poly.h,v 1.3 2007/04/20 02:02:27 kichiki Exp $
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
#ifndef	_CHECK_POLY_H_
#define	_CHECK_POLY_H_


/* compare scalar functions for two routines
 *  scalars_nonewald() and scalars_nonewald_poly()
 * for equal-sphere case
 * INPUT
 *  r       := x_2 - x_1
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_scalars_nonewald_poly (double r,
			     int verbose, double tiny);


/* check the symmetry between (ab) and (ba) for M^inf
 * INPUT
 *  r       := x_2 - x_1
 *  a1, a2  : radius of particle 1 and 2
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_scalars_nonewald_poly_symmetry (double r, double a1, double a2,
				      int verbose, double tiny);

/* compare scalar functions for two routines
 *  scalars_ewald_real() and scalars_ewald_real_poly()
 * for equal-sphere case
 * INPUT
 *  r       := x_2 - x_1
 *  xi      : splitting parameter for ewald summation
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_scalars_ewald_real_poly (double r, double xi,
			       int verbose, double tiny);


#endif /* !_CHECK_POLY_H_ */
