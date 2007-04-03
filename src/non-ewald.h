/* header file for non-ewald.c --
 * utility for non-Ewald routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: non-ewald.h,v 1.3 2007/04/03 02:10:34 kichiki Exp $
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
#ifndef	_NON_EWALD_H_
#define	_NON_EWALD_H_


/* calculate scalar functions under no periodic boundary condition
 * all functions are explicitly shown in Durlofsky-Brady (1987).
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 * OUTPUT
 *  scalar [11]:
 *   0, 1,    : (xa12, ya12) for F version
 *   2,       : (yb12)
 *   3, 4,    : (xc12, yc12) for FT version
 *   5, 6,    : (xg12, yg12)
 *   7,       : (yh12)
 *   8, 9, 10 : (xm12, ym12, zm12) for FTS version
 */
void
scalars_nonewald (int version,
		  double r,
		  double *scalar);

/* calculate scalar functions under no periodic boundary condition
 * radii for two particles are taken into account
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar [11]:
 *   0, 1,    : (xa12, ya12) for F version
 *   2,       : (yb12)
 *   3, 4,    : (xc12, yc12) for FT version
 *   5, 6,    : (xg12, yg12)
 *   7,       : (yh12)
 *   8, 9, 10 : (xm12, ym12, zm12) for FTS version
 */
void
scalars_nonewald_poly (int version,
		       double r,
		       double aa, double ab,
		       double *scalar);

/* ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_nonewald_3all (int n, const double *x, double *y, void *user_data);


/* make plain mobility matrix for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_nonewald_3all (struct stokes *sys, double *mat);


/* ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * through matrix procedure
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_nonewald_3all_matrix (int n, const double *x,
			     double *y, void *user_data);


#endif /* !_NON_EWALD_H_ */
