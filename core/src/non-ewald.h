/* header file for non-ewald.c --
 * utility for non-Ewald routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: non-ewald.h,v 1.4 2007/04/12 03:40:00 kichiki Exp $
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

/* calculate scalar functions for unequal spheres
 * under no periodic boundary condition in dimensional form
 * to convert them in the SD form, use scalars_mob_poly_scale_SD ().
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

/* convert scalar functions for mobility from dimensional to SD form
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  a1      : radius for the particle 1
 *            Note that the scalar functions are for (12)-interaction.
 *  scalar [11]:
 *    0, 1,    : (xa12, ya12) for F version
 *    2,       : (yb12)
 *    3, 4,    : (xc12, yc12) for FT version
 *    5, 6,    : (xg12, yg12)
 *    7,       : (yh12)
 *    8, 9, 10 : (xm12, ym12, zm12) for FTS version
 * OUTPUT
 *  scalar [11]: scaled
 */
void
scalars_mob_poly_scale_SD (int version,
			   double a1,
			   double *scalar);

/* calculate scalar functions for unequal spheres
 * under no periodic boundary condition in dimensional form
 * to convert them in the SD form, use scalars_mob_poly_scale_SD ().
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar [44] : scalar functions in dimensional form!
 *   0, 1, 2, 3 : (xa11, xa12, xa21, xa22)
 *   4, 5, 6, 7 : (ya11, ya12, ya21, ya22)
 *   8, 9,10,11 : (yb11, yb12, yb21, yb22)
 *  12,13,14,15 : (xc11, xc12, xc21, xc22)
 *  16,17,18,19 : (yc11, yc12, yc21, yc22)
 *  20,21,22,23 : (xg11, xg12, xg21, xg22)
 *  24,25,26,27 : (yg11, yg12, yg21, yg22)
 *  28,29,30,31 : (yh11, yh12, yh21, yh22)
 *  32,33,34,35 : (xm11, xm12, xm21, xm22)
 *  36,37,38,39 : (ym11, ym12, ym21, ym22)
 *  40,41,42,43 : (zm11, zm12, zm21, zm22)
 */
void
scalars_nonewald_poly_full (int version,
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


#include "stokes.h" // struct stokes

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
