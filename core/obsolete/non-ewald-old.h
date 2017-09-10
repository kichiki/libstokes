/* header file for non-ewald-old.c --
 * backup of bug fixing for polydisperse systems
 * utility for non-Ewald routines
 * Copyright (C) 2007-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_NON_EWALD_OLD_H_
#define	_NON_EWALD_OLD_H_


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


/* backup of bug fixing for polydisperse systems of
 * ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_nonewald_3all_old
(int n, const double *x, double *y, void *user_data);


#include "stokes.h" // struct stokes

/* backup of bug fixing for polydisperse systems of
 * make plain mobility matrix for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_nonewald_3all_old
(struct stokes *sys, double *mat);


/* backup of bug fixing for polydisperse systems of
 * ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
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
atimes_nonewald_3all_matrix_old
(int n, const double *x,
 double *y, void *user_data);


#endif /* !_NON_EWALD_OLD_H_ */
