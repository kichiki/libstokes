/* header file for ewald.c --
 * utility for Ewald summation calculation
 * Copyright (C) 2006-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_EWALD_H_
#define	_EWALD_H_


/*
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 * OUTPUT
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 */
void
scalars_ewald_real (int version,
		    double xi, double r,
		    double *xa, double *ya,
		    double *yb,
		    double *xc, double *yc,
		    double *xg, double *yg,
		    double *yh,
		    double *xm, double *ym, double *zm);

/* calculate scalar functions of (12)-interaction for unequal spheres
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  r      := x_2 - x_1
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 */
void
scalars_ewald_real_poly (int version,
			 double xi, double r,
			 double aa, double ab,
			 double *xa, double *ya,
			 double *yb,
			 double *xc, double *yc,
			 double *xg, double *yg,
			 double *yh,
			 double *xm, double *ym, double *zm);


/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * this is a wrapper for non-periodic and periodic cases
 * also polydisperse systems for non-periodic
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_3all
(int n, const double *x, double *y, void * user_data);


/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all
(int n, const double *x, double *y, void * user_data);


#endif /* !_EWALD_H_ */
