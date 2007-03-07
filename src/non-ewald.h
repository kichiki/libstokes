/* header file for non-ewald.c --
 * utility for non-Ewald routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: non-ewald.h,v 1.2 2007/03/07 22:32:29 kichiki Exp $
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
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 */
void
scalars_nonewald (int version,
		  double r,
		  double *xa, double *ya,
		  double *yb,
		  double *xc, double *yc,
		  double *xg, double *yg,
		  double *yh,
		  double *xm, double *ym, double *zm);


/* ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_nonewald_3all (int n, const double *x, double *y, void * user_data);


/* make plain mobility matrix for F/FT/FTS versions for non-periodic case
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_nonewald_3all (struct stokes * sys, double * mat);


/* ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
 * through matrix with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_nonewald_3all_matrix (int n, const double *x,
			     double *y, void * user_data);


#endif /* !_NON_EWALD_H_ */
