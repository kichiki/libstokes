/* header file for ewald.c --
 * utility for Ewald summation calculation
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald.h,v 1.4 2007/02/15 03:27:44 kichiki Exp $
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


/** table version **/

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
 * this routine also can handle non-periodic case seamlessly
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all (int n, const double *x, double *y, void * user_data);

/* make ewald-summed mobility matrix for F/FT/FTS versions
 * with the ewald table
 * this routine also can handle non-periodic case seamlessly
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_ewald_3all (struct stokes * sys, double * mat);

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * through matrix with the ewald table
 * this routine also can handle non-periodic case seamlessly
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_matrix (int n, const double *x,
			  double *y, void * user_data);

/** non-table version **/

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_notbl (int n, const double *x,
			 double *y, void * user_data);
/* make ewald-summed mobility matrix for F/FT/FTS versions
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_ewald_3all_notbl (struct stokes * sys, double * mat);
/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * through matrix
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_matrix_notbl (int n, const double *x,
				double *y, void * user_data);

#endif /* !_EWALD_H_ */
