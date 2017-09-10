/* header file for ewald-old.c --
 * backup of bug fixing for polydisperse systems
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
#ifndef	_EWALD_OLD_H_
#define	_EWALD_OLD_H_


/* backup of bug fixing for polydisperse systems of
 * convert scalar functions for mobility from dimensional to SD form
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  a1      : radius for the particle 1
 *            Note that the scalar functions are for (12)-interaction.
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 * OUTPUT
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 */
void
scalars_mob_poly_scale_SD_ (int version,
			    double a1,
			    double *xa, double *ya,
			    double *yb,
			    double *xc, double *yc,
			    double *xg, double *yg,
			    double *yh,
			    double *xm, double *ym, double *zm);


/* backup of bug fixing for polydisperse systems of
 * ATIMES of calc ewald-summed mobility for F/FT/FTS versions
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
atimes_3all_old
(int n, const double *x, double *y, void * user_data);

/* backup of bug fixing for polydisperse systems of
 * make mobility matrix for F/FT/FTS versions
 * this is a wrapper for non-periodic and periodic cases
 * also polydisperse systems for non-periodic
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_3all_old
(struct stokes * sys, double * mat);

/* backup of bug fixing for polydisperse systems of
 * ATIMES of calc ewald-summed mobility for F/FT/FTS versions through matrix
 * for both periodic and non-periodic boundary conditions
 * (now polydisperse can be handled for non-periodic case)
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_3all_matrix_old
(int n, const double *x,
 double *y, void * user_data);


/** table version **/

/* backup of bug fixing for polydisperse systems of
 * ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_old
(int n, const double *x, double *y, void * user_data);

/* backup of bug fixing for polydisperse systems of
 * make ewald-summed mobility matrix for F/FT/FTS versions
 * with the ewald table
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_ewald_3all_old
(struct stokes * sys, double * mat);


/** non-table version **/

/* backup of bug fixing for polydisperse systems of
 * ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_notbl_old
(int n, const double *x,
 double *y, void * user_data);

/* backup of bug fixing for polydisperse systems of
 * make ewald-summed mobility matrix for F/FT/FTS versions
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_ewald_3all_notbl_old
(struct stokes * sys, double * mat);

/* backup of bug fixing for polydisperse systems of
 * ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * through matrix
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_matrix_notbl_old
(int n, const double *x,
 double *y, void * user_data);


#endif /* !_EWALD_OLD_H_ */
