/* header file for ewald.c --
 * utility for Ewald summation calculation
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald.h,v 1.1 2006/10/05 00:36:38 ichiki Exp $
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


/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_table (int n, const double *x, double *y, void * user_data);

/* make ewald-summed mobility matrix for F/FT/FTS versions
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_ewald_3all_table (struct stokes * sys, double * mat);

/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for F/FT/FTS versions
 * INPUT
 *  n := np * 11
 *  x [n * 11] : FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n * 11] : UOE
 */
void
atimes_ewald_3all_table_matrix (int n, const double *x,
				double *y, void * user_data);

#endif /* !_EWALD_H_ */
