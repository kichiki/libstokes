/* header file for check-minv.c --
 * test code for the calculation of minv
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-minv.h,v 1.1 2007/04/12 05:28:46 kichiki Exp $
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
#ifndef	_CHECK_MINV_H_
#define	_CHECK_MINV_H_


/* check matrix_mob_nonewald_3all() in FTS version with e=(0,0,1)
 */
int
check_matrix_mob_nonewald_fts (double r, int verbose, double tiny);


/* test the following functions
 *  make_matrix_mob_nonewald_3all() in ...
 *  scalar_minv_fts() in fts.c,
 *  scalars_minv_fts_poly() in minv-poly.c.
 * 1) calculate scalar functions of (M^infty)^-1 in FTS version directly by
 *    make_matrix_mob_nonewald_3all() inversing by lapack_inv_().
 *    check the symmetry on the inverse matrix.
 * 2) calculate scalar functions by the function
 *    scalar_minv_fts() in fts.c,
 *    and compare with the results in 1).
 * 3) calculate scalar functions by the function
 *    scalars_minv_fts_poly() in minv-poly.c
 *    with a1=a2=a (monodisperse case)
 *    and compare with the results in 1).
 */
int
check_minv_fts (double r, int verbose, double tiny);


#endif /* !_CHECK_MINV_H_ */
