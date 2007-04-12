/* header file for minv-poly.c --
 * calc (M^inf)^-1 for unequal spheres
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: minv-poly.h,v 1.1 2007/04/12 05:16:11 kichiki Exp $
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
#ifndef	_MINV_POLY_H_
#define	_MINV_POLY_H_


/* calc scalar functions of (M^inf)^-1 in F for unequal spheres
 * INPUT
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar_f [4] : scalar functions
 *          0, 1 : (XA11, XA12)
 *          2, 3 : (YA11, YA12)
 */
void
scalars_minv_f_poly (double r, double aa, double ab,
		    double *scalar_f);

/* calc scalar functions of (M^inf)^-1 in FT for unequal spheres
 * INPUT
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar_f [10] : scalar functions
 *           0, 1 : (XA11, XA12)
 *           2, 3 : (YA11, YA12)
 *           4, 5 : (YB11, YB12)
 *           6, 7 : (XC11, XC12)
 *           8  9 : (YC11, YC12)
 */
void
scalars_minv_ft_poly (double r, double aa, double ab,
		     double *scalar_ft);

/* calc scalar functions of (M^inf)^-1 in FTS for unequal spheres
 * INPUT
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar_f [22] : scalar functions
 *           0, 1 : (XA11, XA12)
 *           2, 3 : (YA11, YA12)
 *           4, 5 : (YB11, YB12)
 *           6, 7 : (XC11, XC12)
 *           8  9 : (YC11, YC12)
 *          10,11 : (XG11, XG12)
 *          12,13 : (YG11, YG12)
 *          14,15 : (YH11, YH12)
 *          16,17 : (XM11, XM12)
 *          18,19 : (YM11, YM12)
 *          20,21 : (ZM11, ZM12)
 */
void
scalars_minv_fts_poly (double r, double aa, double ab,
		      double *scalar_fts);


/** analytic version **/

/* calc scalar functions of (M^inf)^-1 in F for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle 1 and 2
 * OUTPUT
 *  scalar_f [4] : scalar functions
 *          0, 1 : (XA11, XA12)
 *          2, 3 : (YA11, YA12)
 */
void
scalars_minv_f_poly_ana (double r, double a1, double a2,
			 double *scalar_f);

/* calc scalar functions of (M^inf)^-1 in FT for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle 1 and 2
 * OUTPUT
 *  scalar_f [10] : scalar functions
 *           0, 1 : (XA11, XA12)
 *           2, 3 : (YA11, YA12)
 *           4, 5 : (YB11, YB12)
 *           6, 7 : (XC11, XC12)
 *           8  9 : (YC11, YC12)
 */
void
scalars_minv_ft_poly_ana (double r, double aa, double ab,
			  double *scalar_ft);


#endif /* !_MINV_POLY_H_ */
