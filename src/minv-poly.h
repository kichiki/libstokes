/* header file for minv-poly.c --
 * calc (M^inf)^-1 for unequal spheres
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: minv-poly.h,v 1.2 2007/04/13 02:06:34 kichiki Exp $
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


/** lubrication functions for polydisperse systems **/

/* calc scalar functions of lubrication correction for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle a and b
 *  n : max order
 *  flag_lub   : 0 to use twobody_far()
 *               1 to use twobody_lub()
 * OUTPUT
 *  lub [22] : scalar functions in dimensional form!
 *      0, 1 : (XA11, XA12)
 *      2, 3 : (YA11, YA12)
 *      4, 5 : (YB11, YB12)
 *      6, 7 : (XC11, XC12)
 *      8  9 : (YC11, YC12)
 *     10,11 : (XG11, XG12)
 *     12,13 : (YG11, YG12)
 *     14,15 : (YH11, YH12)
 *     16,17 : (XM11, XM12)
 *     18,19 : (YM11, YM12)
 *     20,21 : (ZM11, ZM12)
 */
void
scalars_lub_poly (int version,
		  double r, double a1, double a2,
		  int n, int flag_lub,
		  double *lub);

/* calc scalar functions of lubrication correction for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle a and b
 *  n : max order
 *  flag_lub   : 0 to use twobody_far()
 *               1 to use twobody_lub()
 * OUTPUT
 *  lub [44] : scalar functions in dimensional form!
 *    0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *    4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 *    8, 9,10,11 : (YB11, YB12, YB21, YB22)
 *   12,13,14,15 : (XC11, XC12, XC21, XC22)
 *   16,17,18,19 : (YC11, YC12, YC21, YC22)
 *   20,21,22,23 : (XG11, XG12, XG21, XG22)
 *   24,25,26,27 : (YG11, YG12, YG21, YG22)
 *   28,29,30,31 : (YH11, YH12, YH21, YH22)
 *   32,33,34,35 : (XM11, XM12, XM21, XM22)
 *   36,37,38,39 : (YM11, YM12, YM21, YM22)
 *   40,41,42,43 : (ZM11, ZM12, ZM21, ZM22)
 */
void
scalars_lub_poly_full (int version,
		       double r, double a1, double a2,
		       int n, int flag_lub,
		       double *lub);


#endif /* !_MINV_POLY_H_ */
