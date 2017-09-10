/* header file for minv-poly.c --
 * backup of bug fixing for polydisperse systems
 * calc (M^inf)^-1 for unequal spheres
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
#ifndef	_MINV_POLY_OLD_H_
#define	_MINV_POLY_OLD_H_


/* convert scalar functions for resistance from dimensional to SD form
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  a1, a2  : radii for the particles 1 and 2
 *  scalar [44] : scalar functions in dimensional form!
 *   0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *   4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 *   8, 9,10,11 : (YB11, YB12, YB21, YB22)
 *  12,13,14,15 : (XC11, XC12, XC21, XC22)
 *  16,17,18,19 : (YC11, YC12, YC21, YC22)
 *  20,21,22,23 : (XG11, XG12, XG21, XG22)
 *  24,25,26,27 : (YG11, YG12, YG21, YG22)
 *  28,29,30,31 : (YH11, YH12, YH21, YH22)
 *  32,33,34,35 : (XM11, XM12, XM21, XM22)
 *  36,37,38,39 : (YM11, YM12, YM21, YM22)
 *  40,41,42,43 : (ZM11, ZM12, ZM21, ZM22)
 * OUTPUT
 *  scalar [44] : scaled in the SD form.
 */
void
scalars_res_poly_scale_SD (int version,
			   double a1, double a2,
			   double *scalar);


#endif /* !_MINV_POLY_OLD_H_ */
