/* backup of bug fixing for polydisperse systems
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
#include "minv-poly-old.h"


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
			   double *scalar)
{
  // A part
  scalar [0] /= a1;
  scalar [1] /= a1;
  scalar [2] /= a2;
  scalar [3] /= a2;

  scalar [4] /= a1;
  scalar [5] /= a1;
  scalar [6] /= a2;
  scalar [7] /= a2;
  // check point for F version
  if (version == 0) return;

  double a12 = a1*a1;
  double a13 = a12*a1;

  double a22 = a2*a2;
  double a23 = a22*a2;
  // B part
  scalar [8] /= a12;
  scalar [9] /= a12;
  scalar[10] /= a22;
  scalar[11] /= a22;

  // C part
  scalar[12] /= a13;
  scalar[13] /= a13;
  scalar[14] /= a23;
  scalar[15] /= a23;

  scalar[16] /= a13;
  scalar[17] /= a13;
  scalar[18] /= a23;
  scalar[19] /= a23;
  // check point for FT version
  if (version == 1) return;

  // G part
  scalar[20] /= a12;
  scalar[21] /= a12;
  scalar[22] /= a22;
  scalar[23] /= a22;

  scalar[24] /= a12;
  scalar[25] /= a12;
  scalar[26] /= a22;
  scalar[27] /= a22;

  // H part
  scalar[28] /= a13;
  scalar[29] /= a13;

  scalar[30] /= a23;
  scalar[31] /= a23;

  // M part
  scalar[32] /= a13;
  scalar[33] /= a13;
  scalar[34] /= a23;
  scalar[35] /= a23;

  scalar[36] /= a13;
  scalar[37] /= a13;
  scalar[38] /= a23;
  scalar[39] /= a23;

  scalar[40] /= a13;
  scalar[41] /= a13;
  scalar[42] /= a23;
  scalar[43] /= a23;
}
