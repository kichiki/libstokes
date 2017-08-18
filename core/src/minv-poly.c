/* calc (M^inf)^-1 for unequal spheres
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: minv-poly.c,v 1.4 2007/04/26 05:11:12 kichiki Exp $
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
#include <math.h>
#include <stdio.h> // printf()
#include <stdlib.h> // exit()

#include "dgetri_c.h" // lapack_inv_()
#include "non-ewald.h" // scalar_nonewald_poly()

#include "minv-poly.h"


/* calc scalar functions of (M^inf)^-1 in F for unequal spheres
 * INPUT
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar_f [8] : scalar functions in dimensional form!
 *    0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *    4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 */
void
scalars_minv_f_poly (double r, double aa, double ab,
		     double *scalar_f)
{
  double scalar[11];

  // (12)-interaction
  // scalar[] is in the dimensional form (not the SD scaling)
  scalars_nonewald_poly (0, /* F version */
			 r, aa, ab, scalar);
  double xa12, ya12;
  xa12 = scalar [0];
  ya12 = scalar [1];


  // (21)-interaction
  double xa21, ya21;
  /* note:
   * [xy]a12 = [xy]a21 in dimensional form by symmetry
   * therefore, we do not need to calculate (21)-interaction explicitly
   */
  xa21 = xa12;
  ya21 = ya12;

  // self part
  double xa11, ya11;
  xa11 = 1.0 / aa;
  ya11 = 1.0 / aa;

  double xa22, ya22;
  xa22 = 1.0 / ab;
  ya22 = 1.0 / ab;


  double det;

  // XA part
  double XA11, XA12, XA21, XA22;
  det = xa11*xa22 - xa12*xa21;
  XA11 =  xa22 / det;
  XA12 = -xa12 / det;
  XA21 = -xa21 / det;
  XA22 =  xa11 / det;

  // YA part
  double YA11, YA12, YA21, YA22;
  det = ya11*ya22 - ya12*ya21;
  YA11 =  ya22 / det;
  YA12 = -ya12 / det;
  YA21 = -ya21 / det;
  YA22 =  ya11 / det;

  scalar_f [0] = XA11;
  scalar_f [1] = XA12;
  scalar_f [2] = XA21;
  scalar_f [3] = XA22;
  scalar_f [4] = YA11;
  scalar_f [5] = YA12;
  scalar_f [6] = YA21;
  scalar_f [7] = YA22;
}

/* calc scalar functions of (M^inf)^-1 in FT for unequal spheres
 * INPUT
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar_ft [20] : scalar functions in dimensional form!
 *      0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *      4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 *      8, 9,10,11 : (YB11, YB12, YB21, YB22)
 *     12,13,14,15 : (XC11, XC12, XC21, XC22)
 *     16,17,18,19 : (YC11, YC12, YC21, YC22)
 */
void
scalars_minv_ft_poly (double r, double aa, double ab,
		      double *scalar_ft)
{
  double mat[16];
  double scalar[11];

  // (12)-interaction
  // scalar[] is in the dimensional form (not the SD scaling)
  scalars_nonewald_poly (1, // FT version
			 r, aa, ab, scalar);
  double xa12, ya12;
  double yb12;
  double xc12, yc12;

  double aa2 = aa * aa;
  double aa3 = aa2 * aa;
  xa12 = scalar [0];
  ya12 = scalar [1];
  yb12 = scalar [2];
  xc12 = scalar [3];
  yc12 = scalar [4];


  // (21)-interaction
  double xa21, ya21;
  double yb21;
  double xc21, yc21;
  double ab2 = ab * ab;
  double ab3 = ab2 * ab;
  /* note:
   * [xy]a12 = [xy]a21 in dimensional form by symmetry
   * [xy]c12 = [xy]c21 in dimensional form by symmetry
   * yb12 = yb21 in dimensional form for M^infty (lack of a dependence)
   * therefore, we do not need to calculate (21)-interaction explicitly
  */
  xa21 = xa12;
  ya21 = ya12;
  yb21 = yb12;
  xc21 = xc12;
  yc21 = yc12;

  double xa11, ya11;
  xa11 = 1.0 / aa;
  ya11 = 1.0 / aa;

  double xa22, ya22;
  xa22 = 1.0 / ab;
  ya22 = 1.0 / ab;

  double xc11, yc11;
  xc11 = 0.75 / aa3;
  yc11 = 0.75 / aa3;

  double xc22, yc22;
  xc22 = 0.75 / ab3;
  yc22 = 0.75 / ab3;

  double det;

  // XA part
  double XA11, XA12, XA21, XA22;
  det = xa11*xa22 - xa12*xa21;
  XA11 =  xa22 / det;
  XA12 = -xa12 / det;
  XA21 = -xa21 / det;
  XA22 =  xa11 / det;

  // XC part
  double XC11, XC12, XC21, XC22;
  det = xc11*xc22 - xc12*xc21;
  XC11 =  xc22 / det;
  XC12 = -xc12 / det;
  XC21 = -xc21 / det;
  XC22 =  xc11 / det;

  // Y part
  double YA11, YA12, YA21, YA22;
  double YB11, YB12, YB21, YB22;
  double YC11, YC12, YC21, YC22;
  mat [0] = ya11; //  ya11
  mat [1] = ya12; //  ya12
  mat [2] = 0.0;  // -yb11
  mat [3] = yb21; //  yb21

  mat [4] = ya21;  //  ya21
  mat [5] = ya22;  //  ya22
  mat [6] = -yb12; // -yb12
  mat [7] = 0.0;   //  yb22

  mat [8] = 0.0;   // -yb11
  mat [9] = -yb12; // -yb12
  mat[10] = yc11;  //  yc11
  mat[11] = yc12;  //  yc12

  mat[12] = yb21; // yb21
  mat[13] = 0.0;  // yb22
  mat[14] = yc21; // yc21
  mat[15] = yc22; // yc22

  lapack_inv_ (4, mat);
  YA11 =  mat [0];
  YA12 =  mat [1];
  YA21 =  mat [4];
  YA22 =  mat [5];
  YB11 = -mat [8];
  YB12 = -mat [9];
  YB21 =  mat[12];
  YB22 =  mat[13];
  YC11 =  mat[10];
  YC12 =  mat[11];
  YC21 =  mat[14];
  YC22 =  mat[15];

  scalar_ft [0] = XA11;
  scalar_ft [1] = XA12;
  scalar_ft [2] = XA21;
  scalar_ft [3] = XA22;
  scalar_ft [4] = YA11;
  scalar_ft [5] = YA12;
  scalar_ft [6] = YA21;
  scalar_ft [7] = YA22;
  scalar_ft [8] = YB11;
  scalar_ft [9] = YB12;
  scalar_ft[10] = YB21;
  scalar_ft[11] = YB22;
  scalar_ft[12] = XC11;
  scalar_ft[13] = XC12;
  scalar_ft[14] = XC21;
  scalar_ft[15] = XC22;
  scalar_ft[16] = YC11;
  scalar_ft[17] = YC12;
  scalar_ft[18] = YC21;
  scalar_ft[19] = YC22;
}

/* calc scalar functions of (M^inf)^-1 in FTS for unequal spheres
 * INPUT
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar_fts [44] : scalar functions in dimensional form!
 *       0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *       4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 *       8, 9,10,11 : (YB11, YB12, YB21, YB22)
 *      12,13,14,15 : (XC11, XC12, XC21, XC22)
 *      16,17,18,19 : (YC11, YC12, YC21, YC22)
 *      20,21,22,23 : (XG11, XG12, XG21, XG22)
 *      24,25,26,27 : (YG11, YG12, YG21, YG22)
 *      28,29,30,31 : (YH11, YH12, YH21, YH22)
 *      32,33,34,35 : (XM11, XM12, XM21, XM22)
 *      36,37,38,39 : (YM11, YM12, YM21, YM22)
 *      40,41,42,43 : (ZM11, ZM12, ZM21, ZM22)
 */
void
scalars_minv_fts_poly (double r, double aa, double ab,
		       double *scalar_fts)
{
  double scalar [44];
  scalars_nonewald_poly_full (2, // FTS version
			      r, aa, ab, scalar);
  double xa11 = scalar [0];
  double xa12 = scalar [1];
  double xa21 = scalar [2];
  double xa22 = scalar [3];

  double ya11 = scalar [4];
  double ya12 = scalar [5];
  double ya21 = scalar [6];
  double ya22 = scalar [7];

  double yb11 = scalar [8];
  double yb12 = scalar [9];
  double yb21 = scalar[10];
  double yb22 = scalar[11];

  double xc11 = scalar[12];
  double xc12 = scalar[13];
  double xc21 = scalar[14];
  double xc22 = scalar[15];

  double yc11 = scalar[16];
  double yc12 = scalar[17];
  double yc21 = scalar[18];
  double yc22 = scalar[19];

  double xg11 = scalar[20];
  double xg12 = scalar[21];
  double xg21 = scalar[22];
  double xg22 = scalar[23];

  double yg11 = scalar[24];
  double yg12 = scalar[25];
  double yg21 = scalar[26];
  double yg22 = scalar[27];

  double yh11 = scalar[28];
  double yh12 = scalar[29];
  double yh21 = scalar[30];
  double yh22 = scalar[31];

  double xm11 = scalar[32];
  double xm12 = scalar[33];
  double xm21 = scalar[34];
  double xm22 = scalar[35];

  double ym11 = scalar[36];
  double ym12 = scalar[37];
  double ym21 = scalar[38];
  double ym22 = scalar[39];

  double zm11 = scalar[40];
  double zm12 = scalar[41];
  double zm21 = scalar[42];
  double zm22 = scalar[43];


  double mat [36];
  double det;

  // XC part
  double XC11, XC12, XC21, XC22;
  det = xc11*xc22 - xc12*xc21;
  XC11 =  xc22 / det;
  XC12 = -xc12 / det;
  XC21 = -xc21 / det;
  XC22 =  xc11 / det;

  // ZM part
  double ZM11, ZM12, ZM21, ZM22;
  det = zm11*zm22 - zm12*zm21;
  ZM11 =  zm22 / det;
  ZM12 = -zm12 / det;
  ZM21 = -zm21 / det;
  ZM22 =  zm11 / det;

  // Y part
  double YA11, YA12, YA21, YA22;
  double YB11, YB12, YB21, YB22;
  double YC11, YC12, YC21, YC22;
  double YG11, YG12, YG21, YG22;
  double YH11, YH12, YH21, YH22;
  double YM11, YM12, YM21, YM22;
  
  mat [0] = ya11;     //  ya11
  mat [1] = ya12;     //  ya12
  mat [2] = -yb11;    // -yb11
  mat [3] = yb21;     //  yb21
  mat [4] = 2.0*yg11; // 2 yg11
  mat [5] = -2.0*yg21;//-2 yg21

  mat [6] = ya21;     //  ya21
  mat [7] = ya22;     //  ya22
  mat [8] = -yb12;    //  -yb12
  mat [9] = yb22;     //   yb22
  mat[10] = 2.0*yg12; // 2 yg12
  mat[11] = -2.0*yg22;//-2 yg22

  mat[12] = -yb11;    //  -yb11
  mat[13] = -yb12;    //  -yb12
  mat[14] = yc11;     //  yc11
  mat[15] = yc12;     //  yc12
  mat[16] = 2.0*yh11; // 2 yh11
  mat[17] = 2.0*yh21; // 2 yh21

  mat[18] = yb21;     //  yb21
  mat[19] = yb22;     //  yb22
  mat[20] = yc21;     //  yc21
  mat[21] = yc22;     //  yc22
  mat[22] = 2.0*yh12; // 2 yh12
  mat[23] = 2.0*yh22; // 2 yh22

  mat[24] = yg11;    //  yg11
  mat[25] = yg12;    //  yg12
  mat[26] = yh11;    //  yh11
  mat[27] = yh12;    //  yh12
  mat[28] = ym11;    //  ym11
  mat[29] = ym12;    //  ym12

  mat[30] = -yg21;   // -yg21
  mat[31] = -yg22;   // -yg22
  mat[32] = yh21;    //  yh21
  mat[33] = yh22;    //  yh22
  mat[34] = ym21;    //  ym21
  mat[35] = ym22;    //  ym22

  lapack_inv_ (6, mat);
  YA11 =  mat [0];
  YA12 =  mat [1];
  YA21 =  mat [6];
  YA22 =  mat [7];
  YB11 = -mat[12];
  YB12 = -mat[13];
  YB21 =  mat[18];
  YB22 =  mat[19];
  YC11 =  mat[14];
  YC12 =  mat[15];
  YC21 =  mat[20];
  YC22 =  mat[21];
  YG11 =  mat[24];
  YG12 =  mat[25];
  YG21 = -mat[30];
  YG22 = -mat[31];
  YH11 =  mat[26];
  YH12 =  mat[27];
  YH21 =  mat[32];
  YH22 =  mat[33];
  YM11 =  mat[28];
  YM12 =  mat[29];
  YM21 =  mat[34];
  YM22 =  mat[35];


  // X part
  double XA11, XA12, XA21, XA22;
  double XG11, XG12, XG21, XG22;
  double XM11, XM12, XM21, XM22;

  double mp11 = 0.5*(xm11 + zm11); // = 0.9/aa3
  double mp12 = 0.5*(xm12 + zm12);
  double mm11 = 0.5*(xm11 - zm11); // = 0.0
  double mm12 = 0.5*(xm12 - zm12);

  double mp21 = 0.5*(xm21 + zm21); // = 0.9/aa3
  double mp22 = 0.5*(xm22 + zm22);
  double mm21 = 0.5*(xm21 - zm21); // = 0.0
  double mm22 = 0.5*(xm22 - zm22);
  
  mat [0] = xa11;     //  xa11
  mat [1] = xa12;     //  xa12
  mat [2] = -xg11;    // -xg11
  mat [3] = xg21;     //  xg21
  mat [4] = -xg11;    // -xg11
  mat [5] = xg21;     //  xg21

  mat [6] = xa21;     //  xa21
  mat [7] = xa22;     //  xa22
  mat [8] = -xg12;    // -xg12
  mat [9] = xg22;     //  xg22
  mat[10] = -xg12;    // -xg12
  mat[11] = xg22;     //  xg22

  mat[12] = -xg11/3.0;// -xg11/3
  mat[13] = -xg12/3.0;// -xg12/3
  mat[14] = mp11;     // m+11
  mat[15] = mp12;     // m+12
  mat[16] = mm11;     // m-11
  mat[17] = mm12;     // m-12

  mat[18] = xg21/3.0; // xg21/3
  mat[19] = xg22/3.0; // xg22/3
  mat[20] = mp21;     // m+21
  mat[21] = mp22;     // m+22
  mat[22] = mm21;     // m-21
  mat[23] = mm22;     // m-22

  mat[24] = -xg11/3.0;// -xg11/3
  mat[25] = -xg12/3.0;// -xg12/3
  mat[26] = mm11;     // m-11
  mat[27] = mm12;     // m-12
  mat[28] = mp11;     // m+11
  mat[29] = mp12;     // m+12

  mat[30] = xg21/3.0; // xg21/3
  mat[31] = xg22/3.0; // xg22/3
  mat[32] = mm21;     // m-21
  mat[33] = mm22;     // m-22
  mat[34] = mp21;     // m+21
  mat[35] = mp22;     // m+22

  lapack_inv_ (6, mat);
  XA11 = mat [0];
  XA12 = mat [1];
  XA21 = mat [6];
  XA22 = mat [7];
  XG11 = -3.0 * mat[12];
  XG12 = -3.0 * mat[13];
  XG21 = 3.0 * mat[18];
  XG22 = 3.0 * mat[19];
  XM11 = 2.0 * mat[14] - ZM11;
  XM12 = 2.0 * mat[15] - ZM12;
  XM21 = 2.0 * mat[20] - ZM21;
  XM22 = 2.0 * mat[21] - ZM22;

  scalar_fts [0] = XA11;
  scalar_fts [1] = XA12;
  scalar_fts [2] = XA21;
  scalar_fts [3] = XA22;
  scalar_fts [4] = YA11;
  scalar_fts [5] = YA12;
  scalar_fts [6] = YA21;
  scalar_fts [7] = YA22;
  scalar_fts [8] = YB11;
  scalar_fts [9] = YB12;
  scalar_fts[10] = YB21;
  scalar_fts[11] = YB22;
  scalar_fts[12] = XC11;
  scalar_fts[13] = XC12;
  scalar_fts[14] = XC21;
  scalar_fts[15] = XC22;
  scalar_fts[16] = YC11;
  scalar_fts[17] = YC12;
  scalar_fts[18] = YC21;
  scalar_fts[19] = YC22;
  scalar_fts[20] = XG11;
  scalar_fts[21] = XG12;
  scalar_fts[22] = XG21;
  scalar_fts[23] = XG22;
  scalar_fts[24] = YG11;
  scalar_fts[25] = YG12;
  scalar_fts[26] = YG21;
  scalar_fts[27] = YG22;
  scalar_fts[28] = YH11;
  scalar_fts[29] = YH12;
  scalar_fts[30] = YH21;
  scalar_fts[31] = YH22;
  scalar_fts[32] = XM11;
  scalar_fts[33] = XM12;
  scalar_fts[34] = XM21;
  scalar_fts[35] = XM22;
  scalar_fts[36] = YM11;
  scalar_fts[37] = YM12;
  scalar_fts[38] = YM21;
  scalar_fts[39] = YM22;
  scalar_fts[40] = ZM11;
  scalar_fts[41] = ZM12;
  scalar_fts[42] = ZM21;
  scalar_fts[43] = ZM22;
}


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


/** lubrication functions for polydisperse systems **/

/* calc scalar functions of lubrication correction for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle a and b
 *  f12    : (struct twobody_f *).
 *           you can give NULL for them.
 *           then, the coefs are calculated on-the-fly (terribly slow).
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
		  struct twobody_f *f12,
		  int n, int flag_lub,
		  double *lub)
{
  // zero clear
  int i;
  for (i = 0; i < 22; i ++)
    {
      lub[i] = 0.0;
    }

  double res  [22];
  twobody_scalars_res (version,
		       r, a1, a2, f12,
		       n, flag_lub, 1, // dimensional
		       res);

  double minv [44];
  if (version == 0)
    {
      // F version
      scalars_minv_f_poly (r, a1, a2, minv);
    }
  else if (version == 1)
    {
      // FT version
      scalars_minv_ft_poly (r, a1, a2, minv);
    }
  else
    {
      // FTS version
      scalars_minv_fts_poly (r, a1, a2, minv);
    }

  lub [ 0] = res [ 0] - minv [0]; // XA11
  lub [ 1] = res [ 1] - minv [1]; // XA12
  lub [ 2] = res [ 2] - minv [4]; // YA11
  lub [ 3] = res [ 3] - minv [5]; // YA12
  if (version == 0) return;

  lub [ 4] = res [ 4] - minv [8]; // YB11
  lub [ 5] = res [ 5] - minv [9]; // YB12
  lub [ 6] = res [ 6] - minv[12]; // XC11
  lub [ 7] = res [ 7] - minv[13]; // XC12
  lub [ 8] = res [ 8] - minv[16]; // YC11
  lub [ 9] = res [ 9] - minv[17]; // YC12
  if (version == 1) return;

  lub [10] = res [10] - minv[20]; // XG11
  lub [11] = res [11] - minv[21]; // XG12
  lub [12] = res [12] - minv[24]; // YG11
  lub [13] = res [13] - minv[25]; // YG12
  lub [14] = res [14] - minv[28]; // YH11
  lub [15] = res [15] - minv[29]; // YH12
  lub [16] = res [16] - minv[32]; // XM11
  lub [17] = res [17] - minv[33]; // XM12
  lub [18] = res [18] - minv[36]; // YM11
  lub [19] = res [19] - minv[37]; // YM12
  lub [20] = res [20] - minv[40]; // ZM11
  lub [21] = res [21] - minv[41]; // ZM12

}

/* calc scalar functions of lubrication correction for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle a and b
 *  f12,f21: (struct twobody_f *).
 *           you can give NULL for them.
 *           then, the coefs are calculated on-the-fly (terribly slow).
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
		       struct twobody_f *f12, struct twobody_f *f21,
		       int n, int flag_lub,
		       double *lub)
{
  // zero clear
  int i;
  for (i = 0; i < 44; i ++)
    {
      lub[i] = 0.0;
    }

  double res12 [22];
  twobody_scalars_res (version,
		       r, a1, a2, f12,
		       n, flag_lub, 1, // dimensional
		       res12);

  double res21 [22];
  twobody_scalars_res (version,
		       r, a2, a1, f21,
		       n, flag_lub, 1, // dimensional
		       res21);

  double minv [44];
  if (version == 0)
    {
      // F version
      scalars_minv_f_poly (r, a1, a2, minv);
    }
  else if (version == 1)
    {
      // FT version
      scalars_minv_ft_poly (r, a1, a2, minv);
    }
  else
    {
      // FTS version
      scalars_minv_fts_poly (r, a1, a2, minv);
    }

  lub [ 0] = res12 [ 0] - minv [0]; // XA11
  lub [ 1] = res12 [ 1] - minv [1]; // XA12
  lub [ 2] = res21 [ 1] - minv [2]; // XA21
  lub [ 3] = res21 [ 0] - minv [3]; // XA22

  lub [ 4] = res12 [ 2] - minv [4]; // YA11
  lub [ 5] = res12 [ 3] - minv [5]; // YA12
  lub [ 6] = res21 [ 3] - minv [6]; // YA21
  lub [ 7] = res21 [ 2] - minv [7]; // YA22
  if (version == 0) return;

  lub  [8] = res12 [ 4] - minv [8]; // YB11
  lub  [9] = res12 [ 5] - minv [9]; // YB12
  lub [10] = res21 [ 5] - minv[10]; // YB21
  lub [11] = res21 [ 4] - minv[11]; // YB22

  lub [12] = res12 [ 6] - minv[12]; // XC11
  lub [13] = res12 [ 7] - minv[13]; // XC12
  lub [14] = res21 [ 7] - minv[14]; // XC21
  lub [15] = res21 [ 6] - minv[15]; // XC22

  lub [16] = res12 [ 8] - minv[16]; // YC11
  lub [17] = res12 [ 9] - minv[17]; // YC12
  lub [18] = res21 [ 9] - minv[18]; // YC21
  lub [19] = res21 [ 8] - minv[19]; // YC22
  if (version == 1) return;

  lub [20] = res12 [10] - minv[20]; // XG11
  lub [21] = res12 [11] - minv[21]; // XG12
  lub [22] = res21 [11] - minv[22]; // XG21
  lub [23] = res21 [10] - minv[23]; // XG22

  lub [24] = res12 [12] - minv[24]; // YG11
  lub [25] = res12 [13] - minv[25]; // YG12
  lub [26] = res21 [13] - minv[26]; // YG21
  lub [27] = res21 [12] - minv[27]; // YG22

  lub [28] = res12 [14] - minv[28]; // YH11
  lub [29] = res12 [15] - minv[29]; // YH12
  lub [30] = res21 [15] - minv[30]; // YH21
  lub [31] = res21 [14] - minv[31]; // YH22

  lub [32] = res12 [16] - minv[32]; // XM11
  lub [33] = res12 [17] - minv[33]; // XM12
  lub [34] = res21 [17] - minv[34]; // XM21
  lub [35] = res21 [16] - minv[35]; // XM22

  lub [36] = res12 [18] - minv[36]; // YM11
  lub [37] = res12 [19] - minv[37]; // YM12
  lub [38] = res21 [19] - minv[38]; // YM21
  lub [39] = res21 [18] - minv[39]; // YM22

  lub [40] = res12 [20] - minv[40]; // ZM11
  lub [41] = res12 [21] - minv[41]; // ZM12
  lub [42] = res21 [21] - minv[42]; // ZM21
  lub [43] = res21 [20] - minv[43]; // ZM22
}
