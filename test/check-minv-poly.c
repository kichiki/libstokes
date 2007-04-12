/* test code for minv-poly.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-minv-poly.c,v 1.1 2007/04/12 05:31:36 kichiki Exp $
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory-check.h"

#include <stokes.h> // struct stokes
#include <f.h>   // scalar_minv_f()
#include <ft.h>  // scalar_minv_ft()
#include <fts.h> // scalar_minv_fts()
#include <minv-poly.h>

#include "check.h" // compare()


/** analytic version **/

/* calc scalar functions of (M^inf)^-1 in F for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle 1 and 2
 * OUTPUT
 *  scalar_f [8] : scalar functions in dimensional form!
 *    0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *    4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 */
void
scalars_minv_f_poly_ana (double r, double a1, double a2,
			 double *scalar_f)
{
  double r2 = r*r;
  double r3 = r2*r;
  double r6 = r3*r3;
  double a1a2 = a1 * a2;
  double a12a22 = a1*a1 + a2*a2;
  double a12a22_2 = a12a22 * a12a22;

  double dxa =
    (- a1a2*a12a22_2 +r2*
     (+ 6.0*a1a2*a12a22 +r2*
      (- 9.0*a1a2 +r2*
       (4.0))));

  // XA11
  scalar_f [0] = 4.0*a1*r6 / dxa;
  /*
                                            6
                                      4 a1 r
     --------------------------------------------------------------------------
        6            4           3       3      2        5       3   3     5
     4 r  - 9 a1 a2 r  + (6 a1 a2  + 6 a1  a2) r  - a1 a2  - 2 a1  a2  - a1  a2
  */
  // XA12
  scalar_f [1] = -2.0*a1a2*r3*(3.0*r2 - a12a22) / dxa;
  /*
                               5             3       3      3
                      6 a1 a2 r  + (- 2 a1 a2  - 2 a1  a2) r
   - --------------------------------------------------------------------------
        6            4           3       3      2        5       3   3     5
     4 r  - 9 a1 a2 r  + (6 a1 a2  + 6 a1  a2) r  - a1 a2  - 2 a1  a2  - a1  a2
  */
  // XA21
  scalar_f [2] = -2.0*a1a2*r3*(3.0*r2 - a12a22) / dxa; // == XA12
  /*
                               5             3       3      3
                      6 a1 a2 r  + (- 2 a1 a2  - 2 a1  a2) r
   - --------------------------------------------------------------------------
        6            4           3       3      2        5       3   3     5
     4 r  - 9 a1 a2 r  + (6 a1 a2  + 6 a1  a2) r  - a1 a2  - 2 a1  a2  - a1  a2
  */
  // XA22
  scalar_f [3] = 4.0*a2*r6 / dxa; // == (a2/a1) XA11
  /*
                                            6
                                      4 a2 r
     --------------------------------------------------------------------------
        6            4           3       3      2        5       3   3     5
     4 r  - 9 a1 a2 r  + (6 a1 a2  + 6 a1  a2) r  - a1 a2  - 2 a1  a2  - a1  a2
  */



  double dy =
    (-a1a2*a12a22*a12a22 +r2*
     (-6.0*a1a2*a12a22 +r2*
      (- 9.0*a1a2 +r2*
       (16.0))));
  // YA11
  scalar_f [4] = 16.0*a1*r6 / dy;
  /*
                                           6
                                    16 a1 r
  -----------------------------------------------------------------------------
      6            4             3       3      2        5       3   3     5
  16 r  - 9 a1 a2 r  + (- 6 a1 a2  - 6 a1  a2) r  - a1 a2  - 2 a1  a2  - a1  a2
  */
  // YA12
  scalar_f [5] = -4.0*a1a2*r3*(3.0*r2 + a12a22) / dy;
  /*
                    5           3       3      3
(%o55) - (12 a1 a2 r  + (4 a1 a2  + 4 a1  a2) r )
      6            4             3       3      2        5       3   3
/(16 r  - 9 a1 a2 r  + (- 6 a1 a2  - 6 a1  a2) r  - a1 a2  - 2 a1  a2
     5
 - a1  a2)
  */
    ;
  // YA21
  scalar_f [6] = -4.0*a1a2*r3*(3.0*r2 + a12a22) / dy; // == YA12
  /*
                    5           3       3      3
(%o56) - (12 a1 a2 r  + (4 a1 a2  + 4 a1  a2) r )
      6            4             3       3      2        5       3   3
/(16 r  - 9 a1 a2 r  + (- 6 a1 a2  - 6 a1  a2) r  - a1 a2  - 2 a1  a2
     5
 - a1  a2)
  */
  // YA22
  scalar_f [7] = 16.0*a2*r6 / dy; // == (a2/a1)YA11
  /*
                                           6
                                    16 a2 r
  -----------------------------------------------------------------------------
      6            4             3       3      2        5       3   3     5
  16 r  - 9 a1 a2 r  + (- 6 a1 a2  - 6 a1  a2) r  - a1 a2  - 2 a1  a2  - a1  a2
  */
}

/* calc scalar functions of (M^inf)^-1 in FT for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle 1 and 2
 * OUTPUT
 *  scalar_ft [20] : scalar functions in dimensional form!
 *      0, 1, 2, 3 : (XA11, XA12, XA21, XA22)
 *      4, 5, 6, 7 : (YA11, YA12, YA21, YA22)
 *      8, 9,10,11 : (YB11, YB12, YB21, YB22)
 *     12,13,14,15 : (XC11, XC12, XC21, XC22)
 *     16,17,18,19 : (YC11, YC12, YC21, YC22)
 */
void
scalars_minv_ft_poly_ana (double r, double a1, double a2,
			  double *scalar_ft)
{
  double r2 = r*r;
  double r3 = r2*r;
  double r6 = r3*r3;
  double a1a2 = a1 * a2;
  double a12a22 = a1*a1 + a2*a2;
  double a12a22_2 = a12a22 * a12a22;

  double dxa =
    (- a1a2*a12a22_2 +r2*
     (+ 6.0*a1a2*a12a22 +r2*
      (- 9.0*a1a2 +r2*
       (4.0))));

  // XA -- same in F version
  // XA11
  scalar_ft [0] = 4.0*a1*r6 / dxa;
  // XA12
  scalar_ft [1] = -2.0*a1a2*r3*(3.0*r2 - a12a22) / dxa;
  // XA21
  scalar_ft [2] = -2.0*a1a2*r3*(3.0*r2 - a12a22) / dxa; // == XA12
  // XA22
  scalar_ft [3] = 4.0*a2*r6 / dxa; // == (a2/a1) XA11

  double a14 = a1*a1*a1*a1;
  double a24 = a2*a2*a2*a2;
  double a1a2_2 = a1a2 * a1a2;
  double a1a2_3 = a1a2_2 * a1a2;
  double a1a2_4 = a1a2_2 * a1a2_2;
  double dy = 
    (a1a2_4*a12a22_2 +r2*
     (-6.0*a1a2_4*a12a22 +r2*
      (+9.0*a1a2_4 +r2*
       (-4.0*a1a2*(a24  + 6.0*a1a2_2 + a14) +r2*
	(-72.0*a1a2*a12a22 +r2*
	 (-36.0*a1a2 +r2*
	  (+64.0)))))));

  double a12 = a1*a1;
  // YA11
  scalar_ft [4] = 16.0*a1*r6*(4.0*r6 - 3.0*a1a2*a12*r2  - a1a2_3) / dy;
  /*
               12        4     8        4   3  6
(%o65) (64 a1 r   - 48 a1  a2 r  - 16 a1  a2  r )
      12             10              3        3      8
/(64 r   - 36 a1 a2 r   + (- 72 a1 a2  - 72 a1  a2) r
             5        3   3       5      6       4   4  4
 + (- 4 a1 a2  - 24 a1  a2  - 4 a1  a2) r  + 9 a1  a2  r
          4   6       6   4   2     4   8       6   6     8   4
 + (- 6 a1  a2  - 6 a1  a2 ) r  + a1  a2  + 2 a1  a2  + a1  a2 )
   */
  // YA12
  scalar_ft [5] = -4.0*a1a2*r3*
	 (-a1a2_3*a12a22 +r2*
	  (+3.0*a1a2_3 +r2*r2*
	   (+4.0*a12a22 +r2*
	    (12.0))))
	 / dy;
  /*
                    11            3        3      9        4   4  5
(%o67) - (48 a1 a2 r   + (16 a1 a2  + 16 a1  a2) r  + 12 a1  a2  r
          4   6       6   4   3       12             10
 + (- 4 a1  a2  - 4 a1  a2 ) r )/(64 r   - 36 a1 a2 r
              3        3      8             5        3   3       5      6
 + (- 72 a1 a2  - 72 a1  a2) r  + (- 4 a1 a2  - 24 a1  a2  - 4 a1  a2) r
       4   4  4          4   6       6   4   2     4   8       6   6     8   4
 + 9 a1  a2  r  + (- 6 a1  a2  - 6 a1  a2 ) r  + a1  a2  + 2 a1  a2  + a1  a2 )
  */
  // YA21 == YA12
  scalar_ft [6] = scalar_ft [5];
  /*
                    11            3        3      9        4   4  5
(%o68) - (48 a1 a2 r   + (16 a1 a2  + 16 a1  a2) r  + 12 a1  a2  r
          4   6       6   4   3       12             10
 + (- 4 a1  a2  - 4 a1  a2 ) r )/(64 r   - 36 a1 a2 r
              3        3      8             5        3   3       5      6
 + (- 72 a1 a2  - 72 a1  a2) r  + (- 4 a1 a2  - 24 a1  a2  - 4 a1  a2) r
       4   4  4          4   6       6   4   2     4   8       6   6     8   4
 + 9 a1  a2  r  + (- 6 a1  a2  - 6 a1  a2 ) r  + a1  a2  + 2 a1  a2  + a1  a2 )
  */
  double a22 = a2*a2;
  // YA22 != (a2/a1) YA11... (of course, YB21(a2,a1) = YB21(a1,a2), though)
  scalar_ft [7] = 16.0*a2*r6*(4.0*r6 - 3.0*a1a2*a22*r2  - a1a2_3) / dy;
  /*
               12           4  8        3   4  6
(%o66) (64 a2 r   - 48 a1 a2  r  - 16 a1  a2  r )
      12             10              3        3      8
/(64 r   - 36 a1 a2 r   + (- 72 a1 a2  - 72 a1  a2) r
             5        3   3       5      6       4   4  4
 + (- 4 a1 a2  - 24 a1  a2  - 4 a1  a2) r  + 9 a1  a2  r
          4   6       6   4   2     4   8       6   6     8   4
 + (- 6 a1  a2  - 6 a1  a2 ) r  + a1  a2  + 2 a1  a2  + a1  a2 )
   */


  double a13 = a1*a1*a1;
  double r7 = r6*r;
  // YB11
  scalar_ft [8] = -16.0*a13*a1a2*r7*
	 (3.0*r2 +3.0*a22 + a12)
	 / dy;
  /*
               4     9         4   3        6      7
(%o69) - (48 a1  a2 r  + (48 a1  a2  + 16 a1  a2) r )
      12             10              3        3      8
/(64 r   - 36 a1 a2 r   + (- 72 a1 a2  - 72 a1  a2) r
             5        3   3       5      6       4   4  4
 + (- 4 a1 a2  - 24 a1  a2  - 4 a1  a2) r  + 9 a1  a2  r
          4   6       6   4   2     4   8       6   6     8   4
 + (- 6 a1  a2  - 6 a1  a2 ) r  + a1  a2  + 2 a1  a2  + a1  a2 )
   */
  double r4 = r2*r2;
  // YB12
  scalar_ft [9] = 8.0*a1a2*r4*
	 (+a1a2_3*a12a22 +r2*
	  (-3.0*a1a2_3 +r2*r2*
	   (+8.0*a12)))
	 / dy;
  /*
             3     10        4   4  6        4   6       6   4   4
(%o71) (64 a1  a2 r   - 24 a1  a2  r  + (8 a1  a2  + 8 a1  a2 ) r )
      12             10              3        3      8
/(64 r   - 36 a1 a2 r   + (- 72 a1 a2  - 72 a1  a2) r
             5        3   3       5      6       4   4  4
 + (- 4 a1 a2  - 24 a1  a2  - 4 a1  a2) r  + 9 a1  a2  r
          4   6       6   4   2     4   8       6   6     8   4
 + (- 6 a1  a2  - 6 a1  a2 ) r  + a1  a2  + 2 a1  a2  + a1  a2 )
   */
  // YB21 != YB12... (of course, YB21(a2,a1) = YB21(a1,a2), though)
  scalar_ft [10] = 8.0*a1a2*r4*
	 (+a1a2_3*a12a22 +r2*
	  (-3.0*a1a2_3 +r2*r2*
	   (+8.0*a22)))
	 / dy;
  /*
                3  10        4   4  6        4   6       6   4   4
(%o72) (64 a1 a2  r   - 24 a1  a2  r  + (8 a1  a2  + 8 a1  a2 ) r )
      12             10              3        3      8
/(64 r   - 36 a1 a2 r   + (- 72 a1 a2  - 72 a1  a2) r
             5        3   3       5      6       4   4  4
 + (- 4 a1 a2  - 24 a1  a2  - 4 a1  a2) r  + 9 a1  a2  r
          4   6       6   4   2     4   8       6   6     8   4
 + (- 6 a1  a2  - 6 a1  a2 ) r  + a1  a2  + 2 a1  a2  + a1  a2 )
   */
  double a23 = a22*a2;
  // YB22 != (a2/a1) YB11... (of course, YB21(a2,a1) = YB21(a1,a2), though)
  scalar_ft [11] = -16.0*a23*a1a2*r7*
	 (3.0*r2 +3.0*a12 + a22)
	 / dy;
  /*
                  4  9            6        3   4   7
(%o70) - (48 a1 a2  r  + (16 a1 a2  + 48 a1  a2 ) r )
      12             10              3        3      8
/(64 r   - 36 a1 a2 r   + (- 72 a1 a2  - 72 a1  a2) r
             5        3   3       5      6       4   4  4
 + (- 4 a1 a2  - 24 a1  a2  - 4 a1  a2) r  + 9 a1  a2  r
          4   6       6   4   2     4   8       6   6     8   4
 + (- 6 a1  a2  - 6 a1  a2 ) r  + a1  a2  + 2 a1  a2  + a1  a2 )
   */

  double a13a23 = a1a2 * a1a2 * a1a2;
  double dxc = 3.0*(r6 - a13a23);
  // XC11
  scalar_ft [12] = 4.0*a13*r6 / dxc;
  /*
   *          3  6
   *      4 a1  r
   *  ----------------
   *     6       3   3
   *  3 r  - 3 a1  a2
   */

  // XC12
  scalar_ft [13] = -4.0*a13a23*r3 / dxc;
  /*
   *          3   3  3
   *      4 a1  a2  r
   *  - ----------------
   *       6       3   3
   *    3 r  - 3 a1  a2
   */
  // XC21 == XC12
  scalar_ft [14] = scalar_ft [13];
  /*
   *          3   3  3
   *      4 a1  a2  r
   *  - ----------------
   *       6       3   3
   *    3 r  - 3 a1  a2
   */
  // XC22 == (a2/a1)^3 XC11
  scalar_ft [15] = 4.0*a23*r6 / dxc;
  /*
   *          3  6
   *      4 a2  r
   *  ----------------
   *     6       3   3
   *  3 r  - 3 a1  a2
   */


  // YC11
  scalar_ft [16] = 16.0*a13*r6*
    (-a1a2*a12a22_2 +r2*
     (-6.0*a1a2*(3.0*a22 + a12) +r2*
      (-9.0*a1a2 +r2*
       (16.0))))
    / 3.0 / dy;
  /*
              3  12         4     10            4   3        6      8
(%o73) (256 a1  r   - 144 a1  a2 r   + (- 288 a1  a2  - 96 a1  a2) r
           4   5        6   3        8      6
 + (- 16 a1  a2  - 32 a1  a2  - 16 a1  a2) r )
       12              10               3         3      8
/(192 r   - 108 a1 a2 r   + (- 216 a1 a2  - 216 a1  a2) r
              5        3   3        5      6        4   4  4
 + (- 12 a1 a2  - 72 a1  a2  - 12 a1  a2) r  + 27 a1  a2  r
           4   6        6   4   2       4   8       6   6       8   4
 + (- 18 a1  a2  - 18 a1  a2 ) r  + 3 a1  a2  + 6 a1  a2  + 3 a1  a2 )
   */
  // YC12
  scalar_ft [17] = 8.0*r3*a1a2_3*
	 (-a1a2*a12a22_2 +r2*r2*
	  (+9.0*a1a2 +r2*
	   (16.0)))
	 / 3.0 / dy;
  /*
              3   3  9        4   4  7          4   8        6   6       8   4
(%o75) (128 a1  a2  r  + 72 a1  a2  r  + (- 8 a1  a2  - 16 a1  a2  - 8 a1  a2 )
  3        12              10               3         3      8
 r )/(192 r   - 108 a1 a2 r   + (- 216 a1 a2  - 216 a1  a2) r
              5        3   3        5      6        4   4  4
 + (- 12 a1 a2  - 72 a1  a2  - 12 a1  a2) r  + 27 a1  a2  r
           4   6        6   4   2       4   8       6   6       8   4
 + (- 18 a1  a2  - 18 a1  a2 ) r  + 3 a1  a2  + 6 a1  a2  + 3 a1  a2 )
   */
  // YC21 == YC12
  scalar_ft [18] = scalar_ft [17];
  /*
              3   3  9        4   4  7          4   8        6   6       8   4
(%o76) (128 a1  a2  r  + 72 a1  a2  r  + (- 8 a1  a2  - 16 a1  a2  - 8 a1  a2 )
  3        12              10               3         3      8
 r )/(192 r   - 108 a1 a2 r   + (- 216 a1 a2  - 216 a1  a2) r
              5        3   3        5      6        4   4  4
 + (- 12 a1 a2  - 72 a1  a2  - 12 a1  a2) r  + 27 a1  a2  r
           4   6        6   4   2       4   8       6   6       8   4
 + (- 18 a1  a2  - 18 a1  a2 ) r  + 3 a1  a2  + 6 a1  a2  + 3 a1  a2 )
   */
  // YC22 != (a2/a1)YC11 (of course, YC22(a2,a1) = YC11(a1,a2), though)
  scalar_ft [19] = 16.0*a23*r6*
    (-a1a2*a12a22_2 +r2*
     (-6.0*a1a2*(3.0*a12 + a22) +r2*
      (-9.0*a1a2 +r2*
       (16.0))))
    / 3.0 / dy;
  /*
              3  12            4  10              6         3   4   8
(%o74) (256 a2  r   - 144 a1 a2  r   + (- 96 a1 a2  - 288 a1  a2 ) r
              8        3   6        5   4   6
 + (- 16 a1 a2  - 32 a1  a2  - 16 a1  a2 ) r )
       12              10               3         3      8
/(192 r   - 108 a1 a2 r   + (- 216 a1 a2  - 216 a1  a2) r
              5        3   3        5      6        4   4  4
 + (- 12 a1 a2  - 72 a1  a2  - 12 a1  a2) r  + 27 a1  a2  r
           4   6        6   4   2       4   8       6   6       8   4
 + (- 18 a1  a2  - 18 a1  a2 ) r  + 3 a1  a2  + 6 a1  a2  + 3 a1  a2 )
   */

}

/* calc scalar functions of (M^inf)^-1 in FTS for unequal spheres
 * INPUT
 *  r      := x_2 - x_1
 *  a1, a2 : radius of particle 1 and 2
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
scalars_minv_fts_poly_ana (double r, double a1, double a2,
			   double *scalar_fts)
{
  double r2 = r*r;
  double r3 = r2*r;
  double r6 = r3*r3;
  double a1a2 = a1 * a2;
  double a12a22 = a1*a1 + a2*a2;
  double a12a22_2 = a12a22 * a12a22;

  double a13 = a1 * a1 * a1;
  double a13a23 = a1a2 * a1a2 * a1a2;
  double a23 = a2 * a2 * a2;
  double dxc = 3.0*(r6 - a13a23);
  // XC11 (same in FT version)
  scalar_fts [12] = 4.0*a13*r6 / dxc;
  // XC12 (same in FT version)
  scalar_fts [13] = -4.0*a13a23*r3 / dxc;
  // XC21 == XC12
  scalar_fts [14] = scalar_fts [13];
  // XC22 == (a2/a1)^3 XC11
  scalar_fts [15] = 4.0*a23*r6 / dxc;


  // ZM part
  double a1a2_3 = a1a2 * a1a2 * a1a2;
  double r5 = r3*r2;
  double r10 = r5*r5;
  double dz = 9.0*(r10 - a1a2_3*a12a22_2);
  // ZM11
  scalar_fts[40] = 10.0*a13*r10 / dz;
  /*
   *                       3  10
   *                  10 a1  r
   *  ------------------------------------------
   *     10       3   7        5   5       7   3
   *  9 r   - 9 a1  a2  - 18 a1  a2  - 9 a1  a2
   */
  // ZM12
  scalar_fts[41] = -10.0*r5*a1a2_3*a12a22 / dz;
  /*
   *                 3   5        5   3   5
   *           (10 a1  a2  + 10 a1  a2 ) r
   *  - ------------------------------------------
   *       10       3   7        5   5       7   3
   *    9 r   - 9 a1  a2  - 18 a1  a2  - 9 a1  a2
   */
  // ZM21 == ZM12
  scalar_fts[42] = scalar_fts[41];
  /*
   *                 3   5        5   3   5
   *           (10 a1  a2  + 10 a1  a2 ) r
   *  - ------------------------------------------
   *       10       3   7        5   5       7   3
   *    9 r   - 9 a1  a2  - 18 a1  a2  - 9 a1  a2
   */
  // ZM22 == (a2/a1)^3 ZM11
  scalar_fts[43] = 10.0*a23*r10 / dz;
  /*
   *                       3  10
   *                  10 a2  r
   *  ------------------------------------------
   *     10       3   7        5   5       7   3
   *  9 r   - 9 a1  a2  - 18 a1  a2  - 9 a1  a2
   */



  double r4 = r2*r2;
  double r8 = r6*r2;
  double r12 = r6*r6;
  double r14 = r12*r2;
  double r16 = r14*r2;
  double r18 = r16*r2;
  double r20 = r18*r2;
  double r22 = r20*r2;

  double a12 = a1*a1;
  double a14 = a12 * a12;
  double a15 = a14 * a1;
  double a16 = a14 * a12;
  double a17 = a16 * a1;
  double a18 = a16 * a12;
  double a19 = a18 * a1;
  double a110 = a18 * a12;
  double a111 = a19 * a12;

  double a22 = a2*a2;
  double a24 = a22 * a22;
  double a25 = a23 * a22;
  double a26 = a23 * a23;
  double a27 = a26 * a2;
  double a28 = a26 * a22;
  double a29 = a28 * a2;
  double a210 = a28 * a22;
  double a211 = a210 * a2;

  double a1a2_4 = a1a2 * a1a2 * a1a2 * a1a2;
  double a1a2_7 = a1a2_4 * a1a2 * a1a2 * a1a2;

  double dy =
    1600.0*r22
    -900.0*a1a2*r20
    -1800.0*a1a2*a12a22*r18
    -100.0*a1a2*(a24
		 +226.0*a12*a22
		 +a14)*r16
    +a1a2*(-720.0*a26
	   +27600.0*a12*a24
	   +12600.0*a13*a23
	   +27600.0*a14*a22
	   -720.0*a16)*r14
    +a1a2_3*(-25600.0*a24
	     -5100.0*a1*a23
	     -51200.0*a12*a22
	     -5100.0*a13*a2
	     -25600.0*a14)*r12
    +a1a2_4*(-8860.0*a24
	     +32600.0*a12*a22
	     -8860.0*a14)*r10
    +a1a2_4*(9580.0*a26
	     -15020.0*a12*a24
	     -5625.0*a13*a23
	     -15020.0*a14*a22
	     +9580.0*a16)*r8
    +a1a2_4*(100.0*a28
	     +6640.0*a12*a26
	     +1500.0*a13*a25
	     +13144.0*a14*a24
	     +1500.0*a15*a23
	     +6640.0*a16*a22
	     +100.0*a18)*r6
    +a1a2_7*(650.0*a24
	     +700.0*a12*a22
	     +650.0*a14)*r4
    +a1a2_7*(-100.0*a26
	     -220.0*a12*a24
	     -220.0*a14*a22
	     -100.0*a16)*r2
    +a1a2_7*(-25.0*a28
	     -60.0*a12*a26
	     -86.0*a14*a24
	     -60.0*a16*a22
	     -25.0*a18)
    ;
  // YA11
  scalar_fts [4] = r6*
    (
     (6400.0*a17*a210
      +12800.0*a19*a28
      +6400.0*a111*a26) +r2*
     (+(-28800.0*a17*a28
	-22400.0*a19*a26
	+7680.0*a111*a24) +r2*
      (+(28000.0*a17*a26
	 -16800.0*a19*a24) +r2*
       (+(-25600.0*a14*a27
	  -51200.0*a16*a25
	  +12000.0*a17*a24
	  -25600.0*a18*a23) +r2*
	(+(30000.0*a14*a25
	   +29600.0*a16*a23
	   -720.0*a18*a2) +r2*
	 (-22400.0*a14*a23 +r2*
	  (-1200.0*a14*a2 +r2*r2*
	   (+1600.0*a1))))))))
    / dy;
  /*
                 22          4     18           4   3  16
(%o85) (1600 a1 r   - 1200 a1  a2 r   - 22400 a1  a2  r
            4   5           6   3         8      14
 + (30000 a1  a2  + 29600 a1  a2  - 720 a1  a2) r
              4   7           6   5           7   4           8   3   12
 + (- 25600 a1  a2  - 51200 a1  a2  + 12000 a1  a2  - 25600 a1  a2 ) r
            7   6           9   4   10              7   8           9   6
 + (28000 a1  a2  - 16800 a1  a2 ) r   + (- 28800 a1  a2  - 22400 a1  a2
          11   4   8           7   10           9   8          11   6   6
 + 7680 a1   a2 ) r  + (6400 a1  a2   + 12800 a1  a2  + 6400 a1   a2 ) r )
        22              20                3          3      18
/(1600 r   - 900 a1 a2 r   + (- 1800 a1 a2  - 1800 a1  a2) r
               5           3   3         5      16
 + (- 100 a1 a2  - 22600 a1  a2  - 100 a1  a2) r
               7           3   5           4   4           5   3         7
 + (- 720 a1 a2  + 27600 a1  a2  + 12600 a1  a2  + 27600 a1  a2  - 720 a1  a2)
  14              3   7          4   6           5   5          6   4
 r   + (- 25600 a1  a2  - 5100 a1  a2  - 51200 a1  a2  - 5100 a1  a2
           7   3   12             4   8           6   6          8   4   10
 - 25600 a1  a2 ) r   + (- 8860 a1  a2  + 32600 a1  a2  - 8860 a1  a2 ) r
           4   10           6   8          7   7           8   6
 + (9580 a1  a2   - 15020 a1  a2  - 5625 a1  a2  - 15020 a1  a2
          10   4   8          4   12          6   10          7   9
 + 9580 a1   a2 ) r  + (100 a1  a2   + 6640 a1  a2   + 1500 a1  a2
           8   8          9   7          10   6         12   4   6
 + 13144 a1  a2  + 1500 a1  a2  + 6640 a1   a2  + 100 a1   a2 ) r
          7   11         9   9         11   7   4
 + (650 a1  a2   + 700 a1  a2  + 650 a1   a2 ) r
            7   13         9   11         11   9         13   7   2
 + (- 100 a1  a2   - 220 a1  a2   - 220 a1   a2  - 100 a1   a2 ) r
        7   15        9   13        11   11        13   9        15   7
 - 25 a1  a2   - 60 a1  a2   - 86 a1   a2   - 60 a1   a2  - 25 a1   a2 )
   */
  // YA12
  /*
                      21             3         3      19           4   4  15
(%o87) - (1200 a1 a2 r   + (400 a1 a2  + 400 a1  a2) r   - 16200 a1  a2  r
            4   6           6   4   13              4   8           6   6
 + (23200 a1  a2  + 23200 a1  a2 ) r   + (- 14200 a1  a2  - 29200 a1  a2
           8   4   11             4   10          6   8           7   7
 - 14200 a1  a2 ) r   + (- 1600 a1  a2   - 3520 a1  a2  + 15000 a1  a2
          8   6          10   4   9             7   9          9   7   7
 - 3520 a1  a2  - 1600 a1   a2 ) r  + (- 8000 a1  a2  - 8000 a1  a2 ) r
            7   11         9   9         11   7   5
 + (- 200 a1  a2   + 400 a1  a2  - 200 a1   a2 ) r
          7   13         9   11         11   9         13   7   3
 + (400 a1  a2   + 880 a1  a2   + 880 a1   a2  + 400 a1   a2 ) r )
        22              20                3          3      18
/(1600 r   - 900 a1 a2 r   + (- 1800 a1 a2  - 1800 a1  a2) r
               5           3   3         5      16
 + (- 100 a1 a2  - 22600 a1  a2  - 100 a1  a2) r
               7           3   5           4   4           5   3         7
 + (- 720 a1 a2  + 27600 a1  a2  + 12600 a1  a2  + 27600 a1  a2  - 720 a1  a2)
  14              3   7          4   6           5   5          6   4
 r   + (- 25600 a1  a2  - 5100 a1  a2  - 51200 a1  a2  - 5100 a1  a2
           7   3   12             4   8           6   6          8   4   10
 - 25600 a1  a2 ) r   + (- 8860 a1  a2  + 32600 a1  a2  - 8860 a1  a2 ) r
           4   10           6   8          7   7           8   6
 + (9580 a1  a2   - 15020 a1  a2  - 5625 a1  a2  - 15020 a1  a2
          10   4   8          4   12          6   10          7   9
 + 9580 a1   a2 ) r  + (100 a1  a2   + 6640 a1  a2   + 1500 a1  a2
           8   8          9   7          10   6         12   4   6
 + 13144 a1  a2  + 1500 a1  a2  + 6640 a1   a2  + 100 a1   a2 ) r
          7   11         9   9         11   7   4
 + (650 a1  a2   + 700 a1  a2  + 650 a1   a2 ) r
            7   13         9   11         11   9         13   7   2
 + (- 100 a1  a2   - 220 a1  a2   - 220 a1   a2  - 100 a1   a2 ) r
        7   15        9   13        11   11        13   9        15   7
 - 25 a1  a2   - 60 a1  a2   - 86 a1   a2   - 60 a1   a2  - 25 a1   a2 )
   */
  // YA21
  /*
                      21             3         3      19           4   4  15
(%o88) - (1200 a1 a2 r   + (400 a1 a2  + 400 a1  a2) r   - 16200 a1  a2  r
            4   6           6   4   13              4   8           6   6
 + (23200 a1  a2  + 23200 a1  a2 ) r   + (- 14200 a1  a2  - 29200 a1  a2
           8   4   11             4   10          6   8           7   7
 - 14200 a1  a2 ) r   + (- 1600 a1  a2   - 3520 a1  a2  + 15000 a1  a2
          8   6          10   4   9             7   9          9   7   7
 - 3520 a1  a2  - 1600 a1   a2 ) r  + (- 8000 a1  a2  - 8000 a1  a2 ) r
            7   11         9   9         11   7   5
 + (- 200 a1  a2   + 400 a1  a2  - 200 a1   a2 ) r
          7   13         9   11         11   9         13   7   3
 + (400 a1  a2   + 880 a1  a2   + 880 a1   a2  + 400 a1   a2 ) r )
        22              20                3          3      18
/(1600 r   - 900 a1 a2 r   + (- 1800 a1 a2  - 1800 a1  a2) r
               5           3   3         5      16
 + (- 100 a1 a2  - 22600 a1  a2  - 100 a1  a2) r
               7           3   5           4   4           5   3         7
 + (- 720 a1 a2  + 27600 a1  a2  + 12600 a1  a2  + 27600 a1  a2  - 720 a1  a2)
  14              3   7          4   6           5   5          6   4
 r   + (- 25600 a1  a2  - 5100 a1  a2  - 51200 a1  a2  - 5100 a1  a2
           7   3   12             4   8           6   6          8   4   10
 - 25600 a1  a2 ) r   + (- 8860 a1  a2  + 32600 a1  a2  - 8860 a1  a2 ) r
           4   10           6   8          7   7           8   6
 + (9580 a1  a2   - 15020 a1  a2  - 5625 a1  a2  - 15020 a1  a2
          10   4   8          4   12          6   10          7   9
 + 9580 a1   a2 ) r  + (100 a1  a2   + 6640 a1  a2   + 1500 a1  a2
           8   8          9   7          10   6         12   4   6
 + 13144 a1  a2  + 1500 a1  a2  + 6640 a1   a2  + 100 a1   a2 ) r
          7   11         9   9         11   7   4
 + (650 a1  a2   + 700 a1  a2  + 650 a1   a2 ) r
            7   13         9   11         11   9         13   7   2
 + (- 100 a1  a2   - 220 a1  a2   - 220 a1   a2  - 100 a1   a2 ) r
        7   15        9   13        11   11        13   9        15   7
 - 25 a1  a2   - 60 a1  a2   - 86 a1   a2   - 60 a1   a2  - 25 a1   a2 )
   */
  // YA22
  scalar_fts [7] = r6*
    (
     (6400.0*a27*a110
      +12800.0*a29*a18
      +6400.0*a211*a16) +r2*
     (+(-28800.0*a27*a18
	-22400.0*a29*a16
	+7680.0*a211*a14) +r2*
      (+(28000.0*a27*a16
	 -16800.0*a29*a14) +r2*
       (+(-25600.0*a24*a17
	  -51200.0*a26*a15
	  +12000.0*a27*a14
	  -25600.0*a28*a13) +r2*
	(+(30000.0*a24*a15
	   +29600.0*a26*a13
	   -720.0*a28*a1) +r2*
	 (-22400.0*a24*a13 +r2*
	  (-1200.0*a24*a1 +r2*r2*
	   (+1600.0*a2))))))))
    / dy;
  /*
                 22             4  18           3   4  16
(%o86) (1600 a2 r   - 1200 a1 a2  r   - 22400 a1  a2  r
               8           3   6           5   4   14
 + (- 720 a1 a2  + 29600 a1  a2  + 30000 a1  a2 ) r
              3   8           4   7           5   6           7   4   12
 + (- 25600 a1  a2  + 12000 a1  a2  - 51200 a1  a2  - 25600 a1  a2 ) r
            6   7           4   9   10           4   11           6   9
 + (28000 a1  a2  - 16800 a1  a2 ) r   + (7680 a1  a2   - 22400 a1  a2
           8   7   8           6   11           8   9          10   7   6
 - 28800 a1  a2 ) r  + (6400 a1  a2   + 12800 a1  a2  + 6400 a1   a2 ) r )
        22              20                3          3      18
/(1600 r   - 900 a1 a2 r   + (- 1800 a1 a2  - 1800 a1  a2) r
               5           3   3         5      16
 + (- 100 a1 a2  - 22600 a1  a2  - 100 a1  a2) r
               7           3   5           4   4           5   3         7
 + (- 720 a1 a2  + 27600 a1  a2  + 12600 a1  a2  + 27600 a1  a2  - 720 a1  a2)
  14              3   7          4   6           5   5          6   4
 r   + (- 25600 a1  a2  - 5100 a1  a2  - 51200 a1  a2  - 5100 a1  a2
           7   3   12             4   8           6   6          8   4   10
 - 25600 a1  a2 ) r   + (- 8860 a1  a2  + 32600 a1  a2  - 8860 a1  a2 ) r
           4   10           6   8          7   7           8   6
 + (9580 a1  a2   - 15020 a1  a2  - 5625 a1  a2  - 15020 a1  a2
          10   4   8          4   12          6   10          7   9
 + 9580 a1   a2 ) r  + (100 a1  a2   + 6640 a1  a2   + 1500 a1  a2
           8   8          9   7          10   6         12   4   6
 + 13144 a1  a2  + 1500 a1  a2  + 6640 a1   a2  + 100 a1   a2 ) r
          7   11         9   9         11   7   4
 + (650 a1  a2   + 700 a1  a2  + 650 a1   a2 ) r
            7   13         9   11         11   9         13   7   2
 + (- 100 a1  a2   - 220 a1  a2   - 220 a1   a2  - 100 a1   a2 ) r
        7   15        9   13        11   11        13   9        15   7
 - 25 a1  a2   - 60 a1  a2   - 86 a1   a2   - 60 a1   a2  - 25 a1   a2 )
   */
}


/** check routines **/

/* check scalars_minv_f_poly() with scalar_minv_f() for equal sphere
 */
int
check_scalars_minv_f_poly_with_equal (double r, int verbose, double tiny)
{
  int check = 0;

  double mono[4];
  double poly[8];

  scalar_minv_f (r, mono);
  scalars_minv_f_poly (r, 1.0, 1.0, poly);
  // poly gives the results in dimensional form but for a1=1, it's the same

  check += compare (mono [0], poly [0], "check_scalars_minv_f_poly_with_equal: XA11", verbose, tiny);
  check += compare (mono [1], poly [1], "check_scalars_minv_f_poly_with_equal: XA12", verbose, tiny);
  check += compare (mono [2], poly [4], "check_scalars_minv_f_poly_with_equal: YA11", verbose, tiny);
  check += compare (mono [3], poly [5], "check_scalars_minv_f_poly_with_equal: YA12", verbose, tiny);

  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_scalars_minv_f_poly_with_equal : PASSED\n");
    }

  return (check);
}

/* check scalars_minv_ft_poly() with scalar_minv_ft() for equal sphere
 */
int
check_scalars_minv_ft_poly_with_equal (double r, int verbose, double tiny)
{
  int check = 0;

  double mono[10];
  double poly[20];

  scalar_minv_ft (r, mono);
  scalars_minv_ft_poly (r, 1.0, 1.0, poly);
  // poly gives the results in dimensional form but for a1=1, it's the same

  check += compare (mono [0], poly [0], "check_scalars_minv_ft_poly_with_equal: XA11", verbose, tiny);
  check += compare (mono [1], poly [1], "check_scalars_minv_ft_poly_with_equal: XA12", verbose, tiny);
  check += compare (mono [2], poly [4], "check_scalars_minv_ft_poly_with_equal: YA11", verbose, tiny);
  check += compare (mono [3], poly [5], "check_scalars_minv_ft_poly_with_equal: YA12", verbose, tiny);
  check += compare (mono [4], poly [8], "check_scalars_minv_ft_poly_with_equal: YB11", verbose, tiny);
  check += compare (mono [5], poly [9], "check_scalars_minv_ft_poly_with_equal: YB12", verbose, tiny);
  check += compare (mono [6], poly[12], "check_scalars_minv_ft_poly_with_equal: XC11", verbose, tiny);
  check += compare (mono [7], poly[13], "check_scalars_minv_ft_poly_with_equal: XC12", verbose, tiny);
  check += compare (mono [8], poly[16], "check_scalars_minv_ft_poly_with_equal: YC11", verbose, tiny);
  check += compare (mono [9], poly[17], "check_scalars_minv_ft_poly_with_equal: YC12", verbose, tiny);

  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_scalars_minv_ft_poly_with_equal : PASSED\n");
    }

  return (check);
}

/* check scalars_minv_fts_poly() with scalar_minv_fts() for equal sphere
 */
int
check_scalars_minv_fts_poly_with_equal (double r, int verbose, double tiny)
{
  int check = 0;

  double mono[22];
  double poly[44];

  scalar_minv_fts (r, mono);
  scalars_minv_fts_poly (r, 1.0, 1.0, poly);
  // poly gives the results in dimensional form but for a1=1, it's the same

  check += compare (mono [0], poly [0], "check_scalars_minv_fts_poly_with_equal: XA11", verbose, tiny);
  check += compare (mono [1], poly [1], "check_scalars_minv_fts_poly_with_equal: XA12", verbose, tiny);
  check += compare (mono [2], poly [4], "check_scalars_minv_fts_poly_with_equal: YA11", verbose, tiny);
  check += compare (mono [3], poly [5], "check_scalars_minv_fts_poly_with_equal: YA12", verbose, tiny);
  check += compare (mono [4], poly [8], "check_scalars_minv_fts_poly_with_equal: YB11", verbose, tiny);
  check += compare (mono [5], poly [9], "check_scalars_minv_fts_poly_with_equal: YB12", verbose, tiny);
  check += compare (mono [6], poly[12], "check_scalars_minv_fts_poly_with_equal: XC11", verbose, tiny);
  check += compare (mono [7], poly[13], "check_scalars_minv_fts_poly_with_equal: XC12", verbose, tiny);
  check += compare (mono [8], poly[16], "check_scalars_minv_fts_poly_with_equal: YC11", verbose, tiny);
  check += compare (mono [9], poly[17], "check_scalars_minv_fts_poly_with_equal: YC12", verbose, tiny);
  check += compare (mono[10], poly[20], "check_scalars_minv_fts_poly_with_equal: XG11", verbose, tiny);
  check += compare (mono[11], poly[21], "check_scalars_minv_fts_poly_with_equal: XG12", verbose, tiny);
  check += compare (mono[12], poly[24], "check_scalars_minv_fts_poly_with_equal: YG11", verbose, tiny);
  check += compare (mono[13], poly[25], "check_scalars_minv_fts_poly_with_equal: YG12", verbose, tiny);
  check += compare (mono[14], poly[28], "check_scalars_minv_fts_poly_with_equal: YH11", verbose, tiny);
  check += compare (mono[15], poly[29], "check_scalars_minv_fts_poly_with_equal: YH12", verbose, tiny);
  check += compare (mono[16], poly[32], "check_scalars_minv_fts_poly_with_equal: XM11", verbose, tiny);
  check += compare (mono[17], poly[33], "check_scalars_minv_fts_poly_with_equal: XM12", verbose, tiny);
  check += compare (mono[18], poly[36], "check_scalars_minv_fts_poly_with_equal: YM11", verbose, tiny);
  check += compare (mono[19], poly[37], "check_scalars_minv_fts_poly_with_equal: YM12", verbose, tiny);
  check += compare (mono[20], poly[40], "check_scalars_minv_fts_poly_with_equal: ZM11", verbose, tiny);
  check += compare (mono[21], poly[41], "check_scalars_minv_fts_poly_with_equal: ZM12", verbose, tiny);

  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_scalars_minv_fts_poly_with_equal : PASSED\n");
    }

  return (check);
}


/* check scalars_minv_f_poly() with scalars_minv_f_poly_ana()
 * for non-equal spheres
 */
int
check_scalars_minv_f_poly_ana (double r, double a1, double a2,
			       int verbose, double tiny)
{
  int check = 0;

  double poly[8];
  double ana [8];

  scalars_minv_f_poly     (r, a1, a2, poly);
  scalars_minv_f_poly_ana (r, a1, a2, ana);
  // poly's give the results in dimensional form

  check += compare (ana [0], poly [0], "check_scalars_minv_f_poly_ana: XA11", verbose, tiny);
  check += compare (ana [1], poly [1], "check_scalars_minv_f_poly_ana: XA12", verbose, tiny);
  check += compare (ana [2], poly [2], "check_scalars_minv_f_poly_ana: XA21", verbose, tiny);
  check += compare (ana [3], poly [3], "check_scalars_minv_f_poly_ana: XA22", verbose, tiny);
  check += compare (ana [4], poly [4], "check_scalars_minv_f_poly_ana: YA11", verbose, tiny);
  check += compare (ana [5], poly [5], "check_scalars_minv_f_poly_ana: YA12", verbose, tiny);
  check += compare (ana [6], poly [6], "check_scalars_minv_f_poly_ana: YA21", verbose, tiny);
  check += compare (ana [7], poly [7], "check_scalars_minv_f_poly_ana: YA22", verbose, tiny);

  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_scalars_minv_f_poly_ana : PASSED\n");
    }

  return (check);
}

/* check scalars_minv_ft_poly() with scalars_minv_ft_poly_ana()
 * for non-equal spheres
 */
int
check_scalars_minv_ft_poly_ana (double r, double a1, double a2,
				int verbose, double tiny)
{
  int check = 0;

  double poly[20];
  double ana [20];

  scalars_minv_ft_poly     (r, a1, a2, poly);
  scalars_minv_ft_poly_ana (r, a1, a2, ana);
  // poly's give the results in dimensional form

  check += compare (ana [0], poly [0], "check_scalars_minv_ft_poly_ana: XA11", verbose, tiny);
  check += compare (ana [1], poly [1], "check_scalars_minv_ft_poly_ana: XA12", verbose, tiny);
  check += compare (ana [2], poly [2], "check_scalars_minv_ft_poly_ana: XA21", verbose, tiny);
  check += compare (ana [3], poly [3], "check_scalars_minv_ft_poly_ana: XA22", verbose, tiny);
  check += compare (ana [4], poly [4], "check_scalars_minv_ft_poly_ana: YA11", verbose, tiny);
  check += compare (ana [5], poly [5], "check_scalars_minv_ft_poly_ana: YA12", verbose, tiny);
  check += compare (ana [6], poly [6], "check_scalars_minv_ft_poly_ana: YA21", verbose, tiny);
  check += compare (ana [7], poly [7], "check_scalars_minv_ft_poly_ana: YA22", verbose, tiny);
  check += compare (ana [8], poly [8], "check_scalars_minv_ft_poly_ana: YB11", verbose, tiny);
  check += compare (ana [9], poly [9], "check_scalars_minv_ft_poly_ana: YB12", verbose, tiny);
  check += compare (ana[10], poly[10], "check_scalars_minv_ft_poly_ana: YB21", verbose, tiny);
  check += compare (ana[11], poly[11], "check_scalars_minv_ft_poly_ana: YB22", verbose, tiny);
  check += compare (ana[12], poly[12], "check_scalars_minv_ft_poly_ana: XC11", verbose, tiny);
  check += compare (ana[13], poly[13], "check_scalars_minv_ft_poly_ana: XC12", verbose, tiny);
  check += compare (ana[14], poly[14], "check_scalars_minv_ft_poly_ana: XC21", verbose, tiny);
  check += compare (ana[15], poly[15], "check_scalars_minv_ft_poly_ana: XC22", verbose, tiny);
  check += compare (ana[16], poly[16], "check_scalars_minv_ft_poly_ana: YC11", verbose, tiny);
  check += compare (ana[17], poly[17], "check_scalars_minv_ft_poly_ana: YC12", verbose, tiny);
  check += compare (ana[18], poly[18], "check_scalars_minv_ft_poly_ana: YC21", verbose, tiny);
  check += compare (ana[19], poly[19], "check_scalars_minv_ft_poly_ana: YC22", verbose, tiny);

  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_scalars_minv_ft_poly_ana : PASSED\n");
    }

  return (check);
}

/* check scalars_minv_fts_poly() with scalars_minv_fts_poly_ana()
 * for non-equal spheres
 */
int
check_scalars_minv_fts_poly_ana (double r, double a1, double a2,
				 int verbose, double tiny)
{
  int check = 0;

  double poly[44];
  double ana [44];

  scalars_minv_fts_poly     (r, a1, a2, poly);
  scalars_minv_fts_poly_ana (r, a1, a2, ana);
  // poly's give the results in dimensional form

  //check += compare (ana [0], poly [0], "check_scalars_minv_fts_poly_ana: XA11", verbose, tiny);
  //check += compare (ana [1], poly [1], "check_scalars_minv_fts_poly_ana: XA12", verbose, tiny);
  check += compare (ana [4], poly [4], "check_scalars_minv_fts_poly_ana: YA11", verbose, tiny);
  //check += compare (ana [5], poly [5], "check_scalars_minv_fts_poly_ana: YA12", verbose, tiny);
  //check += compare (ana [6], poly [6], "check_scalars_minv_fts_poly_ana: YA21", verbose, tiny);
  check += compare (ana [7], poly [7], "check_scalars_minv_fts_poly_ana: YA22", verbose, tiny);
  //check += compare (ana [8], poly [8], "check_scalars_minv_fts_poly_ana: YB11", verbose, tiny);
  //check += compare (ana [9], poly [9], "check_scalars_minv_fts_poly_ana: YB12", verbose, tiny);
  check += compare (ana[12], poly[12], "check_scalars_minv_fts_poly_ana: XC11", verbose, tiny);
  check += compare (ana[13], poly[13], "check_scalars_minv_fts_poly_ana: XC12", verbose, tiny);
  check += compare (ana[14], poly[14], "check_scalars_minv_fts_poly_ana: XC21", verbose, tiny);
  check += compare (ana[15], poly[15], "check_scalars_minv_fts_poly_ana: XC22", verbose, tiny);
  //check += compare (ana[16], poly[16], "check_scalars_minv_fts_poly_ana: YC11", verbose, tiny);
  //check += compare (ana[17], poly[17], "check_scalars_minv_fts_poly_ana: YC12", verbose, tiny);
  //check += compare (ana[20], poly[20], "check_scalars_minv_fts_poly_ana: XG11", verbose, tiny);
  //check += compare (ana[21], poly[21], "check_scalars_minv_fts_poly_ana: XG12", verbose, tiny);
  //check += compare (ana[24], poly[24], "check_scalars_minv_fts_poly_ana: YG11", verbose, tiny);
  //check += compare (ana[25], poly[25], "check_scalars_minv_fts_poly_ana: YG12", verbose, tiny);
  //check += compare (ana[28], poly[28], "check_scalars_minv_fts_poly_ana: YH11", verbose, tiny);
  //check += compare (ana[29], poly[29], "check_scalars_minv_fts_poly_ana: YH12", verbose, tiny);
  //check += compare (ana[32], poly[32], "check_scalars_minv_fts_poly_ana: XM11", verbose, tiny);
  //check += compare (ana[33], poly[33], "check_scalars_minv_fts_poly_ana: XM12", verbose, tiny);
  //check += compare (ana[36], poly[36], "check_scalars_minv_fts_poly_ana: YM11", verbose, tiny);
  //check += compare (ana[37], poly[37], "check_scalars_minv_fts_poly_ana: YM12", verbose, tiny);
  check += compare (ana[40], poly[40], "check_scalars_minv_fts_poly_ana: ZM11", verbose, tiny);
  check += compare (ana[41], poly[41], "check_scalars_minv_fts_poly_ana: ZM12", verbose, tiny);
  check += compare (ana[42], poly[42], "check_scalars_minv_fts_poly_ana: ZM21", verbose, tiny);
  check += compare (ana[43], poly[43], "check_scalars_minv_fts_poly_ana: ZM22", verbose, tiny);

  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_scalars_minv_fts_poly_ana : PASSED\n");
    }

  return (check);
}
