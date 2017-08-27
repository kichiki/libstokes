/* bug fixing for polydisperse systems
 * subroutine for the procedure of FT version
 * Copyright (C) 2000-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include <math.h> // sqrt ()

#include "stokes.h" /* struct stokeks */
#include "minv-poly.h" // scalars_lub_poly_full()
#include "twobody.h" // struct twobody_f

#include "ft.h"

#include "ft-new.h"


/* calculate ft by uo for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 * INPUT
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin2      : square of min distance for lub calculation.
 *         sys->twobody_nmax : max order in twobody.
 *         sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   uo1 [6] : velocity, angular velocity
 *   uo2 [6] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 * OUTPUT
 *   ft1 [6] : force, torque
 *   ft2 [6] :
 */
void
calc_lub_ft_2b_poly_new
(struct stokes *sys,
 const double *uo1, const double *uo2,
 const double *x1, const double *x2,
 int i1, int i2,
 double *ft1, double *ft2)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1 = sys->a[i1];
  double a2 = sys->a[i2];
  double rs;
  rs = a1 + a2;
  rs *= rs; // = (a1 + a2)^2
  rs *= 0.25; // = (a1 + a2)^2 / 4
  double s2 = r2 / rs;
  if (s2 < sys->lubmin2)
    {
      s2 = sys->lubmin2;
      r2 = rs * s2;
    }

  double rr = sqrt (r2);
  double ex = xx / rr;
  double ey = yy / rr;
  double ez = zz / rr;

  /* calc scalar functions of lubrication */
  struct twobody_f *f12
    = sys->twobody_f_list->f[sys->poly_table[i1*sys->np+i2]];
  struct twobody_f *f21
    = sys->twobody_f_list->f[sys->poly_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_poly_full (1, // FT version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  /*
  scalars_res_poly_scale_SD (1, // FT version
			     a1, a2, lub);
  */

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  double yb11, yb12, yb21, yb22;
  double xc11, xc12, xc21, xc22;
  double yc11, yc12, yc21, yc22;
  xa11 = lub [0];
  xa12 = lub [1];
  xa21 = lub [2];
  xa22 = lub [3];
  ya11 = lub [4];
  ya12 = lub [5];
  ya21 = lub [6];
  ya22 = lub [7];
  yb11 = lub [8];
  yb12 = lub [9];
  yb21 = lub[10];
  yb22 = lub[11];
  xc11 = lub[12];
  xc12 = lub[13];
  xc21 = lub[14];
  xc22 = lub[15];
  yc11 = lub[16];
  yc12 = lub[17];
  yc21 = lub[18];
  yc22 = lub[19];

  matrix_ft_self_atimes (uo1, ft1,
			 ex, ey, ez,
			 xa11, ya11,
			 yb11,
			 xc11, yc11);
  matrix_ft_atimes (uo2, ft1,
		    ex, ey, ez,
		    xa12, ya12,
		    yb12,
		    xc12, yc12);

  matrix_ft_self_atimes (uo2, ft2,
			 -ex, -ey, -ez,
			 xa22, ya22,
			 yb22,
			 xc22, yc22);
  matrix_ft_atimes (uo1, ft2,
		    -ex, -ey, -ez,
		    xa21, ya21,
		    yb21,
		    xc21, yc21);
}
