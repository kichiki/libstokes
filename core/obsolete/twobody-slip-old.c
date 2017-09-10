/* backup of bug fixing for polydisperse systems
 * twobody solutions for slip particles
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h> // sqrt()
#include <stokes.h> // struct stokes
#include <twobody.h> // twobody_scale(), twobody_scale_SD()
#include <minv-poly.h> // scalars_minv_f_poly()
#include "memory-check.h" // macro CHECK_MALLOC

#include "f.h"
#include "ft.h"
#include "fts.h"

#include "twobody-slip-old.h"


/* calculate f by u for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 * INPUT
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin2      : square of min distance for lub calculation.
 *         sys->twobody_nmax : max order in twobody.
 *         sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   u1 [3] : velocity of particle 1
 *   u2 [3] : velocity of particle 2
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 * OUTPUT
 *   f1 [3] : force of particle 1
 *   f2 [3] : force of particle 2
 */
void
calc_lub_f_2b_slip_old
(struct stokes *sys,
 const double *u1, const double *u2,
 const double *x1, const double *x2,
 int i1, int i2,
 double *f1, double *f2)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
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
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (0, // F version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (0, // F version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  xa11 = lub [0];
  xa12 = lub [1];
  xa21 = lub [2];
  xa22 = lub [3];
  ya11 = lub [4];
  ya12 = lub [5];
  ya21 = lub [6];
  ya22 = lub [7];

  matrix_f_atimes (u1, f1,
		   ex, ey, ez,
		   xa11, ya11);
  matrix_f_atimes (u2, f1,
		   ex, ey, ez,
		   xa12, ya12);

  matrix_f_atimes (u2, f2,
		   -ex, -ey, -ez,
		   xa22, ya22);
  matrix_f_atimes (u1, f2,
		   -ex, -ey, -ez,
		   xa21, ya21);
}

/* calculate lub-matrix in F version for pair of unequal spheres 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ matrix_lub_f_2b(i,j); }}
 * INPUT
 *   sys    : system parameters. the followings are referred:
 *            sys->lubmin2      : square of min distance for lub calculation.
 *            sys->twobody_nmax : max order in twobody.
 *            sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   i      : particle index for '1'
 *   j      : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *   n      : dimension of matrix 'mat' (must be np*3)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_f_2b_slip_old
(struct stokes *sys,
 int i, int j,
 const double *x1, const double *x2,
 int i1, int i2,
 int n, double *mat)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
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
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (0, // F version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (0, // F version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  xa11 = lub [0];
  xa12 = lub [1];
  xa21 = lub [2];
  xa22 = lub [3];
  ya11 = lub [4];
  ya12 = lub [5];
  ya21 = lub [6];
  ya22 = lub [7];

  matrix_f_ij (i, i,
	       ex, ey, ez,
	       xa11, ya11,
	       n, mat);
  matrix_f_ij (i, j,
	       ex, ey, ez,
	       xa12, ya12,
	       n, mat);

  matrix_f_ij (j, j,
	       -ex, -ey, -ez,
	       xa22, ya22,
	       n, mat);
  matrix_f_ij (j, i,
	       -ex, -ey, -ez,
	       xa21, ya21,
	       n, mat);
}


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
calc_lub_ft_2b_slip_old
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

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
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
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (1, // FT version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (1, // FT version
			     a1, a2, lub);

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

/* calculate lub-matrix in FT version for pair of unequal spheres 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ matrix_lub_f_2b(i,j); }}
 * INPUT
 *   sys    : system parameters. the followings are referred:
 *            sys->lubmin2      : square of min distance for lub calculation.
 *            sys->twobody_nmax : max order in twobody.
 *            sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   i      : particle index for '1'
 *   j      : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *   n      : dimension of matrix 'mat' (must be np*6)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_ft_2b_slip_old
(struct stokes *sys,
 int i, int j,
 const double *x1, const double *x2,
 int i1, int i2,
 int n, double *mat)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
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
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (1, // FT version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (1, // FT version
			     a1, a2, lub);

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

  matrix_ft_ij (i, i,
		ex, ey, ez,
		xa11, ya11,
		yb11,
		xc11, yc11,
		n, mat);
  matrix_ft_ij (i, j,
		ex, ey, ez,
		xa12, ya12,
		yb12,
		xc12, yc12,
		n, mat);

  matrix_ft_ij (j, j,
		-ex, -ey, -ez,
		xa22, ya22,
		yb22,
		xc22, yc22,
		n, mat);
  matrix_ft_ij (j, i,
		-ex, -ey, -ez,
		xa21, ya21,
		yb21,
		xc21, yc21,
		n, mat);
}


/* calculate fts by uoe for pair of particles 1 and 2 for unequal spheres
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 *   sys : system parameters. the followings are referred:
 *         sys->lubmin2      : square of min distance for lub calculation.
 *         sys->twobody_nmax : max order in twobody.
 *         sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   uoe1 [11] : velocity, angular velocity, strain
 *   uoe2 [11] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 * OUTPUT
 *   fts1 [11] : force, torque, stresslet
 *   fts2 [11] :
 */
void
calc_lub_fts_2b_slip_old
(struct stokes *sys,
 const double *uoe1, const double *uoe2,
 const double *x1, const double *x2,
 int i1, int i2,
 double *fts1, double *fts2)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
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
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (2, // FTS version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (2, // FTS version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  double yb11, yb12, yb21, yb22;
  double xc11, xc12, xc21, xc22;
  double yc11, yc12, yc21, yc22;
  double xg11, xg12, xg21, xg22;
  double yg11, yg12, yg21, yg22;
  double yh11, yh12, yh21, yh22;
  double xm11, xm12, xm21, xm22;
  double ym11, ym12, ym21, ym22;
  double zm11, zm12, zm21, zm22;
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
  xg11 = lub[20];
  xg12 = lub[21];
  xg21 = lub[22];
  xg22 = lub[23];
  yg11 = lub[24];
  yg12 = lub[25];
  yg21 = lub[26];
  yg22 = lub[27];
  yh11 = lub[28];
  yh12 = lub[29];
  yh21 = lub[30];
  yh22 = lub[31];
  xm11 = lub[32];
  xm12 = lub[33];
  xm21 = lub[34];
  xm22 = lub[35];
  ym11 = lub[36];
  ym12 = lub[37];
  ym21 = lub[38];
  ym22 = lub[39];
  zm11 = lub[40];
  zm12 = lub[41];
  zm21 = lub[42];
  zm22 = lub[43];

  matrix_fts_self_atimes (uoe1, fts1,
			  ex, ey, ez,
			  xa11, ya11,
			  yb11,
			  xc11, yc11,
			  xg11, yg11,
			  yh11,
			  xm11, ym11, zm11);
  matrix_fts_atimes (uoe2, fts1,
		     ex, ey, ez,
		     xa12, ya12,
		     yb12,
		     xc12, yc12,
		     xg12, yg12,
		     yh12,
		     xm12, ym12, zm12);

  matrix_fts_self_atimes (uoe2, fts2,
			  -ex, -ey, -ez,
			  xa22, ya22,
			  yb22,
			  xc22, yc22,
			  xg22, yg22,
			  yh22,
			  xm22, ym22, zm22);
  matrix_fts_atimes (uoe1, fts2,
		     -ex, -ey, -ez,
		     xa21, ya21,
		     yb21,
		     xc21, yc21,
		     xg21, yg21,
		     yh21,
		     xm21, ym21, zm21);
}

/* calculate lub-matrix in FTS version for pair of unequal spheres 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ matrix_lub_f_2b(i,j); }}
 * INPUT
 *   sys    : system parameters. the followings are referred:
 *            sys->lubmin2      : square of min distance for lub calculation.
 *            sys->twobody_nmax : max order in twobody.
 *            sys->twobody_lub  : 0 for far form, 1 for lub form in twobody.
 *   i      : particle index for '1'
 *   j      : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   i1, i2 : particle index for particles 1 and 2
 *   n      : dimension of matrix 'mat' (must be np*11)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_fts_2b_slip_old
(struct stokes *sys,
 int i, int j,
 const double *x1, const double *x2,
 int i1, int i2,
 int n, double *mat)
{
  /* r := x[j] - x[i] for (j -> i) interaction */
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;

  double a1;
  double a2;
  if (sys->a == NULL)
    {
      a1 = 1.0;
      a2 = 1.0;
    }
  else
    {
      a1 = sys->a[i1];
      a2 = sys->a[i2];
    }
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
  struct twobody_slip_f *f12
    = sys->twobody_slip_f_list->f[sys->slip_table[i1*sys->np+i2]];
  struct twobody_slip_f *f21
    = sys->twobody_slip_f_list->f[sys->slip_table[i2*sys->np+i1]];

  double lub [44];
  scalars_lub_slip_full (2, // FTS version
			 rr, a1, a2,
			 f12, f21,
			 sys->twobody_nmax,
			 sys->twobody_lub,
			 lub);
  scalars_res_poly_scale_SD (2, // FTS version
			     a1, a2, lub);

  double xa11, xa12, xa21, xa22;
  double ya11, ya12, ya21, ya22;
  double yb11, yb12, yb21, yb22;
  double xc11, xc12, xc21, xc22;
  double yc11, yc12, yc21, yc22;
  double xg11, xg12, xg21, xg22;
  double yg11, yg12, yg21, yg22;
  double yh11, yh12, yh21, yh22;
  double xm11, xm12, xm21, xm22;
  double ym11, ym12, ym21, ym22;
  double zm11, zm12, zm21, zm22;
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
  xg11 = lub[20];
  xg12 = lub[21];
  xg21 = lub[22];
  xg22 = lub[23];
  yg11 = lub[24];
  yg12 = lub[25];
  yg21 = lub[26];
  yg22 = lub[27];
  yh11 = lub[28];
  yh12 = lub[29];
  yh21 = lub[30];
  yh22 = lub[31];
  xm11 = lub[32];
  xm12 = lub[33];
  xm21 = lub[34];
  xm22 = lub[35];
  ym11 = lub[36];
  ym12 = lub[37];
  ym21 = lub[38];
  ym22 = lub[39];
  zm11 = lub[40];
  zm12 = lub[41];
  zm21 = lub[42];
  zm22 = lub[43];

  matrix_fts_ij (i, i,
		 ex, ey, ez,
		 xa11, ya11,
		 yb11,
		 xc11, yc11,
		 xg11, yg11,
		 yh11,
		 xm11, ym11, zm11,
		 n, mat);
  matrix_fts_ij (i, j,
		 ex, ey, ez,
		 xa12, ya12,
		 yb12,
		 xc12, yc12,
		 xg12, yg12,
		 yh12,
		 xm12, ym12, zm12,
		 n, mat);

  matrix_fts_ij (j, j,
		 -ex, -ey, -ez,
		 xa22, ya22,
		 yb22,
		 xc22, yc22,
		 xg22, yg22,
		 yh22,
		 xm22, ym22, zm22,
		 n, mat);
  matrix_fts_ij (j, i,
		 -ex, -ey, -ez,
		 xa21, ya21,
		 yb21,
		 xc21, yc21,
		 xg21, yg21,
		 yh21,
		 xm21, ym21, zm21,
		 n, mat);
}