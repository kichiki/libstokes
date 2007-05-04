/* subroutine for the procedure of F version
 * Copyright (C) 2001-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: f.c,v 2.11 2007/05/04 01:12:14 kichiki Exp $
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
#include <stdio.h> // fprintf ()
#include <stdlib.h> // malloc ()
#include <math.h> // sqrt ()

#include "two-body-res.h" /* scalar_two_body_res () */
#include "stokes.h" /* struct stokeks */
#include "minv-poly.h" // scalars_lub_poly_full()
#include "twobody.h" // struct twobody_f

#include "f.h"

/* store matrix in F format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   i, j : particle index
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xa, ya, ... : scalar functions
 *   n3 : dimension of the matrix mat []
 * OUTPUT
 *   mat [n3 * n3] :
 */
void
matrix_f_ij (int i, int j,
	     double ex, double ey, double ez,
	     double xa, double ya,
	     int n3, double *mat)
{
  int i1;
  int j1;

  i1  = i * 3;

  j1  = j * 3;

  matrix_ij_A (n3, & mat [i1 * n3 + j1],
	       ex, ey, ez,
	       xa, ya);
}

/* store matrix in A part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xa, ya : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2) * n + (0,1,2)] : added (not cleared!)
 */
void
matrix_ij_A (int n, double *mat,
	     double ex, double ey, double ez,
	     double xa, double ya)
{
  double a1, a2;

  double exx, eyy, ezz, exy, eyz, exz;

  a1 = ya;
  a2 = xa - ya;

  exx = ex * ex;
  eyy = ey * ey;
  ezz = ez * ez;
  exy = ex * ey;
  eyz = ey * ez;
  exz = ez * ex;

  /* A part */
  mat [0 * n + 0] += a1 + a2 * exx;
  mat [1 * n + 1] += a1 + a2 * eyy;
  mat [2 * n + 2] += a1 + a2 * ezz;
  
  mat [0 * n + 1] += a2 * exy;
  mat [1 * n + 0] += a2 * exy;
  mat [1 * n + 2] += a2 * eyz;
  mat [2 * n + 1] += a2 * eyz;
  mat [2 * n + 0] += a2 * exz;
  mat [0 * n + 2] += a2 * exz;
}

/* ATIMES version (for O(N^2) scheme) of
 * store matrix in F format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * NOTE that only 'alpha(i) <- beta(j)' interaction is stored.
 * INPUT
 *   x [3] : F of particle 'i'
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for y[] and 'j' is for x[].
 *   xa, ya, ... : scalar functions
 * OUTPUT
 *   y [3] : U of particle 'j'
 */
void
matrix_f_atimes (const double *x,
		 double *y,
		 double ex, double ey, double ez,
		 double xa, double ya)
{
  double a1, a2;

  double exx, eyy, ezz, exy, eyz, exz;


  a1 = ya;
  a2 = xa - ya;

  exx = ex * ex;
  eyy = ey * ey;
  ezz = ez * ez;
  exy = ex * ey;
  eyz = ey * ez;
  exz = ez * ex;

  /* A part */
  y [ 0] += x [ 0] * (a1 + a2 * exx);
  y [ 1] += x [ 1] * (a1 + a2 * eyy);
  y [ 2] += x [ 2] * (a1 + a2 * ezz);
  
  y [ 0] += x [ 1] * (a2 * exy);
  y [ 1] += x [ 0] * (a2 * exy);
  y [ 1] += x [ 2] * (a2 * eyz);
  y [ 2] += x [ 1] * (a2 * eyz);
  y [ 2] += x [ 0] * (a2 * exz);
  y [ 0] += x [ 2] * (a2 * exz);
}

/* convert f[] to f[] (this is applicable for U)
 * INPUT
 *  n : # particles
 *  f2 [n * 3] :
 * OUTPUT
 *  f1 [n * 3] := f2 [n * 3]
 */
void
set_F_by_f (int n,
	    double *f1,
	    const double *f2)
{
  int i;
  int j;
  int i3;


  for (i = 0; i < n; i ++)
    {
      i3 = i * 3;
      for (j = 0; j < 3; j ++)
	{
	  f1 [i3 + j] = f2 [i3 + j];
	}
    }
}

/* calc scalar functions of (M^inf)^-1 in F
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *  scalar_f [4] :
 */
void
scalar_minv_f (double s, double * scalar_f)
{
  double xa11, xa12;
  double ya11, ya12;

  double s2, s3, s6;
  double dx, dy;


  s2 = s * s;
  s3 = s * s2;
  s6 = s3 * s3;

  dx = 4.0 * s6 - (3.0 * s2 - 2.0) * (3.0 * s2 - 2.0);
  dy = (((16.0 * s2 - 9.0) * s2 - 12.0 ) * s2 - 4.0);

  xa11 =
    4.0 * s6
    / dx;
  xa12 =
    - 2.0 * s3 * (3.0 * s2 - 2.0)
    / dx;
  ya11 =
    16.0 * s6
    / dy;
  ya12 =
    - 4.0 * s3 * (3.0 * s2 + 2.0)
    / dy;

  scalar_f [0] = xa11;
  scalar_f [1] = xa12;
  scalar_f [2] = ya11;
  scalar_f [3] = ya12;
}

/* calculate f by u for pair of particles 1 and 2
 * Note that this take care of both (12)- and (21)-interactions,
 * so that this is called in the loop
 *   for(i=0;i<n;i++){ for(j=i+1;j<n;j++){ calc_lub_f_2b(i,j); }}
 * INPUT
 *   sys : system parameters
 *         sys->lubmin2 is used.
 *   u1 [3] : velocity
 *   u2 [3] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 * OUTPUT
 *   f1 [3] : force of particle 1
 *   f2 [3] : force of particle 2
 */
void
calc_lub_f_2b (struct stokes * sys,
	       const double * u1, const double * u2,
	       const double * x1, const double * x2,
	       double * f1, double * f2)
{
  // r := x[j] - x[i] for (j -> i) interaction
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;
  if (r2 < sys->lubmin2)
    {
      r2 = sys->lubmin2;
    }

  double rr = sqrt (r2);
  double ex = xx / rr;
  double ey = yy / rr;
  double ez = zz / rr;

  // calc scalar functions of lubrication
  double res2b [22];
  double resinf[22];
  scalar_two_body_res (rr, res2b);
  /*
  // checking
  twobody_scalars_res (0, // F
		       rr, 1.0, 1.0, NULL, 100,
		       1, 1, res2b);
  // for a=1 in dimensional form is equivalent to SD form in mono
  */
  scalar_minv_f (rr, resinf);

  double xa11, ya11;
  double xa12, ya12;
  xa11 = res2b [0] - resinf [0];
  xa12 = res2b [1] - resinf [1];
  ya11 = res2b [2] - resinf [2];
  ya12 = res2b [3] - resinf [3];

  matrix_f_atimes (u1, f1,
		   ex, ey, ez,
		   xa11, ya11);
  matrix_f_atimes (u2, f1,
		   ex, ey, ez,
		   xa12, ya12);

  matrix_f_atimes (u2, f2,
		   -ex, -ey, -ez,
		   xa11, ya11);
  matrix_f_atimes (u1, f2,
		   -ex, -ey, -ez,
		   xa12, ya12);
}

/* calculate lub-matrix in F version for pair of particles 1 and 2
 * INPUT
 *   sys : system parameters
 *         sys->lubmin2 is used.
 *   i : particle index for '1'
 *   j : particle index for '2'
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 *   n : dimension of matrix 'mat' (must be np*3)
 * OUTPUT
 *   mat [n * n] : add for (i,j)- and (j,i)-pairs.
 */
void
matrix_lub_f_2b (struct stokes * sys,
		 int i, int j,
		 const double *x1, const double *x2,
		 int n, double * mat)
{
  // r := x[j] - x[i] for (j -> i) interaction
  double xx = x2 [0] - x1 [0];
  double yy = x2 [1] - x1 [1];
  double zz = x2 [2] - x1 [2];
  double r2 = xx * xx + yy * yy + zz * zz;
  if (r2 < sys->lubmin2)
    {
      r2 = sys->lubmin2;
    }

  double rr = sqrt (r2);
  double ex = xx / rr;
  double ey = yy / rr;
  double ez = zz / rr;

  // calc scalar functions of lubrication
  double res2b [22];
  double resinf[22];
  scalar_two_body_res (rr, res2b);
  /*
  // checking
  twobody_scalars_res (0, // F
		       rr, 1.0, 1.0, NULL, 100,
		       1, 1, res2b);
  // for a=1 in dimensional form is equivalent to SD form in mono
  */
  scalar_minv_f (rr, resinf);

  double xa11, ya11;
  double xa12, ya12;
  xa11 = res2b [ 0] - resinf [ 0];
  xa12 = res2b [ 1] - resinf [ 1];
  ya11 = res2b [ 2] - resinf [ 2];
  ya12 = res2b [ 3] - resinf [ 3];

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
	       xa11, ya11,
	       n, mat);
  matrix_f_ij (j, i,
	       -ex, -ey, -ez,
	       xa12, ya12,
	       n, mat);
}

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
calc_lub_f_2b_poly (struct stokes *sys,
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
  scalars_lub_poly_full (0, // F version
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
matrix_lub_f_2b_poly (struct stokes *sys,
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
  scalars_lub_poly_full (0, // F version
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

/* pre-process for imposed flow shifting, that is, converting U
 * from the labo frame
 *    u(x) is given by the imposed flow field as |x|-> infty
 * to the fluid-rest frame
 *    u(x) = 0 as |x|-> infty
 * INPUT
 *  sys     : struct stokes
 *  np      : number of particles to shift (and defined in u[])
 *  u[np*3] : velocity in the fluid-rest frame
 *            (data is preserved)
 * OUTPUT
 *  u0[np*3] : velocity in the labo frame
 */
void
shift_labo_to_rest_U (struct stokes * sys,
		      int np, const double *u,
		      double *u0)
{
  int i, ix, iy, iz;
  double uOx, uOy, uOz;
  double uEx, uEy, uEz;
  double Ezz;
  Ezz = - sys->Ei[0] - sys->Ei[4];

  for (i = 0; i < np; i ++)
    {
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;

      // O\times x
      uOx = sys->Oi[1] * sys->pos[iz] - sys->Oi[2] * sys->pos[iy];
      uOy = sys->Oi[2] * sys->pos[ix] - sys->Oi[0] * sys->pos[iz];
      uOz = sys->Oi[0] * sys->pos[iy] - sys->Oi[1] * sys->pos[ix];

      // E.x
      uEx = sys->Ei[0] * sys->pos[ix] // xx . x
	+   sys->Ei[1] * sys->pos[iy] // xy . y
	+   sys->Ei[2] * sys->pos[iz];// xz . z
      uEy = sys->Ei[1] * sys->pos[ix] // yx . x
	+   sys->Ei[4] * sys->pos[iy] // yy . y
	+   sys->Ei[3] * sys->pos[iz];// yz . z
      uEz = sys->Ei[2] * sys->pos[ix] // zx . x
	+   sys->Ei[3] * sys->pos[iy] // zy . y
	+   Ezz        * sys->pos[iz];// zz . z

      u0 [ix] = u[ix] - (sys->Ui[0] + uOx + uEx);
      u0 [iy] = u[iy] - (sys->Ui[1] + uOy + uEy);
      u0 [iz] = u[iz] - (sys->Ui[2] + uOz + uEz);
    }
}

/* post-process for imposed flow shifting, that is, converting U
 * from the fluid-rest frame
 *    u(x) = 0 as |x|-> infty
 * to the labo frame
 *    u(x) is given by the imposed flow field as |x|-> infty
 * INPUT
 *  sys     : struct stokes
 *  np      : number of particles to shift (and defined in u[])
 *  u[np*3] : velocity in the fluid-rest frame
 *            (data is overwritten after the process)
 * OUTPUT
 *  u[np*3] : velocity in the labo frame
 */
void
shift_rest_to_labo_U (struct stokes * sys,
		      int np, double *u)
{
  int i, ix, iy, iz;
  double uOx, uOy, uOz;
  double uEx, uEy, uEz;
  double Ezz;
  Ezz = - sys->Ei[0] - sys->Ei[4];

  for (i = 0; i < np; i ++)
    {
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;

      // O\times x
      uOx = sys->Oi[1] * sys->pos[iz] - sys->Oi[2] * sys->pos[iy];
      uOy = sys->Oi[2] * sys->pos[ix] - sys->Oi[0] * sys->pos[iz];
      uOz = sys->Oi[0] * sys->pos[iy] - sys->Oi[1] * sys->pos[ix];

      // E.x
      uEx = sys->Ei[0] * sys->pos[ix] // xx . x
	+   sys->Ei[1] * sys->pos[iy] // xy . y
	+   sys->Ei[2] * sys->pos[iz];// xz . z
      uEy = sys->Ei[1] * sys->pos[ix] // yx . x
	+   sys->Ei[4] * sys->pos[iy] // yy . y
	+   sys->Ei[3] * sys->pos[iz];// yz . z
      uEz = sys->Ei[2] * sys->pos[ix] // zx . x
	+   sys->Ei[3] * sys->pos[iy] // zy . y
	+   Ezz        * sys->pos[iz];// zz . z

      u [ix] += sys->Ui[0] + uOx + uEx;
      u [iy] += sys->Ui[1] + uOy + uEy;
      u [iz] += sys->Ui[2] + uOz + uEz;
    }
}
