/* subroutine for the procedure of FT version
 * Copyright (C) 2000-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ft.c,v 2.2 2006/09/27 00:09:58 ichiki Exp $
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
#include "f.h" // matrix_ij_A ()

#include "ft.h"

/* store matrix in FT format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   i, j : particle index
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xa, ya, ... : scalar functions
 *   n6 : dimension of the matrix mat []
 * OUTPUT
 *   mat [n6 * n6] :
 */
void
matrix_ft_ij (int i, int j,
	      double ex, double ey, double ez,
	      double xa, double ya,
	      double yb,
	      double xc, double yc,
	      int n6, double *mat)
{
  int i1, i4;
  int j1, j4;

  i1  = i * 6;
  i4  = i1 + 3;

  j1  = j * 6;
  j4  = j1 + 3;

  matrix_ij_A (n6, & mat [i1 * n6 + j1],
	       ex, ey, ez,
	       xa, ya);
  matrix_ij_B (n6, & mat [i4 * n6 + j1],
	       ex, ey, ez,
	       yb);
  matrix_ij_Bt (n6, & mat [i1 * n6 + j4],
		ex, ey, ez,
		yb);
  matrix_ij_C (n6, & mat [i4 * n6 + j4],
	       ex, ey, ez,
	       xc, yc);
}

/* store matrix in B part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   yb : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2) * n + (0,1,2)] : added (not cleared!)
 */
void
matrix_ij_B (int n, double *mat,
	     double ex, double ey, double ez,
	     double yb)
{
  double b1;
  double b1x, b1y, b1z;

  b1 = yb;

  b1x = b1 * ex;
  b1y = b1 * ey;
  b1z = b1 * ez;

  /* B part */
  /* eps_ijk e_k */
  mat [0 * n + 1] += + b1z ; /* x,y */
  mat [0 * n + 2] += - b1y ; /* x,z */
  mat [1 * n + 0] += - b1z ; /* y,x */
  mat [1 * n + 2] += + b1x ; /* y,z */
  mat [2 * n + 0] += + b1y ; /* z,x */
  mat [2 * n + 1] += - b1x ; /* z,y */
}

/* store matrix in B-tilde part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   yb : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2) * n + (0,1,2)] : added (not cleared!)
 */
void
matrix_ij_Bt (int n, double *mat,
	      double ex, double ey, double ez,
	      double yb)
{
  double b1;
  double b1x, b1y, b1z;

  b1 = yb;

  b1x = b1 * ex;
  b1y = b1 * ey;
  b1z = b1 * ez;

  /* BT part */
  mat [1 * n + 0] -= + b1z ; /* y,x */
  mat [2 * n + 0] -= - b1y ; /* z,x */
  mat [0 * n + 1] -= - b1z ; /* x,y */
  mat [2 * n + 1] -= + b1x ; /* z,y */
  mat [0 * n + 2] -= + b1y ; /* x,z */
  mat [1 * n + 2] -= - b1x ; /* y,z */
}

/* store matrix in C part with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xc, yc : scalar functions
 *   n : # dimension of matrix mat[] (should be more tha 'np * 3')
 * OUTPUT
 *   mat [(0,1,2) * n + (0,1,2)] : added (not cleared!)
 */
void
matrix_ij_C (int n, double *mat,
	     double ex, double ey, double ez,
	     double xc, double yc)
{
  double c1, c2;

  double exx, eyy, ezz, exy, eyz, exz;

  c1 = yc;
  c2 = xc - yc;

  exx = ex * ex;
  eyy = ey * ey;
  ezz = ez * ez;
  exy = ex * ey;
  eyz = ey * ez;
  exz = ez * ex;

  /* C part */
  mat [0 * n + 0] += c1 + c2 * exx;
  mat [1 * n + 1] += c1 + c2 * eyy;
  mat [2 * n + 2] += c1 + c2 * ezz;

  mat [0 * n + 1] += c2 * exy;
  mat [1 * n + 0] += c2 * exy;
  mat [1 * n + 2] += c2 * eyz;
  mat [2 * n + 1] += c2 * eyz;
  mat [2 * n + 0] += c2 * exz;
  mat [0 * n + 2] += c2 * exz;
}

/* ATIMES version (for O(N^2) scheme) of
 * store matrix in FT format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * NOTE that only 'alpha(i) <- beta(j)' interaction is stored.
 * INPUT
 *   x [6] : FTS of particle 'i'
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for y[] and 'j' is for x[].
 *   xa, ya, ... : scalar functions
 * OUTPUT
 *   y [6] : UOE of particle 'j'
 */
void
matrix_ft_atimes (const double *x,
		  double *y,
		  double ex, double ey, double ez,
		  double xa, double ya,
		  double yb,
		  double xc, double yc)
{
  double a1, a2;
  double b1;
  double c1, c2;

  double b1x, b1y, b1z;

  double exx, eyy, ezz, exy, eyz, exz;


  a1 = ya;
  a2 = xa - ya;

  b1 = yb;

  c1 = yc;
  c2 = xc - yc;

  exx = ex * ex;
  eyy = ey * ey;
  ezz = ez * ez;
  exy = ex * ey;
  eyz = ey * ez;
  exz = ez * ex;

  b1x = b1 * ex;
  b1y = b1 * ey;
  b1z = b1 * ez;

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

  /* B part */
  /* eps_ijk e_k */
  y [ 3] += x [ 1] * (+ b1z ); /* x,y */
  y [ 3] += x [ 2] * (- b1y ); /* x,z */
  y [ 4] += x [ 0] * (- b1z ); /* y,x */
  y [ 4] += x [ 2] * (+ b1x ); /* y,z */
  y [ 5] += x [ 0] * (+ b1y ); /* z,x */
  y [ 5] += x [ 1] * (- b1x ); /* z,y */

  /* BT part */
  y [ 1] -= x [ 3] * (+ b1z ); /* y,x */
  y [ 2] -= x [ 3] * (- b1y ); /* z,x */
  y [ 0] -= x [ 4] * (- b1z ); /* x,y */
  y [ 2] -= x [ 4] * (+ b1x ); /* z,y */
  y [ 0] -= x [ 5] * (+ b1y ); /* x,z */
  y [ 1] -= x [ 5] * (- b1x ); /* y,z */

  /* C part */
  y [ 3] += x [ 3] * (c1 + c2 * exx);
  y [ 4] += x [ 4] * (c1 + c2 * eyy);
  y [ 5] += x [ 5] * (c1 + c2 * ezz);

  y [ 3] += x [ 4] * (c2 * exy);
  y [ 4] += x [ 3] * (c2 * exy);
  y [ 4] += x [ 5] * (c2 * eyz);
  y [ 5] += x [ 4] * (c2 * eyz);
  y [ 5] += x [ 3] * (c2 * exz);
  y [ 3] += x [ 5] * (c2 * exz);
}

/* convert ft[] to f[], t[] (this is applicable for UO)
 * INPUT
 *  n : # particles
 *  ft [n * 6] :
 * OUTPUT
 *  f[n * 3] :
 *  t[n * 3] :
 */
void
set_FT_by_ft (int n,
	      double *f, double *t,
	      const double *ft)
{
  int i;
  int j;
  int i3, i6;


  for (i = 0; i < n; i ++)
    {
      i3 = i * 3;
      i6 = i * 6;
      for (j = 0; j < 3; j ++)
	{
	  f [i3 + j] = ft [i6 + j];
	  t [i3 + j] = ft [i6 + 3 + j];
	}
    }
}

/* convert ft[] to f[], t[] (this is applicable for UO)
 * INPUT
 *  n : # particles
 *  f[n * 3] :
 *  t[n * 3] :
 * OUTPUT
 *  ft [n * 6] :
 */
void
set_ft_by_FT (int n,
	      double *ft,
	      const double *f, const double *t)
{
  int i;
  int j;
  int i3, i6;


  for (i = 0; i < n; i ++)
    {
      i3 = i * 3;
      i6 = i * 6;
      for (j = 0; j < 3; j ++)
	{
	  ft [i6 + j] = f [i3 + j];
	  ft [i6 + 3 + j] = t [i3 + j];
	}
    }
}

/* calc scalar functions of (M^inf)^-1 in FT
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *  scalar_ft [10] :
 */
void
scalar_minv_ft (double s, double * scalar_ft)
{
  double xa11, xa12;
  double ya11, ya12;
  double yb11, yb12;
  double xc11, xc12;
  double yc11, yc12;

  double s2, s3, s4, s6;
  double dx, dy;


  s2 = s * s;
  s3 = s * s2;
  s4 = s2 * s2;
  s6 = s3 * s3;

  dx = 4.0 * s6 - (3.0 * s2 - 2.0) * (3.0 * s2 - 2.0);
  dy = 4.0 + s2 *
    (- 12.0 + s2 *
     (9.0 + s2 *
      (- 32.0 + s2 *
       (- 144.0 + s2 *
	(- 36.0 + s2 *
	 (64.0))))));

  xa11 =
    4.0 * s6
    / dx;
  xa12 =
    - 2.0 * s3 * (3.0 * s2 - 2.0)
    / dx;
  ya11 =
    16.0 * s6 * ((4.0 * s4 - 3.0) * s2 - 1.0)
    / dy;
  ya12 =
    - 4.0 * s3 * (((12.0 * s2 + 8.0) * s4 + 3.0) * s2 - 2.0)
    / dy;
  yb11 =
    - 16.0 * s6 * s * (3.0 * s2 + 4.0)
    / dy;
  yb12 =
    8.0 * s4 * ((8.0 * s4 - 3.0) * s2 + 2.0)
    / dy;
  xc11 =
    4.0 * s6 / 3.0 / (s6 - 1.0);
  xc12 =
    - 4.0 * s3 / 3.0 / (s6 - 1.0);
  yc11 =
    16.0 * s6 * (((16.0 * s2 - 9.0) * s2 - 24.0) * s2 - 4.0)
    / 3.0 / dy;
  yc12 =
    8.0 * s3 * ((16.0 * s2 + 9.0) * s4 - 4.0)
    / 3.0 / dy;

  scalar_ft [0] = xa11;
  scalar_ft [1] = xa12;
  scalar_ft [2] = ya11;
  scalar_ft [3] = ya12;
  scalar_ft [4] = yb11;
  scalar_ft [5] = yb12;
  scalar_ft [6] = xc11;
  scalar_ft [7] = xc12;
  scalar_ft [8] = yc11;
  scalar_ft [9] = yc12;
}

/* calculate lubrication ft by uoe for all particles
 * INPUT
 *  sys : system parameters
 *   uo [np * 6] : velocity, angular velocity
 * OUTPUT
 *   ft [np * 6] : force, torque
 */
void
calc_lub_3ft (struct stokes * sys,
	      const double * uo,
	      double * ft)
{
  int np; 

  int i, j;
  int i3, i6;
  int j3, j6;


  np = sys->np;

  /* clear ft [np * 6] */
  for (i = 0; i < np * 6; ++i)
    {
      ft [i] = 0.0;
    }

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i6 = i * 6;
      for (j = i + 1; j < np; ++j)
	{
	  j3 = j * 3;
	  j6 = j * 6;
	  calc_lub_ft_2b (uo + i6, uo + j6,
			  sys->pos + i3, sys->pos + j3,
			  ft + i6, ft + j6);
	  
	}
    }
}

/* calculate ft by uoe for pair of particles 1 and 2
 * INPUT
 *   (global) p : order of expansion
 *   uo1 [6] : velocity, angular velocity, strain
 *   uo2 [6] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 * OUTPUT
 *   ft1 [6] : force, torque, stresslet
 *   ft2 [6] :
 */
void
calc_lub_ft_2b (const double *uo1, const double *uo2,
		const double *x1, const double *x2,
		double *ft1, double *ft2)
{
  double *res2b, *resinf;

  double xx, yy, zz, rr;
  double ex, ey, ez;

  double xa11, ya11;
  double xa12, ya12;
  double yb11, yb12;
  double xc11, yc11;
  double xc12, yc12;


  res2b = malloc (sizeof (double) * 44);
  if (res2b == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_2b ().\n");
      exit (1);
    }
  resinf = res2b + 22;

  /* r := x[j] - x[i] for (j -> i) interaction */
  xx = x2 [0] - x1 [0];
  yy = x2 [1] - x1 [1];
  zz = x2 [2] - x1 [2];
  rr = sqrt (xx * xx + yy * yy + zz * zz);

  if (rr <= 2.0)
    rr = 2.0 + 1.0e-12;

  ex = xx / rr;
  ey = yy / rr;
  ez = zz / rr;

  /* calc scalar functions of lubrication */
  scalar_two_body_res (rr, res2b);
  scalar_minv_ft (rr, resinf);

  xa11 = res2b [ 0] - resinf [ 0];
  xa12 = res2b [ 1] - resinf [ 1];
  ya11 = res2b [ 2] - resinf [ 2];
  ya12 = res2b [ 3] - resinf [ 3];
  yb11 = res2b [ 4] - resinf [ 4];
  yb12 = res2b [ 5] - resinf [ 5];
  xc11 = res2b [ 6] - resinf [ 6];
  xc12 = res2b [ 7] - resinf [ 7];
  yc11 = res2b [ 8] - resinf [ 8];
  yc12 = res2b [ 9] - resinf [ 9];

  matrix_ft_atimes (uo1, ft1,
		    ex, ey, ez,
		    xa11, ya11,
		    yb11,
		    xc11, yc11);
  matrix_ft_atimes (uo2, ft1,
		    ex, ey, ez,
		    xa12, ya12,
		    yb12,
		    xc12, yc12);

  matrix_ft_atimes (uo2, ft2,
		    - ex, - ey, - ez,
		    xa11, ya11,
		    yb11,
		    xc11, yc11);
  matrix_ft_atimes (uo1, ft2,
		    - ex, - ey, - ez,
		    xa12, ya12,
		    yb12,
		    xc12, yc12);

  free (res2b);
}
