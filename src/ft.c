/* subroutine for the procedure of FT version
 * Copyright (C) 2000-2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ft.c,v 1.6 2001/02/03 12:50:34 ichiki Exp $
 */
#include "ft.h"

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
matrix_ft_atimes (double *x, double *y,
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
 *  fts [n * 6] :
 * OUTPUT
 *  f[n * 3] :
 *  t[n * 3] :
 */
void
set_FT_by_ft (int n,
	      double *f, double *t,
	      double *ft)
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
	      double *f, double *t)
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

/*
 * OUTPUT
 *  scalar_ft [10] :
 */
void
scalar_minv_FT (double s, double *scalar_ft)
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
    - 4.0 * s6
    / dx;
  xa12 =
    + 2.0 * s3 * (3.0 * s2 - 2.0)
    / dx;
  ya11 =
    - 16.0 * s6 * ((4.0 * s4 - 3.0) * s2 - 1.0)
    / dy;
  ya12 =
    + 4.0 * s3 * (((12.0 * s2 + 8.0) * s4 + 3.0) * s2 - 2.0)
    / dy;
  yb11 =
    + 16.0 * s6 * s * (3.0 * s2 + 4.0)
    / dy;
  yb12 =
    - 8.0 * s4 * ((8.0 * s4 - 3.0) * s2 + 2.0)
    / dy;
  xc11 =
    - 4.0 * s6 / 3.0 / (s6 - 1.0);
  xc12 =
    + 4.0 * s3 / 3.0 / (s6 - 1.0);
  yc11 =
    - 16.0 * s6 * (((16.0 * s2 - 9.0) * s2 - 24.0) * s2 - 4.0)
    / 3.0 / dy;
  yc12 =
    - 8.0 * s3 * ((16.0 * s2 + 9.0) * s4 - 4.0)
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
