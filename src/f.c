/* subroutine for the procedure of F version
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: f.c,v 1.1 2001/02/03 12:51:51 ichiki Exp $
 */
#include "f.h"

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
matrix_f_atimes (double *x, double *y,
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
	    double *f2)
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

/*
 * OUTPUT
 *  scalar_f [4] :
 */
void
scalar_minv_F (double s, double *scalar_f)
{
  double xa11, xa12;
  double ya11, ya12;

  double s2, s3, s6;
  double dx, dy;


  s2 = s * s;
  s3 = s * s2;
  s6 = s3 * s3;

  dx = 4.0 * s6 - (3.0 * s2 - 2.0) * (3.0 * s2 - 2.0);
  dy = 16.0 + s2 *
    (- 9.0 + s2 *
     (- 12.0 + s2 *
      (- 4.0)));

  xa11 =
    - 4.0 * s6
    / dx;
  xa12 =
    + 2.0 * s3 * (3.0 * s2 - 2.0)
    / dx;
  ya11 =
    - 16.0 * s6
    / dy;
  ya12 =
    + 4.0 * s3 * (3.0 * s2 + 2.0)
    / dy;

  scalar_f [0] = xa11;
  scalar_f [1] = xa12;
  scalar_f [2] = ya11;
  scalar_f [3] = ya12;
}
