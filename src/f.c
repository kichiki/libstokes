/* subroutine for the procedure of F version
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: f.c,v 1.5 2001/02/12 08:41:11 ichiki Exp $
 */
#include <stdio.h> // fprintf ()
#include <stdlib.h> // malloc ()
#include <math.h> // sqrt ()
#include "../FINITE/two-body-res.h" /* scalar_two_body_res () */

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

/* calc scalar functions of (M^inf)^-1 in FT
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

/* calculate lubrication f by u for all particles
 * INPUT
 *   (global) pos [np * 3] : position of particles
 *   np : # particles
 *   u [np * 3] : velocity
 * OUTPUT
 *   f [np * 3] : force
 */
void
calc_lub_3f (int np, double * u, double * f)
{
  extern double * pos;

  int i, j;
  int i3;
  int j3;


  /* clear f [np * 3] */
  for (i = 0; i < np * 3; ++i)
    f [i] = 0.0;

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i + 1; j < np; ++j)
	{
	  j3 = j * 3;
	  calc_lub_f_2b (u + i3, u + j3,
			 pos + i3, pos + j3,
			 f + i3, f + j3);
	  
	}
    }
}

/* calculate f by u for pair of particles 1 and 2
 * INPUT
 *   (global) p : order of expansion
 *   u1 [3] : velocity
 *   u2 [3] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 * OUTPUT
 *   f1 [3] : force
 *   f2 [3] :
 */
void
calc_lub_f_2b (double * u1, double * u2,
	       double * x1, double * x2,
	       double * f1, double * f2)
{
  double * res2b, * resinf;

  double xx, yy, zz, rr;
  double ex, ey, ez;

  double xa11, ya11;
  double xa12, ya12;


  res2b = malloc (sizeof (double) * 44);
  if (res2b == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_f_2b ().\n");
      exit (1);
    }
  resinf = res2b + 22;

  /* r := x[j] - x[i] for (j -> i) interaction */
  xx = x2 [0] - x1 [0];
  yy = x2 [1] - x1 [1];
  zz = x2 [2] - x1 [2];
  rr = sqrt (xx * xx + yy * yy + zz * zz);

  if (rr <= 2.0)
    //rr = 2.0 + 1.0e-12;
    rr = 2.0 + 1.0e-6;

  ex = xx / rr;
  ey = yy / rr;
  ez = zz / rr;

  /* calc scalar functions of lubrication */
  scalar_two_body_res (rr, res2b);
  scalar_minv_f (rr, resinf);

  xa11 = res2b [0] - resinf [0];
  xa12 = res2b [1] - resinf [1];
  ya11 = res2b [2] - resinf [2];
  ya12 = res2b [3] - resinf [3];

  /* for check
  fprintf (stdout, "%e %e %e %e %e\n",
	   rr,
	   xa11, xa12, ya11, ya12); */

  matrix_f_atimes (u1, f1,
		   ex, ey, ez,
		   xa11, ya11);
  matrix_f_atimes (u2, f1,
		   ex, ey, ez,
		   xa12, ya12);

  matrix_f_atimes (u2, f2,
		   - ex, - ey, - ez,
		   xa11, ya11);
  matrix_f_atimes (u1, f2,
		   - ex, - ey, - ez,
		   xa12, ya12);

  free (res2b);
}
