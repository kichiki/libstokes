/* subroutine for the procedure of FT version
 * Copyright (C) 2000-2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ft.c,v 1.5 2001/02/03 09:33:17 ichiki Exp $
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
	      double *fts,
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
