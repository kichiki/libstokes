/* subroutine for the procedure of FT version
 * Copyright (C) 2000 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ft.c,v 1.3 2001/01/22 02:59:43 ichiki Exp $
 */
#include "ft.h"

/* store matrix in FT format with scalar functions
 * r := x [alpha(i)] - x [beta(j)]
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   i, j : particle index
 *   ex, ey, ez := (x[i] - x[j]) / r
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
  int i1, i2, i3, i4, i5, i6;
  int j1, j2, j3, j4, j5, j6;

  double a1, a2;
  double b1;
  double c1, c2;

  double b1x, b1y, b1z;

  double exx, eyy, ezz, exy, eyz, exz;


  /* I DON'T KNOW WHY I SHOULD DO THE FOLLOWING */
  ex = - ex;
  ey = - ey;
  ez = - ez;


  a1 = ya;
  a2 = xa - ya;

  b1 = yb;

  c1 = yc;
  c2 = xc - yc;

  i1  = i * 6;
  i2  = i1 + 1;
  i3  = i1 + 2;
  i4  = i1 + 3;
  i5  = i1 + 4;
  i6  = i1 + 5;

  j1  = j * 6;
  j2  = j1 + 1;
  j3  = j1 + 2;
  j4  = j1 + 3;
  j5  = j1 + 4;
  j6  = j1 + 5;

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
  mat [i1 * n6 + j1] += a1 + a2 * exx;
  mat [i2 * n6 + j2] += a1 + a2 * eyy;
  mat [i3 * n6 + j3] += a1 + a2 * ezz;
  
  mat [i1 * n6 + j2] += a2 * exy;
  mat [i2 * n6 + j1] += a2 * exy;
  mat [i2 * n6 + j3] += a2 * eyz;
  mat [i3 * n6 + j2] += a2 * eyz;
  mat [i3 * n6 + j1] += a2 * exz;
  mat [i1 * n6 + j3] += a2 * exz;

  /* B part */
  /* eps_ijk e_k */
  mat [i4 * n6 + j2] += + b1z ; /* x,y */
  mat [i4 * n6 + j3] += - b1y ; /* x,z */
  mat [i5 * n6 + j1] += - b1z ; /* y,x */
  mat [i5 * n6 + j3] += + b1x ; /* y,z */
  mat [i6 * n6 + j1] += + b1y ; /* z,x */
  mat [i6 * n6 + j2] += - b1x ; /* z,y */

  /* BT part */
  mat [i2 * n6 + j4] -= + b1z ; /* y,x */
  mat [i3 * n6 + j4] -= - b1y ; /* z,x */
  mat [i1 * n6 + j5] -= - b1z ; /* x,y */
  mat [i3 * n6 + j5] -= + b1x ; /* z,y */
  mat [i1 * n6 + j6] -= + b1y ; /* x,z */
  mat [i2 * n6 + j6] -= - b1x ; /* y,z */

  /* C part */
  mat [i4 * n6 + j4] += c1 + c2 * exx;
  mat [i5 * n6 + j5] += c1 + c2 * eyy;
  mat [i6 * n6 + j6] += c1 + c2 * ezz;

  mat [i4 * n6 + j5] += c2 * exy;
  mat [i5 * n6 + j4] += c2 * exy;
  mat [i5 * n6 + j6] += c2 * eyz;
  mat [i6 * n6 + j5] += c2 * eyz;
  mat [i6 * n6 + j4] += c2 * exz;
  mat [i4 * n6 + j6] += c2 * exz;
}

/* ATIMES version (for O(N^2) scheme) of
 * store matrix in FT format with scalar functions
 * r := pos [alpha(i)] - pos [beta(j)]
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   x [6] : related to 'i'
 *   ex, ey, ez := (x[i] - x[j]) / r
 *   xa, ya, ... : scalar functions
 * OUTPUT
 *   y [n6] : related to 'j'
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


  /* I DON'T KNOW WHY I SHOULD DO THE FOLLOWING */
  ex = - ex;
  ey = - ey;
  ez = - ez;


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
