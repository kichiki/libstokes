/* subroutine for the procedure of FTS version
 * Copyright (C) 2000-2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: fts.c,v 1.10 2001/02/04 11:46:33 ichiki Exp $
 */
#include <stdio.h> // fprintf ()
#include <stdlib.h> // malloc ()
#include <math.h> // sqrt ()
#include "../FINITE/two-body-res.h" /* scalar_two_body_res () */

#include "fts.h"


/* ATIMES version (for O(N^2) scheme) of
 * store matrix in FTS format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * NOTE that only 'alpha(i) <- beta(j)' interaction is stored.
 * INPUT
 *   x [11] : FTS of particle 'i' (extracted form)
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for y[] and 'j' is for x[].
 *   xa, ya, ... : scalar functions
 * OUTPUT
 *   y [11] : UOE of particle 'j' (extracted form)
 */
void
matrix_fts_atimes (double *x, double *y,
		   double ex, double ey, double ez,
		   double xa, double ya,
		   double yb,
		   double xc, double yc,
		   double xg, double yg,
		   double yh,
		   double xm, double ym, double zm)
{
  double a1, a2;
  double b1;
  double c1, c2;
  double g1, g2, g3;
  double h1;
  double mm1, mm2, mm3, mm4, mm5;

  double b1x, b1y, b1z;
  double g1x, g1y, g1z;
  double g2x, g2y, g2z;
  double g3xxx, g3xxy, g3xxz, g3xyy, g3xyz, g3yyy, g3yyz, g3yzz, g3xzz, g3zzz;
  double h1xx, h1xy, h1xz, h1yy, h1yz, h1zz;
  double m1xxxx, m1xxxy, m1xxxz, m1xxyy, m1xxyz, m1xxzz;
  double m1xyyy, m1xyyz, m1xyzz, m1xzzz;
  double m1yyyy, m1yyyz, m1yyzz, m1yzzz, m1zzzz;
  double m2xx, m2xy, m2xz, m2yy, m2yz, m2zz;
  double m4xx, m4xy, m4xz, m4yy, m4yz, m4zz;

  double exx, eyy, ezz, exy, eyz, exz;
  double exxx, exxy, exxz, exyy, exyz, exzz, eyyy, eyyz, eyzz, ezzz;
  double exxxx, exxxy, exxxz, exxyy, exxyz, exxzz;
  double exyyy, exyyz, exyzz, exzzz;
  double eyyyy, eyyyz, eyyzz, eyzzz, ezzzz;

  int i;
  double *z; // modefied for extracted matrix elements


  z = malloc (sizeof (double) * 11);
  if (z == NULL)
    {
      fprintf (stderr, "allocation error in matrix_fts_atimes ().\n");
      exit (1);
    }

  for (i = 0; i < 6; i ++)
    {
      z [i] = x [i];
    }
  z [6] = 2.0 * x [6] + x [10];
  z [7] = 2.0 * x [7];
  z [8] = 2.0 * x [8];
  z [9] = 2.0 * x [9];
  z [10] = 2.0 * x [10] + x [6];


  a1 = ya;
  a2 = xa - ya;

  b1 = yb;

  c1 = yc;
  c2 = xc - yc;

  g1 = - 1.0 / 3.0 * xg;
  g2 = yg;
  g3 = xg - 2.0 * yg;

  h1 = yh;

  mm1 =   1.5     * xm - 2.0 * ym + 0.5 * zm;
  mm2 = - 0.5     * xm            + 0.5 * zm;
  mm3 = 1.0 / 6.0 * xm            - 0.5 * zm;
  mm4 =                  0.5 * ym - 0.5 * zm;
  mm5 =                             0.5 * zm;

  exx = ex * ex;
  eyy = ey * ey;
  ezz = ez * ez;
  exy = ex * ey;
  eyz = ey * ez;
  exz = ez * ex;
  exxx = exx * ex;
  exxy = exx * ey;
  exxz = exx * ez;
  exyy = exy * ey;
  exyz = exy * ez;
  exzz = exz * ez;
  eyyy = eyy * ey;
  eyyz = eyy * ez;
  eyzz = eyz * ez;
  ezzz = ezz * ez;
  exxxx = exx * exx;
  exxxy = exx * exy;
  exxxz = exx * exz;
  exxyy = exx * eyy;
  exxyz = exx * eyz;
  exxzz = exx * ezz;
  exyyy = exy * eyy;
  exyyz = exy * eyz;
  exyzz = exy * ezz;
  exzzz = exz * ezz;
  eyyyy = eyy * eyy;
  eyyyz = eyy * eyz;
  eyyzz = eyy * ezz;
  eyzzz = eyz * ezz;
  ezzzz = ezz * ezz;

  b1x = b1 * ex;
  b1y = b1 * ey;
  b1z = b1 * ez;

  g1x = g1 * ex;
  g1y = g1 * ey;
  g1z = g1 * ez;

  g2x = g2 * ex;
  g2y = g2 * ey;
  g2z = g2 * ez;

  g3xxx = g3 * exxx;
  g3xxy = g3 * exxy;
  g3xxz = g3 * exxz;
  g3xyy = g3 * exyy;
  g3xyz = g3 * exyz;
  g3xzz = g3 * exzz;
  g3yyy = g3 * eyyy;
  g3yyz = g3 * eyyz;
  g3yzz = g3 * eyzz;
  g3zzz = g3 * ezzz;
	      
  h1xx = h1 * exx;
  h1xy = h1 * exy;
  h1xz = h1 * exz;
  h1yy = h1 * eyy;
  h1yz = h1 * eyz;
  h1zz = h1 * ezz;
	      
  m1xxxx = mm1 * exxxx;
  m1xxxy = mm1 * exxxy;
  m1xxxz = mm1 * exxxz;
  m1xxyy = mm1 * exxyy;
  m1xxyz = mm1 * exxyz;
  m1xxzz = mm1 * exxzz;
  m1xyyy = mm1 * exyyy;
  m1xyyz = mm1 * exyyz;
  m1xyzz = mm1 * exyzz;
  m1xzzz = mm1 * exzzz;
  m1yyyy = mm1 * eyyyy;
  m1yyyz = mm1 * eyyyz;
  m1yyzz = mm1 * eyyzz;
  m1yzzz = mm1 * eyzzz;
  m1zzzz = mm1 * ezzzz;

  m2xx = mm2 * exx;
  m2xy = mm2 * exy;
  m2xz = mm2 * exz;
  m2yy = mm2 * eyy;
  m2yz = mm2 * eyz;
  m2zz = mm2 * ezz;

  m4xx = mm4 * exx;
  m4xy = mm4 * exy;
  m4xz = mm4 * exz;
  m4yy = mm4 * eyy;
  m4yz = mm4 * eyz;
  m4zz = mm4 * ezz;

  /* A part */
  y [ 0] += z [ 0] * (a1 + a2 * exx);
  y [ 1] += z [ 1] * (a1 + a2 * eyy);
  y [ 2] += z [ 2] * (a1 + a2 * ezz);
  
  y [ 0] += z [ 1] * (a2 * exy);
  y [ 1] += z [ 0] * (a2 * exy);
  y [ 1] += z [ 2] * (a2 * eyz);
  y [ 2] += z [ 1] * (a2 * eyz);
  y [ 2] += z [ 0] * (a2 * exz);
  y [ 0] += z [ 2] * (a2 * exz);

  /* B part */
  /* eps_ijk e_k */
  y [ 3] += z [ 1] * (+ b1z ); /* x,y */
  y [ 3] += z [ 2] * (- b1y ); /* x,z */
  y [ 4] += z [ 0] * (- b1z ); /* y,x */
  y [ 4] += z [ 2] * (+ b1x ); /* y,z */
  y [ 5] += z [ 0] * (+ b1y ); /* z,x */
  y [ 5] += z [ 1] * (- b1x ); /* z,y */

  /* BT part */
  y [ 1] -= z [ 3] * (+ b1z ); /* y,x */
  y [ 2] -= z [ 3] * (- b1y ); /* z,x */
  y [ 0] -= z [ 4] * (- b1z ); /* x,y */
  y [ 2] -= z [ 4] * (+ b1x ); /* z,y */
  y [ 0] -= z [ 5] * (+ b1y ); /* x,z */
  y [ 1] -= z [ 5] * (- b1x ); /* y,z */

  /* C part */
  y [ 3] += z [ 3] * (c1 + c2 * exx);
  y [ 4] += z [ 4] * (c1 + c2 * eyy);
  y [ 5] += z [ 5] * (c1 + c2 * ezz);

  y [ 3] += z [ 4] * (c2 * exy);
  y [ 4] += z [ 3] * (c2 * exy);
  y [ 4] += z [ 5] * (c2 * eyz);
  y [ 5] += z [ 4] * (c2 * eyz);
  y [ 5] += z [ 3] * (c2 * exz);
  y [ 3] += z [ 5] * (c2 * exz);
  
  /* G part */
  /* g1: d_ij e_k part */
  y [ 6] += z [ 0] * (g1x); /* xx,x */
  y [ 6] += z [ 1] * (g1y); /* xx,y */
  y [ 6] += z [ 2] * (g1z); /* xx,z */
  y [10] += z [ 0] * (g1x); /* yy,x */
  y [10] += z [ 1] * (g1y); /* yy,y */
  y [10] += z [ 2] * (g1z); /* yy,z */

  /* g2: (e_i d_jk + e_j d_ik) part */
  y [ 6] += z [ 0] * (2.0 * g2x); /* xx,x */
  y [ 7] += z [ 0] * (g2y);       /* xy,x */
  y [ 7] += z [ 1] * (g2x);       /* xy,y */
  y [ 8] += z [ 0] * (g2z);       /* xz,x */
  y [ 8] += z [ 2] * (g2x);       /* xz,z */
  y [ 9] += z [ 1] * (g2z);       /* yz,y */
  y [ 9] += z [ 2] * (g2y);       /* yz,z */
  y [10] += z [ 1] * (2.0 * g2y); /* yy,y */

  /* g3: e_ijk part */
  y [ 6] += z [ 0] * (g3xxx); /* xx,x */
  y [ 6] += z [ 1] * (g3xxy); /* xx,y */
  y [ 6] += z [ 2] * (g3xxz); /* xx,z */
  y [ 7] += z [ 0] * (g3xxy); /* xy,x */
  y [ 7] += z [ 1] * (g3xyy); /* xy,y */
  y [ 7] += z [ 2] * (g3xyz); /* xy,z */
  y [ 8] += z [ 0] * (g3xxz); /* xz,x */
  y [ 8] += z [ 1] * (g3xyz); /* xz,y */
  y [ 8] += z [ 2] * (g3xzz); /* xz,z */
  y [ 9] += z [ 0] * (g3xyz); /* yz,x */
  y [ 9] += z [ 1] * (g3yyz); /* yz,y */
  y [ 9] += z [ 2] * (g3yzz); /* yz,z */
  y [10] += z [ 0] * (g3xyy); /* yy,x */
  y [10] += z [ 1] * (g3yyy); /* yy,y */
  y [10] += z [ 2] * (g3yyz); /* yy,z */

  /* G21  =  (GT12)t  =  (-GT21)t */
  /* g1: d_ij e_k part */
  y [ 0] -= z [ 6] * g1x; /* x,xx */
  y [ 1] -= z [ 6] * g1y; /* y,xx */
  y [ 2] -= z [ 6] * g1z; /* z,xx */
  y [ 0] -= z [10] * g1x; /* x,yy */
  y [ 1] -= z [10] * g1y; /* y,yy */
  y [ 2] -= z [10] * g1z; /* z,yy */

  /* g2: (e_i d_jk + e_j d_ik) part */
  y [ 0] -= z [ 6] * (2.0 * g2x); /* x,xx */
  y [ 0] -= z [ 7] * (g2y);       /* x,xy */
  y [ 1] -= z [ 7] * (g2x);       /* y,xy */
  y [ 0] -= z [ 8] * (g2z);       /* x,xz */
  y [ 2] -= z [ 8] * (g2x);       /* z,xz */
  y [ 1] -= z [ 9] * (g2z);       /* y,yz */
  y [ 2] -= z [ 9] * (g2y);       /* z,yz */
  y [ 1] -= z [10] * (2.0 * g2y); /* y,yy */

  /* g3: e_ijk part */
  y [ 0] -= z [ 6] * g3xxx; /* x,xx */
  y [ 1] -= z [ 6] * g3xxy; /* y,xx */
  y [ 2] -= z [ 6] * g3xxz; /* z,xx */
  y [ 0] -= z [ 7] * g3xxy; /* x,xy */
  y [ 1] -= z [ 7] * g3xyy; /* y,xy */
  y [ 2] -= z [ 7] * g3xyz; /* z,xy */
  y [ 0] -= z [ 8] * g3xxz; /* x,xz */
  y [ 1] -= z [ 8] * g3xyz; /* y,xz */
  y [ 2] -= z [ 8] * g3xzz; /* z,xz */
  y [ 0] -= z [ 9] * g3xyz; /* x,yz */
  y [ 1] -= z [ 9] * g3yyz; /* y,yz */
  y [ 2] -= z [ 9] * g3yzz; /* z,yz */
  y [ 0] -= z [10] * g3xyy; /* x,yy */
  y [ 1] -= z [10] * g3yyy; /* y,yy */
  y [ 2] -= z [10] * g3yyz; /* z,yy */

  /* H part */
  /* e_i eps_jkl e_l + e_j eps_ikl e_l */
  y [ 6] += z [ 4] * (+ 2.0 * h1xz);  /* xx,y */
  y [ 6] += z [ 5] * (- 2.0 * h1xy);  /* xx,z */
  y [ 7] += z [ 3] * (- h1xz);        /* xy,x */
  y [ 7] += z [ 4] * (+ h1yz);        /* xy,y */
  y [ 7] += z [ 5] * (+ h1xx - h1yy); /* xy,z */
  y [ 8] += z [ 3] * (+ h1xy);        /* xz,x */
  y [ 8] += z [ 4] * (- h1xx + h1zz); /* xz,y */
  y [ 8] += z [ 5] * (- h1yz);        /* xz,z */
  y [ 9] += z [ 3] * (+ h1yy - h1zz); /* yz,x */
  y [ 9] += z [ 4] * (- h1xy);        /* yz,y */
  y [ 9] += z [ 5] * (+ h1xz);        /* yz,z */
  y [10] += z [ 3] * (- 2.0 * h1yz);  /* yy,x */
  y [10] += z [ 5] * (+ 2.0 * h1xy);  /* yy,z */
  
  /* HT part */
  /* H12  =  (HT21)t  =  (HT12)t */
  y [ 4] += z [ 6] * (+ 2.0 * h1xz);  /* y,xx */
  y [ 5] += z [ 6] * (- 2.0 * h1xy);  /* z,xx */
  y [ 3] += z [ 7] * (- h1xz);        /* x,xy */
  y [ 4] += z [ 7] * (+ h1yz);        /* y,xy */
  y [ 5] += z [ 7] * (+ h1xx - h1yy); /* z,xy */
  y [ 3] += z [ 8] * (+ h1xy);        /* x,xz */
  y [ 4] += z [ 8] * (- h1xx + h1zz); /* y,xz */
  y [ 5] += z [ 8] * (- h1yz);        /* z,xz */
  y [ 3] += z [ 9] * (+ h1yy - h1zz); /* x,yz */
  y [ 4] += z [ 9] * (- h1xy);        /* y,yz */
  y [ 5] += z [ 9] * (+ h1xz);        /* z,yz */
  y [ 3] += z [10] * (- 2.0 * h1yz);  /* x,yy */
  y [ 5] += z [10] * (+ 2.0 * h1xy);  /* z,yy */
  
  /* M part */
  /* e_ijkl part */
  y [ 6] += z [ 6] * (m1xxxx); /* xx,xx */
  y [ 6] += z [ 7] * (m1xxxy); /* xx,xy */
  y [ 6] += z [ 8] * (m1xxxz); /* xx,xz */
  y [ 6] += z [ 9] * (m1xxyz); /* xx,yz */
  y [ 6] += z [10] * (m1xxyy); /* xx,yy */
  y [ 7] += z [ 6] * (m1xxxy); /* xy,xx */
  y [ 7] += z [ 7] * (m1xxyy); /* xy,xy */
  y [ 7] += z [ 8] * (m1xxyz); /* xy,xz */
  y [ 7] += z [ 9] * (m1xyyz); /* xy,yz */
  y [ 7] += z [10] * (m1xyyy); /* xy,yy */
  y [ 8] += z [ 6] * (m1xxxz); /* xz,xx */
  y [ 8] += z [ 7] * (m1xxyz); /* xz,xy */
  y [ 8] += z [ 8] * (m1xxzz); /* xz,xz */
  y [ 8] += z [ 9] * (m1xyzz); /* xz,yz */
  y [ 8] += z [10] * (m1xyyz); /* xz,yy */
  y [ 9] += z [ 6] * (m1xxyz); /* yz,xx */
  y [ 9] += z [ 7] * (m1xyyz); /* yz,xy */
  y [ 9] += z [ 8] * (m1xyzz); /* yz,xz */
  y [ 9] += z [ 9] * (m1yyzz); /* yz,yz */
  y [ 9] += z [10] * (m1yyyz); /* yz,yy */
  y [10] += z [ 6] * (m1xxyy); /* yy,xx */
  y [10] += z [ 7] * (m1xyyy); /* yy,xy */
  y [10] += z [ 8] * (m1xyyz); /* yy,xz */
  y [10] += z [ 9] * (m1yyyz); /* yy,yz */
  y [10] += z [10] * (m1yyyy); /* yy,yy */

  /* (d_ij e_kl + e_ij d_kl)part */
  y [ 6] += z [ 6] * (2.0 * m2xx);  /* xx,xx */
  y [ 6] += z [ 7] * (m2xy);        /* xx,xy */
  y [ 6] += z [ 8] * (m2xz);        /* xx,xz */
  y [ 6] += z [ 9] * (m2yz);        /* xx,yz */
  y [ 6] += z [10] * (m2yy + m2xx); /* xx,yy */
  y [ 7] += z [ 6] * (m2xy);        /* xy,xx */
  y [ 7] += z [10] * (m2xy);        /* xy,yy */
  y [ 8] += z [ 6] * (m2xz);        /* xz,xx */
  y [ 8] += z [10] * (m2xz);        /* xz,yy */
  y [ 9] += z [ 6] * (m2yz);        /* yz,xx */
  y [ 9] += z [10] * (m2yz);        /* yz,yy */
  y [10] += z [ 6] * (m2xx + m2yy); /* yy,xx */
  y [10] += z [ 7] * (m2xy);        /* yy,xy */
  y [10] += z [ 8] * (m2xz);        /* yy,xz */
  y [10] += z [ 9] * (m2yz);        /* yy,yz */
  y [10] += z [10] * (2.0 * m2yy);  /* yy,yy */

  /* d_ij d_kl part */
  y [ 6] += z [ 6] * (mm3); /* xx,xx */
  y [ 6] += z [10] * (mm3); /* xx,yy */
  y [10] += z [ 6] * (mm3); /* yy,xx */
  y [10] += z [10] * (mm3); /* yy,yy */

  /* e_i d_jk e_l + e_i d_jl e_k + */
  /* e_j d_ik e_l + e_j d_il e_k part */
  y [ 6] += z [ 6] * (4.0 * m4xx);  /* xx,xx */
  y [ 6] += z [ 7] * (2.0 * m4xy);  /* xx,xy */
  y [ 6] += z [ 8] * (2.0 * m4xz);  /* xx,xz */
  y [ 7] += z [ 6] * (2.0 * m4xy);  /* xy,xx */
  y [ 7] += z [ 7] * (m4xx + m4yy); /* xy,xy */
  y [ 7] += z [ 8] * (m4yz);        /* xy,xz */
  y [ 7] += z [ 9] * (m4xz);        /* xy,yz */
  y [ 7] += z [10] * (2.0 * m4xy);  /* xy,yy */
  y [ 8] += z [ 6] * (2.0 * m4xz);  /* xz,xx */
  y [ 8] += z [ 7] * (m4yz);        /* xz,xy */
  y [ 8] += z [ 8] * (m4xx + m4zz); /* xz,xz */
  y [ 8] += z [ 9] * (m4xy);        /* xz,yz */
  y [ 9] += z [ 7] * (m4xz);        /* yz,xy */
  y [ 9] += z [ 8] * (m4xy);        /* yz,xz */
  y [ 9] += z [ 9] * (m4yy + m4zz); /* yz,yz */
  y [ 9] += z [10] * (2.0 * m4yz);  /* yz,yy */
  y [10] += z [ 7] * (2.0 * m4xy);  /* yy,xy */
  y [10] += z [ 9] * (2.0 * m4yz);  /* yy,yz */
  y [10] += z [10] * (4.0 * m4yy);  /* yy,yy */

  /* (d_ik d_jl + d_jk d_il) part */
  y [ 6] += z [ 6] * (2.0 * mm5); /* xx,xx */
  y [ 7] += z [ 7] * (mm5);       /* xy,xy */
  y [ 8] += z [ 8] * (mm5);       /* xz,xz */
  y [ 9] += z [ 9] * (mm5);       /* yz,yz */
  y [10] += z [10] * (2.0 * mm5); /* yy,yy */

  free (z);
}

/* convert fts[] to f[], t[], s[] (this is applicable for UOE)
 * INPUT
 *  n : # particles
 *  fts [n * 11] :
 * OUTPUT
 *  f[n * 3] :
 *  t[n * 3] :
 *  s[n * 5] :
 */
void
set_FTS_by_fts (int n,
		double *f, double *t, double *s,
		double *fts)
{
  int i;
  int j;
  int i3, i5, i11;


  for (i = 0; i < n; i ++)
    {
      i3 = i * 3;
      i5 = i * 5;
      i11 = i * 11;
      for (j = 0; j < 3; j ++)
	{
	  f [i3 + j] = fts [i11 + j];
	  t [i3 + j] = fts [i11 + 3 + j];
	}
      for (j = 0; j < 5; j ++)
	{
	  s [i5 + j] = fts [i11 + 6 + j];
	}
    }
}

/* convert fts[] to f[], t[], s[] (this is applicable for UOE)
 * INPUT
 *  n : # particles
 *  f[n * 3] :
 *  t[n * 3] :
 *  s[n * 5] :
 * OUTPUT
 *  fts [n * 11] :
 */
void
set_fts_by_FTS (int n,
		double *fts,
		double *f, double *t, double *s)
{
  int i;
  int j;
  int i3, i5, i11;


  for (i = 0; i < n; i ++)
    {
      i3 = i * 3;
      i5 = i * 5;
      i11 = i * 11;
      for (j = 0; j < 3; j ++)
	{
	  fts [i11 + j] = f [i3 + j];
	  fts [i11 + 3 + j] = t [i3 + j];
	}
      for (j = 0; j < 5; j ++)
	{
	  fts [i11 + 6 + j] = s [i5 + j];
	}
    }
}

/* calc scalar functions of (M^inf)^-1 in FTS
 * INPUT
 *   s : distance of particles
 * OUTPUT
 *   lub [22] : scalar functions
 */
void
scalar_minv_fts (double s,  double * scalar_fts)
{
  double xa11, xa12;
  double ya11, ya12;
  double yb11, yb12;
  double xc11, xc12;
  double yc11, yc12;
  double xg11, xg12;
  double yg11, yg12;
  double yh11, yh12;
  double xm11, xm12;
  double ym11, ym12;
  double zm11, zm12;

  double s2, s3, s4, s5, s6, s7, s8, s10;
  double gx, gy;


  s2 = s * s;
  s3 = s * s2;
  s4 = s2 * s2;
  s5 = s2 * s3;
  s6 = s2 * s4;
  s7 = s * s6;
  s8 = s4 * s4;
  s10 = s5 * s5;
  
  gx = 
    (2304.0 + s2 *
     (- 21120.0 + s2 *
      (55600.0 + s2 *
       (- 90600.0 + s2 *
	(45945.0 + s2 *
	 (- 800.0 + s2 *
	  (- 1800.0 + s2 *
	   (- 900.0 + s2 *
	    (400.0)))))))));
  gy =
    (256.0 + s2 *
     (640.0 + s2 *
      (- 2000.0 + s2 *
       (- 29624.0 + s2 *
	(16505.0 + s2 *
	 (- 14880.0 + s2 *
	  (112600.0 + s2 *
	   (- 66360 + s2 *
	    (22800.0 + s2 *
	     (3600.0 + s2 *
	      (900.0 + s2 *
	       (- 1600.0))))))))))));
  xa11 =
    20.0 * s6 *
    (- 2880.0 + s2 *
     (2208.0 + s2 *
      (- 260.0 + s2 *
       (- 75.0 + s2  * s2*
	(20.0))))) / gx;
  xa12 =
    20.0 * s3 *
    (- 576.0 + s2 *
     (2880.0 + s2 *
      (- 2000.0 + s2 *
       (375.0 + s2  * s2*
	(20.0 + s2 *
	 (- 30.0)))))) / gx;
  ya11 =
    - 80.0 * s6 *
    (320.0 + s2 *
     (- 544.0 + s2 *
      (140.0 + s2 *
       (- 1130.0 + s2 *
	(736.0 + s2 *
	 (- 280.0 + s2 *
	  (- 15.0 + s2  * s2*
	   (20.0)))))))) / gy;
  ya12 =
    40.0 * s3 *
    (64.0 + s2  * s2*
     (- 400.0 + s2 *
      (119.0 + s2 *
       (- 1440.0 + s2 *
	(1160.0 + s2 *
	 (- 405.0 + s2  * s2*
	  (20.0 + s2 *
	   (30.0)))))))) / gy;
  yb11 =
    - 80.0 * s7 *
    (384.0 + s2 *
     (- 248.0 + s2 *
      (- 190.0 + s2 *
       (150.0 + s2 *
	(- 80.0 + s2 *
	 (- 20.0 + s2 *
	  (- 15.0))))))) / gy;
  yb12 =
    - 20.0 * s4 *
    (- 128.0 + s2 *
     (- 80.0 + s2 *
      (700.0 + s2 *
       (- 2935.0 + s2 *
	(2464.0 + s2 *
	 (- 540.0 + s2 *
	  (- 30.0 + s2  * s2*
	   (80.0)))))))) / gy;
  xc11 =
    4.0 * s6
    / 3.0 / (s6 - 1.0);
  xc12 =
    - 4.0 * s3
    / 3.0 / (s6 - 1.0);
  yc11 =
    - 16.0 * s6 *
    (256.0 + s2 *
     (7840.0 + s2 *
      (1960.0 + s2 *
       (- 31525.0 + s2 *
	(15690.0 + s2 *
	 (- 4100.0 + s2 *
	  (- 600.0 + s2 *
	   (- 225.0 + s2 *
	    (400.0)))))))))
    / 3.0 / gy;
  yc12 =
    - 4.0 * s3 *
    (512.0 + s2 *
     (3680.0 + s2 *
      (200.0 + s2 *
       (- 66950.0 + s2 *
	(80505.0 + s2 *
	 (- 10600.0 + s2  * s2*
	  (450.0 + s2 *
	   (800.0))))))))
    / 3.0 / gy;
  xg11 =
    100.0 * s7 *
    (192.0 + s2 *
     (- 184.0 + s2 *
      (16.0 + s2 *
       (15.0)))) / gx;
  xg12 =
    10.0 * s4 *
    (384.0 + s2 *
     (- 2000.0 + s2 *
      (1700.0 + s2 *
       (- 375.0 + s2 *
	(160.0 + s2 *
	 (- 100.0)))))) / gx;
  yg11 =
    - 100.0 * s7 *
    (128.0 + s2 *
     (320.0 + s2 *
      (- 48.0 + s2 *
       (377.0 + s2 *
	(- 128.0 + s2 *
	 (108.0))))))
    / 3.0 / gy;
  yg12 =
    - 10.0 * s4 *
    (128.0 + s2 *
     (- 80.0 + s2 *
      (- 900.0 + s2 *
       (613.0 + s2 *
	(- 5280.0 + s2 *
	 (2580.0 + s2 *
	  (- 450.0 + s2 *
	   (- 640.0))))))))
    / 3.0 / gy;
  yh11 =
    - 20.0 * s8 *
    (736.0 + s2 *
     (- 1240.0 + s2 *
      (- 1020.0 + s2 *
       (3875.0 + s2 *
	(- 480.0)))))
    / 3.0 / gy;
  yh12 =
    - 10.0 * s5 *
    (96.0 + s2 *
     (- 120.0 + s2 *
      (- 750.0 + s2 *
       (6885.0 + s2 *
	(- 5160.0 + s2 *
	 (- 720.0 + s2 *
	  (- 450.0 + s2 *
	   (800.0))))))))
    / 3.0 / gy;
  xm11 =
    200.0 * s8 *
    (- 192.0 + s2 *
     (220.0 + s2 *
      (- 15.0 + s2 *
       (- 45.0 + s2 *
	(20.0)))))
    / 9.0 / gx;
  xm12 =
    100.0 * s5 *
    (96.0 + s2 *
     (- 584.0 + s2 *
      (810.0 + s2 *
       (- 705.0 + s2 *
	(200.0)))))
    / 9.0 / gx;
  ym11 =
    - 200.0 * s8 *
    (64.0 + s2 *
     (- 208.0 + s2 *
      (75.0 + s2 *
       (- 76.0 + s2 *
	(- 340.0 + s2 *
	 (- 180.0 + s2 *
	  (- 45.0 + s2 *
	   (80.0))))))))
    / 9.0 / gy;
  ym12 =
    - 50.0 * s5 *
    (32.0 + s2 *
     (- 8.0 + s2 *
      (- 210.0 + s2 *
       (- 543.0 + s2 *
	(- 2072.0 + s2 *
	 (360.0 + s2 *
	  (3010.0 + s2 *
	   (- 800.0))))))))
    / 9.0 / gy;
  zm11 =
    10.0 * s10
    / 9.0 / (s10 - 4.0);
  zm12 =
    - 20.0 * s5
    / 9.0 / (s10 - 4.0);

  scalar_fts [ 0] = xa11;
  scalar_fts [ 1] = xa12;
  scalar_fts [ 2] = ya11;
  scalar_fts [ 3] = ya12;
  scalar_fts [ 4] = yb11;
  scalar_fts [ 5] = yb12;
  scalar_fts [ 6] = xc11;
  scalar_fts [ 7] = xc12;
  scalar_fts [ 8] = yc11;
  scalar_fts [ 9] = yc12;
  scalar_fts [10] = xg11;
  scalar_fts [11] = xg12;
  scalar_fts [12] = yg11;
  scalar_fts [13] = yg12;
  scalar_fts [14] = yh11;
  scalar_fts [15] = yh12;
  scalar_fts [16] = xm11;
  scalar_fts [17] = xm12;
  scalar_fts [18] = ym11;
  scalar_fts [19] = ym12;
  scalar_fts [20] = zm11;
  scalar_fts [21] = zm12;
}

/* calculate lubrication fts by uoe for all particles
 * INPUT
 *   (global) pos [np * 3] : position of particles
 *   np : # particles
 *   uoe [np * 11] : velocity, angular velocity, strain
 * OUTPUT
 *   fts [np * 11] : force, torque, stresslet
 */
void
calc_lub_3fts (int np, double * uoe, double * fts)
{
  extern double * pos;

  int i, j;
  int i3, i11;
  int j3, j11;


  /* clear fts [np * 11] */
  for (i = 0; i < np * 11; ++i)
    fts [i] = 0.0;

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i11 = i * 11;
      for (j = i + 1; j < np; ++j)
	{
	  j3 = j * 3;
	  j11 = j * 11;
	  calc_lub_fts_2b (uoe + i11, uoe + j11,
			   pos + i3, pos + j3,
			   fts + i11, fts + j11);
	  
	}
    }
}

/* calculate fts by uoe for pair of particles 1 and 2
 * INPUT
 *   (global) p : order of expansion
 *   uoe1 [11] : velocity, angular velocity, strain
 *   uoe2 [11] :
 *   x1 [3] : position of particle 1
 *   x2 [3] : position of particle 2
 * OUTPUT
 *   fts1 [11] : force, torque, stresslet
 *   fts2 [11] :
 */
void
calc_lub_fts_2b (double *uoe1, double *uoe2,
		 double *x1, double *x2,
		 double *fts1, double *fts2)
{
  double *res2b, *resinf;

  double xx, yy, zz, rr;
  double ex, ey, ez;

  double xa11, ya11;
  double xa12, ya12;
  double yb11, yb12;
  double xc11, yc11;
  double xc12, yc12;
  double xg11, xg12, yg11, yg12;
  double yh11, yh12;
  double xm11, xm12, ym11, ym12, zm11, zm12;


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
  scalar_minv_fts (rr, resinf);

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
  xg11 = res2b [10] - resinf [10];
  xg12 = res2b [11] - resinf [11];
  yg11 = res2b [12] - resinf [12];
  yg12 = res2b [13] - resinf [13];
  yh11 = res2b [14] - resinf [14];
  yh12 = res2b [15] - resinf [15];
  xm11 = res2b [16] - resinf [16];
  xm12 = res2b [17] - resinf [17];
  ym11 = res2b [18] - resinf [18];
  ym12 = res2b [19] - resinf [19];
  zm11 = res2b [20] - resinf [20];
  zm12 = res2b [21] - resinf [21];

  matrix_fts_atimes (uoe1, fts1,
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

  matrix_fts_atimes (uoe2, fts2,
		     - ex, - ey, - ez,
		     xa11, ya11,
		     yb11,
		     xc11, yc11,
		     xg11, yg11,
		     yh11,
		     xm11, ym11, zm11);
  matrix_fts_atimes (uoe1, fts2,
		     - ex, - ey, - ez,
		     xa12, ya12,
		     yb12,
		     xc12, yc12,
		     xg12, yg12,
		     yh12,
		     xm12, ym12, zm12);

  free (res2b);
}
