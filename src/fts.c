/* subroutine for the procedure of FTS version
 * Copyright (C) 2000-2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: fts.c,v 1.5 2001/01/31 09:02:56 ichiki Exp $
 */
#include <stdio.h> // fprintf ()
#include <stdlib.h> // malloc ()
#include "fts.h"

/* store matrix in FTS format with scalar functions
 * r := pos [beta(j)] - pos [alpha(i)] (NOTE THE SIGN!)
 * only m [alpha(i), beta(j)] is stored.
 * INPUT
 *   i, j : particle index
 *   ex, ey, ez := (pos[j] - pos[i]) / r,
 *                 where 'i' is for UOE (y[]) and 'j' is for FTS(x[]).
 *   xa, ya, ... : scalar functions
 *   n11 : dimension of the matrix mat []
 * OUTPUT
 *   mat [n11 * n11] :
 */
void
matrix_fts_ij (int i, int j,
	       double ex, double ey, double ez,
	       double xa, double ya,
	       double yb,
	       double xc, double yc,
	       double xg, double yg,
	       double yh,
	       double xm, double ym, double zm,
	       int n11, double *mat)
{
  int i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11;
  int j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11;

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

  i1  = i * 11;
  i2  = i1 + 1;
  i3  = i1 + 2;
  i4  = i1 + 3;
  i5  = i1 + 4;
  i6  = i1 + 5;
  i7  = i1 + 6;
  i8  = i1 + 7;
  i9  = i1 + 8;
  i10 = i1 + 9;
  i11 = i1 + 10;

  j1  = j * 11;
  j2  = j1 + 1;
  j3  = j1 + 2;
  j4  = j1 + 3;
  j5  = j1 + 4;
  j6  = j1 + 5;
  j7  = j1 + 6;
  j8  = j1 + 7;
  j9  = j1 + 8;
  j10 = j1 + 9;
  j11 = j1 + 10;

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
  mat [i1 * n11 + j1] += a1 + a2 * exx;
  mat [i2 * n11 + j2] += a1 + a2 * eyy;
  mat [i3 * n11 + j3] += a1 + a2 * ezz;
  
  mat [i1 * n11 + j2] += a2 * exy;
  mat [i2 * n11 + j1] += a2 * exy;
  mat [i2 * n11 + j3] += a2 * eyz;
  mat [i3 * n11 + j2] += a2 * eyz;
  mat [i3 * n11 + j1] += a2 * exz;
  mat [i1 * n11 + j3] += a2 * exz;

  /* B part */
  /* eps_ijk e_k */
  mat [i4 * n11 + j2] += + b1z ; /* x,y */
  mat [i4 * n11 + j3] += - b1y ; /* x,z */
  mat [i5 * n11 + j1] += - b1z ; /* y,x */
  mat [i5 * n11 + j3] += + b1x ; /* y,z */
  mat [i6 * n11 + j1] += + b1y ; /* z,x */
  mat [i6 * n11 + j2] += - b1x ; /* z,y */

  /* BT part */
  mat [i2 * n11 + j4] -= + b1z ; /* y,x */
  mat [i3 * n11 + j4] -= - b1y ; /* z,x */
  mat [i1 * n11 + j5] -= - b1z ; /* x,y */
  mat [i3 * n11 + j5] -= + b1x ; /* z,y */
  mat [i1 * n11 + j6] -= + b1y ; /* x,z */
  mat [i2 * n11 + j6] -= - b1x ; /* y,z */

  /* C part */
  mat [i4 * n11 + j4] += c1 + c2 * exx;
  mat [i5 * n11 + j5] += c1 + c2 * eyy;
  mat [i6 * n11 + j6] += c1 + c2 * ezz;

  mat [i4 * n11 + j5] += c2 * exy;
  mat [i5 * n11 + j4] += c2 * exy;
  mat [i5 * n11 + j6] += c2 * eyz;
  mat [i6 * n11 + j5] += c2 * eyz;
  mat [i6 * n11 + j4] += c2 * exz;
  mat [i4 * n11 + j6] += c2 * exz;
  
  /* G part */
  /* g1: d_ij e_k part */
  mat [i7  * n11 + j1] += g1x; /* xx,x */
  mat [i7  * n11 + j2] += g1y; /* xx,y */
  mat [i7  * n11 + j3] += g1z; /* xx,z */
  mat [i11 * n11 + j1] += g1x; /* yy,x */
  mat [i11 * n11 + j2] += g1y; /* yy,y */
  mat [i11 * n11 + j3] += g1z; /* yy,z */

  /* g2: (e_i d_jk + e_j d_ik) part */
  mat [i7  * n11 + j1] += 2.0 * g2x; /* xx,x */
  mat [i8  * n11 + j1] += g2y;       /* xy,x */
  mat [i8  * n11 + j2] += g2x;       /* xy,y */
  mat [i9  * n11 + j1] += g2z;       /* xz,x */
  mat [i9  * n11 + j3] += g2x;       /* xz,z */
  mat [i10 * n11 + j2] += g2z;       /* yz,y */
  mat [i10 * n11 + j3] += g2y;       /* yz,z */
  mat [i11 * n11 + j2] += 2.0 * g2y; /* yy,y */

  /* g3: e_ijk part */
  mat [i7  * n11 + j1] += g3xxx; /* xx,x */
  mat [i7  * n11 + j2] += g3xxy; /* xx,y */
  mat [i7  * n11 + j3] += g3xxz; /* xx,z */
  mat [i8  * n11 + j1] += g3xxy; /* xy,x */
  mat [i8  * n11 + j2] += g3xyy; /* xy,y */
  mat [i8  * n11 + j3] += g3xyz; /* xy,z */
  mat [i9  * n11 + j1] += g3xxz; /* xz,x */
  mat [i9  * n11 + j2] += g3xyz; /* xz,y */
  mat [i9  * n11 + j3] += g3xzz; /* xz,z */
  mat [i10 * n11 + j1] += g3xyz; /* yz,x */
  mat [i10 * n11 + j2] += g3yyz; /* yz,y */
  mat [i10 * n11 + j3] += g3yzz; /* yz,z */
  mat [i11 * n11 + j1] += g3xyy; /* yy,x */
  mat [i11 * n11 + j2] += g3yyy; /* yy,y */
  mat [i11 * n11 + j3] += g3yyz; /* yy,z */

  /* G21  =  (GT12)t  =  (-GT21)t */
  /* g1: d_ij e_k part */
  mat [i1 * n11 + j7 ] -= g1x; /* x,xx */
  mat [i2 * n11 + j7 ] -= g1y; /* y,xx */
  mat [i3 * n11 + j7 ] -= g1z; /* z,xx */
  mat [i1 * n11 + j11] -= g1x; /* x,yy */
  mat [i2 * n11 + j11] -= g1y; /* y,yy */
  mat [i3 * n11 + j11] -= g1z; /* z,yy */

  /* g2: (e_i d_jk + e_j d_ik) part */
  mat [i1 * n11 + j7 ] -= 2.0 * g2x; /* x,xx */
  mat [i1 * n11 + j8 ] -= g2y;       /* x,xy */
  mat [i2 * n11 + j8 ] -= g2x;       /* y,xy */
  mat [i1 * n11 + j9 ] -= g2z;       /* x,xz */
  mat [i3 * n11 + j9 ] -= g2x;       /* z,xz */
  mat [i2 * n11 + j10] -= g2z;       /* y,yz */
  mat [i3 * n11 + j10] -= g2y;       /* z,yz */
  mat [i2 * n11 + j11] -= 2.0 * g2y; /* y,yy */

  /* g3: e_ijk part */
  mat [i1 * n11 + j7 ] -= g3xxx; /* x,xx */
  mat [i2 * n11 + j7 ] -= g3xxy; /* y,xx */
  mat [i3 * n11 + j7 ] -= g3xxz; /* z,xx */
  mat [i1 * n11 + j8 ] -= g3xxy; /* x,xy */
  mat [i2 * n11 + j8 ] -= g3xyy; /* y,xy */
  mat [i3 * n11 + j8 ] -= g3xyz; /* z,xy */
  mat [i1 * n11 + j9 ] -= g3xxz; /* x,xz */
  mat [i2 * n11 + j9 ] -= g3xyz; /* y,xz */
  mat [i3 * n11 + j9 ] -= g3xzz; /* z,xz */
  mat [i1 * n11 + j10] -= g3xyz; /* x,yz */
  mat [i2 * n11 + j10] -= g3yyz; /* y,yz */
  mat [i3 * n11 + j10] -= g3yzz; /* z,yz */
  mat [i1 * n11 + j11] -= g3xyy; /* x,yy */
  mat [i2 * n11 + j11] -= g3yyy; /* y,yy */
  mat [i3 * n11 + j11] -= g3yyz; /* z,yy */

  /* H part */
  /* e_i eps_jkl e_l + e_j eps_ikl e_l */
  mat [i7  * n11 + j5] += + 2.0 * h1xz;  /* xx,y */
  mat [i7  * n11 + j6] += - 2.0 * h1xy;  /* xx,z */
  mat [i8  * n11 + j4] += - h1xz;        /* xy,x */
  mat [i8  * n11 + j5] += + h1yz;        /* xy,y */
  mat [i8  * n11 + j6] += + h1xx - h1yy; /* xy,z */
  mat [i9  * n11 + j4] += + h1xy;        /* xz,x */
  mat [i9  * n11 + j5] += - h1xx + h1zz; /* xz,y */
  mat [i9  * n11 + j6] += - h1yz;        /* xz,z */
  mat [i10 * n11 + j4] += + h1yy - h1zz; /* yz,x */
  mat [i10 * n11 + j5] += - h1xy;        /* yz,y */
  mat [i10 * n11 + j6] += + h1xz;        /* yz,z */
  mat [i11 * n11 + j4] += - 2.0 * h1yz;  /* yy,x */
  mat [i11 * n11 + j6] += + 2.0 * h1xy;  /* yy,z */
  
  /* HT part */
  /* H12  =  (HT21)t  =  (HT12)t */
  mat [i5 * n11 + j7 ] += + 2.0 * h1xz;  /* y,xx */
  mat [i6 * n11 + j7 ] += - 2.0 * h1xy;  /* z,xx */
  mat [i4 * n11 + j8 ] += - h1xz;        /* x,xy */
  mat [i5 * n11 + j8 ] += + h1yz;        /* y,xy */
  mat [i6 * n11 + j8 ] += + h1xx - h1yy; /* z,xy */
  mat [i4 * n11 + j9 ] += + h1xy;        /* x,xz */
  mat [i5 * n11 + j9 ] += - h1xx + h1zz; /* y,xz */
  mat [i6 * n11 + j9 ] += - h1yz;        /* z,xz */
  mat [i4 * n11 + j10] += + h1yy - h1zz; /* x,yz */
  mat [i5 * n11 + j10] += - h1xy;        /* y,yz */
  mat [i6 * n11 + j10] += + h1xz;        /* z,yz */
  mat [i4 * n11 + j11] += - 2.0 * h1yz;  /* x,yy */
  mat [i6 * n11 + j11] += + 2.0 * h1xy;  /* z,yy */
  
  /* M part */
  /* e_ijkl part */
  mat [i7  * n11 + j7 ] += m1xxxx; /* xx,xx */
  mat [i7  * n11 + j8 ] += m1xxxy; /* xx,xy */
  mat [i7  * n11 + j9 ] += m1xxxz; /* xx,xz */
  mat [i7  * n11 + j10] += m1xxyz; /* xx,yz */
  mat [i7  * n11 + j11] += m1xxyy; /* xx,yy */
  mat [i8  * n11 + j7 ] += m1xxxy; /* xy,xx */
  mat [i8  * n11 + j8 ] += m1xxyy; /* xy,xy */
  mat [i8  * n11 + j9 ] += m1xxyz; /* xy,xz */
  mat [i8  * n11 + j10] += m1xyyz; /* xy,yz */
  mat [i8  * n11 + j11] += m1xyyy; /* xy,yy */
  mat [i9  * n11 + j7 ] += m1xxxz; /* xz,xx */
  mat [i9  * n11 + j8 ] += m1xxyz; /* xz,xy */
  mat [i9  * n11 + j9 ] += m1xxzz; /* xz,xz */
  mat [i9  * n11 + j10] += m1xyzz; /* xz,yz */
  mat [i9  * n11 + j11] += m1xyyz; /* xz,yy */
  mat [i10 * n11 + j7 ] += m1xxyz; /* yz,xx */
  mat [i10 * n11 + j8 ] += m1xyyz; /* yz,xy */
  mat [i10 * n11 + j9 ] += m1xyzz; /* yz,xz */
  mat [i10 * n11 + j10] += m1yyzz; /* yz,yz */
  mat [i10 * n11 + j11] += m1yyyz; /* yz,yy */
  mat [i11 * n11 + j7 ] += m1xxyy; /* yy,xx */
  mat [i11 * n11 + j8 ] += m1xyyy; /* yy,xy */
  mat [i11 * n11 + j9 ] += m1xyyz; /* yy,xz */
  mat [i11 * n11 + j10] += m1yyyz; /* yy,yz */
  mat [i11 * n11 + j11] += m1yyyy; /* yy,yy */

  /* (d_ij e_kl + e_ij d_kl)part */
  mat [i7  * n11 + j7 ] += 2.0 * m2xx;  /* xx,xx */
  mat [i7  * n11 + j8 ] += m2xy;        /* xx,xy */
  mat [i7  * n11 + j9 ] += m2xz;        /* xx,xz */
  mat [i7  * n11 + j10] += m2yz;        /* xx,yz */
  mat [i7  * n11 + j11] += m2yy + m2xx; /* xx,yy */
  mat [i8  * n11 + j7 ] += m2xy;        /* xy,xx */
  mat [i8  * n11 + j11] += m2xy;        /* xy,yy */
  mat [i9  * n11 + j7 ] += m2xz;        /* xz,xx */
  mat [i9  * n11 + j11] += m2xz;        /* xz,yy */
  mat [i10 * n11 + j7 ] += m2yz;        /* yz,xx */
  mat [i10 * n11 + j11] += m2yz;        /* yz,yy */
  mat [i11 * n11 + j7 ] += m2xx + m2yy; /* yy,xx */
  mat [i11 * n11 + j8 ] += m2xy;        /* yy,xy */
  mat [i11 * n11 + j9 ] += m2xz;        /* yy,xz */
  mat [i11 * n11 + j10] += m2yz;        /* yy,yz */
  mat [i11 * n11 + j11] += 2.0 * m2yy;  /* yy,yy */

  /* d_ij d_kl part */
  mat [i7  * n11 + j7 ] += mm3; /* xx,xx */
  mat [i7  * n11 + j11] += mm3; /* xx,yy */
  mat [i11 * n11 + j7 ] += mm3; /* yy,xx */
  mat [i11 * n11 + j11] += mm3; /* yy,yy */

  /* e_i d_jk e_l + e_i d_jl e_k + */
  /* e_j d_ik e_l + e_j d_il e_k part */
  mat [i7  * n11 + j7 ] += 4.0 * m4xx;  /* xx,xx */
  mat [i7  * n11 + j8 ] += 2.0 * m4xy;  /* xx,xy */
  mat [i7  * n11 + j9 ] += 2.0 * m4xz;  /* xx,xz */
  mat [i8  * n11 + j7 ] += 2.0 * m4xy;  /* xy,xx */
  mat [i8  * n11 + j8 ] += m4xx + m4yy; /* xy,xy */
  mat [i8  * n11 + j9 ] += m4yz;        /* xy,xz */
  mat [i8  * n11 + j10] += m4xz;        /* xy,yz */
  mat [i8  * n11 + j11] += 2.0 * m4xy;  /* xy,yy */
  mat [i9  * n11 + j7 ] += 2.0 * m4xz;  /* xz,xx */
  mat [i9  * n11 + j8 ] += m4yz;        /* xz,xy */
  mat [i9  * n11 + j9 ] += m4xx + m4zz; /* xz,xz */
  mat [i9  * n11 + j10] += m4xy;        /* xz,yz */
  mat [i10 * n11 + j8 ] += m4xz;        /* yz,xy */
  mat [i10 * n11 + j9 ] += m4xy;        /* yz,xz */
  mat [i10 * n11 + j10] += m4yy + m4zz; /* yz,yz */
  mat [i10 * n11 + j11] += 2.0 * m4yz;  /* yz,yy */
  mat [i11 * n11 + j8 ] += 2.0 * m4xy;  /* yy,xy */
  mat [i11 * n11 + j10] += 2.0 * m4yz;  /* yy,yz */
  mat [i11 * n11 + j11] += 4.0 * m4yy;  /* yy,yy */

  /* (d_ik d_jl + d_jk d_il) part */
  mat [i7  * n11 + j7 ] += 2.0 * mm5; /* xx,xx */
  mat [i8  * n11 + j8 ] += mm5;       /* xy,xy */
  mat [i9  * n11 + j9 ] += mm5;       /* xz,xz */
  mat [i10 * n11 + j10] += mm5;       /* yz,yz */
  mat [i11 * n11 + j11] += 2.0 * mm5; /* yy,yy */
}

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
