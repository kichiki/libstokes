/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-1996,1999-2000 Kengo Ichiki
 *               <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-3fts.c,v 2.4 2000/12/08 14:10:42 ichiki Exp $
 *
 * 3 dimensional hydrodynamics, 3D configuration
 * periodic boundary condition in 3 direction,
 * FTS version
 * non-dimension formulation
 */
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */

#ifdef ZETA
#include <time.h> /* clock() */
#else /* not ZETA */
#include "cholesky.h" /* cholesky() */
#endif /* ZETA */

#include "ewald-3fts.h"

/* calc ewald-summed mobility for FT version
 * INPUT
 *  n : # of all (mobile and fixed) particles
 *  x [n * 3] : position of particles
 * OUTPUT
 *  mat [n * 11 * n * 11] : mobility matrix
 */
void
calc_mob_ewald_3fts (int n, double *x, double *mat)
{
  extern int pcellx, pcelly, pcellz;
  extern int kmaxx, kmaxy, kmaxz;

  extern double zeta, zeta2, zaspi, za2;
  extern double pi2;
  extern double pi6vol;
  extern double lx, ly, lz; /* cell size */

#ifdef ZETA
  extern double cpu1, cpu2, cpu3;
  clock_t ctmp1, ctmp2, ctmp3;
#endif /* ZETA */

  double a1, a2;
  double axx, ayy, azz, axy, ayz, azx;
  double b1;
  double bx, by, bz;
  double c1, c2;
  double cxx, cyy, czz, cxy, cyz, czx;
  double g1, g2, g3;
  double h1;
  double mm1, mm2, mm3, mm4, mm5;
  double g1x, g1y, g1z;
  double g2x, g2y, g2z;
  double g3xxx, g3xxy, g3xxz, g3xyy, g3xyz, g3yyy, g3yyz, g3yzz, g3xzz, g3zzz;
  double h1xx, h1xy, h1xz, h1yy, h1yz, h1zz;
  double m1xxxx, m1xxxy, m1xxxz, m1xxyy, m1xxyz, m1xxzz;
  double m1xyyy, m1xyyz, m1xyzz, m1xzzz;
  double m1yyyy, m1yyyz, m1yyzz, m1yzzz, m1zzzz;
  double m2xx, m2xy, m2xz, m2yy, m2yz, m2zz;
  double m4xx, m4xy, m4xz, m4yy, m4yz, m4zz;

  double ex, ey, ez;
  double exx, eyy, ezz, exy, eyz, exz;
  double exxx, exxy, exxz, exyy, exyz, exzz, eyyy, eyyz, eyzz, ezzz;
  double exxxx, exxxy, exxxz, exxyy, exxyz, exxzz;
  double exyyy, exyyz, exyzz, exzzz;
  double eyyyy, eyyyz, eyyzz, eyzzz, ezzzz;

  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;
  double rlx, rly, rlz;

  int n11;
  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11;
  int j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double k;
  double cfactor;
  double mob;
  double pi8vol;


  /* temporary setting : we should be modefy GLOBAL pi6vol */
  pi8vol = pi6vol * 8.0 / 6.0; /* 8 pi / Vol */

  n11 = n * 11;
  /* clear matrix */
  for (i = 0; i < n11 * n11; i++)
    {
      mat [i] = 0.0;
    }
  /* diagonal part ( self part ) */
  for (i = 0; i < n; i++)
    {
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

      /* A part */
      mat [i1 * n11 + i1] = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
      mat [i2 * n11 + i2] = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
      mat [i3 * n11 + i3] = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);

      /* C part */
      mat [i4 * n11 + i4] = 0.75 - zaspi * 10.0 * za2;
      mat [i5 * n11 + i5] = 0.75 - zaspi * 10.0 * za2;
      mat [i6 * n11 + i6] = 0.75 - zaspi * 10.0 * za2;

      /* M part */
      mm5 = 9.0 / 20.0 - 0.3 * zaspi * za2 * (20.0 - 252.0 / 5.0 * za2);
      mat [i7  * n11 + i7 ] = 4.0 / 3.0 * mm5; /* xx,xx */
      mat [i7  * n11 + i11] = - 2.0 / 3.0 * mm5; /* xx,yy */
      mat [i8  * n11 + i8 ] = mm5; /* xy,xy */
      mat [i9  * n11 + i9 ] = mm5; /* xz,xz */
      mat [i10 * n11 + i10] = mm5; /* yz,yz */
      mat [i11 * n11 + i7 ] = - 2.0 / 3.0 * mm5; /* yy,xx */
      mat [i11 * n11 + i11] = 4.0 / 3.0 * mm5; /* yy,yy */
    }

#ifdef ZETA
  ctmp1 = clock ();
#endif /* ZETA */

  /* first Ewald part ( real space ) */
  for (i = 0; i < n; i++)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
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
      for (j = 0; j < n; j++)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;
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

	  for (m1 = - pcellx; m1 <= pcellx; m1++)
	    {
	      rlx = lx * (double) m1;
	      for (m2 = - pcelly; m2 <= pcelly; m2++)
		{
		  rly = ly * (double) m2;
		  for (m3 = - pcellz; m3 <= pcellz; m3++)
		    {
		      rlz = lz * (double) m3;
  
		      xx = x [jx] - x [ix] + rlx;
		      yy = x [jy] - x [iy] + rly;
		      zz = x [jz] - x [iz] + rlz;
		      rr = sqrt (xx * xx + yy * yy + zz * zz);

		      if (rr > 0.0)
			{
			  zr = zeta * rr;
			  zr2 = zr * zr;
			  s  = rr;
			  s2 = s * s;

			  ex = xx / rr;
			  ey = yy / rr;
			  ez = zz / rr;
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

			  a1 = (0.75 + 0.5 / s2) / s * erfc (zr)
			    + ((1.0 + zr2 *
				(14.0 + 4.0 * zr2 *
				 (- 5.0 + zr2))) / s2
			       - 4.5 + 3.0 * zr2) * zaspi * exp (- zr2);
			  a2 = (0.75 - 1.5 / s2) / s * erfc (zr)
			    + ((- 3.0 + zr2 *
				(- 2.0 + 4.0 * zr2 *
				 (4.0 - zr2))) / s2
			       + 1.5 - 3.0 * zr2) * zaspi * exp (- zr2);

			  b1 = - 0.75 / s2 * erfc (zr)
			    - 1.5 * (1.0 + zr2 *
				     (- 6.0 + zr2 *
				      2.0)) / s * zaspi * exp (- zr2);

			  c1 = - 3.0 / 8.0 / s2 / s * erfc (zr)
			    - 0.75 * (1.0 + zr2 *
				      (14.0 + zr2 *
				       (-20.0 + zr2 *
					4.0))) / s2 * zaspi * exp (- zr2);
			  c2 = 9.0 / 8.0 / s2 / s * erfc (zr)
			    - 0.75 * (- 3.0 + zr2 *
				      (- 2.0 + zr2 *
				       (16.0 + zr2 *
					- 4.0))) / s2 * zaspi * exp (- zr2);

			  g1 = (- 0.75 + 1.2 / s2) / s2 * erfc (zr)
			    + (- 1.5 * (1.0 + zr2 *
					(- 2.0))
			       - 0.8 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (- 16.0 + zr2 *
					  (- 4.0 )))) / s2)
			    / s * zaspi * exp (- zr2);
			  g2 = 1.2 / s2 / s2 * erfc (zr)
			    + (- 3.0 * ( zr2 *
					 (2.0 + zr2 *
					  (- 1.0)))
			       - 0.8 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (- 26.0 + zr2 *
					  (+ 26.0 + zr2 *
					   (- 4.0))))) / s2)
			    / s * zaspi * exp (- zr2);
			  g3 = (2.25 - 6.0 / s2) / s2 * erfc (zr)
			    + (- 1.5 * (- 3.0 +  zr2 *
					(- 2.0 + zr2 *
					 (+ 4.0)))
			       - 0.8 * (15.0 + zr2 *
					(10.0 + zr2 *
					 (4.0 + zr2 *
					  (- 40.0 + zr2 *
					   (8.0))))) / s2)
			    / s * zaspi * exp (- zr2);
      
			  h1 = - 9.0 / 8.0 / s2 / s * erfc (zr)
			    + 1.5 * (- 1.5 + zr2 *
				     (- 1.0 + zr2 *
				      (+ 8.0 + zr2 *
				       (- 2.0))))
			    / s2 * zaspi * exp (- zr2);
      
			  mm1 = (- 45.0 / 4.0 + 31.5 / s2) / s / s2 * erfc (zr)
			    + (- 1.5 * (15.0 +  zr2 *
					(+ 10.0 + zr2 *
					 (+ 4.0 + zr2 *
					  (- 8.0))))
			       - 0.6 * (- 105.0 + zr2 *
					(- 70.0 + zr2 *
					 (- 28.0 + zr2 *
					  (- 8.0 + zr2 *
					   (+ 96.0 + zr2 *
					    (- 16.0)))))) / s2)
			    / s2 * zaspi * exp (- zr2);
			  mm2 = (9.0 / 4.0 - 4.5 / s2) / s / s2 * erfc (zr)
			    + (- 1.5 * (- 3.0 +  zr2 *
					(- 2.0 + zr2 *
					 (+ 4.0)))
			       - 0.6 * (15.0 + zr2 *
					(10.0 + zr2 *
					 (4.0 + zr2 *
					  (- 40.0 + zr2 *
					   (8.0))))) / s2)
			    / s2 * zaspi * exp (- zr2);
			  mm3 = (- 0.75 + 0.9 / s2) / s / s2 * erfc (zr)
			    + (- 1.5 * (1.0 +  zr2 *
					(2.0))
			       - 0.6 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (+ 16.0 + zr2 *
					  (- 4.0)))) / s2)
			    / s2 * zaspi * exp (- zr2);
			  mm4 = (9.0 / 8.0 - 4.5 / s2) / s / s2 * erfc (zr)
			    + (- 1.5 * (- 3.0 +  zr2 *
					(- 2.0 + zr2 *
					 (- 8.0 + zr2 *
					  (4.0))))
			       - 0.6 * (15.0 + zr2 *
					(10.0 + zr2 *
					 (4.0 + zr2 *
					  (32.0 + zr2 *
					   (- 30.0 + zr2 *
					    (4.0)))))) / s2)
			    / s2 * zaspi * exp (- zr2);
			  mm5 = (0.9 / s2) / s / s2 * erfc (zr)
			    + (- 1.5 * (0.0 +  zr2 *
					(+ 4.0 + zr2 *
					 (- 2.0)))
			       - 0.6 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (- 26.0 + zr2 *
					  (+ 26.0 + zr2 *
					   (- 4.0))))) / s2)
			    / s2 * zaspi * exp (- zr2);
      
			  axx = (a1 + a2 * exx);
			  ayy = (a1 + a2 * eyy);
			  azz = (a1 + a2 * ezz);
			  axy = a2 * exy;
			  ayz = a2 * eyz;
			  azx = a2 * exz;

			  bx = b1 * ex;
			  by = b1 * ey;
			  bz = b1 * ez;

			  cxx = (c1 + c2 * exx);
			  cyy = (c1 + c2 * eyy);
			  czz = (c1 + c2 * ezz);
			  cxy = c2 * exy;
			  cyz = c2 * eyz;
			  czx = c2 * exz;

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

			  /*m3 = ;*/

			  m4xx = mm4 * exx;
			  m4xy = mm4 * exy;
			  m4xz = mm4 * exz;
			  m4yy = mm4 * eyy;
			  m4yz = mm4 * eyz;
			  m4zz = mm4 * ezz;

			  /*m5 = ;*/
      
			  /* A part */
			  mat [i1 * n11 + j1] += axx;
			  mat [i2 * n11 + j2] += ayy;
			  mat [i3 * n11 + j3] += azz;

			  mat [i1 * n11 + j2] += axy;
			  mat [i2 * n11 + j1] += axy;
			  mat [i2 * n11 + j3] += ayz;
			  mat [i3 * n11 + j2] += ayz;
			  mat [i3 * n11 + j1] += azx;
			  mat [i1 * n11 + j3] += azx;

			  /* B part */
			  /* note: r = x_i-x_j so that jj = 11 and ji = 12 */
			  /* B21 = -B12*/
			  mat [i4 * n11 + j2] += - bz;
			  mat [i4 * n11 + j3] +=   by;
			  mat [i5 * n11 + j1] +=   bz;
			  mat [i5 * n11 + j3] += - bx;
			  mat [i6 * n11 + j1] += - by;
			  mat [i6 * n11 + j2] +=   bx;

			  /* B_ part */
			  /* BT21 = -B12 */
			  mat [i1 * n11 + j5] += - bz;
			  mat [i1 * n11 + j6] +=   by;
			  mat [i2 * n11 + j4] +=   bz;
			  mat [i2 * n11 + j6] += - bx;
			  mat [i3 * n11 + j4] += - by;
			  mat [i3 * n11 + j5] +=   bx;

			  /* C part */
			  mat [i4 * n11 + j4] += cxx;
			  mat [i5 * n11 + j5] += cyy;
			  mat [i6 * n11 + j6] += czz;

			  mat [i4 * n11 + j5] += cxy;
			  mat [i5 * n11 + j4] += cxy;
			  mat [i5 * n11 + j6] += cyz;
			  mat [i6 * n11 + j5] += cyz;
			  mat [i6 * n11 + j4] += czx;
			  mat [i4 * n11 + j6] += czx;

			  /* G part */
			  mat [i7  * n11 + j1]
			    += g1x + 2.0 * g2x + g3xxx; /* xx,x */
			  mat [i7  * n11 + j2]
			    += g1y             + g3xxy; /* xx,y */
			  mat [i7  * n11 + j3]
			    += g1z             + g3xxz; /* xx,z */
			  mat [i8  * n11 + j1]
			    +=       g2y       + g3xxy; /* xy,x */
			  mat [i8  * n11 + j2]
			    +=       g2x       + g3xyy; /* xy,y */
			  mat [i8  * n11 + j3]
			    +=                   g3xyz; /* xy,z */
			  mat [i9  * n11 + j1]
			    +=       g2z       + g3xxz; /* xz,x */
			  mat [i9  * n11 + j2]
			    +=                   g3xyz; /* xz,y */
			  mat [i9  * n11 + j3]
			    +=       g2x       + g3xzz; /* xz,z */
			  mat [i10 * n11 + j1]
			    +=                   g3xyz; /* yz,x */
			  mat [i10 * n11 + j2]
			    +=      g2z        + g3yyz; /* yz,y */
			  mat [i10 * n11 + j3]
			    +=      g2y        + g3yzz; /* yz,z */
			  mat [i11 * n11 + j1]
			    += g1x             + g3xyy; /* yy,x */
			  mat [i11 * n11 + j2]
			    += g1y + 2.0 * g2y + g3yyy; /* yy,y */
			  mat [i11 * n11 + j3]
			    += g1z             + g3yyz; /* yy,z */

			  /* G21  =  (GT12)t  =  (-GT21)t */
			  mat [i1 * n11 + j7 ]
			    -= g1x + 2.0 * g2x + g3xxx; /* x,xx */
			  mat [i2 * n11 + j7 ]
			    -= g1y             + g3xxy; /* y,xx */
			  mat [i3 * n11 + j7 ]
			    -= g1z             + g3xxz; /* z,xx */
			  mat [i1 * n11 + j8 ]
			    -=       g2y       + g3xxy; /* x,xy */
			  mat [i2 * n11 + j8 ]
			    -=       g2x       + g3xyy; /* y,xy */
			  mat [i3 * n11 + j8 ]
			    -=                   g3xyz; /* z,xy */
			  mat [i1 * n11 + j9 ]
			    -=       g2z       + g3xxz; /* x,xz */
			  mat [i2 * n11 + j9 ]
			    -=                   g3xyz; /* y,xz */
			  mat [i3 * n11 + j9 ]
			    -=       g2x       + g3xzz; /* z,xz */
			  mat [i1 * n11 + j10]
			    -=                   g3xyz; /* x,yz */
			  mat [i2 * n11 + j10]
			    -=      g2z        + g3yyz; /* y,yz */
			  mat [i3 * n11 + j10]
			    -=      g2y        + g3yzz; /* z,yz */
			  mat [i1 * n11 + j11]
			    -= g1x             + g3xyy; /* x,yy */
			  mat [i2 * n11 + j11]
			    -= g1y + 2.0 * g2y + g3yyy; /* y,yy */
			  mat [i3 * n11 + j11]
			    -= g1z             + g3yyz; /* z,yy */

			  /* H part */
			  /*mat [i7  * n11 + j4] +=;*/           /* xx,x */
			  mat [i7  * n11 + j5] += 2.0 * h1xz;    /* xx,y */
			  mat [i7  * n11 + j6] += - 2.0 * h1xy;  /* xx,z */
			  mat [i8  * n11 + j4] += - h1xz;        /* xy,x */
			  mat [i8  * n11 + j5] += h1yz;          /* xy,y */
			  mat [i8  * n11 + j6] += - h1yy + h1xx; /* xy,z */
			  mat [i9  * n11 + j4] += - h1xy;        /* xz,x */
			  mat [i9  * n11 + j5] += h1zz - h1xx;   /* xz,y */
			  mat [i9  * n11 + j6] += - h1yz;        /* xz,z */
			  mat [i10 * n11 + j4] += - h1zz + h1yy; /* yz,x */
			  mat [i10 * n11 + j5] += - h1xy;        /* yz,y */
			  mat [i10 * n11 + j6] += h1xz;          /* yz,z */
			  mat [i11 * n11 + j4] += - 2.0 * h1yz;  /* yy,x */
			  /*mat [i11 * n11 + j5] +=;*/           /* yy,y */
			  mat [i11 * n11 + j6] += 2.0 * h1xy;    /* yy,z */

			  /* HT part */
			  /* H12  =  (HT21)t  =  (HT12)t */
			  /*mat [i4 * n11 + j7 ] +=;*/           /* x,xx */
			  mat [i5 * n11 + j7 ] += 2.0 * h1xz;    /* y,xx */
			  mat [i6 * n11 + j7 ] += - 2.0 * h1xy;  /* z,xx */
			  mat [i4 * n11 + j8 ] += - h1xz;        /* x,xy */
			  mat [i5 * n11 + j8 ] += h1yz;          /* y,xy */
			  mat [i6 * n11 + j8 ] += - h1yy + h1xx; /* z,xy */
			  mat [i4 * n11 + j9 ] += - h1xy;        /* x,xz */
			  mat [i5 * n11 + j9 ] += h1zz - h1xx;   /* y,xz */
			  mat [i6 * n11 + j9 ] += - h1yz;        /* z,xz */
			  mat [i4 * n11 + j10] += - h1zz + h1yy; /* x,yz */
			  mat [i5 * n11 + j10] += - h1xy;        /* y,yz */
			  mat [i6 * n11 + j10] += h1xz;          /* z,yz */
			  mat [i4 * n11 + j11] += - 2.0 * h1yz;  /* x,yy */
			  /*mat [i5 * n11 + j11] += ;*/           /* y,yy */
			  mat [i6 * n11 + j11] += 2.0 * h1xy;    /* z,yy */

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
			  mat [i7  * n11 + j7 ] += 2.0 * m2xx; /* xx,xx */
			  mat [i7  * n11 + j8 ] += m2xy; /* xx,xy */
			  mat [i7  * n11 + j9 ] += m2xz; /* xx,xz */
			  mat [i7  * n11 + j10] += m2yz; /* xx,yz */
			  mat [i7  * n11 + j11] += m2yy + m2xx; /* xx,yy */
			  mat [i8  * n11 + j7 ] += m2xy; /* xy,xx */
			  /*mat [i8  * n11 + j8 ] += ;*/ /* xy,xy */
			  /*mat [i8  * n11 + j9 ] += ;*/ /* xy,xz */
			  /*mat [i8  * n11 + j10] += ;*/ /* xy,yz */
			  mat [i8  * n11 + j11] += m2xy; /* xy,yy */
			  mat [i9  * n11 + j7 ] += m2xz; /* xz,xx */
			  /*mat [i9  * n11 + j8 ] += ;*/ /* xz,xy */
			  /*mat [i9  * n11 + j9 ] += ;*/ /* xz,xz */
			  /*mat [i9  * n11 + j10] += ;*/ /* xz,yz */
			  mat [i9  * n11 + j11] += m2xz; /* xz,yy */
			  mat [i10 * n11 + j7 ] += m2yz; /* yz,xx */
			  /*mat [i10 * n11 + j8 ] += ;*/ /* yz,xy */
			  /*mat [i10 * n11 + j9 ] += ;*/ /* yz,xz */
			  /*mat [i10 * n11 + j10] += ;*/ /* yz,yz */
			  mat [i10 * n11 + j11] += m2yz; /* yz,yy */
			  mat [i11 * n11 + j7 ] += m2xx + m2yy; /* yy,xx */
			  mat [i11 * n11 + j8 ] += m2xy; /* yy,xy */
			  mat [i11 * n11 + j9 ] += m2xz; /* yy,xz */
			  mat [i11 * n11 + j10] += m2yz; /* yy,yz */
			  mat [i11 * n11 + j11] += 2.0 * m2yy; /* yy,yy */

			  /* d_ij d_kl part */
			  mat [i7  * n11 + j7 ] += mm3; /* xx,xx */
			  /*mat [i7  * n11 + j8 ] +=;*/ /* xx,xy */
			  /*mat [i7  * n11 + j9 ] +=;*/ /* xx,xz */
			  /*mat [i7  * n11 + j10] +=;*/ /* xx,yz */
			  mat [i7  * n11 + j11] += mm3; /* xx,yy */
			  /*mat [i8  * n11 + j7 ] +=;*/ /* xy,xx */
			  /*mat [i8  * n11 + j8 ] +=;*/ /* xy,xy */
			  /*mat [i8  * n11 + j9 ] +=;*/ /* xy,xz */
			  /*mat [i8  * n11 + j10] +=;*/ /* xy,yz */
			  /*mat [i8  * n11 + j11] +=;*/ /* xy,yy */
			  /*mat [i9  * n11 + j7 ] +=;*/ /* xz,xx */
			  /*mat [i9  * n11 + j8 ] +=;*/ /* xz,xy */
			  /*mat [i9  * n11 + j9 ] +=;*/ /* xz,xz */
			  /*mat [i9  * n11 + j10] +=;*/ /* xz,yz */
			  /*mat [i9  * n11 + j11] +=;*/ /* xz,yy */
			  /*mat [i10 * n11 + j7 ] +=;*/ /* yz,xx */
			  /*mat [i10 * n11 + j8 ] +=;*/ /* yz,xy */
			  /*mat [i10 * n11 + j9 ] +=;*/ /* yz,xz */
			  /*mat [i10 * n11 + j10] +=;*/ /* yz,yz */
			  /*mat [i10 * n11 + j11] +=;*/ /* yz,yy */
			  mat [i11 * n11 + j7 ] += mm3; /* yy,xx */
			  /*mat [i11 * n11 + j8 ] +=;*/ /* yy,xy */
			  /*mat [i11 * n11 + j9 ] +=;*/ /* yy,xz */
			  /*mat [i11 * n11 + j10] +=;*/ /* yy,yz */
			  mat [i11 * n11 + j11] += mm3; /* yy,yy */

			  /* e_i d_jk e_l + e_i d_jl e_k + */
			  /* e_j d_ik e_l + e_j d_il e_k part */
			  mat [i7  * n11 + j7 ] += 4.0 * m4xx; /* xx,xx */
			  mat [i7  * n11 + j8 ] += 2.0 * m4xy; /* xx,xy */
			  mat [i7  * n11 + j9 ] += 2.0 * m4xz; /* xx,xz */
			  /*mat [i7  * n11 + j10] +=;*/ /* xx,yz */
			  /*mat [i7  * n11 + j11] +=;*/ /* xx,yy */
			  mat [i8  * n11 + j7 ] += 2.0 * m4xy; /* xy,xx */
			  mat [i8  * n11 + j8 ] += m4xx + m4yy; /* xy,xy */
			  mat [i8  * n11 + j9 ] += m4yz; /* xy,xz */
			  mat [i8  * n11 + j10] += m4xz; /* xy,yz */
			  mat [i8  * n11 + j11] += 2.0 * m4xy; /* xy,yy */
			  mat [i9  * n11 + j7 ] += 2.0 * m4xz; /* xz,xx */
			  mat [i9  * n11 + j8 ] += m4yz; /* xz,xy */
			  mat [i9  * n11 + j9 ] += m4xx + m4zz; /* xz,xz */
			  mat [i9  * n11 + j10] += m4xy; /* xz,yz */
			  /*mat [i9  * n11 + j11] +=;*/ /* xz,yy */
			  /*mat [i10 * n11 + j7 ] +=;*/ /* yz,xx */
			  mat [i10 * n11 + j8 ] += m4xz; /* yz,xy */
			  mat [i10 * n11 + j9 ] += m4xy; /* yz,xz */
			  mat [i10 * n11 + j10] += m4yy + m4zz; /* yz,yz */
			  mat [i10 * n11 + j11] += 2.0 * m4yz; /* yz,yy */
			  /*mat [i11 * n11 + j7 ] +=;*/ /* yy,xx */
			  mat [i11 * n11 + j8 ] += 2.0 * m4xy; /* yy,xy */
			  /*mat [i11 * n11 + j9 ] +=;*/ /* yy,xz */
			  mat [i11 * n11 + j10] += 2.0 * m4yz; /* yy,yz */
			  mat [i11 * n11 + j11] += 4.0 * m4yy; /* yy,yy */

			  /* (d_ik d_jl + d_jk d_il) part */
			  mat [i7  * n11 + j7 ] += 2.0 * mm5; /* xx,xx */
			  /* mat [i7  * n11 + j8 ] += ;*/ /* xx,xy */
			  /* mat [i7  * n11 + j9 ] += ;*/ /* xx,xz */
			  /* mat [i7  * n11 + j10] += ;*/ /* xx,yz */
			  /* mat [i7  * n11 + j11] += ;*/ /* xx,yy */
			  /* mat [i8  * n11 + j7 ] += ;*/ /* xy,xx */
			  mat [i8  * n11 + j8 ] += mm5; /* xy,xy */
			  /* mat [i8  * n11 + j9 ] += ;*/ /* xy,xz */
			  /* mat [i8  * n11 + j10] += ;*/ /* xy,yz */
			  /* mat [i8  * n11 + j11] += ;*/ /* xy,yy */
			  /* mat [i9  * n11 + j7 ] += ;*/ /* xz,xx */
			  /* mat [i9  * n11 + j8 ] += ;*/ /* xz,xy */
			  mat [i9  * n11 + j9 ] += mm5; /* xz,xz */
			  /* mat [i9  * n11 + j10] += ;*/ /* xz,yz */
			  /* mat [i9  * n11 + j11] += ;*/ /* xz,yy */
			  /* mat [i10 * n11 + j7 ] += ;*/ /* yz,xx */
			  /* mat [i10 * n11 + j8 ] += ;*/ /* yz,xy */
			  /* mat [i10 * n11 + j9 ] += ;*/ /* yz,xz */
			  mat [i10 * n11 + j10] += mm5; /* yz,yz */
			  /* mat [i10 * n11 + j11] += ;*/ /* yz,yy */
			  /* mat [i11 * n11 + j7 ] += ;*/ /* yy,xx */
			  /* mat [i11 * n11 + j8 ] += ;*/ /* yy,xy */
			  /* mat [i11 * n11 + j9 ] += ;*/ /* yy,xz */
			  /* mat [i11 * n11 + j10] += ;*/ /* yy,yz */
			  mat [i11 * n11 + j11] += 2.0 * mm5; /* yy,yy */
			}
		    }
		}
	    }
	}
    }

#ifdef ZETA
  ctmp2 = clock ();
#endif /* ZETA */

  /* Second Ewald part ( reciprocal space ) */
  for (m1 = - kmaxx; m1 <= kmaxx; m1++)
    {
      k1 = pi2 * (double) m1 / lx;
      for (m2 = - kmaxy; m2 <= kmaxy; m2++)
	{
	  k2 = pi2 * (double) m2 / ly;
	  for (m3 = - kmaxz; m3 <= kmaxz; m3++)
	    {
	      k3 = pi2 * (double) m3 / lz;
	      if (m1 != 0 || m2 != 0 || m3 != 0)
		{
		  kk = k1 * k1 + k2 * k2 + k3 * k3;
		  k = sqrt (kk);
		  k4z = kk / 4.0 / zeta2;
		  mob = pi8vol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  ex = k1 / k;
		  ey = k2 / k;
		  ez = k3 / k;
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

		  a1 = 0.75 * (1.0 - kk / 3.0) * mob;
		  b1 = - 3.0 / 8.0 * k * mob;
		  c1 = 3.0 / 16.0 * kk * mob;

		  g1 = 0.0;
		  g2 = - 3.0 / 8.0 * (1.0 - 4.0 / 15.0 * kk) * k * mob;
		  g3 = - 2.0 * g2;

		  h1 = 3.0 / 16.0 * kk * mob;
      
		  mm4 = 3.0 / 16.0 * (1.0 - kk / 5.0) * kk * mob;
		  mm1 = - 4.0 * mm4;
		  mm2 = 0.0;
		  mm3 = 0.0;
		  mm5 = 0.0;
      
		  axx = (1.0 - k1 * k1 / kk) * a1;
		  ayy = (1.0 - k2 * k2 / kk) * a1;
		  azz = (1.0 - k3 * k3 / kk) * a1;
		  axy = - k1 * k2 / kk * a1;
		  ayz = - k2 * k3 / kk * a1;
		  azx = - k3 * k1 / kk * a1;

		  bx = b1 * k1 / k;
		  by = b1 * k2 / k;
		  bz = b1 * k3 / k;

		  cxx = (1.0 - k1 * k1 / kk) * c1;
		  cyy = (1.0 - k2 * k2 / kk) * c1;
		  czz = (1.0 - k3 * k3 / kk) * c1;
		  cxy = - k1 * k2 / kk * c1;
		  cyz = - k2 * k3 / kk * c1;
		  czx = - k3 * k1 / kk * c1;

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

		  /*m3 = ;*/

		  m4xx = mm4 * exx;
		  m4xy = mm4 * exy;
		  m4xz = mm4 * exz;
		  m4yy = mm4 * eyy;
		  m4yz = mm4 * eyz;
		  m4zz = mm4 * ezz;

		  /*m5 = ;*/

		  for (i = 0; i < n; i++)
		    {
		      ix = i * 3;
		      iy = ix + 1;
		      iz = ix + 2;
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
		      for (j = 0; j < n; j++)
			{
			  jx = j * 3;
			  jy = jx + 1;
			  jz = jx + 2;
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

			  xx = x [jx] - x [ix];
			  yy = x [jy] - x [iy];
			  zz = x [jz] - x [iz];

			  cfactor = cos (+ k1 * xx
					 + k2 * yy
					 + k3 * zz);

			  /* A part */
			  mat [i1 * n11 + j1] += cfactor * axx;
			  mat [i2 * n11 + j2] += cfactor * ayy;
			  mat [i3 * n11 + j3] += cfactor * azz;

			  mat [i1 * n11 + j2] += cfactor * axy;
			  mat [i2 * n11 + j1] += cfactor * axy;
			  mat [i2 * n11 + j3] += cfactor * ayz;
			  mat [i3 * n11 + j2] += cfactor * ayz;
			  mat [i3 * n11 + j1] += cfactor * azx;
			  mat [i1 * n11 + j3] += cfactor * azx;

			  /* B part */
			  /* note: r = x_i-x_j so that jj = 11 and ji = 12 */
			  /* B21 = -B12*/
			  mat [i4 * n11 + j2] += - cfactor * bz;
			  mat [i4 * n11 + j3] +=   cfactor * by;
			  mat [i5 * n11 + j1] +=   cfactor * bz;
			  mat [i5 * n11 + j3] += - cfactor * bx;
			  mat [i6 * n11 + j1] += - cfactor * by;
			  mat [i6 * n11 + j2] +=   cfactor * bx;

			  /* B_ part */
			  /* BT21 = -B12 */
			  mat [i1 * n11 + j5] += - cfactor * bz;
			  mat [i1 * n11 + j6] +=   cfactor * by;
			  mat [i2 * n11 + j4] +=   cfactor * bz;
			  mat [i2 * n11 + j6] += - cfactor * bx;
			  mat [i3 * n11 + j4] += - cfactor * by;
			  mat [i3 * n11 + j5] +=   cfactor * bx;

			  /* C part */
			  mat [i4 * n11 + j4] += cfactor * cxx;
			  mat [i5 * n11 + j5] += cfactor * cyy;
			  mat [i6 * n11 + j6] += cfactor * czz;

			  mat [i4 * n11 + j5] += cfactor * cxy;
			  mat [i5 * n11 + j4] += cfactor * cxy;
			  mat [i5 * n11 + j6] += cfactor * cyz;
			  mat [i6 * n11 + j5] += cfactor * cyz;
			  mat [i6 * n11 + j4] += cfactor * czx;
			  mat [i4 * n11 + j6] += cfactor * czx;

			  /* G part */
			  mat [i7  * n11 + j1]
			    += cfactor * (g1x + 2.0 * g2x + g3xxx); /* xx,x */
			  mat [i7  * n11 + j2]
			    += cfactor * (g1y             + g3xxy); /* xx,y */
			  mat [i7  * n11 + j3]
			    += cfactor * (g1z             + g3xxz); /* xx,z */
			  mat [i8  * n11 + j1]
			    += cfactor * (      g2y       + g3xxy); /* xy,x */
			  mat [i8  * n11 + j2]
			    += cfactor * (      g2x       + g3xyy); /* xy,y */
			  mat [i8  * n11 + j3]
			    += cfactor * (                  g3xyz); /* xy,z */
			  mat [i9  * n11 + j1]
			    += cfactor * (      g2z       + g3xxz); /* xz,x */
			  mat [i9  * n11 + j2]
			    += cfactor * (                  g3xyz); /* xz,y */
			  mat [i9  * n11 + j3]
			    += cfactor * (      g2x       + g3xzz); /* xz,z */
			  mat [i10 * n11 + j1]
			    += cfactor * (                  g3xyz); /* yz,x */
			  mat [i10 * n11 + j2]
			    += cfactor * (     g2z        + g3yyz); /* yz,y */
			  mat [i10 * n11 + j3]
			    += cfactor * (     g2y        + g3yzz); /* yz,z */
			  mat [i11 * n11 + j1]
			    += cfactor * (g1x             + g3xyy); /* yy,x */
			  mat [i11 * n11 + j2]
			    += cfactor * (g1y + 2.0 * g2y + g3yyy); /* yy,y */
			  mat [i11 * n11 + j3]
			    += cfactor * (g1z             + g3yyz); /* yy,z */

			  /* G21  =  (GT12)t  =  (-GT21)t */
			  mat [i1 * n11 + j7 ]
			    -= cfactor * (g1x + 2.0 * g2x + g3xxx); /* x,xx */
			  mat [i2 * n11 + j7 ]
			    -= cfactor * (g1y             + g3xxy); /* y,xx */
			  mat [i3 * n11 + j7 ]
			    -= cfactor * (g1z             + g3xxz); /* z,xx */
			  mat [i1 * n11 + j8 ]
			    -= cfactor * (      g2y       + g3xxy); /* x,xy */
			  mat [i2 * n11 + j8 ]
			    -= cfactor * (      g2x       + g3xyy); /* y,xy */
			  mat [i3 * n11 + j8 ]
			    -= cfactor * (                  g3xyz); /* z,xy */
			  mat [i1 * n11 + j9 ]
			    -= cfactor * (      g2z       + g3xxz); /* x,xz */
			  mat [i2 * n11 + j9 ]
			    -= cfactor * (                  g3xyz); /* y,xz */
			  mat [i3 * n11 + j9 ]
			    -= cfactor * (      g2x       + g3xzz); /* z,xz */
			  mat [i1 * n11 + j10]
			    -= cfactor * (                  g3xyz); /* x,yz */
			  mat [i2 * n11 + j10]
			    -= cfactor * (     g2z        + g3yyz); /* y,yz */
			  mat [i3 * n11 + j10]
			    -= cfactor * (     g2y        + g3yzz); /* z,yz */
			  mat [i1 * n11 + j11]
			    -= cfactor * (g1x             + g3xyy); /* x,yy */
			  mat [i2 * n11 + j11]
			    -= cfactor * (g1y + 2.0 * g2y + g3yyy); /* y,yy */
			  mat [i3 * n11 + j11]
			    -= cfactor * (g1z             + g3yyz); /* z,yy */

			  /* H part */
			  /*mat [i7  * n11 + j4] +=;*/           /* xx,x */
			  mat [i7  * n11 + j5] += cfactor * (2.0 * h1xz);    /* xx,y */
			  mat [i7  * n11 + j6] += cfactor * (- 2.0 * h1xy);  /* xx,z */
			  mat [i8  * n11 + j4] += cfactor * (- h1xz);        /* xy,x */
			  mat [i8  * n11 + j5] += cfactor * (h1yz);          /* xy,y */
			  mat [i8  * n11 + j6] += cfactor * (- h1yy + h1xx); /* xy,z */
			  mat [i9  * n11 + j4] += cfactor * (- h1xy);        /* xz,x */
			  mat [i9  * n11 + j5] += cfactor * (h1zz - h1xx);   /* xz,y */
			  mat [i9  * n11 + j6] += cfactor * (- h1yz);        /* xz,z */
			  mat [i10 * n11 + j4] += cfactor * (- h1zz + h1yy); /* yz,x */
			  mat [i10 * n11 + j5] += cfactor * (- h1xy);        /* yz,y */
			  mat [i10 * n11 + j6] += cfactor * (h1xz);          /* yz,z */
			  mat [i11 * n11 + j4] += cfactor * (- 2.0 * h1yz);  /* yy,x */
			  /*mat [i11 * n11 + j5] +=;*/           /* yy,y */
			  mat [i11 * n11 + j6] += cfactor * (2.0 * h1xy);    /* yy,z */

			  /* HT part */
			  /* H12  =  (HT21)t  =  (HT12)t */
			  /*mat [i4 * n11 + j7 ] +=;*/           /* x,xx */
			  mat [i5 * n11 + j7 ] += cfactor * (2.0 * h1xz);    /* y,xx */
			  mat [i6 * n11 + j7 ] += cfactor * (- 2.0 * h1xy);  /* z,xx */
			  mat [i4 * n11 + j8 ] += cfactor * (- h1xz);        /* x,xy */
			  mat [i5 * n11 + j8 ] += cfactor * (h1yz);          /* y,xy */
			  mat [i6 * n11 + j8 ] += cfactor * (- h1yy + h1xx); /* z,xy */
			  mat [i4 * n11 + j9 ] += cfactor * (- h1xy);        /* x,xz */
			  mat [i5 * n11 + j9 ] += cfactor * (h1zz - h1xx);   /* y,xz */
			  mat [i6 * n11 + j9 ] += cfactor * (- h1yz);        /* z,xz */
			  mat [i4 * n11 + j10] += cfactor * (- h1zz + h1yy); /* x,yz */
			  mat [i5 * n11 + j10] += cfactor * (- h1xy);        /* y,yz */
			  mat [i6 * n11 + j10] += cfactor * (h1xz);          /* z,yz */
			  mat [i4 * n11 + j11] += cfactor * (- 2.0 * h1yz);  /* x,yy */
			  /*mat [i5 * n11 + j11] +=;*/           /* y,yy */
			  mat [i6 * n11 + j11] += cfactor * (2.0 * h1xy);    /* z,yy */

			  /* M part */
			  /* e_ijkl part */
			  mat [i7  * n11 + j7 ] += cfactor * m1xxxx; /* xx,xx */
			  mat [i7  * n11 + j8 ] += cfactor * m1xxxy; /* xx,xy */
			  mat [i7  * n11 + j9 ] += cfactor * m1xxxz; /* xx,xz */
			  mat [i7  * n11 + j10] += cfactor * m1xxyz; /* xx,yz */
			  mat [i7  * n11 + j11] += cfactor * m1xxyy; /* xx,yy */
			  mat [i8  * n11 + j7 ] += cfactor * m1xxxy; /* xy,xx */
			  mat [i8  * n11 + j8 ] += cfactor * m1xxyy; /* xy,xy */
			  mat [i8  * n11 + j9 ] += cfactor * m1xxyz; /* xy,xz */
			  mat [i8  * n11 + j10] += cfactor * m1xyyz; /* xy,yz */
			  mat [i8  * n11 + j11] += cfactor * m1xyyy; /* xy,yy */
			  mat [i9  * n11 + j7 ] += cfactor * m1xxxz; /* xz,xx */
			  mat [i9  * n11 + j8 ] += cfactor * m1xxyz; /* xz,xy */
			  mat [i9  * n11 + j9 ] += cfactor * m1xxzz; /* xz,xz */
			  mat [i9  * n11 + j10] += cfactor * m1xyzz; /* xz,yz */
			  mat [i9  * n11 + j11] += cfactor * m1xyyz; /* xz,yy */
			  mat [i10 * n11 + j7 ] += cfactor * m1xxyz; /* yz,xx */
			  mat [i10 * n11 + j8 ] += cfactor * m1xyyz; /* yz,xy */
			  mat [i10 * n11 + j9 ] += cfactor * m1xyzz; /* yz,xz */
			  mat [i10 * n11 + j10] += cfactor * m1yyzz; /* yz,yz */
			  mat [i10 * n11 + j11] += cfactor * m1yyyz; /* yz,yy */
			  mat [i11 * n11 + j7 ] += cfactor * m1xxyy; /* yy,xx */
			  mat [i11 * n11 + j8 ] += cfactor * m1xyyy; /* yy,xy */
			  mat [i11 * n11 + j9 ] += cfactor * m1xyyz; /* yy,xz */
			  mat [i11 * n11 + j10] += cfactor * m1yyyz; /* yy,yz */
			  mat [i11 * n11 + j11] += cfactor * m1yyyy; /* yy,yy */

			  /* (d_ij e_kl + e_ij d_kl)part */
			  mat [i7  * n11 + j7 ] += cfactor * 2.0 * m2xx; /* xx,xx */
			  mat [i7  * n11 + j8 ] += cfactor * m2xy; /* xx,xy */
			  mat [i7  * n11 + j9 ] += cfactor * m2xz; /* xx,xz */
			  mat [i7  * n11 + j10] += cfactor * m2yz; /* xx,yz */
			  mat [i7  * n11 + j11] += cfactor * (m2yy + m2xx); /* xx,yy */
			  mat [i8  * n11 + j7 ] += cfactor * m2xy; /* xy,xx */
			  /*mat [i8  * n11 + j8 ] +=;*/ /* xy,xy */
			  /*mat [i8  * n11 + j9 ] +=;*/ /* xy,xz */
			  /*mat [i8  * n11 + j10] +=;*/ /* xy,yz */
			  mat [i8  * n11 + j11] += cfactor * m2xy; /* xy,yy */
			  mat [i9  * n11 + j7 ] += cfactor * m2xz; /* xz,xx */
			  /*mat [i9  * n11 + j8 ] +=;*/ /* xz,xy */
			  /*mat [i9  * n11 + j9 ] +=;*/ /* xz,xz */
			  /*mat [i9  * n11 + j10] +=;*/ /* xz,yz */
			  mat [i9  * n11 + j11] += cfactor * m2xz; /* xz,yy */
			  mat [i10 * n11 + j7 ] += cfactor * m2yz; /* yz,xx */
			  /*mat [i10 * n11 + j8 ] +=;*/ /* yz,xy */
			  /*mat [i10 * n11 + j9 ] +=;*/ /* yz,xz */
			  /*mat [i10 * n11 + j10] +=;*/ /* yz,yz */
			  mat [i10 * n11 + j11] += cfactor * m2yz; /* yz,yy */
			  mat [i11 * n11 + j7 ] += cfactor * (m2xx + m2yy); /* yy,xx */
			  mat [i11 * n11 + j8 ] += cfactor * m2xy; /* yy,xy */
			  mat [i11 * n11 + j9 ] += cfactor * m2xz; /* yy,xz */
			  mat [i11 * n11 + j10] += cfactor * m2yz; /* yy,yz */
			  mat [i11 * n11 + j11] += cfactor * 2.0 * m2yy; /* yy,yy */

			  /* d_ij d_kl part */
			  mat [i7  * n11 + j7 ] += cfactor * m3; /* xx,xx */
			  /*mat [i7  * n11 + j8 ] +=;*/ /* xx,xy */
			  /*mat [i7  * n11 + j9 ] +=;*/ /* xx,xz */
			  /*mat [i7  * n11 + j10] +=;*/ /* xx,yz */
			  mat [i7  * n11 + j11] += cfactor * m3; /* xx,yy */
			  /*mat [i8  * n11 + j7 ] +=;*/ /* xy,xx */
			  /*mat [i8  * n11 + j8 ] +=;*/ /* xy,xy */
			  /*mat [i8  * n11 + j9 ] +=;*/ /* xy,xz */
			  /*mat [i8  * n11 + j10] +=;*/ /* xy,yz */
			  /*mat [i8  * n11 + j11] +=;*/ /* xy,yy */
			  /*mat [i9  * n11 + j7 ] +=;*/ /* xz,xx */
			  /*mat [i9  * n11 + j8 ] +=;*/ /* xz,xy */
			  /*mat [i9  * n11 + j9 ] +=;*/ /* xz,xz */
			  /*mat [i9  * n11 + j10] +=;*/ /* xz,yz */
			  /*mat [i9  * n11 + j11] +=;*/ /* xz,yy */
			  /*mat [i10 * n11 + j7 ] +=;*/ /* yz,xx */
			  /*mat [i10 * n11 + j8 ] +=;*/ /* yz,xy */
			  /*mat [i10 * n11 + j9 ] +=;*/ /* yz,xz */
			  /*mat [i10 * n11 + j10] +=;*/ /* yz,yz */
			  /*mat [i10 * n11 + j11] +=;*/ /* yz,yy */
			  mat [i11 * n11 + j7 ] += cfactor * m3; /* yy,xx */
			  /*mat [i11 * n11 + j8 ] +=;*/ /* yy,xy */
			  /*mat [i11 * n11 + j9 ] +=;*/ /* yy,xz */
			  /*mat [i11 * n11 + j10] +=;*/ /* yy,yz */
			  mat [i11 * n11 + j11] += cfactor * m3; /* yy,yy */

			  /* e_i d_jk e_l + e_i d_jl e_k + */
			  /* e_j d_ik e_l + e_j d_il e_k part */
			  mat [i7  * n11 + j7 ] += cfactor * 4.0 * m4xx; /* xx,xx */
			  mat [i7  * n11 + j8 ] += cfactor * 2.0 * m4xy; /* xx,xy */
			  mat [i7  * n11 + j9 ] += cfactor * 2.0 * m4xz; /* xx,xz */
			  /*mat [i7  * n11 + j10] +=;*/ /* xx,yz */
			  /*mat [i7  * n11 + j11] +=;*/ /* xx,yy */
			  mat [i8  * n11 + j7 ] += cfactor * 2.0 * m4xy; /* xy,xx */
			  mat [i8  * n11 + j8 ] += cfactor * (m4xx + m4yy); /* xy,xy */
			  mat [i8  * n11 + j9 ] += cfactor * m4yz; /* xy,xz */
			  mat [i8  * n11 + j10] += cfactor * m4xz; /* xy,yz */
			  mat [i8  * n11 + j11] += cfactor * 2.0 * m4xy; /* xy,yy */
			  mat [i9  * n11 + j7 ] += cfactor * 2.0 * m4xz; /* xz,xx */
			  mat [i9  * n11 + j8 ] += cfactor * m4yz; /* xz,xy */
			  mat [i9  * n11 + j9 ] += cfactor * (m4xx + m4zz); /* xz,xz */
			  mat [i9  * n11 + j10] += cfactor * m4xy; /* xz,yz */
			  /*mat [i9  * n11 + j11] +=;*/ /* xz,yy */
			  /*mat [i10 * n11 + j7 ] +=;*/ /* yz,xx */
			  mat [i10 * n11 + j8 ] += cfactor * m4xz; /* yz,xy */
			  mat [i10 * n11 + j9 ] += cfactor * m4xy; /* yz,xz */
			  mat [i10 * n11 + j10] += cfactor * (m4yy + m4zz); /* yz,yz */
			  mat [i10 * n11 + j11] += cfactor * 2.0 * m4yz; /* yz,yy */
			  /*mat [i11 * n11 + j7 ] +=;*/ /* yy,xx */
			  mat [i11 * n11 + j8 ] += cfactor * 2.0 * m4xy; /* yy,xy */
			  /*mat [i11 * n11 + j9 ] +=;*/ /* yy,xz */
			  mat [i11 * n11 + j10] += cfactor * 2.0 * m4yz; /* yy,yz */
			  mat [i11 * n11 + j11] += cfactor * 4.0 * m4yy; /* yy,yy */

			  /* (d_ik d_jl + d_jk d_il) part */
			  mat [i7  * n11 + j7 ] += cfactor * 2.0 * mm5; /* xx,xx */
			  /* mat [i7  * n11 + j8 ] +=;*/ /* xx,xy */
			  /* mat [i7  * n11 + j9 ] +=;*/ /* xx,xz */
			  /* mat [i7  * n11 + j10] +=;*/ /* xx,yz */
			  /* mat [i7  * n11 + j11] +=;*/ /* xx,yy */
			  /* mat [i8  * n11 + j7 ] +=;*/ /* xy,xx */
			  mat [i8  * n11 + j8 ] += cfactor * mm5; /* xy,xy */
			  /* mat [i8  * n11 + j9 ] +=;*/ /* xy,xz */
			  /* mat [i8  * n11 + j10] +=;*/ /* xy,yz */
			  /* mat [i8  * n11 + j11] +=;*/ /* xy,yy */
			  /* mat [i9  * n11 + j7 ] +=;*/ /* xz,xx */
			  /* mat [i9  * n11 + j8 ] +=;*/ /* xz,xy */
			  mat [i9  * n11 + j9 ] += cfactor * mm5; /* xz,xz */
			  /* mat [i9  * n11 + j10] +=;*/ /* xz,yz */
			  /* mat [i9  * n11 + j11] +=;*/ /* xz,yy */
			  /* mat [i10 * n11 + j7 ] +=;*/ /* yz,xx */
			  /* mat [i10 * n11 + j8 ] +=;*/ /* yz,xy */
			  /* mat [i10 * n11 + j9 ] +=;*/ /* yz,xz */
			  mat [i10 * n11 + j10] += cfactor * mm5; /* yz,yz */
			  /* mat [i10 * n11 + j11] +=;*/ /* yz,yy */
			  /* mat [i11 * n11 + j7 ] +=;*/ /* yy,xx */
			  /* mat [i11 * n11 + j8 ] +=;*/ /* yy,xy */
			  /* mat [i11 * n11 + j9 ] +=;*/ /* yy,xz */
			  /* mat [i11 * n11 + j10] +=;*/ /* yy,yz */
			  mat [i11 * n11 + j11] += cfactor * 2.0 * mm5; /* yy,yy */
			}
		    }
		}
	    }
	}
    }

#ifdef ZETA
  ctmp3 = clock ();

  cpu1 = (double) (ctmp3 - ctmp1);
  cpu2 = (double) (ctmp2 - ctmp1);
  cpu3 = (double) (ctmp3 - ctmp2);
#else
  cholesky (mat, n11);
#endif /* ZETA */
}
