/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-1996,1999-2000 Kengo Ichiki
 *               <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-3fts.c,v 2.5 2000/12/09 10:11:17 ichiki Exp $
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
/*#else
#include "cholesky.h" /* cholesky() */
#endif /* ZETA */

#include "fts.h"
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

  double xa, ya; 
  double yb;
  double xc, yc;
  double xg, yg;
  double yh;
  double xm, ym, zm;

  double ex, ey, ez;

  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;
  double rlx, rly, rlz;

  int n11;
  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double k;
  double cf, sf;
  double mob;
  double pi8vol;

  double mm5; /* for self part */
  double a2, c2;


  /* temporary setting : we should be modefy GLOBAL pi6vol */
  pi8vol = pi6vol * 8.0 / 6.0; /* 8 pi / Vol */

  n11 = n * 11;
  /* clear matrix */
  for (i = 0; i < n11 * n11; i++)
    mat [i] = 0.0;

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
      for (j = 0; j < n; j++)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;

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

			  ya = (0.75 + 0.5 / s2) / s * erfc (zr)
			    + ((1.0 + zr2 *
				(14.0 + 4.0 * zr2 *
				 (- 5.0 + zr2))) / s2
			       - 4.5 + 3.0 * zr2) * zaspi * exp (- zr2);
			  a2 = (0.75 - 1.5 / s2) / s * erfc (zr)
			    + ((- 3.0 + zr2 *
				(- 2.0 + 4.0 * zr2 *
				 (4.0 - zr2))) / s2
			       + 1.5 - 3.0 * zr2) * zaspi * exp (- zr2);
			  xa = a2 + ya;
	      
			  yb = - 0.75 / s2 * erfc (zr)
			    - 1.5 * (1.0 + zr2 *
				     (- 6.0 + zr2 *
				      2.0)) / s * zaspi * exp (- zr2);

			  yc = - 3.0 / 8.0 / s2 / s * erfc (zr)
			    - 0.75 * (1.0 + zr2 *
				      (14.0 + zr2 *
				       (-20.0 + zr2 *
					4.0))) / s2 * zaspi * exp (- zr2);
			  c2 = 9.0 / 8.0 / s2 / s * erfc (zr)
			    - 0.75 * (- 3.0 + zr2 *
				      (- 2.0 + zr2 *
				       (16.0 + zr2 *
					- 4.0))) / s2 * zaspi * exp (- zr2);
			  xc = c2 + yc;
	      
			  xg = (2.25 - 3.6 / s2) / s2 * erfc (zr)
			    + (- 1.5 * (- 3.0 + zr2 *
					6.0)
			       - 0.8 * (9.0 + zr2 *
					(6.0 + zr2 *
					 (- 48.0 + zr2 *
					  12.0))) / s2)
			    / s * zaspi * exp (- zr2);
			  yg = 1.2 / s2 / s2 * erfc (zr)
			    + (- 3.0 * ( zr2 *
					 (2.0 + zr2 *
					  (- 1.0)))
			       - 0.8 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (- 26.0 + zr2 *
					  (+ 26.0 + zr2 *
					   (- 4.0))))) / s2)
			    / s * zaspi * exp (- zr2);

			  yh = - 9.0 / 8.0 / s2 / s * erfc (zr)
			    + 1.5 * (- 1.5 + zr2 *
				     (- 1.0 + zr2 *
				      (+ 8.0 + zr2 *
				       (- 2.0))))
			    / s2 * zaspi * exp (- zr2);

			  xm = (- 4.5 + 10.8 / s2) / s / s2 * erfc (zr)
			    + (1.5 * (- 6.0 +  zr2 *
				      (- 12.0 + zr2 *
				       12.0))
			       + 0.6 * (18.0 + zr2 *
					(12.0 + zr2 *
					 (30.0 + zr2 *
					  (- 66.0 + zr2 *
					   12.0)))) / s2)
			    / s2 * zaspi * exp (- zr2);
			  ym = (4.25 - 7.2 / s2) / s / s2 * erfc (zr)
			    + (- 1.5 * (- 3.0 +  zr2 *
					(6.0 + zr2 *
					 (- 12.0 + zr2 *
					  4.0)))
			       - 0.6 * (12.0 + zr2 *
					(8.0 + zr2 *
					 (- 22.0 + zr2 *
					  (58.0 + zr2 *
					   (- 34.0 + zr2 *
					    4.0))))) / s2)
			    / s2 * zaspi * exp (- zr2);
			  zm = 4.5 / s2 / s / s2 * erfc (zr)
			    + (- 1.5 * (0.0 +  zr2 *
					(8.0 + zr2 *
					 (- 4.0)))
			       - 0.6 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (- 26.0 + zr2 *
					  (+ 26.0 + zr2 *
					   (- 4.0))))) / s2)
			    / s2 * zaspi * exp (- zr2);
	      
			  matrix_fts_ij (j, i,
					 ex, ey, ez,
					 xa, ya,
					 yb,
					 xc, yc,
					 xg, yg,
					 yh,
					 xm, ym, zm,
					 n11, mat);
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

		  xa = 0.0;
		  ya = 0.75 * (1.0 - kk / 3.0) * mob;

		  yb = 3.0 / 8.0 * k * mob;

		  xc = 0.0;
		  yc = 3.0 / 16.0 * kk * mob;

		  xg = 0.0;
		  yg = 3.0 / 8.0 * (1.0 - 4.0 / 15.0 * kk) * k * mob;

		  yh = 3.0 / 16.0 * kk * mob;
      
		  xm = 0.0;
		  ym = 3.0 / 8.0 * (1.0 - kk / 5.0) * kk * mob;
		  zm = 0.0;
      
		  for (i = 0; i < n; i++)
		    {
		      ix = i * 3;
		      iy = ix + 1;
		      iz = ix + 2;
		      for (j = 0; j < n; j++)
			{
			  jx = j * 3;
			  jy = jx + 1;
			  jz = jx + 2;

			  xx = x [jx] - x [ix];
			  yy = x [jy] - x [iy];
			  zz = x [jz] - x [iz];

			  cf = cos (+ k1 * xx
				    + k2 * yy
				    + k3 * zz);

			  sf = - sin (+ k1 * xx
				      + k2 * yy
				      + k3 * zz);

			  matrix_fts_ij (j, i,
					 ex, ey, ez,
					 cf * xa, cf * ya,
					 sf * yb,
					 cf * xc, cf * yc,
					 sf * xg, sf * yg,
					 cf * yh,
					 cf * xm, cf * ym,
					 cf * zm,
					 n11, mat);
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
  /*#else
    cholesky (mat, n11);*/
#endif /* ZETA */
}
