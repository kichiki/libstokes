/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-1996,1999-2000 Kengo Ichiki
 *               <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-3fts.c,v 2.2 2000/12/08 06:13:40 ichiki Exp $
 *
 * 3 dimensional hydrodynamics, 3D configuration
 * periodic boundary condition in 3 direction,
 * FT version
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

#include "ewald-3ft.h"

/* calc ewald-summed mobility for FT version
 * INPUT
 *  n : # of all (mobile and fixed) particles
 *  x [n * 3] : position of particles
 * OUTPUT
 *  mat [n * 6 * n * 6] : mobility matrix
 */
void
calc_mob_ewald_3ft (int n, double *x, double *mat)
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
  double ex, ey, ez;
  double exx, eyy, ezz, exy, eyz, ezx;
  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;
  double rlx, rly, rlz;

  int n6;
  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int i1, i2, i3, i4, i5, i6;
  int j1, j2, j3, j4, j5, j6;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double k;
  double cfactor;
  double mob;
  double pi8vol;


  /* temporary setting : we should be modefy GLOBAL pi6vol */
  pi8vol = pi6vol * 8.0 / 6.0; /* 8 pi / Vol */

  n6 = n * 6;
  /* clear matrix */
  for (i = 0; i < n6 * n6; i++)
    {
      mat [i] = 0.0;
    }
  /* diagonal part ( self part ) */
  for (i = 0; i < n; i++)
    {
      i1 = i * 6;
      i2 = i1 + 1;
      i3 = i1 + 2;
      i4 = i1 + 3;
      i5 = i1 + 4;
      i6 = i1 + 5;

      /* A part */
      mat [i1 * n6 + i1] = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
      mat [i2 * n6 + i2] = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
      mat [i3 * n6 + i3] = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);

      /* C part */
      mat [i4 * n6 + i4] = 0.75 - zaspi * 10.0 * za2;
      mat [i5 * n6 + i5] = 0.75 - zaspi * 10.0 * za2;
      mat [i6 * n6 + i6] = 0.75 - zaspi * 10.0 * za2;
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
      i1 = i * 6;
      i2 = i1 + 1;
      i3 = i1 + 2;
      i4 = i1 + 3;
      i5 = i1 + 4;
      i6 = i1 + 5;
      for (j = 0; j < n; j++)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;
	  j1 = j * 6;
	  j2 = j1 + 1;
	  j3 = j1 + 2;
	  j4 = j1 + 3;
	  j5 = j1 + 4;
	  j6 = j1 + 5;

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
			  ezx = ez * ex;

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

			  axx = (a1 + a2 * exx);
			  ayy = (a1 + a2 * eyy);
			  azz = (a1 + a2 * ezz);
			  axy = a2 * exy;
			  ayz = a2 * eyz;
			  azx = a2 * ezx;

			  bx = b1 * ex;
			  by = b1 * ey;
			  bz = b1 * ez;

			  cxx = (c1 + c2 * exx);
			  cyy = (c1 + c2 * eyy);
			  czz = (c1 + c2 * ezz);
			  cxy = c2 * exy;
			  cyz = c2 * eyz;
			  czx = c2 * ezx;

			  /* A part */
			  mat [i1 * n6 + j1] += axx;
			  mat [i2 * n6 + j2] += ayy;
			  mat [i3 * n6 + j3] += azz;

			  mat [i1 * n6 + j2] += axy;
			  mat [i2 * n6 + j1] += axy;
			  mat [i2 * n6 + j3] += ayz;
			  mat [i3 * n6 + j2] += ayz;
			  mat [i3 * n6 + j1] += azx;
			  mat [i1 * n6 + j3] += azx;

			  /* B part */
			  /* note: r = x_i-x_j so that jj = 11 and ji = 12 */
			  /* B21 = -B12*/
			  mat [i4 * n6 + j2] += - bz;
			  mat [i4 * n6 + j3] +=   by;
			  mat [i5 * n6 + j1] +=   bz;
			  mat [i5 * n6 + j3] += - bx;
			  mat [i6 * n6 + j1] += - by;
			  mat [i6 * n6 + j2] +=   bx;

			  /* B_ part */
			  /* BT21 = -B12 */
			  mat [i1 * n6 + j5] += - bz;
			  mat [i1 * n6 + j6] +=   by;
			  mat [i2 * n6 + j4] +=   bz;
			  mat [i2 * n6 + j6] += - bx;
			  mat [i3 * n6 + j4] += - by;
			  mat [i3 * n6 + j5] +=   bx;

			  /* C part */
			  mat [i4 * n6 + j4] += cxx;
			  mat [i5 * n6 + j5] += cyy;
			  mat [i6 * n6 + j6] += czz;

			  mat [i4 * n6 + j5] += cxy;
			  mat [i5 * n6 + j4] += cxy;
			  mat [i5 * n6 + j6] += cyz;
			  mat [i6 * n6 + j5] += cyz;
			  mat [i6 * n6 + j4] += czx;
			  mat [i4 * n6 + j6] += czx;
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

		  a1 = 0.75 * (1.0 - kk / 3.0) * mob;
		  b1 = - 3.0 / 8.0 * k * mob;
		  c1 = 3.0 / 16.0 * kk * mob;

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

		  for (i = 0; i < n; i++)
		    {
		      ix = i * 3;
		      iy = ix + 1;
		      iz = ix + 2;
		      i1 = i * 6;
		      i2 = i1 + 1;
		      i3 = i1 + 2;
		      i4 = i1 + 3;
		      i5 = i1 + 4;
		      i6 = i1 + 5;
		      for (j = 0; j < n; j++)
			{
			  jx = j * 3;
			  jy = jx + 1;
			  jz = jx + 2;
			  j1 = j * 6;
			  j2 = j1 + 1;
			  j3 = j1 + 2;
			  j4 = j1 + 3;
			  j5 = j1 + 4;
			  j6 = j1 + 5;

			  xx = x [jx] - x [ix];
			  yy = x [jy] - x [iy];
			  zz = x [jz] - x [iz];

			  cfactor = cos (+ k1 * xx
					 + k2 * yy
					 + k3 * zz);

			  /* A part */
			  mat [i1 * n6 + j1] += cfactor * axx;
			  mat [i2 * n6 + j2] += cfactor * ayy;
			  mat [i3 * n6 + j3] += cfactor * azz;

			  mat [i1 * n6 + j2] += cfactor * axy;
			  mat [i2 * n6 + j1] += cfactor * axy;
			  mat [i2 * n6 + j3] += cfactor * ayz;
			  mat [i3 * n6 + j2] += cfactor * ayz;
			  mat [i3 * n6 + j1] += cfactor * azx;
			  mat [i1 * n6 + j3] += cfactor * azx;

			  /* B part */
			  /* note: r = x_i-x_j so that jj = 11 and ji = 12 */
			  /* B21 = -B12*/
			  mat [i4 * n6 + j2] += - cfactor * bz;
			  mat [i4 * n6 + j3] +=   cfactor * by;
			  mat [i5 * n6 + j1] +=   cfactor * bz;
			  mat [i5 * n6 + j3] += - cfactor * bx;
			  mat [i6 * n6 + j1] += - cfactor * by;
			  mat [i6 * n6 + j2] +=   cfactor * bx;

			  /* B_ part */
			  /* BT21 = -B12 */
			  mat [i1 * n6 + j5] += - cfactor * bz;
			  mat [i1 * n6 + j6] +=   cfactor * by;
			  mat [i2 * n6 + j4] +=   cfactor * bz;
			  mat [i2 * n6 + j6] += - cfactor * bx;
			  mat [i3 * n6 + j4] += - cfactor * by;
			  mat [i3 * n6 + j5] +=   cfactor * bx;

			  /* C part */
			  mat [i4 * n6 + j4] += cfactor * cxx;
			  mat [i5 * n6 + j5] += cfactor * cyy;
			  mat [i6 * n6 + j6] += cfactor * czz;

			  mat [i4 * n6 + j5] += cfactor * cxy;
			  mat [i5 * n6 + j4] += cfactor * cxy;
			  mat [i5 * n6 + j6] += cfactor * cyz;
			  mat [i6 * n6 + j5] += cfactor * cyz;
			  mat [i6 * n6 + j4] += cfactor * czx;
			  mat [i4 * n6 + j6] += cfactor * czx;
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
  cholesky (mat, n6);
#endif /* ZETA */
}
