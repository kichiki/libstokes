/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-1996,1999-2000 Kengo Ichiki
 *               <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-3f.c,v 2.2 2000/12/11 06:28:03 ichiki Exp $
 *
 * 3 dimensional hydrodynamics, 3D configuration
 * periodic boundary condition in 3 direction,
 * F version
 * non-dimension formulation
 */
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */

#ifdef ZETA
#include <time.h> /* clock() */
/*#else
  #include "cholesky.h"*/ /* cholesky() */
#endif /* ZETA */

#include "ewald-3f.h"

/* calc ewald-summed mobility for F version
 * INPUT
 *  n : # of all (mobile and fixed) particles
 *  x [n * 3] : position of particles
 * OUTPUT
 *  mat [n * 3 * n * 3] : mobility matrix
 */
void
calc_mob_ewald_3f (int n, double *x, double *mat)
{
  extern int pcellx, pcelly, pcellz;
  extern int kmaxx, kmaxy, kmaxz;

  extern double zeta, zeta2, zaspi, za2;
  extern double pi2;
  extern double pivol;
  extern double lx, ly, lz; /* cell size */

#ifdef ZETA
  extern double cpu1, cpu2, cpu3;
  clock_t ctmp1, ctmp2, ctmp3;
#endif /* ZETA */

  double ya;

  double ex, ey, ez;
  double exx, eyy, ezz, exy, eyz, ezx;

  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;
  double rlx, rly, rlz;

  int n3;
  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double cf;
  double kexp;

  double mob11, mob22, mob33, mob12, mob23, mob31;

  double erfczr;
  double expzr2;

  double a2;


  n3 = n * 3;
  /* clear matrix */
  for (i = 0; i < n3 * n3; i++)
    mat [i] = 0.0;

  /* diagonal part ( self part ) */
  ya = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
  for (i = 0; i < n3; i++)
    mat [i * n3 + i] = ya;

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

			  erfczr = erfc (zr);
			  expzr2 = zaspi * exp (- zr2);

			  ex = xx / rr;
			  ey = yy / rr;
			  ez = zz / rr;

			  exx = ex * ex;
			  eyy = ey * ey;
			  ezz = ez * ez;
			  exy = ex * ey;
			  eyz = ey * ez;
			  ezx = ez * ex;

			  ya = (0.75 + 0.5 / s2) / s * erfczr
			    + ((1.0 + zr2 *
				(14.0 + 4.0 * zr2 *
				 (- 5.0 + zr2))) / s2
			       - 4.5 + 3.0 * zr2)
			    * expzr2;
			  a2 = (0.75 - 1.5 / s2) / s * erfczr
			    + ((- 3.0 + zr2 *
				(- 2.0 + 4.0 * zr2 *
				 (4.0 - zr2))) / s2
			       + 1.5 - 3.0 * zr2)
			    * expzr2;
			  /*xa = a2 + ya;*/

			  mat [ix * n3 + jx] += ya + a2 * exx;
			  mat [iy * n3 + jy] += ya + a2 * eyy;
			  mat [iz * n3 + jz] += ya + a2 * ezz;

			  mat [ix * n3 + jy] += a2 * exy;
			  mat [iy * n3 + jx] += a2 * exy;
			  mat [iy * n3 + jz] += a2 * eyz;
			  mat [iz * n3 + jy] += a2 * eyz;
			  mat [iz * n3 + jx] += a2 * ezx;
			  mat [ix * n3 + jz] += a2 * ezx;
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
		  k4z = kk / 4.0 / zeta2;
		  kexp = pivol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  ya = 6.0 * (1.0 - kk / 3.0) * kexp;

		  mob11 = (1.0 - k1 * k1 / kk) * ya;
		  mob22 = (1.0 - k2 * k2 / kk) * ya;
		  mob33 = (1.0 - k3 * k3 / kk) * ya;
		  mob12 = - k1 * k2 / kk * ya;
		  mob23 = - k2 * k3 / kk * ya;
		  mob31 = - k3 * k1 / kk * ya;

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

			  mat [ix * n3 + jx] += cf * mob11;
			  mat [iy * n3 + jy] += cf * mob22;
			  mat [iz * n3 + jz] += cf * mob33;

			  mat [ix * n3 + jy] += cf * mob12;
			  mat [iy * n3 + jx] += cf * mob12;
			  mat [iy * n3 + jz] += cf * mob23;
			  mat [iz * n3 + jy] += cf * mob23;
			  mat [iz * n3 + jx] += cf * mob31;
			  mat [ix * n3 + jz] += cf * mob31;
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
    cholesky (mat, n3);*/
#endif /* ZETA */
}
