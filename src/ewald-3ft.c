/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-1996,1999-2000 Kengo Ichiki
 *               <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-3ft.c,v 2.1 2000/12/08 05:04:14 ichiki Exp $
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
#else /* not ZETA */
#include "cholesky.h" /* cholesky() */
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
  extern double pi6vol;
  extern double lx, ly, lz; /* cell size */

#ifdef ZETA
  extern double cpu1, cpu2, cpu3;
  clock_t ctmp1, ctmp2, ctmp3;
#endif /* ZETA */

  double m1i, m1n;
  double m1xx, m1yy, m1zz, m1xy, m1yz, m1zx;
  double ex, ey, ez;
  double exx, eyy, ezz, exy, eyz, ezx;
  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;
  double rlx, rly, rlz;

  int n3;
  int i, j;
  int i1, i2, i3, j1, j2, j3;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double cfactor;
  double mob, mob11, mob22, mob33, mob12, mob23, mob31;


  n3 = n * 3;
  /* clear matrix */
  for (i = 0; i < n3 * n3; i++)
    {
      mat [i] = 0.0;
    }
  /* diagonal part ( self part ) */
  for (i = 0; i < n3; i++)
    {
      mat [i * n3 + i] = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
    }

#ifdef ZETA
  ctmp1 = clock ();
#endif /* ZETA */

  /* first Ewald part ( real space ) */
  for (i = 0; i < n; i++)
    {
      i1 = i * 3;
      i2 = i * 3 + 1;
      i3 = i * 3 + 2;
      for (j = 0; j < n; j++)
	{
	  j1 = j * 3;
	  j2 = j * 3 + 1;
	  j3 = j * 3 + 2;

	  for (m1 = - pcellx; m1 <= pcellx; m1++)
	    {
	      rlx = lx * (double) m1;
	      for (m2 = - pcelly; m2 <= pcelly; m2++)
		{
		  rly = ly * (double) m2;
		  for (m3 = - pcellz; m3 <= pcellz; m3++)
		    {
		      rlz = lz * (double) m3;
  
		      xx = x [j1] - x [i1] + rlx;
		      yy = x [j2] - x [i2] + rly;
		      zz = x [j3] - x [i3] + rlz;
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

			  m1i = (0.75 + 0.5 / s2) / s * erfc (zr)
			    + ((1.0 + zr2 *
				(14.0 + 4.0 * zr2 *
				 (- 5.0 + zr2))) / s2
			       - 4.5 + 3.0 * zr2) * zaspi * exp (- zr2);
			  m1n = (0.75 - 1.5 / s2) / s * erfc (zr)
			    + ((- 3.0 + zr2 *
				(- 2.0 + 4.0 * zr2 *
				 (4.0 - zr2))) / s2
			       + 1.5 - 3.0 * zr2) * zaspi * exp (- zr2);

			  m1xx = (m1i + m1n * exx);
			  m1yy = (m1i + m1n * eyy);
			  m1zz = (m1i + m1n * ezz);
			  m1xy = m1n * exy;
			  m1yz = m1n * eyz;
			  m1zx = m1n * ezx;

			  mat[i1 * n3 + j1] += m1xx;
			  mat[i2 * n3 + j2] += m1yy;
			  mat[i3 * n3 + j3] += m1zz;

			  mat[i1 * n3 + j2] += m1xy;
			  mat[i2 * n3 + j1] += m1xy;
			  mat[i2 * n3 + j3] += m1yz;
			  mat[i3 * n3 + j2] += m1yz;
			  mat[i3 * n3 + j1] += m1zx;
			  mat[i1 * n3 + j3] += m1zx;
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
		  mob = pi6vol
		    * (1.0 - kk / 3.0)
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  mob11 = (1.0 - k1 * k1 / kk) * mob;
		  mob22 = (1.0 - k2 * k2 / kk) * mob;
		  mob33 = (1.0 - k3 * k3 / kk) * mob;
		  mob12 = - k1 * k2 / kk * mob;
		  mob23 = - k2 * k3 / kk * mob;
		  mob31 = - k3 * k1 / kk * mob;

		  for (i = 0; i < n; i++)
		    {
		      i1 = i * 3;
		      i2 = i * 3 + 1;
		      i3 = i * 3 + 2;
		      for (j = 0; j < n; j++)
			{
			  j1 = j * 3;
			  j2 = j * 3 + 1;
			  j3 = j * 3 + 2;

			  xx = x [j1] - x [i1];
			  yy = x [j2] - x [i2];
			  zz = x [j3] - x [i3];

			  cfactor = cos (+ k1 * xx
					 + k2 * yy
					 + k3 * zz);

			  mat [i1 * n3 + j1] += cfactor * mob11;
			  mat [i2 * n3 + j2] += cfactor * mob22;
			  mat [i3 * n3 + j3] += cfactor * mob33;

			  mat [i1 * n3 + j2] += cfactor * mob12;
			  mat [i2 * n3 + j1] += cfactor * mob12;
			  mat [i2 * n3 + j3] += cfactor * mob23;
			  mat [i3 * n3 + j2] += cfactor * mob23;
			  mat [i3 * n3 + j1] += cfactor * mob31;
			  mat [i1 * n3 + j3] += cfactor * mob31;
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
  cholesky (mat, n3);
#endif /* ZETA */
}
