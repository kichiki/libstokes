/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-1996,1999-2001 Kengo Ichiki
 *               <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-3fts.c,v 3.2 2001/01/24 08:14:23 ichiki Exp $
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
#endif /* ZETA */

#include <libiter.h> /* solve_iter_stab (), gpb () */

#include "fts.h"
#include "ewald-3fts.h"


/** function prototypes for local routines **/
/* utility routines for calc_mob_ewald_3fts () */
static void
calc_b_mob_ewald_3fts (int n,
		       double *f, double *t, double *e,
		       double *b);
static void
atimes_mob_ewald_3fts (int n, double *x, double *y);



/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for FTS version
 * INPUT
 *  (global) pos [] : position of particles
 *  n := np * 11
 *  x [n * 11] : FTS
 * OUTPUT
 *  y [n * 11] : UOE
 */
void
atimes_ewald_3fts (int n, double *x, double *y)
{
  extern int pcellx, pcelly, pcellz;
  extern int kmaxx, kmaxy, kmaxz;

  extern double zeta, zeta2, zaspi, za2;
  extern double pi2;
  extern double pivol;
  extern double lx, ly, lz; /* cell size */

  extern double *pos;

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

  int np;
  int i, j;
  int i11, j11;
  int ix, iy, iz;
  int jx, jy, jz;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double k;
  double cf, sf;
  double kexp;

  double erfczr;
  double expzr2;

  double a2, c2;


  np = n / 11;

  /* clear result */
  for (i = 0; i < n; i ++)
    y [i] = 0.0;

  /* diagonal part ( self part ) */
  xa = ya = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
  xc = yc = 0.75 - zaspi * za2 * 10.0;
  xm = ym = zm = 0.9 - zaspi * za2 * (12.0 - 30.24 * za2);
  /*xm = ym = zm = 0.9 - 3.0 * zaspi * za2 * (2.0 - 5.04 * za2);*/
  for (i = 0; i < np; i++)
    {
      i11 = i * 11;
      matrix_fts_atimes (x + i11, y + i11,
			 0.0, 0.0, 0.0,
			 xa, ya,
			 0.0,
			 xc, yc,
			 0.0, 0.0,
			 0.0,
			 xm, ym, zm);
    }

#ifdef ZETA
  ctmp1 = clock ();
#endif /* ZETA */

  /* first Ewald part ( real space ) */
  for (i = 0; i < np; i++)
    {
      i11 = i * 11;
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = 0; j < np; j++)
	{
	  j11 = j * 11;
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
  
		      xx = pos [jx] - pos [ix] + rlx;
		      yy = pos [jy] - pos [iy] + rly;
		      zz = pos [jz] - pos [iz] + rlz;
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
			  xa = a2 + ya;
	      
			  yb = - 0.75 / s2 * erfczr
			    - 1.5 * (+ 1.0 + zr2 *
				     (- 6.0 + zr2 *
				      (+ 2.0)))
			    / s * expzr2;

			  yc = - 3.0 / 8.0 / s2 / s * erfczr
			    - 0.75 * (+ 1.0 + zr2 *
				      (+ 14.0 + zr2 *
				       (-20.0 + zr2 *
					( + 4.0))))
			    / s2 * expzr2;
			  c2 = 9.0 / 8.0 / s2 / s * erfczr
			    - 0.75 * (- 3.0 + zr2 *
				      (- 2.0 + zr2 *
				       (+ 16.0 + zr2 *
					(- 4.0))))
			    / s2 * expzr2;
			  xc = c2 + yc;
	      
			  xg = (2.25 - 3.6 / s2) / s2 * erfczr
			    + (- 1.5 * (- 3.0 + zr2 *
					(+ 6.0))
			       - 0.8 * (+ 9.0 + zr2 *
					(+ 6.0 + zr2 *
					 (- 48.0 + zr2 *
					  (+ 12.0)))) / s2)
			    / s * expzr2;
			  yg = 1.2 / s2 / s2 * erfczr
			    + (- 3.0 * ( zr2 *
					 (2.0 + zr2 *
					  (- 1.0)))
			       - 0.8 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (- 26.0 + zr2 *
					  (+ 26.0 + zr2 *
					   (- 4.0))))) / s2)
			    / s * expzr2;

			  yh = - 9.0 / 8.0 / s2 / s * erfczr
			    + 1.5 * (- 1.5 + zr2 *
				     (- 1.0 + zr2 *
				      (+ 8.0 + zr2 *
				       (- 2.0))))
			    / s2 * expzr2;

			  xm = (- 4.5 + 10.8 / s2) / s / s2 * erfczr
			    + (+ 1.5 * (- 6.0 +  zr2 *
					(- 12.0 + zr2 *
					 (+ 12.0)))
			       + 1.2 * (+ 18.0 + zr2 *
					(+ 12.0 + zr2 *
					 (+ 30.0 + zr2 *
					  (- 66.0 + zr2 *
					   (+ 12.0))))) / s2)
			    / s2 * expzr2;
			  ym = (+ 2.25 - 7.2 / s2) / s / s2 * erfczr
			    + (- 1.5 * (- 3.0 +  zr2 *
					(+ 6.0 + zr2 *
					 (- 12.0 + zr2 *
					  (+ 4.0))))
			       - 1.2 * (+ 12.0 + zr2 *
					(+ 8.0 + zr2 *
					 (- 22.0 + zr2 *
					  (+ 58.0 + zr2 *
					   (- 34.0 + zr2 *
					    (+ 4.0)))))) / s2)
			    / s2 * expzr2;
			  zm = + 1.8 / s2 / s / s2 * erfczr
			    + (- 1.5 * (+ 0.0 +  zr2 *
					(+ 8.0 + zr2 *
					 (- 4.0)))
			       - 1.2 * (- 3.0 + zr2 *
					(- 2.0 + zr2 *
					 (- 26.0 + zr2 *
					  (+ 26.0 + zr2 *
					   (- 4.0))))) / s2)
			    / s2 * expzr2;
	      
			  matrix_fts_atimes (x + i11, y + j11,
					     ex, ey, ez,
					     xa, ya,
					     yb,
					     xc, yc,
					     xg, yg,
					     yh,
					     xm, ym, zm);
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
		  kexp = pivol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  ex = k1 / k;
		  ey = k2 / k;
		  ez = k3 / k;

		  ya = 6.0 * (1.0 - kk / 3.0) * kexp;
		  yb = 3.0 * k * kexp;
		  yc = 3.0 / 2.0 * kk * kexp;
		  yg = 3.0 * (1.0 - 4.0 / 15.0 * kk) * k * kexp;
		  yh = 3.0 / 2.0 * kk * kexp;
		  ym = 3.0 * (1.0 - kk / 5.0) * kk * kexp;
      
		  for (i = 0; i < np; i++)
		    {
		      i11 = i * 11;
		      ix = i * 3;
		      iy = ix + 1;
		      iz = ix + 2;
		      for (j = 0; j < np; j++)
			{
			  j11 = j * 11;
			  jx = j * 3;
			  jy = jx + 1;
			  jz = jx + 2;

			  xx = pos [jx] - pos [ix];
			  yy = pos [jy] - pos [iy];
			  zz = pos [jz] - pos [iz];

			  cf = cos (+ k1 * xx
				    + k2 * yy
				    + k3 * zz);

			  sf = - sin (+ k1 * xx
				      + k2 * yy
				      + k3 * zz);

			  matrix_fts_atimes (x + i11, y + j11,
					     ex, ey, ez,
					     0.0, cf * ya,
					     sf * yb,
					     0.0, cf * yc,
					     0.0, sf * yg,
					     cf * yh,
					     0.0, cf * ym, 0.0);
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
#endif /* ZETA */
}


/** natural mobility problem **/
/* solve natural mobility problem under Ewald sum
 * INPUT
 *  np : # particles
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
calc_mob_ewald_3fts (int np,
		     double *f, double *t, double *e,
		     double *u, double *o, double *s)
{
  int i;
  int n11;

  double *b;
  double *x;


  n11 = np * 11;
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3fts ().\n");
      exit (1);
    }

  calc_b_mob_ewald_3fts (np, f, t, e, b);

  /* first guess */
  for (i = 0; i < n11; ++i)
    x [i] = 0.0;

  solve_iter_stab (n11, b, x, atimes_mob_ewald_3fts,
		   gpb, 2000, -6);

  set_FTS_by_fts (np, u, o, s, x);

  free (b);
  free (x);
}

/* calc b-term (constant term) of (natural) mobility problem under Ewald sum
 * where b := -(0,0,e) + M.(f,t,0).
 * INPUT
 *  np : # particles (not # elements in b[]!)
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] :
 * OUTPUT
 *  b [np * 11] : constant vector
 */
static void
calc_b_mob_ewald_3fts (int np,
		       double *f, double *t, double *e,
		       double *b)
{
  int i;
  int i5, i11;
  int j;
  int n3, n5, n11;

  double *x;
  double *v3_0;
  double *v5_0;


  n3 = np * 3;
  n5 = np * 5;
  n11 = np * 11;

  x = malloc (sizeof (double) * n11);
  v3_0 = malloc (sizeof (double) * n3);
  v5_0 = malloc (sizeof (double) * n5);
  if (x == NULL
      || v3_0 == NULL
      || v5_0 == NULL)
    {
      fprintf (stderr, "allocation error in calc_b_mob_ewald_3fts ().\n");
      exit (1);
    }

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }
  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  set_fts_by_FTS (np, x, f, t, v5_0);
  atimes_ewald_3fts (n11, x, b);

  for (i = 0; i < np; ++i)
    {
      i5 = i * 5;
      i11 = i * 11;
      for (j = 0; j < 5; j ++)
	{
	  b [i11 + 6 + j] -= e [i5 + j];
	}
    }

  free (x);
  free (v3_0);
  free (v5_0);
}

/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := (u,o,0) - M.(0,0,s).
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_ewald_3fts (int n, double *x, double *y)
{
  int i;
  int n3, n5, n11;

  double *z;
  double *v3_0;
  double *v5_0;
  double *u;
  double *o;
  double *s;


  n3 = n * 3;
  n5 = n * 5;
  n11 = n * 11;

  z = malloc (sizeof (double) * n11);
  v3_0 = malloc (sizeof (double) * n3);
  v5_0 = malloc (sizeof (double) * n5);
  u = malloc (sizeof (double) * n3);
  o = malloc (sizeof (double) * n3);
  s = malloc (sizeof (double) * n5);
  if (z == NULL
      || v3_0 == NULL
      || v5_0 == NULL
      || u == NULL
      || o == NULL
      || s == NULL)
    {
      fprintf (stderr, "allocation error in atimes_mob_ewald_3fts ().\n");
      exit (1);
    }


  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }
  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  set_FTS_by_fts (n, u, o, s, x);

  set_fts_by_FTS (n, y, v3_0, v3_0, s);
  atimes_ewald_3fts (n11, y, z);

  set_fts_by_FTS (n, y, u, 0, v5_0);

  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
    }

  free (z);
  free (v3_0);
  free (v5_0);
  free (u);
  free (o);
  free (s);
}
