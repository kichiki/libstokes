/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3fts.c,v 5.3 2006/09/28 04:49:16 kichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 3D configuration
 * periodic boundary condition in 3 direction
 * FTS version
 * non-dimension formulation
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include <libiter.h> /* solve_iter() */

#include "bench.h" // ptime_ms_d()
#include "stokes.h" /* struct stokeks */
#include "fts.h"

#include "ewald-3fts.h"


/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for FTS version
 * INPUT
 *  n := np * 11
 *  x [n * 11] : FTS
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n * 11] : UOE
 */
void
atimes_ewald_3fts (int n, const double *x, double *y, void * user_data)
{
  struct stokes * sys;
  int pcellx, pcelly, pcellz;
  int kmaxx, kmaxy, kmaxz;
  double zeta, zeta2, zaspi, za2;
  double pi2;
  double pivol;
  double lx, ly, lz;
  double *pos;
  int np_sys; /* for check */

  double cpu0, cpu; /* for ptime_ms_d() */

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


  sys = (struct stokes *) user_data;
  pcellx = sys->pcellx;
  pcelly = sys->pcelly;
  pcellz = sys->pcellz;
  kmaxx  = sys->kmaxx;
  kmaxy  = sys->kmaxy;
  kmaxz  = sys->kmaxz;
  zeta   = sys->zeta;
  zeta2  = sys->zeta2;
  zaspi  = sys->zaspi;
  za2    = sys->za2;
  pi2    = sys->pi2;
  pivol  = sys->pivol;
  lx     = sys->lx;
  ly     = sys->ly;
  lz     = sys->lz;
  pos    = sys->pos;
  np_sys = sys->np;

  np = n / 11;
  if (np_sys != np)
    {
      fprintf (stderr, "wrong n %d <-> np %d\n", n, np_sys);
      exit (1);
    }

  /* clear result */
  for (i = 0; i < n; i ++)
    {
      y [i] = 0.0;
    }

  /* diagonal part ( self part ) */
  xa = ya = 1.0 - zaspi * (6.0 - 40.0 / 3.0 * za2);
  xc = yc = 0.75 - zaspi * za2 * 10.0;
  xm = ym = zm = 0.9 - zaspi * za2 * (12.0 - 30.24 * za2);

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

  /* for zeta code to measure CPU times */
  cpu0 = ptime_ms_d ();

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

  /* for zeta code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu2 = cpu - cpu0;
  cpu0 = cpu;

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

  /* for zeta code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu3 = cpu - cpu0;
  sys->cpu1 = sys->cpu2 + sys->cpu3;
}


/** natural resistance problem **/
/* solve natural resistance problem in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_ewald_3fts (struct stokes * sys,
		     const double *u, const double *o, const double *e,
		     double *f, double *t, double *s)
{
  int np;
  int i;
  int n11;

  double *b;
  double *x;


  np = sys->np;

  n11 = np * 11;
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_3fts ().\n");
      exit (1);
    }

  set_fts_by_FTS (np, b, u, o, e);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_ewald_3fts, (void *) sys,
	      sys->it);
  /*
  solve_iter_stab (n11, b, x, atimes_ewald_3fts, (void *) sys,
		   gpb,
		   sys->iter_max,
		   sys->iter_log10_eps);
  */

  set_FTS_by_fts (np, f, t, s, x);

  free (b);
  free (x);
}

/** natural mobility problem **/
/* calc b-term (constant term) of (natural) mobility problem under Ewald sum
 * where b := -(0,0,e) + M.(f,t,0).
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 *  t [np * 3] :
 *  e [np * 5] :
 * OUTPUT
 *  b [np * 11] : constant vector
 */
static void
calc_b_mob_ewald_3fts (struct stokes * sys,
		       const double *f, const double *t, const double *e,
		       double *b)
{
  int np;
  int i;
  int i5, i11;
  int j;
  int n3, n5, n11;

  double *x;
  double *v5_0;


  np = sys->np;

  n3 = np * 3;
  n5 = np * 5;
  n11 = np * 11;

  x = malloc (sizeof (double) * n11);
  v5_0 = malloc (sizeof (double) * n5);
  if (x == NULL
      || v5_0 == NULL)
    {
      fprintf (stderr, "allocation error in calc_b_mob_ewald_3fts ().\n");
      exit (1);
    }

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  set_fts_by_FTS (np, x, f, t, v5_0);
  atimes_ewald_3fts (n11, x, b, (void *) sys);

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
  free (v5_0);
}
/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := (u,o,0) - M.(0,0,s).
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_ewald_3fts (int n, const double *x, double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int n3, n5;

  double *z;
  double *v5_0;
  double *u;
  double *o;
  double *s;


  sys = (struct stokes *) user_data;
  np = sys->np;

  n3 = np * 3;
  n5 = np * 5;

  z = malloc (sizeof (double) * n);
  v5_0 = malloc (sizeof (double) * n5);
  u = malloc (sizeof (double) * n3);
  o = malloc (sizeof (double) * n3);
  s = malloc (sizeof (double) * n5);
  if (z == NULL
      || v5_0 == NULL
      || u == NULL
      || o == NULL
      || s == NULL)
    {
      fprintf (stderr, "allocation error in atimes_mob_ewald_3fts ().\n");
      exit (1);
    }


  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  set_FTS_by_fts (np, u, o, s, x);

  set_fts_by_FTS (np, y, v5_0, v5_0, s);
  atimes_ewald_3fts (n, y, z, (void *) sys);

  set_fts_by_FTS (np, y, u, o, v5_0);

  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
    }

  free (z);
  free (v5_0);
  free (u);
  free (o);
  free (s);
}
/* solve natural mobility problem in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
calc_mob_ewald_3fts (struct stokes * sys,
		     const double *f, const double *t, const double *e,
		     double *u, double *o, double *s)
{
  int np;
  int i;
  int n11;

  double *b;
  double *x;


  np = sys->np;
  n11 = np * 11;
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3fts ().\n");
      exit (1);
    }

  calc_b_mob_ewald_3fts (sys, f, t, e, b);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_mob_ewald_3fts, (void *) sys,
	      sys->it);
  /*
  solve_iter_stab (n11, b, x, atimes_mob_ewald_3fts, (void *) sys,
		   gpb,
		   sys->iter_max,
		   sys->iter_log10_eps);
  */

  set_FTS_by_fts (np, u, o, s, x);

  free (b);
  free (x);
}

/** natural mobility problem with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem under Ewald sum
 * where b := -(0,0,e) + M.(f,t,0).
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [nm * 5] :
 *  uf [nf * 3] :
 *  of [nf * 3] :
 *  ef [nf * 5] :
 * OUTPUT
 *  b [np * 11] : constant vector
 */
static void
calc_b_mob_fix_ewald_3fts (struct stokes * sys,
			   const double *f, const double *t, const double *e,
			   const double *uf, const double *of,
			   const double *ef,
			   double *b)
{
  int np, nm;
  int i;
  int i3, i5, i11;
  int j;
  int nf;
  int n3, n5, n11;
  int nm11;

  double *x;
  double *v5_0;


  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  n5 = np * 5;
  n11 = np * 11;
  nm11 = nm * 11;

  x = malloc (sizeof (double) * n11);
  v5_0 = malloc (sizeof (double) * n5);
  if (x == NULL
      || v5_0 == NULL)
    {
      fprintf (stderr, "allocation error in calc_b_mob_ewald_3fts ().\n");
      exit (1);
    }

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  /* set x := [(F,T,0)_m,(0,0,0)_f] */
  set_fts_by_FTS (nm, x, f, t, v5_0);
  set_fts_by_FTS (nf, x + nm11, v5_0, v5_0, v5_0);
  atimes_ewald_3fts (n11, x, b, (void *) sys);

  /* set b := M.x - [(0,0,E)_m,(U,O,E)_f] */
  for (i = 0; i < nm; ++i)
    {
      i5 = i * 5;
      i11 = i * 11;
      for (j = 0; j < 5; ++j)
	{
	  b [i11 + 6 + j] -= e [i5 + j];
	}
    }
  for (i = 0; i < nf; ++i)
    {
      i3 = i * 3;
      i5 = i * 5;
      i11 = (i + nm) * 11;
      for (j = 0; j < 3; ++j)
	{
	  b [i11 + j] -= uf [i3 + j];
	  b [i11 + 3 + j] -= of [i3 + j];
	}
      for (j = 0; j < 5; ++j)
	{
	  b [i11 + 6 + j] -= ef [i5 + j];
	}
    }

  free (x);
  free (v5_0);
}
/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := [(u,o,0)_m,(0,0,0)_f] - M.[(0,0,s)_m,(f,t,s)_f].
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_fix_ewald_3fts (int n, const double *x, double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int n5;
  int nm, nf;
  int nf3, nf5;
  int nm3, nm5, nm11;

  double *z;
  double *v5_0;
  double *u;
  double *o;
  double *s;
  double *ff;
  double *tf;
  double *sf;


  sys = (struct stokes *) user_data;
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n5 = np * 5;
  nf3 = nf * 3;
  nf5 = nf * 5;
  nm3 = nm * 3;
  nm5 = nm * 5;
  nm11 = nm * 11;

  z = malloc (sizeof (double) * n);
  v5_0 = malloc (sizeof (double) * n5);
  u = malloc (sizeof (double) * nm3);
  o = malloc (sizeof (double) * nm3);
  s = malloc (sizeof (double) * nm5);
  ff = malloc (sizeof (double) * nf3);
  tf = malloc (sizeof (double) * nf3);
  sf = malloc (sizeof (double) * nf5);
  if (z == NULL
      || v5_0 == NULL
      || u == NULL
      || o == NULL
      || s == NULL
      || ff == NULL
      || tf == NULL
      || sf == NULL)
    {
      fprintf (stderr, "allocation error in atimes_mob_ewald_3fts ().\n");
      exit (1);
    }

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  /* set (U,O,S)_mobile,(F,T,S)_fixed by x[] */
  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  /* set y := [(0,0,S)_mobile,(F,T,S)_fixed] */
  set_fts_by_FTS (nm, y, v5_0, v5_0, s);
  set_fts_by_FTS (nf, y + nm11, ff, tf, sf);
  atimes_ewald_3fts (n, y, z, (void *) sys);

  /* set y := [(U,O,0)_mobile,(0,0,0)_fixed] */
  set_fts_by_FTS (nm, y, u, o, v5_0);
  set_fts_by_FTS (nf, y + nm11, v5_0, v5_0, v5_0);

  /* set y := [(U,O,0)_m,(0,0,0)_f] - M.[(0,0,S)_m,(F,T,S)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
    }

  free (z);
  free (v5_0);
  free (u);
  free (o);
  free (s);
  free (ff);
  free (tf);
  free (sf);
}
/* solve natural mobility problem with fixed particles in FTS version
 * under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   e [nm * 5] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 *   ef [nf * 5] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
calc_mob_fix_ewald_3fts (struct stokes * sys,
			 const double *f, const double *t, const double *e,
			 const double *uf, const double *of, const double *ef,
			 double *u, double *o, double *s,
			 double *ff, double *tf, double *sf)
{
  int np, nm;
  int i;
  int n11;
  int nf;
  int nm11;

  double *b;
  double *x;


  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n11 = np * 11;
  nm11 = nm * 11;

  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3fts ().\n");
      exit (1);
    }

  calc_b_mob_fix_ewald_3fts (sys, f, t, e, uf, of, ef, b);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_mob_fix_ewald_3fts, (void *) sys,
	      sys->it);
  /*
  solve_iter_stab (n11, b, x, atimes_mob_fix_ewald_3fts, (void *) sys,
		   gpb,
		   sys->iter_max,
		   sys->iter_log10_eps);
  */

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  free (b);
  free (x);
}

/** natural resistance problem with lubrication **/
/* condition for lubrication
 * INPUT
 *  x1 [3], x2 [3] : position
 * OUTPUT (return value)
 *  0 : r != 0 and r < 3.0
 *  1 : otherwise
 */
static int
cond_lub (const double * x1, const double * x2)
{
  double x, y, z;
  double r2;


  x = x1 [0] - x2 [0];
  y = x1 [1] - x2 [1];
  z = x1 [2] - x2 [2];

  r2 = x * x
    + y * y
    + z * z;

  if (r2 != 0.0
      && r2 < 9.0) // r = 3.0 is the critical separation for lubrication now.
    {
      return 0;
    }
  else
    {
      return 1;
    }
}
/* calculate lubrication fts by uoe for all particles
 * under the periodic boundary condition
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   uoe [np * 11] : velocity, angular velocity, strain
 * OUTPUT
 *   fts [np * 11] : force, torque, stresslet
 */
static void
calc_lub_ewald_3fts (struct stokes * sys,
		     const double * uoe, double * fts)
{
  int np;
  int i, j, k;
  int i3, i11;
  int j3, j11;

  double * tmp_pos;


  np = sys->np;

  tmp_pos = malloc (sizeof (double) * 3);
  if (tmp_pos == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_ewald_3fts().\n");
      exit (1);
    }

  /* clear fts [np * 11] */
  for (i = 0; i < np * 11; ++i)
    {
      fts [i] = 0.0;
    }

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i11 = i * 11;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  j11 = j * 11;
	  /* all image cells */
	  for (k = 0; k < 27; ++k)
	    {
	      tmp_pos [0] = sys->pos [j3 + 0] + sys->llx [k];
	      tmp_pos [1] = sys->pos [j3 + 1] + sys->lly [k];
	      tmp_pos [2] = sys->pos [j3 + 2] + sys->llz [k];
	      if (cond_lub (sys->pos + i3, tmp_pos) == 0)
		{
		  calc_lub_fts_2b (sys,
				   uoe + i11, uoe + j11,
				   sys->pos + i3, tmp_pos,
				   fts + i11, fts + j11);
		}
	    }
	}
    }

  free (tmp_pos);
}
/* solve natural resistance problem with lubrication
 * in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_lub_ewald_3fts (struct stokes * sys,
			 const double *u, const double *o, const double *e,
			 double *f, double *t, double *s)
{
  int np;
  int i;
  int n11;

  double *b;
  double *x;
  double *lub;


  np = sys->np;
  n11 = np * 11;
  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  lub = malloc (sizeof (double) * n11);
  if (b == NULL
      || x == NULL
      || lub == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_lub_ewald_3fts ().\n");
      exit (1);
    }

  set_fts_by_FTS (np, b, u, o, e);
  calc_lub_ewald_3fts (sys, b, lub);
  atimes_ewald_3fts (n11, lub, x, (void *) sys); /* x[] is used temporaly */
  for (i = 0; i < n11; ++i)
    {
      b [i] += x [i];
    }

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_ewald_3fts, (void *) sys,
	      sys->it);
  /*
  solve_iter_stab (n11, b, x, atimes_ewald_3fts, (void *) sys,
		   gpb,
		   sys->iter_max,
		   sys->iter_log10_eps);
  */

  set_FTS_by_fts (np, f, t, s, x);

  free (b);
  free (x);
  free (lub);
}

/** natural mobility problem with lubrication with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem
 * with lubrication under Ewald sum
 * where b := -(0,0,e) + M.(f,t,0).
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [nm * 5] :
 *  uf [nf * 3] :
 *  of [nf * 3] :
 *  ef [nf * 5] :
 * OUTPUT
 *  b [np * 11] : constant vector
 */
static void
calc_b_mob_lub_fix_ewald_3fts (struct stokes * sys,
			       const double *f, const double *t,
			       const double *e,
			       const double *uf, const double *of,
			       const double *ef,
			       double *b)
{
  int np, nm;
  int i;
  int nf;
  int n3, n5, n11;
  int nm11;

  double *x;
  double *y;
  double *v5_0;


  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  n5 = np * 5;
  n11 = np * 11;
  nm11 = nm * 11;

  x = malloc (sizeof (double) * n11);
  y = malloc (sizeof (double) * n11);
  v5_0 = malloc (sizeof (double) * n5);
  if (x == NULL
      || y == NULL
      || v5_0 == NULL)
    {
      fprintf (stderr, "allocation error in calc_b_mob_lub_ewald_3fts ().\n");
      exit (1);
    }

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  /* set b := [(0,0,E)_m,(U,O,E)_f] */
  set_fts_by_FTS (nm, b, v5_0, v5_0, e);
  set_fts_by_FTS (nf, b + nm11, uf, of, ef);

  /* set y := L.[(0,0,E)_m,(U,O,E)_f] */
  calc_lub_ewald_3fts (sys, b, y);

  /* set x := [(F,T,0)_m,(0,0,0)_f] */
  set_fts_by_FTS (nm, x, f, t, v5_0);
  set_fts_by_FTS (nf, x + nm11, v5_0, v5_0, v5_0);

  /* set x := x - y */
  for (i = 0; i < n11; ++i)
    {
      x [i] -= y [i];
    }

  atimes_ewald_3fts (n11, x, y, (void *) sys);

  /* set b := - (I + M.L).[(0,0,E)_m,(U,O,E)_f] + M.[(F,T,0)_m,(0,0,0)_f] */
  for (i = 0; i < n11; ++i)
    {
      b [i] = - b[i] + y [i];
    }

  free (x);
  free (y);
  free (v5_0);
}
/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := [(u,o,0)_m,(0,0,0)_f] - M.[(0,0,s)_m,(f,t,s)_f].
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_lub_fix_ewald_3fts (int n, const double *x,
			       double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int n5;
  int nm, nf;
  int nf3, nf5;
  int nm3, nm5, nm11;

  double *w;
  double *z;
  double *v5_0;
  double *u;
  double *o;
  double *s;
  double *ff;
  double *tf;
  double *sf;


  sys = (struct stokes *) user_data;
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n5 = np * 5;
  nf3 = nf * 3;
  nf5 = nf * 5;
  nm3 = nm * 3;
  nm5 = nm * 5;
  nm11 = nm * 11;

  w = malloc (sizeof (double) * n);
  z = malloc (sizeof (double) * n);
  v5_0 = malloc (sizeof (double) * n5);
  u = malloc (sizeof (double) * nm3);
  o = malloc (sizeof (double) * nm3);
  s = malloc (sizeof (double) * nm5);
  ff = malloc (sizeof (double) * nf3);
  tf = malloc (sizeof (double) * nf3);
  sf = malloc (sizeof (double) * nf5);
  if (z == NULL
      || v5_0 == NULL
      || u == NULL
      || o == NULL
      || s == NULL
      || ff == NULL
      || tf == NULL
      || sf == NULL)
    {
      fprintf (stderr, "allocation error in atimes_mob_lub_ewald_3fts ().\n");
      exit (1);
    }

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  /* set (U,O,S)_mobile,(F,T,S)_fixed by x[] */
  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  /* set y := [(U,O,0)_mobile,(0,0,0)_fixed] */
  set_fts_by_FTS (nm, y, u, o, v5_0);
  set_fts_by_FTS (nf, y + nm11, v5_0, v5_0, v5_0);

  /* set w := L.[(U,O,0)_mobile,(0,0,0)_fixed] */
  calc_lub_ewald_3fts (sys, y, w);

  /* set z := [(0,0,S)_mobile,(F,T,S)_fixed] */
  set_fts_by_FTS (nm, z, v5_0, v5_0, s);
  set_fts_by_FTS (nf, z + nm11, ff, tf, sf);

  for (i = 0; i < n; ++i)
    {
      w [i] -= z [i];
    }

  atimes_ewald_3fts (n, w, z, (void *) sys);

  /* set y := (I + M.L).[(U,O,0)_m,(0,0,0)_f] - M.[(0,0,S)_m,(F,T,S)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] += z [i];
    }

  free (w);
  free (z);
  free (v5_0);
  free (u);
  free (o);
  free (s);
  free (ff);
  free (tf);
  free (sf);
}
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   e [nm * 5] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 *   ef [nf * 5] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
calc_mob_lub_fix_ewald_3fts (struct stokes * sys,
			     const double *f, const double *t, const double *e,
			     const double *uf, const double *of,
			     const double *ef,
			     double *u, double *o, double *s,
			     double *ff, double *tf, double *sf)
{
  int np, nm;

  int i;
  int n11;
  int nf;
  int nm11;

  double *b;
  double *x;


  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n11 = np * 11;
  nm11 = nm * 11;

  b = malloc (sizeof (double) * n11);
  x = malloc (sizeof (double) * n11);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3fts ().\n");
      exit (1);
    }

  calc_b_mob_lub_fix_ewald_3fts (sys, f, t, e, uf, of, ef, b);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_mob_lub_fix_ewald_3fts, (void *) sys,
	      sys->it);
  /*
  solve_iter_stab (n11, b, x, atimes_mob_lub_fix_ewald_3fts, (void *) sys,
		   //gpb_chk,
		   //st2_chk,
		   sta2,
		   sys->iter_max,
		   sys->iter_log10_eps);
  solve_iter_gmres (n11, b, x, atimes_mob_lub_fix_ewald_3fts, (void *) sys,
		    sys->iter_max,
		    20,
		    sys->iter_eps);
  solve_iter_otmk (n11, b, x, atimes_mob_lub_fix_ewald_3fts, (void *) sys,
		   otmk,
		   sys->iter_max,
		   sys->iter_log10_eps,
		   10);
  */

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  free (b);
  free (x);
}
