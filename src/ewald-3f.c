/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3f.c,v 4.3 2006/09/28 04:47:32 kichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 3D configuration
 * periodic boundary condition in 3 direction
 * F version
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
#include "f.h"

#include "ewald-3f.h"


/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for F version
 * INPUT
 *  n := np * 11
 *  x [n * 3] : F
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n * 3] : U
 */
void
atimes_ewald_3f (int n, const double *x, double *y, void * user_data)
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

  double ex, ey, ez;

  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;
  double rlx, rly, rlz;

  int np;
  int i, j;
  int i3, j3;
  int ix, iy, iz;
  int jx, jy, jz;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double k;
  double cf, sf;
  double kexp;

  double erfczr;
  double expzr2;

  double a2;


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

  np = n / 3;
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

  for (i = 0; i < np; i++)
    {
      i3 = i * 3;
      matrix_f_atimes (x + i3, y + i3,
		       0.0, 0.0, 0.0,
		       xa, ya);
    }

  /* for zeta code to measure CPU times */
  cpu0 = ptime_ms_d ();

  /* first Ewald part ( real space ) */
  for (i = 0; i < np; i++)
    {
      i3 = i * 3;
      ix = i3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = 0; j < np; j++)
	{
	  j3 = j * 3;
	  jx = j3;
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
	      
			  matrix_f_atimes (x + i3, y + j3,
					   ex, ey, ez,
					   xa, ya);
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
      
		  for (i = 0; i < np; i++)
		    {
		      i3 = i * 3;
		      ix = i * 3;
		      iy = ix + 1;
		      iz = ix + 2;
		      for (j = 0; j < np; j++)
			{
			  j3 = j * 3;
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

			  matrix_f_atimes (x + i3, y + j3,
					   ex, ey, ez,
					   0.0, cf * ya);
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
/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  u [np * 3] :
 * OUTPUT
 *  f [np * 3] :
 */
void
calc_res_ewald_3f (struct stokes * sys,
		   const double *u,
		   double *f)
{
  int np;
  int i;
  int n3;


  np = sys->np;
  n3 = np * 3;

  /* first guess */
  for (i = 0; i < n3; ++i)
    {
      f [i] = 0.0;
    }

  solve_iter (n3, u, f,
	      atimes_ewald_3f, (void *) sys,
	      sys->it);
}

/** natural mobility problem **/
/* solve natural mobility problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 * OUTPUT
 *  u [np * 3] :
 */
void
calc_mob_ewald_3f (struct stokes * sys,
		   const double *f,
		   double *u)
{
  atimes_ewald_3f (sys->np * 3, f, u, (void *) sys);
}


/** natural mobility problem with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem under Ewald sum
 * where b := - [(0,0)_m,(u,o)_f] + M.[(f,t)_m,(0,0)_f].
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *  b [np * 3] : constant vector
 */
static void
calc_b_mob_fix_ewald_3f (struct stokes * sys,
			 const double *f,
			 const double *uf,
			 double *b)
{
  int np, nm;
  int i;
  int i3;
  int im3;
  int j;
  int nf;
  int n3;
  int nm3;

  double *x;
  double *v3_0;


  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  nm3 = nm * 3;

  x = malloc (sizeof (double) * n3);
  v3_0 = malloc (sizeof (double) * n3);
  if (x == NULL
      || v3_0 == NULL)
    {
      fprintf (stderr, "allocation error in calc_b_mob_ewald_3f ().\n");
      exit (1);
    }

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set x := [(F)_m,(0)_f] */
  set_F_by_f (nm, x, f);
  set_F_by_f (nf, x + nm3, v3_0);
  atimes_ewald_3f (n3, x, b, (void *) sys);

  /* set b := M.x - [(0)_m,(U)_f] */
  for (i = 0; i < nf; ++i)
    {
      i3 = i * 3;
      im3 = (i + nm) * 3;
      for (j = 0; j < 3; ++j)
	{
	  b [im3 + j] -= uf [i3 + j];
	}
    }

  free (x);
  free (v3_0);
}
/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := [(u)_m,(0)_f] - M.[(0)_m,(f)_f].
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_fix_ewald_3f (int n, const double *x, double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int np3;
  int nm, nf;
  int nf3;
  int nm3;

  double *z;
  double *v3_0;
  double *u;
  double *ff;


  sys = (struct stokes *) user_data;
  np = sys->np;
  nm = sys->nm;

  np3 = np * 3;
  nf = np - nm;
  nf3 = nf * 3;
  nm3 = nm * 3;

  z = malloc (sizeof (double) * n);
  v3_0 = malloc (sizeof (double) * np3);
  u = malloc (sizeof (double) * nm3);
  ff = malloc (sizeof (double) * nf3);
  if (z == NULL
      || v3_0 == NULL
      || u == NULL
      || ff == NULL)
    {
      fprintf (stderr, "allocation error in atimes_mob_ewald_3f ().\n");
      exit (1);
    }

  for (i = 0; i < np3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set (U)_mobile,(F)_fixed by x[] */
  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

  /* set y := [(0)_mobile,(F)_fixed] */
  set_F_by_f (nm, y, v3_0);
  set_F_by_f (nf, y + nm3, ff);
  atimes_ewald_3f (n, y, z, (void *) sys);

  /* set y := [(U)_mobile,(0)_fixed] */
  set_F_by_f (nm, y, u);
  set_F_by_f (nf, y + nm3, v3_0);

  /* set y := [(U)_m,(0)_f] - M.[(0)_m,(F)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
    }

  free (z);
  free (v3_0);
  free (u);
  free (ff);
}
/* solve natural mobility problem with fixed particles in F version
 * under Ewald sum
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
calc_mob_fix_ewald_3f (struct stokes * sys,
		       const double *f,
		       const double *uf,
		       double *u,
		       double *ff)
{
  int np, nm;
  int i;
  int n3;
  int nf;
  int nm3;

  double *b;
  double *x;


  np = sys->np;
  nm = sys->nm;
  if (np == nm)
    {
      atimes_ewald_3f (np * 3, f, u, (void *) sys);
      return;
    }

  nf = np - nm;
  n3 = np * 3;
  nm3 = nm * 3;

  b = malloc (sizeof (double) * n3);
  x = malloc (sizeof (double) * n3);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3f ().\n");
      exit (1);
    }

  calc_b_mob_fix_ewald_3f (sys, f, uf, b);

  /* first guess */
  for (i = 0; i < n3; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n3, b, x,
	      atimes_mob_fix_ewald_3f, (void *) sys,
	      sys->it);
  /*
  solve_iter_stab (n3, b, x, atimes_mob_fix_ewald_3f, (void *) sys,
		   sta2,
		   sys->iter_max,
		   sys->iter_log10_eps);
  solve_iter_gmres (n3, b, x, atimes_mob_fix_ewald_3f, (void *) sys,
		    sys->iter_max,
		    20,
		    sys->iter_eps);
  */

  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

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
/* calculate lubrication f by uoe for all particles
 * under the periodic boundary condition
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   u [np * 3] : velocity, angular velocity, strain
 * OUTPUT
 *   f [np * 3] : force, torque, stresslet
 */
static void
calc_lub_ewald_3f (struct stokes * sys,
		   const double * u,
		   double * f)
{
  int np;
  int i, j, k;
  int i3;
  int j3;

  double * tmp_pos;


  np = sys->np;

  tmp_pos = malloc (sizeof (double) * 3);
  if (tmp_pos == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_ewald_3f().\n");
      exit (1);
    }

  /* clear f [np * 3] */
  for (i = 0; i < np * 3; ++i)
    {
      f [i] = 0.0;
    }

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  /* all image cells */
	  for (k = 0; k < 27; ++k)
	    {
	      tmp_pos [0] = sys->pos [j3 + 0] + sys->llx [k];
	      tmp_pos [1] = sys->pos [j3 + 1] + sys->lly [k];
	      tmp_pos [2] = sys->pos [j3 + 2] + sys->llz [k];
	      if (cond_lub (sys->pos + i3, tmp_pos) == 0)
		{
		  calc_lub_f_2b (sys,
				 u + i3, u + j3,
				 sys->pos + i3, tmp_pos,
				 f + i3, f + j3);
		}
	    }
	}
    }

  free (tmp_pos);
}

/* solve natural resistance problem with lubrication
 * in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
calc_res_lub_ewald_3f (struct stokes * sys,
		       const double *u,
		       double *f)
{
  int np;
  int i;
  int n3;

  double *lub;


  np = sys->np;
  n3 = np * 3;
  lub = malloc (sizeof (double) * n3);
  if (lub == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_lub_ewald_3f ().\n");
      exit (1);
    }

  calc_lub_ewald_3f (sys, u, lub);
  atimes_ewald_3f (n3, lub, f, (void *) sys); /* f[] is used temporaly */
  for (i = 0; i < n3; ++i)
    {
      lub [i] = u [i] + f [i];
    }

  /* first guess */
  for (i = 0; i < n3; ++i)
    {
      f [i] = 0.0;
    }

  solve_iter (n3, lub, f,
	      atimes_ewald_3f, (void *) sys,
	      sys->it);
  /*
  solve_iter_stab (n3, lub, f, atimes_ewald_3f, (void *) sys,
		   sta2,
		   sys->iter_max,
		   sys->iter_log10_eps);
  */

  free (lub);
}

/** natural mobility problem with lubrication with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem
 * with lubrication under Ewald sum
 * where b := -(0) + M.(f).
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *  b [np * 3] : constant vector
 */
static void
calc_b_mob_lub_fix_ewald_3f (struct stokes * sys,
			     const double *f,
			     const double *uf,
			     double *b)
{
  int np, nm;
  int i;
  int nf;
  int n3;
  int nm3;

  double *x;
  double *y;
  double *v3_0;


  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  nm3 = nm * 3;

  x = malloc (sizeof (double) * n3);
  y = malloc (sizeof (double) * n3);
  v3_0 = malloc (sizeof (double) * n3);
  if (x == NULL
      || y == NULL
      || v3_0 == NULL)
    {
      fprintf (stderr, "allocation error in calc_b_mob_lub_ewald_3f ().\n");
      exit (1);
    }

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set b := [(0)_m,(U)_f] */
  set_F_by_f (nm, b, v3_0);
  set_F_by_f (nf, b + nm3, uf);

  /* set y := L.[(0)_m,(U)_f] */
  calc_lub_ewald_3f (sys, b, y);

  /* set x := [(F)_m,(0)_f] */
  set_F_by_f (nm, x, f);
  set_F_by_f (nf, x + nm3, v3_0);

  /* set x := x - y */
  for (i = 0; i < n3; ++i)
    {
      x [i] -= y [i];
    }

  atimes_ewald_3f (n3, x, y, (void *) sys);

  /* set b := - (I + M.L).[(0)_m,(U)_f] + M.[(F)_m,(0)_f] */
  for (i = 0; i < n3; ++i)
    {
      b [i] = - b [i] + y [i];
    }

  free (x);
  free (y);
  free (v3_0);
}

/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := [(u,o)_m,(0,0)_f] - M.[(0,0)_m,(f,t)_f].
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_lub_fix_ewald_3f (int n, const double *x,
			     double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int np3;
  int nm, nf;
  int nf3;
  int nm3;

  double *w;
  double *z;
  double *v3_0;
  double *u;
  double *ff;


  sys = (struct stokes *) user_data;
  np = sys->np;
  nm = sys->nm;

  np3 = np * 3;
  nf = np - nm;
  nf3 = nf * 3;
  nm3 = nm * 3;

  w = malloc (sizeof (double) * n);
  z = malloc (sizeof (double) * n);
  v3_0 = malloc (sizeof (double) * np3);
  u = malloc (sizeof (double) * nm3);
  ff = malloc (sizeof (double) * nf3);
  if (z == NULL
      || v3_0 == NULL
      || u == NULL
      || ff == NULL)
    {
      fprintf (stderr, "allocation error in atimes_mob_lub_ewald_3f ().\n");
      exit (1);
    }

  for (i = 0; i < np3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set (U)_mobile,(F)_fixed by x[] */
  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

  /* set y := [(U)_mobile,(0)_fixed] */
  set_F_by_f (nm, y, u);
  set_F_by_f (nf, y + nm3, v3_0);

  /* set w := L.[(U,O,0)_mobile,(0,0,0)_fixed] */
  calc_lub_ewald_3f (sys, y, w);

  /* set z := [(0)_mobile,(F)_fixed] */
  set_F_by_f (nm, z, v3_0);
  set_F_by_f (nf, z + nm3, ff);

  for (i = 0; i < n; ++i)
    {
      w [i] -= z [i];
    }

  atimes_ewald_3f (n, w, z, (void *) sys);

  /* set y := (I + M.L).[(U)_m,(0)_f] - M.[(0)_m,(F)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] += z [i];
    }

  free (w);
  free (z);
  free (v3_0);
  free (u);
  free (ff);
}
/* solve natural mobility problem with lubrication
 * with fixed particles in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
calc_mob_lub_fix_ewald_3f (struct stokes * sys,
			   const double *f,
			   const double *uf,
			   double *u,
			   double *ff)
{
  int np, nm;

  int i;
  int n3;
  int nf;
  int nm3;

  double * b;
  double * x;


  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  nm3 = nm * 3;

  b = malloc (sizeof (double) * n3);
  x = malloc (sizeof (double) * n3);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3f ().\n");
      exit (1);
    }

  calc_b_mob_lub_fix_ewald_3f (sys, f, uf, b);

  /* first guess */
  for (i = 0; i < n3; ++i)
    {
      x [i] = 0.0;
    }

  //NUMB_mobile_particles = nm;
  solve_iter (n3, b, x,
	      atimes_mob_lub_fix_ewald_3f, (void *) sys,
	      sys->it);
  /*
  solve_iter_stab (n3, b, x, atimes_mob_lub_fix_ewald_3f, (void *) sys,
		   //gpb,
		   //gpb_chk,
		   sys->iter_max,
		   sys->iter_log10_eps);
  solve_iter_gmres (n3, b, x, atimes_mob_lub_fix_ewald_3f, (void *) sys,
		    sys->iter_max,
		    20,
		    sys->iter_eps);
  solve_iter_otmk (n3, b, x, atimes_mob_lub_fix_ewald_3f, (void *) sys,
		   otmk,
		   sys->iter_max,
		   sys->iter_log10_eps,
		   10);
  */

  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

  free (b);
  free (x);
}
