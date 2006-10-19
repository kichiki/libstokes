/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3f.c,v 4.7 2006/10/19 04:11:28 ichiki Exp $
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
#include "ewald.h" // atimes_ewald_3all()
#include "lub.h" // calc_lub_ewald_3f()

#include "ewald-3f.h"


/** natural resistance problem **/
/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  u [np * 3] :
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_ewald_3f (struct stokes * sys,
		    const double *u,
		    double *f)
{
  int np;
  int i;
  int n3;


  sys->version = 0; // F version
  np = sys->np;
  n3 = np * 3;

  /* first guess */
  for (i = 0; i < n3; ++i)
    {
      f [i] = 0.0;
    }

  solve_iter (n3, u, f,
	      atimes_ewald_3all, (void *) sys,
	      sys->it);
  // for atimes_ewald_3all(), sys->version is 0 (F)
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
solve_mob_ewald_3f (struct stokes * sys,
		    const double *f,
		    double *u)
{
  sys->version = 0; // F version
  atimes_ewald_3all (sys->np * 3, f, u, (void *) sys);
  // sys->version is 0 (F)
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
calc_b_mix_ewald_3f (struct stokes * sys,
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


  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to F.\n");
      sys->version = 0;
    }
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  nm3 = nm * 3;

  x = (double *) malloc (sizeof (double) * n3);
  v3_0 = (double *) malloc (sizeof (double) * n3);

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set x := [(F)_m,(0)_f] */
  set_F_by_f (nm, x, f);
  set_F_by_f (nf, x + nm3, v3_0);
  atimes_ewald_3all (n3, x, b, (void *) sys); // sys->version is 0 (F)

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
atimes_mix_ewald_3f (int n, const double *x, double *y, void * user_data)
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
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to F.\n");
      sys->version = 0;
    }
  np = sys->np;
  nm = sys->nm;

  np3 = np * 3;
  nf = np - nm;
  nf3 = nf * 3;
  nm3 = nm * 3;

  z = (double *) malloc (sizeof (double) * n);
  v3_0 = (double *) malloc (sizeof (double) * np3);
  u = (double *) malloc (sizeof (double) * nm3);
  ff = (double *) malloc (sizeof (double) * nf3);

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
  atimes_ewald_3all (n, y, z, (void *) sys); // sys->version is 0 (F)

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
solve_mix_ewald_3f (struct stokes * sys,
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


  sys->version = 0; // F version
  np = sys->np;
  nm = sys->nm;
  if (np == nm)
    {
      atimes_ewald_3all (np * 3, f, u, (void *) sys); // sys->version is 0 (F)
      return;
    }

  nf = np - nm;
  n3 = np * 3;
  nm3 = nm * 3;

  double *uf0;
  uf0 = (double *) malloc (sizeof (double) * nf * 3);
  b = (double *) malloc (sizeof (double) * n3);
  x = (double *) malloc (sizeof (double) * n3);
  if (uf0 == NULL ||
      b == NULL ||
      x == NULL)
    {
      fprintf (stderr, "libstokes: allocation error"
	       " at solve_mix_ewald_3f()\n");
      exit (1);
    }

  int ix, iy, iz;
  for (i = 0; i < nf; i ++)
    {
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;
      uf0 [ix] = uf[ix] - sys->Ui[0];
      uf0 [iy] = uf[iy] - sys->Ui[1];
      uf0 [iz] = uf[iz] - sys->Ui[2];
    }

  //calc_b_mix_ewald_3f (sys, f, uf, b);
  calc_b_mix_ewald_3f (sys, f, uf0, b);
  free (uf0);

  /* first guess */
  for (i = 0; i < n3; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n3, b, x,
	      atimes_mix_ewald_3f, (void *) sys,
	      sys->it);

  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

  for (i = 0; i < nm; i ++)
    {
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;
      u [ix] += sys->Ui[0];
      u [iy] += sys->Ui[1];
      u [iz] += sys->Ui[2];
    }

  free (b);
  free (x);
}



/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *  u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
void
solve_res_lub_ewald_3f (struct stokes * sys,
			const double *u,
			double *f)
{
  int np;
  int i;
  int n3;

  double *lub;


  sys->version = 0; // F version
  np = sys->np;
  n3 = np * 3;
  lub = (double *) malloc (sizeof (double) * n3);

  calc_lub_ewald_3f (sys, u, lub);
  atimes_ewald_3all (n3, lub, f, (void *) sys);
  // sys->version is 1 (FT)
  // f[] is used temporarily
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
	      atimes_ewald_3all, (void *) sys,
	      sys->it);
  // for atimes_ewald_3all(), sys->version is 0 (F)

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
calc_b_mix_lub_ewald_3f (struct stokes * sys,
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


  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to F.\n");
      sys->version = 0;
    }
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  nm3 = nm * 3;

  x = (double *) malloc (sizeof (double) * n3);
  y = (double *) malloc (sizeof (double) * n3);
  v3_0 = (double *) malloc (sizeof (double) * n3);

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

  atimes_ewald_3all (n3, x, y, (void *) sys); // sys->version is 0 (F)

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
atimes_mix_lub_ewald_3f (int n, const double *x,
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
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to F.\n");
      sys->version = 0;
    }
  np = sys->np;
  nm = sys->nm;

  np3 = np * 3;
  nf = np - nm;
  nf3 = nf * 3;
  nm3 = nm * 3;

  w = (double *) malloc (sizeof (double) * n);
  z = (double *) malloc (sizeof (double) * n);
  v3_0 = (double *) malloc (sizeof (double) * np3);
  u = (double *) malloc (sizeof (double) * nm3);
  ff = (double *) malloc (sizeof (double) * nf3);

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

  atimes_ewald_3all (n, w, z, (void *) sys); // sys->version is 0 (F)

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
solve_mix_lub_ewald_3f (struct stokes * sys,
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


  sys->version = 0; // F version
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  nm3 = nm * 3;

  double *uf0;
  uf0 = (double *) malloc (sizeof (double) * nf * 3);
  b = (double *) malloc (sizeof (double) * n3);
  x = (double *) malloc (sizeof (double) * n3);
  if (uf0 == NULL ||
      b == NULL ||
      x == NULL)
    {
      fprintf (stderr, "libstokes: allocation error"
	       " at solve_mix_lub_ewald_3f()\n");
      exit (1);
    }

  int ix, iy, iz;
  for (i = 0; i < nf; i ++)
    {
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;
      uf0 [ix] = uf[ix] - sys->Ui[0];
      uf0 [iy] = uf[iy] - sys->Ui[1];
      uf0 [iz] = uf[iz] - sys->Ui[2];
    }

  //calc_b_mix_lub_ewald_3f (sys, f, uf, b);
  calc_b_mix_lub_ewald_3f (sys, f, uf0, b);
  free (uf0);

  /* first guess */
  for (i = 0; i < n3; ++i)
    {
      x [i] = 0.0;
    }

  //NUMB_mobile_particles = nm;
  solve_iter (n3, b, x,
	      atimes_mix_lub_ewald_3f, (void *) sys,
	      sys->it);

  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

  for (i = 0; i < nm; i ++)
    {
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;
      u [ix] += sys->Ui[0];
      u [iy] += sys->Ui[1];
      u [iz] += sys->Ui[2];
    }

  free (b);
  free (x);
}
