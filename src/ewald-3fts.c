/* Solvers for 3 dimensional FTS version problems
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3fts.c,v 5.13 2007/03/18 23:48:37 kichiki Exp $
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
#include "ft.h"
#include "f.h"
#include "ewald.h" // atimes_3all()
#include "lub.h" // calc_lub_3fts()

#include "memory-check.h" // macro CHECK_MALLOC

#include "ewald-3fts.h"


/** natural resistance problem **/
/* solve natural resistance problem in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_res_3fts (struct stokes * sys,
		const double *u, const double *o, const double *e,
		double *f, double *t, double *s)
{
  int np;
  int i;
  int n11;


  sys->version = 2; // FTS version
  np = sys->np;
  n11 = np * 11;

  double *u0 = NULL;
  u0 = (double *) malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (u0, "solve_res_3fts");
  double *o0 = NULL;
  o0 = (double *) malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (o0, "solve_res_3fts");
  double *e0 = NULL;
  e0 = (double *) malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (e0, "solve_res_3fts");
  double *b = NULL;
  b = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (b, "solve_res_3fts");
  double *x = NULL;
  x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (x, "solve_res_3fts");


  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  shift_labo_to_rest_E (sys, np, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  set_fts_by_FTS (np, b, u0, o0, e0);
  free (u0);
  free (o0);
  free (e0);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_3all, (void *) sys,
	      sys->it);
  // for atimes_3all(), sys->version is 2 (FTS)

  set_FTS_by_fts (np, f, t, s, x);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}

/** natural mobility problem **/
/* calc b-term (constant term) of (natural) mobility problem
 * for both periodic and non-periodic boundary conditions
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
calc_b_mob_3fts (struct stokes * sys,
		 const double *f, const double *t, const double *e,
		 double *b)
{
  int np;
  int i;
  int i5, i11;
  int j;
  int n3, n5, n11;


  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }
  np = sys->np;

  n3 = np * 3;
  n5 = np * 5;
  n11 = np * 11;

  double *x = NULL;
  x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (x, "calc_b_mob_3fts");
  double *v5_0 = NULL;
  v5_0 = (double *) malloc (sizeof (double) * n5);
  CHECK_MALLOC (v5_0, "calc_b_mob_3fts");

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  set_fts_by_FTS (np, x, f, t, v5_0);
  atimes_3all (n11, x, b, (void *) sys); // sys->version is 2 (FTS)

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
/* calc atimes of (natural) mobility problem
 * for both periodic and non-periodic boundary conditions
 * where A.x := (u,o,0) - M.(0,0,s).
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_3fts (int n, const double *x, double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int n3, n5;


  sys = (struct stokes *) user_data;
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }
  np = sys->np;

  n3 = np * 3;
  n5 = np * 5;

  double *z = NULL;
  z = (double *) malloc (sizeof (double) * n);
  CHECK_MALLOC (z, "atimes_mob_3fts");
  double *v5_0 = NULL;
  v5_0 = (double *) malloc (sizeof (double) * n5);
  CHECK_MALLOC (v5_0, "atimes_mob_3fts");
  double *u = NULL;
  u = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (u, "atimes_mob_3fts");
  double *o = NULL;
  o = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (o, "atimes_mob_3fts");
  double *s = NULL;
  s = (double *) malloc (sizeof (double) * n5);
  CHECK_MALLOC (s, "atimes_mob_3fts");

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  set_FTS_by_fts (np, u, o, s, x);

  set_fts_by_FTS (np, y, v5_0, v5_0, s);
  atimes_3all (n, y, z, (void *) sys); // sys->version is 2 (FTS)

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
/* solve natural mobility problem in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mob_3fts (struct stokes * sys,
		const double *f, const double *t, const double *e,
		double *u, double *o, double *s)
{
  int np;
  int i;
  int n11;


  sys->version = 2; // FTS version
  np = sys->np;
  n11 = np * 11;

  double *e0 = NULL;
  e0 = (double *) malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (e0, "solve_mob_3fts");
  double *b = NULL;
  b = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (b, "solve_mob_3fts");
  double *x = NULL;
  x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (x, "solve_mob_3fts");


  shift_labo_to_rest_E (sys, np, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mob_3fts (sys, f, t, e0, b);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_mob_3fts, (void *) sys,
	      sys->it);

  set_FTS_by_fts (np, u, o, s, x);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
  shift_rest_to_labo_O (sys, sys->np, o);
}

/** natural mobility problem with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem
 * for both periodic and non-periodic boundary conditions
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
calc_b_mix_3fts (struct stokes * sys,
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


  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  n5 = np * 5;
  n11 = np * 11;
  nm11 = nm * 11;

  double *x = NULL;
  x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (x, "calc_b_mix_3fts");
  double *v5_0 = NULL;
  v5_0 = (double *) malloc (sizeof (double) * n5);
  CHECK_MALLOC (v5_0, "calc_b_mix_3fts");

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  /* set x := [(F,T,0)_m,(0,0,0)_f] */
  set_fts_by_FTS (nm, x, f, t, v5_0);
  set_fts_by_FTS (nf, x + nm11, v5_0, v5_0, v5_0);
  atimes_3all (n11, x, b, (void *) sys); // sys->version is 2 (FTS)

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
/* calc atimes of (natural) mobility problem
 * for both periodic and non-periodic boundary conditions
 * where A.x := [(u,o,0)_m,(0,0,0)_f] - M.[(0,0,s)_m,(f,t,s)_f].
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mix_3fts (int n, const double *x, double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int n5;
  int nm, nf;
  int nf3, nf5;
  int nm3, nm5, nm11;


  sys = (struct stokes *) user_data;
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n5 = np * 5;
  nf3 = nf * 3;
  nf5 = nf * 5;
  nm3 = nm * 3;
  nm5 = nm * 5;
  nm11 = nm * 11;

  double *z = NULL;
  z = (double *) malloc (sizeof (double) * n);
  CHECK_MALLOC (z, "atimes_mix_3fts");
  double *v5_0 = NULL;
  v5_0 = (double *) malloc (sizeof (double) * n5);
  CHECK_MALLOC (v5_0, "atimes_mix_3fts");
  double *u = NULL;
  u = (double *) malloc (sizeof (double) * nm3);
  CHECK_MALLOC (u, "atimes_mix_3fts");
  double *o = NULL;
  o = (double *) malloc (sizeof (double) * nm3);
  CHECK_MALLOC (o, "atimes_mix_3fts");
  double *s = NULL;
  s = (double *) malloc (sizeof (double) * nm5);
  CHECK_MALLOC (s, "atimes_mix_3fts");
  double *ff = NULL;
  ff = (double *) malloc (sizeof (double) * nf3);
  CHECK_MALLOC (ff, "atimes_mix_3fts");
  double *tf = NULL;
  tf = (double *) malloc (sizeof (double) * nf3);
  CHECK_MALLOC (tf, "atimes_mix_3fts");
  double *sf = NULL;
  sf = (double *) malloc (sizeof (double) * nf5);
  CHECK_MALLOC (sf, "atimes_mix_3fts");

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
  atimes_3all (n, y, z, (void *) sys); // sys->version is 2 (FTS)

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
 * for both periodic and non-periodic boundary conditions
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
solve_mix_3fts (struct stokes * sys,
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


  sys->version = 2; // FTS version
  np = sys->np;
  nm = sys->nm;
  if (np == nm)
    {
      solve_mob_3fts (sys,
		      f, t, e,
		      u, o, s);
      return;
    }

  nf = np - nm;
  n11 = np * 11;
  nm11 = nm * 11;

  double *uf0 = NULL;
  uf0 = (double *) malloc (sizeof (double) * nf * 3);
  CHECK_MALLOC (uf0, "solve_mix_3fts");
  double *of0 = NULL;
  of0 = (double *) malloc (sizeof (double) * nf * 3);
  CHECK_MALLOC (of0, "solve_mix_3fts");
  double *ef0 = NULL;
  ef0 = (double *) malloc (sizeof (double) * nf * 5);
  CHECK_MALLOC (ef0, "solve_mix_3fts");
  double *e0 = NULL;
  e0  = (double *) malloc (sizeof (double) * nm * 5);
  CHECK_MALLOC (e0, "solve_mix_3fts");
  double *b = NULL;
  b = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (b, "solve_mix_3fts");
  double *x = NULL;
  x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (x, "solve_mix_3fts");


  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  shift_labo_to_rest_E (sys, nf, ef, ef0);
  shift_labo_to_rest_E (sys, nm, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mix_3fts (sys, f, t, e0, uf0, of0, ef0, b);
  free (uf0);
  free (of0);
  free (ef0);
  free (e0);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_mix_3fts, (void *) sys,
	      sys->it);

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_res_lub_3fts (struct stokes * sys,
		    const double *u, const double *o, const double *e,
		    double *f, double *t, double *s)
{
  int np;
  int i;
  int n11;


  sys->version = 2; // FTS version
  np = sys->np;
  n11 = np * 11;

  double *u0 = NULL;
  u0 = (double *) malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (u0, "solve_res_lub_3fts");
  double *o0 = NULL;
  o0 = (double *) malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (o0, "solve_res_lub_3fts");
  double *e0 = NULL;
  e0 = (double *) malloc (sizeof (double) * np * 5);
  CHECK_MALLOC (e0, "solve_res_lub_3fts");
  double *b = NULL;
  b = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (b, "solve_res_lub_3fts");
  double *x = NULL;
  x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (x, "solve_res_lub_3fts");
  double *lub = NULL;
  lub = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (lub, "solve_res_lub_3fts");


  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  shift_labo_to_rest_E (sys, np, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  set_fts_by_FTS (np, b, u0, o0, e0);
  free (u0);
  free (o0);
  free (e0);

  calc_lub_3fts (sys, b, lub);
  atimes_3all (n11, lub, x, (void *) sys);
  // sys->version is 2 (FTS)
  // x[] is used temporarily
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
	      atimes_3all, (void *) sys,
	      sys->it);
  // for atimes_3all(), sys->version is 2 (FTS)

  set_FTS_by_fts (np, f, t, s, x);

  free (b);
  free (x);
  free (lub);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}


/** mob_lub_3fts **/
/* b := M.(f,t,0)-(I+M.L).(0,0,e)
 */
static void
calc_b_mob_lub_3fts (struct stokes * sys,
		     const double *f, const double *t, const double *e,
		     double *b)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }

  int n5  = sys->np * 5;
  int n11 = sys->np * 11;

  double *v5_0 = NULL;
  v5_0 = (double *)malloc (sizeof (double) * n5);
  CHECK_MALLOC (v5_0, "calc_b_mob_lub_3ft");
  int i;
  for (i = 0; i < n5; i ++)
    {
      v5_0 [i] = 0.0;
    }

  double *ft = NULL;
  ft = (double *)malloc (sizeof (double) * n11);
  CHECK_MALLOC (ft, "atimes_mob_lub_3ft");
  double *tmp = NULL;
  tmp = (double *)malloc (sizeof (double) * n11);
  CHECK_MALLOC (tmp, "atimes_mob_lub_3ft");

  // set ft = (f,t,0)
  set_fts_by_FTS (sys->np, ft, f, t, v5_0);
  // set b = M.(f,t,0)
  atimes_3all (n11, ft, b, (void *) sys); // sys->version is 2 (FTS)

  // set ft = (0,0,e)
  set_fts_by_FTS (sys->np, ft, v5_0, v5_0, e);
  // at this point, b [i] = M.(f,t,0) - I.(0.0,e)
  for (i = 0; i < n11; i ++)
    {
      b [i] -= ft [i];
    }

  // set tmp = L.(0,0,e)
  calc_lub_3fts (sys, ft, tmp);
  // set ft = M.L.(0,0,e)
  atimes_3all (n11, tmp, ft, (void *) sys);

  // at this point, b [i] = M.(f,t,0) - (I + ML).(0.0,e)
  for (i = 0; i < n11; i ++)
    {
      b [i] -= ft [i];
    }

  free (v5_0);
  free (ft);
  free (tmp);
}

/* A.x := -M.(0,0,s)+(I+M.L).(u,o,0)
 * where x = (u,o,s)
 */
static void
atimes_mob_lub_3fts (int n, const double *x, double *y, void * user_data)
{
  struct stokes *sys;
  sys = (struct stokes *) user_data;
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }

  int n3  = sys->np * 3;
  int n5  = sys->np * 5;
  int n11 = sys->np * 11;
  if (n != n11)
    {
      fprintf (stderr, "n(%d) != np(%d) * 11\n", n, sys->np);
      exit (1);
    }

  double *v5_0 = NULL;
  v5_0 = (double *)malloc (sizeof (double) * n5);
  CHECK_MALLOC (v5_0, "atimes_mob_lub_3ft");
  int i;
  for (i = 0; i < n5; i ++)
    {
      v5_0 [i] = 0.0;
    }

  double *uo = NULL;
  uo = (double *)malloc (sizeof (double) * n11);
  CHECK_MALLOC (uo, "atimes_mob_lub_3ft");
  double *z = NULL;
  z = (double *)malloc (sizeof (double) * n11);
  CHECK_MALLOC (z, "atimes_mob_lub_3ft");
  double *tmp = NULL;
  tmp = (double *)malloc (sizeof (double) * n11);
  CHECK_MALLOC (tmp, "atimes_mob_lub_3ft");
  double *u;
  double *o;
  double *s;
  u = tmp;
  o = tmp + n3;
  s = tmp + n3 + n3;

  // set (u,o,s) by x
  set_FTS_by_fts (sys->np, u, o, s, x);

  // set uo = (u,o,0)
  set_fts_by_FTS (sys->np, uo, u, o, v5_0);
  // at this point, y [i] = I.(u.o,0)
  for (i = 0; i < n11; i ++)
    {
      y [i] = uo [i];
    }

  // set z = (0,0,s)
  set_fts_by_FTS (sys->np, z, v5_0, v5_0, s); // v5_0 is used as v3_0
  // set tmp = M.(0,0,s)
  atimes_3all (n11, z, tmp, (void *) sys); // sys->version is 2 (FTS)

  // at this point, y [i] = -M.(0,0,s) + I.(u.o,0)
  for (i = 0; i < n11; i ++)
    {
      y [i] -= tmp [i];
    }

  // set z = L.(u,o,0)
  calc_lub_3fts (sys, uo, z);
  // set tmp = M.L.(u,o,0)
  atimes_3all (n11, z, tmp, (void *) sys);

  // finally, b [i] = -M.(0,0,s) + (I + M.L).(u.o,0)
  for (i = 0; i < n11; i ++)
    {
      y [i] += tmp [i];
    }

  free (v5_0);
  free (uo);
  free (z);
  free (tmp);
}


/* solve natural mobility problem with lubrication in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mob_lub_3fts (struct stokes * sys,
		    const double *f, const double *t, const double *e,
		    double *u, double *o, double *s)
{
  sys->version = 2; // FTS version

  int np;
  np = sys->np;
  int n11;
  n11 = np * 11;
  int n5;
  n5 = np * 5;

  double *b = NULL;
  b = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (b, "solve_mob_lub_3fts");
  double *x = NULL;
  x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (x, "solve_mob_lub_3fts");

  double *e0 = NULL;
  e0 = (double *) malloc (sizeof (double) * n5);
  CHECK_MALLOC (e0, "solve_mob_lub_3fts");
  shift_labo_to_rest_E (sys, np, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mob_lub_3fts (sys, f, t, e0, b);
  free (e0);

  /* first guess */
  int i;
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_mob_lub_3fts, (void *) sys,
	      sys->it);

  set_FTS_by_fts (np, u, o, s, x);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
  shift_rest_to_labo_O (sys, sys->np, o);
}


/** natural mobility problem with lubrication with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem with lubrication
 * for both periodic and non-periodic boundary conditions
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
calc_b_mix_lub_3fts (struct stokes * sys,
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


  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  n5 = np * 5;
  n11 = np * 11;
  nm11 = nm * 11;

  double *x = NULL;
  x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (x, "calc_b_mix_lub_3fts");
  double *y = NULL;
  y = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (y, "calc_b_mix_lub_3fts");
  double *v5_0 = NULL;
  v5_0 = (double *) malloc (sizeof (double) * n5);
  CHECK_MALLOC (v5_0, "calc_b_mix_lub_3fts");

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  /* set b := [(0,0,E)_m,(U,O,E)_f] */
  set_fts_by_FTS (nm, b, v5_0, v5_0, e);
  set_fts_by_FTS (nf, b + nm11, uf, of, ef);

  /* set y := L.[(0,0,E)_m,(U,O,E)_f] */
  calc_lub_3fts (sys, b, y);

  /* set x := [(F,T,0)_m,(0,0,0)_f] */
  set_fts_by_FTS (nm, x, f, t, v5_0);
  set_fts_by_FTS (nf, x + nm11, v5_0, v5_0, v5_0);

  /* set x := x - y */
  for (i = 0; i < n11; ++i)
    {
      x [i] -= y [i];
    }

  atimes_3all (n11, x, y, (void *) sys); // sys->version is 2 (FTS)

  /* set b := - (I + M.L).[(0,0,E)_m,(U,O,E)_f] + M.[(F,T,0)_m,(0,0,0)_f] */
  for (i = 0; i < n11; ++i)
    {
      b [i] = - b[i] + y [i];
    }

  free (x);
  free (y);
  free (v5_0);
}
/* calc atimes of (natural) mobility problem
 * for both periodic and non-periodic boundary conditions
 * where A.x := [(u,o,0)_m,(0,0,0)_f] - M.[(0,0,s)_m,(f,t,s)_f].
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mix_lub_3fts (int n, const double *x,
		     double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int n5;
  int nm, nf;
  int nf3, nf5;
  int nm3, nm5, nm11;


  sys = (struct stokes *) user_data;
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n5 = np * 5;
  nf3 = nf * 3;
  nf5 = nf * 5;
  nm3 = nm * 3;
  nm5 = nm * 5;
  nm11 = nm * 11;

  double *w = NULL;
  w = (double *) malloc (sizeof (double) * n);
  CHECK_MALLOC (w, "atimes_mix_lub_3fts");
  double *z = NULL;
  z = (double *) malloc (sizeof (double) * n);
  CHECK_MALLOC (z, "atimes_mix_lub_3fts");
  double *v5_0 = NULL;
  v5_0 = (double *) malloc (sizeof (double) * n5);
  CHECK_MALLOC (v5_0, "atimes_mix_lub_3fts");
  double *u = NULL;
  u = (double *) malloc (sizeof (double) * nm3);
  CHECK_MALLOC (u, "atimes_mix_lub_3fts");
  double *o = NULL;
  o = (double *) malloc (sizeof (double) * nm3);
  CHECK_MALLOC (o, "atimes_mix_lub_3fts");
  double *s = NULL;
  s = (double *) malloc (sizeof (double) * nm5);
  CHECK_MALLOC (s, "atimes_mix_lub_3fts");
  double *ff = NULL;
  ff = (double *) malloc (sizeof (double) * nf3);
  CHECK_MALLOC (ff, "atimes_mix_lub_3fts");
  double *tf = NULL;
  tf = (double *) malloc (sizeof (double) * nf3);
  CHECK_MALLOC (tf, "atimes_mix_lub_3fts");
  double *sf = NULL;
  sf = (double *) malloc (sizeof (double) * nf5);
  CHECK_MALLOC (sf, "atimes_mix_lub_3fts");

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
  calc_lub_3fts (sys, y, w);

  /* set z := [(0,0,S)_mobile,(F,T,S)_fixed] */
  set_fts_by_FTS (nm, z, v5_0, v5_0, s);
  set_fts_by_FTS (nf, z + nm11, ff, tf, sf);

  for (i = 0; i < n; ++i)
    {
      w [i] -= z [i];
    }

  atimes_3all (n, w, z, (void *) sys); // sys->version is 2 (FTS)

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
 * with fixed particles in FTS version
 * for both periodic and non-periodic boundary conditions
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
solve_mix_lub_3fts (struct stokes * sys,
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


  sys->version = 2; // FTS version
  np = sys->np;
  nm = sys->nm;
  if (np == nm)
    {
      solve_mob_lub_3fts (sys, f, t, e, u, o, s);
      return;
    }

  nf = np - nm;
  n11 = np * 11;
  nm11 = nm * 11;

  double *uf0 = NULL;
  uf0 = (double *) malloc (sizeof (double) * nf * 3);
  CHECK_MALLOC (uf0, "solve_mix_lub_3fts");
  double *of0 = NULL;
  of0 = (double *) malloc (sizeof (double) * nf * 3);
  CHECK_MALLOC (of0, "solve_mix_lub_3fts");
  double *ef0 = NULL;
  ef0 = (double *) malloc (sizeof (double) * nf * 5);
  CHECK_MALLOC (ef0, "solve_mix_lub_3fts");
  double *e0 = NULL;
  e0  = (double *) malloc (sizeof (double) * nm * 5);
  CHECK_MALLOC (e0, "solve_mix_lub_3fts");
  double *b = NULL;
  b = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (b, "solve_mix_lub_3fts");
  double *x = NULL;
  x = (double *) malloc (sizeof (double) * n11);
  CHECK_MALLOC (x, "solve_mix_lub_3fts");


  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  shift_labo_to_rest_E (sys, nf, ef, ef0);
  shift_labo_to_rest_E (sys, nm, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mix_lub_3fts (sys, f, t, e0, uf0, of0, ef0, b);
  free (uf0);
  free (of0);
  free (ef0);
  free (e0);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_mix_lub_3fts, (void *) sys,
	      sys->it);

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}
