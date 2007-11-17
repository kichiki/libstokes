/* Solvers for 3 dimensional FT version problems
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3ft.c,v 4.16 2007/11/17 23:31:48 kichiki Exp $
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
#include "memory-check.h" // macro CHECK_MALLOC

#include "bench.h" // ptime_ms_d()
#include "stokes.h" /* struct stokeks */
#include "ft.h"
#include "f.h"
#include "ewald.h" // atimes_3all()
#include "lub.h" // calc_lub_3ft()

#include "ewald-3ft.h"


/** natural resistance problem **/
/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 *  o [np * 3] : angular  velocity in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_3ft (struct stokes * sys,
	       const double *u, const double *o,
	       double *f, double *t)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_res_3ft :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  double *u0 = (double *) malloc (sizeof (double) * np * 3);
  double *o0 = (double *) malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (u0, "solve_res_3ft");
  CHECK_MALLOC (o0, "solve_res_3ft");

  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_3ft_0 (sys,
		   u0, o0,
		   f, t);

  free (u0);
  free (o0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}

/* solve natural resistance problem in FT version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 *  t [np * 3] :
 */
void
solve_res_3ft_0 (struct stokes * sys,
		 const double *u, const double *o,
		 double *f, double *t)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_res_3ft_0 :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int n6 = np * 6;

  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (b, "solve_res_3ft_0");
  CHECK_MALLOC (x, "solve_res_3ft_0");

  set_ft_by_FT (np, b, u, o);

  solve_iter (n6, b, x,
	      atimes_3all, (void *)sys,
	      sys->it);
  // for atimes_3all(), sys->version is 1 (FT)

  set_FT_by_ft (np, f, t, x);

  free (b);
  free (x);
}

/** natural mobility problem **/
/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_3ft (struct stokes * sys,
	       const double *f, const double *t,
	       double *u, double *o)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_mob_3ft :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int n6 = np * 6;

  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (b, "solve_mob_3ft");
  CHECK_MALLOC (x, "solve_mob_3ft");

  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  set_ft_by_FT (np, x, f, t);
  atimes_3all (n6, x, b, (void *) sys); // sys->version is 1 (FT)
  set_FT_by_ft (np, u, o, b);

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
 * where b := - [(0,0)_m,(u,o)_f] + M.[(f,t)_m,(0,0)_f].
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  uf [nf * 3] :
 *  of [nf * 3] :
 * OUTPUT
 *  b [np * 6] : constant vector
 */
static void
calc_b_mix_3ft (struct stokes * sys,
		const double *f, const double *t,
		const double *uf, const double *of,
		double *b)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FT.\n");
      sys->version = 1;
    }
  int np = sys->np;
  int nm = sys->nm;

  int nf = np - nm;
  int n3 = np * 3;
  int n6 = np * 6;
  int nm6 = nm * 6;

  double *x = (double *) malloc (sizeof (double) * n6);
  double *v3_0 = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (x, "calc_b_mix_3ft");
  CHECK_MALLOC (v3_0, "calc_b_mix_3ft");

  int i;
  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set x := [(F,T)_m,(0,0)_f] */
  set_ft_by_FT (nm, x, f, t);
  set_ft_by_FT (nf, x + nm6, v3_0, v3_0);
  atimes_3all (n6, x, b, (void *) sys); // sys->version is 1 (FT)

  /* set b := M.x - [(0,0)_m,(U,O)_f] */
  for (i = 0; i < nf; ++i)
    {
      int i3 = i * 3;
      int i6 = (i + nm) * 6;
      int j;
      for (j = 0; j < 3; ++j)
	{
	  b [i6 + j] -= uf [i3 + j];
	  b [i6 + 3 + j] -= of [i3 + j];
	}
    }

  free (x);
  free (v3_0);
}
/* calc atimes of (natural) mobility problem
 * for both periodic and non-periodic boundary conditions
 * where A.x := [(u,o)_m,(0,0)_f] - M.[(0,0)_m,(f,t)_f].
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mix_3ft (int n, const double *x, double *y, void * user_data)
{
  struct stokes *sys = (struct stokes *) user_data;
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FT.\n");
      sys->version = 1;
    }
  int np = sys->np;
  int nm = sys->nm;

  int np3 = np * 3;
  int nf = np - nm;
  int nf3 = nf * 3;
  int nm3 = nm * 3;
  int nm6 = nm * 6;

  double *z = (double *) malloc (sizeof (double) * n);
  double *v3_0 = (double *) malloc (sizeof (double) * np3);
  double *u = (double *) malloc (sizeof (double) * nm3);
  double *o = (double *) malloc (sizeof (double) * nm3);
  double *ff = (double *) malloc (sizeof (double) * nf3);
  double *tf = (double *) malloc (sizeof (double) * nf3);
  CHECK_MALLOC (z, "atimes_mix_3ft");
  CHECK_MALLOC (v3_0, "atimes_mix_3ft");
  CHECK_MALLOC (u, "atimes_mix_3ft");
  CHECK_MALLOC (o, "atimes_mix_3ft");
  CHECK_MALLOC (ff, "atimes_mix_3ft");
  CHECK_MALLOC (tf, "atimes_mix_3ft");

  int i;
  for (i = 0; i < np3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set (U,O)_mobile,(F,T)_fixed by x[] */
  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  /* set y := [(0,0)_mobile,(F,T)_fixed] */
  set_ft_by_FT (nm, y, v3_0, v3_0);
  set_ft_by_FT (nf, y + nm6, ff, tf);
  atimes_3all (n, y, z, (void *) sys); // sys->version is 1 (FT)

  /* set y := [(U,O,0)_mobile,(0,0,0)_fixed] */
  set_ft_by_FT (nm, y, u, o);
  set_ft_by_FT (nf, y + nm6, v3_0, v3_0);

  /* set y := [(U,O)_m,(0,0)_f] - M.[(0,0)_m,(F,T)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
    }

  free (z);
  free (v3_0);
  free (u);
  free (o);
  free (ff);
  free (tf);
}
/* solve natural mobility problem with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
solve_mix_3ft (struct stokes * sys,
	       const double *f, const double *t,
	       const double *uf, const double *of,
	       double *u, double *o,
	       double *ff, double *tf)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_mix_3ft :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_3ft (sys, f, t,
		     u, o);
      return;
    }

  int nf = np - nm;
  int n6 = np * 6;
  int nm6 = nm * 6;

  double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
  double *of0 = (double *) malloc (sizeof (double) * nf * 3);
  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (uf0, "solve_mix_3ft");
  CHECK_MALLOC (of0, "solve_mix_3ft");
  CHECK_MALLOC (b, "solve_mix_3ft");
  CHECK_MALLOC (x, "solve_mix_3ft");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mix_3ft (sys, f, t, uf0, of0, b);
  free (uf0);
  free (of0);

  solve_iter (n6, b, x,
	      atimes_mix_3ft, (void *)sys,
	      sys->it);

  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_lub_3ft (struct stokes * sys,
		   const double *u, const double *o,
		   double *f, double *t)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_res_lub_3ft :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  double *u0 = (double *) malloc (sizeof (double) * np * 3);
  double *o0 = (double *) malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (u0, "solve_res_lub_3ft");
  CHECK_MALLOC (o0, "solve_res_lub_3ft");

  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_lub_3ft_0 (sys,
		       u0, o0,
		       f, t);

  free (u0);
  free (o0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}

/* solve natural resistance problem with lubrication in FT version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 *  o [np * 3] : = O - O^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_lub_3ft_0 (struct stokes * sys,
		     const double *u, const double *o,
		     double *f, double *t)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_res_lub_3ft_0 :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int n6 = np * 6;

  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  double *lub = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (b, "solve_res_lub_3ft_0");
  CHECK_MALLOC (x, "solve_res_lub_3ft_0");
  CHECK_MALLOC (lub, "solve_res_lub_3ft_0");

  set_ft_by_FT (np, b, u, o);

  calc_lub_3ft (sys, b, lub);
  atimes_3all (n6, lub, x, (void *) sys);
  // sys->version is 1 (FT)
  // x[] is used temporarily
  int i;
  for (i = 0; i < n6; ++i)
    {
      b [i] += x [i];
    }

  solve_iter (n6, b, x,
	      atimes_3all, (void *)sys,
	      sys->it);
  // for atimes_3all(), sys->version is 1 (FT)

  set_FT_by_ft (np, f, t, x);

  free (b);
  free (x);
  free (lub);
}


/** mob_lub_3ft **/
static void
atimes_mob_lub_3ft (int n, const double *x,
		    double *y, void *user_data)
{
  struct stokes *sys = (struct stokes *) user_data;
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FT.\n");
      sys->version = 1;
    }

  double *lub = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (lub, "atimes_mob_lub_3ft");

  calc_lub_3ft (sys, x, y);
  atimes_3all (n, y, lub, (void *) sys);
  // lub = M.L.(UO)

  int i;
  for (i = 0; i < n; i ++)
    {
      y [i] = x [i] + lub [i];
    }
  // y = (I + M.L).(UO)

  free (lub);
}

/* solve natural mobility problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_lub_3ft (struct stokes * sys,
		   const double *f, const double *t,
		   double *u, double *o)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_mob_lub_3ft :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int n6 = np * 6;

  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (b, "solve_mob_lub_3ft");
  CHECK_MALLOC (x, "solve_mob_lub_3ft");

  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  set_ft_by_FT (np, x, f, t);

  atimes_3all (n6, x, b, (void *) sys); // sys->version is 1 (FT)
  // b = M.(FT)

  solve_iter (n6, b, x,
	      atimes_mob_lub_3ft, (void *)sys,
	      sys->it);
  // x = (I + M.L)^-1 . b

  set_FT_by_ft (np, u, o, x);

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
 *  uf [nf * 3] :
 *  of [nf * 3] :
 * OUTPUT
 *  b [np * 6] : constant vector
 */
static void
calc_b_mix_lub_3ft (struct stokes * sys,
		    const double *f, const double *t,
		    const double *uf, const double *of,
		    double *b)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FT.\n");
      sys->version = 1;
    }
  int np = sys->np;
  int nm = sys->nm;

  int nf = np - nm;
  int n3 = np * 3;
  int n6 = np * 6;
  int nm6 = nm * 6;

  double *x = (double *) malloc (sizeof (double) * n6);
  double *y = (double *) malloc (sizeof (double) * n6);
  double *v3_0 = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (x, "calc_b_mix_lub_3ft");
  CHECK_MALLOC (y, "calc_b_mix_lub_3ft");
  CHECK_MALLOC (v3_0, "calc_b_mix_lub_3ft");

  int i;
  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set b := [(0,0)_m,(U,O)_f] */
  set_ft_by_FT (nm, b, v3_0, v3_0);
  set_ft_by_FT (nf, b + nm6, uf, of);

  /* set y := L.[(0,0)_m,(U,O)_f] */
  calc_lub_3ft (sys, b, y);

  /* set x := [(F,T)_m,(0,0)_f] */
  set_ft_by_FT (nm, x, f, t);
  set_ft_by_FT (nf, x + nm6, v3_0, v3_0);

  /* set x := x - y */
  for (i = 0; i < n6; ++i)
    {
      x [i] -= y [i];
    }

  atimes_3all (n6, x, y, (void *) sys); // sys->version is 1 (FT)

  /* set b := - (I + M.L).[(0,0)_m,(U,O)_f] + M.[(F,T)_m,(0,0)_f] */
  for (i = 0; i < n6; ++i)
    {
      b [i] = - b [i] + y [i];
    }

  free (x);
  free (y);
  free (v3_0);
}
/* calc atimes of (natural) mobility problem
 * for both periodic and non-periodic boundary conditions
 * where A.x := [(u,o)_m,(0,0)_f] - M.[(0,0)_m,(f,t)_f].
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mix_lub_3ft (int n, const double *x,
		    double *y, void * user_data)
{
  struct stokes *sys = (struct stokes *) user_data;
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FT.\n");
      sys->version = 1;
    }
  int np = sys->np;
  int nm = sys->nm;

  int np3 = np * 3;
  int nf = np - nm;
  int nf3 = nf * 3;
  int nm3 = nm * 3;
  int nm6 = nm * 6;

  double *w = (double *) malloc (sizeof (double) * n);
  double *z = (double *) malloc (sizeof (double) * n);
  double *v3_0 = (double *) malloc (sizeof (double) * np3);
  double *u = (double *) malloc (sizeof (double) * nm3);
  double *o = (double *) malloc (sizeof (double) * nm3);
  double *ff = (double *) malloc (sizeof (double) * nf3);
  double *tf = (double *) malloc (sizeof (double) * nf3);
  CHECK_MALLOC (w, "atimes_mix_lub_3ft");
  CHECK_MALLOC (z, "atimes_mix_lub_3ft");
  CHECK_MALLOC (v3_0, "atimes_mix_lub_3ft");
  CHECK_MALLOC (u, "atimes_mix_lub_3ft");
  CHECK_MALLOC (o, "atimes_mix_lub_3ft");
  CHECK_MALLOC (ff, "atimes_mix_lub_3ft");
  CHECK_MALLOC (tf, "atimes_mix_lub_3ft");

  int i;
  for (i = 0; i < np3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set (U,O)_mobile,(F,T)_fixed by x[] */
  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  /* set y := [(U,O)_mobile,(0,0)_fixed] */
  set_ft_by_FT (nm, y, u, o);
  set_ft_by_FT (nf, y + nm6, v3_0, v3_0);

  /* set w := L.[(U,O,0)_mobile,(0,0,0)_fixed] */
  calc_lub_3ft (sys, y, w);

  /* set z := [(0,0)_mobile,(F,T)_fixed] */
  set_ft_by_FT (nm, z, v3_0, v3_0);
  set_ft_by_FT (nf, z + nm6, ff, tf);

  for (i = 0; i < n; ++i)
    {
      w [i] -= z [i];
    }

  atimes_3all (n, w, z, (void *) sys); // sys->version is 1 (FT)

  /* set y := (I + M.L).[(U,O)_m,(0,0)_f] - M.[(0,0)_m,(F,T)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] += z [i];
    }

  free (w);
  free (z);
  free (v3_0);
  free (u);
  free (o);
  free (ff);
  free (tf);
}
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [nm * 3] :
 *   t [nm * 3] :
 *   uf [nf * 3] :
 *   of [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
solve_mix_lub_3ft (struct stokes * sys,
		   const double *f, const double *t,
		   const double *uf, const double *of,
		   double *u, double *o,
		   double *ff, double *tf)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_mix_lub_3ft :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_lub_3ft (sys, f, t, u, o);
      return;
    }

  int nf = np - nm;
  int n6 = np * 6;
  int nm6 = nm * 6;

  double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
  double *of0 = (double *) malloc (sizeof (double) * nf * 3);
  double *b = (double *) malloc (sizeof (double) * n6);
  double *x = (double *) malloc (sizeof (double) * n6);
  CHECK_MALLOC (uf0, "solve_mix_lub_3ft");
  CHECK_MALLOC (of0, "solve_mix_lub_3ft");
  CHECK_MALLOC (b, "solve_mix_lub_3ft");
  CHECK_MALLOC (x, "solve_mix_lub_3ft");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mix_lub_3ft (sys, f, t, uf0, of0, b);
  free (uf0);
  free (of0);

  solve_iter (n6, b, x,
	       atimes_mix_lub_3ft, (void *)sys,
	       sys->it);

  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}
