/* backup of bug fixing for polydisperse systems
 * Solvers for 3 dimensional F version problems
 * Copyright (C) 1993-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include "f.h"
#include "ewald-old.h" // atimes_3all()
#include "lub-old.h" // calc_lub_3f()

#include "ewald-3f-old.h"


/** natural resistance problem **/
/* solve natural resistance problem in F version in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_3f_0_old
(struct stokes *sys,
 const double *u,
 double *f)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_res_3f_0_old :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int n3 = sys->np * 3;
  solve_iter (n3, u, f,
	      atimes_3all_old, (void *)sys,
	      sys->it);
  // for atimes_3all(), sys->version is 0 (F)
}
/* solve natural resistance problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_3f_old
(struct stokes *sys,
 const double *u,
 double *f)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_res_3f_old :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;
  int n3 = np * 3;
  double *u0 = (double *)malloc (sizeof (double) * n3);
  CHECK_MALLOC (u0, "solve_res_3f");

  shift_labo_to_rest_U (sys, np, u, u0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_3f_0_old (sys, u0, f);

  free (u0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}


/** natural mobility problem **/
/* solve natural mobility problem in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [np * 3] :
 * OUTPUT
 *  u [np * 3] :
 */
void
solve_mob_3f_old
(struct stokes * sys,
 const double *f,
 double *u)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_mob_3f_old :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  atimes_3all_old (sys->np * 3, f, u, (void *) sys);
  // sys->version is 0 (F)

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
}


/** natural mobility problem with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem
 * for both periodic and non-periodic boundary conditions
 * where b := - [(0,0)_m,(u,o)_f] + M.[(f,t)_m,(0,0)_f].
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *  b [np * 3] : constant vector
 */
static void
calc_b_mix_3f_old
(struct stokes * sys,
 const double *f,
 const double *uf,
 double *b)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to F.\n");
      sys->version = 0;
    }
  int np = sys->np;
  int nm = sys->nm;

  int nf = np - nm;
  int n3 = np * 3;
  int nm3 = nm * 3;

  double *x = (double *)malloc (sizeof (double) * n3);
  double *v3_0 = (double *)malloc (sizeof (double) * n3);
  CHECK_MALLOC (x, "calc_b_mix_3f_old");
  CHECK_MALLOC (v3_0, "calc_b_mix_3f_old");

  int i;
  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set x := [(F)_m,(0)_f] */
  set_F_by_f (nm, x, f);
  set_F_by_f (nf, x + nm3, v3_0);
  atimes_3all_old (n3, x, b, (void *) sys); // sys->version is 0 (F)

  /* set b := M.x - [(0)_m,(U)_f] */
  for (i = 0; i < nf; ++i)
    {
      int i3 = i * 3;
      int im3 = (i + nm) * 3;
      int j;
      for (j = 0; j < 3; ++j)
	{
	  b [im3 + j] -= uf [i3 + j];
	}
    }

  free (x);
  free (v3_0);
}
/* calc atimes of (natural) mobility problem
 * for both periodic and non-periodic boundary conditions
 * where A.x := [(u)_m,(0)_f] - M.[(0)_m,(f)_f].
 * INPUT
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 *  user_data = (struct stokes *) sys : system parameters
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mix_3f_old
(int n, const double *x, double *y, void * user_data)
{
  struct stokes *sys = (struct stokes *) user_data;
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to F.\n");
      sys->version = 0;
    }
  int np = sys->np;
  int nm = sys->nm;

  int np3 = np * 3;
  int nf = np - nm;
  int nf3 = nf * 3;
  int nm3 = nm * 3;

  double *z = (double *)malloc (sizeof (double) * n);
  double *v3_0 = (double *)malloc (sizeof (double) * np3);
  double *u = (double *)malloc (sizeof (double) * nm3);
  double *ff = (double *)malloc (sizeof (double) * nf3);
  CHECK_MALLOC (z, "atimes_mix_3f_old");
  CHECK_MALLOC (v3_0, "atimes_mix_3f_old");
  CHECK_MALLOC (u, "atimes_mix_3f_old");
  CHECK_MALLOC (ff, "atimes_mix_3f_old");

  int i;
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
  atimes_3all_old (n, y, z, (void *) sys); // sys->version is 0 (F)

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
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
solve_mix_3f_old
(struct stokes * sys,
 const double *f,
 const double *uf,
 double *u,
 double *ff)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_mix_3f_old :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_3f_old (sys, f, u);
      return;
    }

  int nf = np - nm;
  int n3 = np * 3;
  int nm3 = nm * 3;

  double *uf0 = (double *)malloc (sizeof (double) * nf * 3);
  double *b = (double *)malloc (sizeof (double) * n3);
  double *x = (double *)malloc (sizeof (double) * n3);
  CHECK_MALLOC (uf0, "solve_mix_3f_old");
  CHECK_MALLOC (b, "solve_mix_3f_old");
  CHECK_MALLOC (x, "solve_mix_3f_old");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mix_3f_old (sys, f, uf0, b);
  free (uf0);

  solve_iter (n3, b, x,
	      atimes_mix_3f_old, (void *)sys,
	      sys->it);

  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
}



/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in F version
 * in the fluid-rest frame
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : = U - u^inf, the velocity in the fluid-rest frame
 * OUTPUT
 *  f [np * 3] :
 */
void
solve_res_lub_3f_0_old
(struct stokes * sys,
 const double *u,
 double *f)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_res_lub_3f_0_old :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int n3 = sys->np * 3;

  double *lub = (double *)malloc (sizeof (double) * n3);
  CHECK_MALLOC (lub, "solve_res_lub_3f");

  calc_lub_3f_old (sys, u, lub);
  atimes_3all_old (n3, lub, f, (void *) sys);
  // sys->version is 0 (F)
  // lub[] is used temporarily
  int i;
  for (i = 0; i < n3; ++i)
    {
      lub [i] = u [i] + f [i];
    }

  solve_iter (n3, lub, f,
	      atimes_3all_old, (void *)sys,
	      sys->it);
  // for atimes_3all(), sys->version is 0 (F)

  free (lub);
}
/* solve natural resistance problem with lubrication in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  u [np * 3] : particle velocity in the labo frame.
 * OUTPUT
 *   f [np * 3] :
 */
void
solve_res_lub_3f_old
(struct stokes * sys,
 const double *u,
 double *f)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_res_lub_3f_old :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;
  int n3 = np * 3;
  double *u0 = (double *)malloc (sizeof (double) * n3);
  CHECK_MALLOC (u0, "solve_res_lub_3f_old");

  shift_labo_to_rest_U (sys, np, u, u0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  solve_res_lub_3f_0_old (sys, u0, f);

  free (u0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}



/** mob_lub_3f **/
void
atimes_mob_lub_3f_old
(int n, const double *x,
 double *y, void *user_data)
{
  struct stokes *sys = (struct stokes *) user_data;
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to F.\n");
      sys->version = 0;
    }

  double *lub = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (lub, "atimes_mob_lub_3f_old");

  calc_lub_3f_old (sys, x, y);
  atimes_3all_old (n, y, lub, (void *) sys);
  // lub = M.L.U

  int i;
  for (i = 0; i < n; i ++)
    {
      y [i] = x [i] + lub [i];
    }
  // y = (I + M.L).U

  free (lub);
}

/* solve natural mobility problem with lubrication in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 */
void
solve_mob_lub_3f_old
(struct stokes * sys,
 const double *f,
 double *u)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_mob_lub_3f_old :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;
  int n3 = np * 3;

  double *b = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (b, "solve_mob_lub_3f_old");

  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  // f itself is F
  atimes_3all_old (n3, f, b, (void *) sys); // sys->version is 1 (FT)
  // b = M.F

  solve_iter (n3, b, u,
	      atimes_mob_lub_3f_old, (void *)sys,
	      sys->it);
  // u = (I + M.L)^-1 . b

  free (b);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
}


/** natural mobility problem with lubrication with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem with lubrication
 * for both periodic and non-periodic boundary conditions
 * where b := -(0) + M.(f).
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *  b [np * 3] : constant vector
 */
static void
calc_b_mix_lub_3f_old
(struct stokes * sys,
 const double *f,
 const double *uf,
 double *b)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to F.\n");
      sys->version = 0;
    }
  int np = sys->np;
  int nm = sys->nm;

  int nf = np - nm;
  int n3 = np * 3;
  int nm3 = nm * 3;

  double *x = (double *) malloc (sizeof (double) * n3);
  double *y = (double *) malloc (sizeof (double) * n3);
  double *v3_0 = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (x, "calc_b_mix_lub_3f_old");
  CHECK_MALLOC (y, "calc_b_mix_lub_3f_old");
  CHECK_MALLOC (v3_0, "calc_b_mix_lub_3f_old");

  int i;
  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set b := [(0)_m,(U)_f] */
  set_F_by_f (nm, b, v3_0);
  set_F_by_f (nf, b + nm3, uf);

  /* set y := L.[(0)_m,(U)_f] */
  calc_lub_3f_old (sys, b, y);

  /* set x := [(F)_m,(0)_f] */
  set_F_by_f (nm, x, f);
  set_F_by_f (nf, x + nm3, v3_0);

  /* set x := x - y */
  for (i = 0; i < n3; ++i)
    {
      x [i] -= y [i];
    }

  atimes_3all_old (n3, x, y, (void *) sys); // sys->version is 0 (F)

  /* set b := - (I + M.L).[(0)_m,(U)_f] + M.[(F)_m,(0)_f] */
  for (i = 0; i < n3; ++i)
    {
      b [i] = - b [i] + y [i];
    }

  free (x);
  free (y);
  free (v3_0);
}

/* calc atimes of (natural) mobility problem with lubrication
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
atimes_mix_lub_3f_old
(int n, const double *x,
 double *y, void * user_data)
{
  struct stokes *sys = (struct stokes *) user_data;
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to F.\n");
      sys->version = 0;
    }
  int np = sys->np;
  int nm = sys->nm;

  int np3 = np * 3;
  int nf = np - nm;
  int nf3 = nf * 3;
  int nm3 = nm * 3;

  double *w = (double *) malloc (sizeof (double) * n);
  double *z = (double *) malloc (sizeof (double) * n);
  double *v3_0 = (double *) malloc (sizeof (double) * np3);
  double *u = (double *) malloc (sizeof (double) * nm3);
  double *ff = (double *) malloc (sizeof (double) * nf3);
  CHECK_MALLOC (w, "atimes_mix_lub_3f_old");
  CHECK_MALLOC (z, "atimes_mix_lub_3f_old");
  CHECK_MALLOC (v3_0, "atimes_mix_lub_3f_old");
  CHECK_MALLOC (u, "atimes_mix_lub_3f_old");
  CHECK_MALLOC (ff, "atimes_mix_lub_3f_old");

  int i;
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
  calc_lub_3f_old (sys, y, w);

  /* set z := [(0)_mobile,(F)_fixed] */
  set_F_by_f (nm, z, v3_0);
  set_F_by_f (nf, z + nm3, ff);

  for (i = 0; i < n; ++i)
    {
      w [i] -= z [i];
    }

  atimes_3all_old (n, w, z, (void *) sys); // sys->version is 0 (F)

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
 * with fixed particles in F version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
void
solve_mix_lub_3f_old
(struct stokes * sys,
 const double *f,
 const double *uf,
 double *u,
 double *ff)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_mix_lub_3f_old :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;
  int nm = sys->nm;
  if (np == nm)
    {
      solve_mob_lub_3f_old (sys, f, u);
      /*
      fprintf (stderr, "solve_mix_lub_3f:"
	       " no fixed particle. mob_lub is not implemented yet."
	       " call plain mob solver\n");
      solve_mob_3f (sys, f, u);
      */
      return;
    }

  int nf = np - nm;
  int n3 = np * 3;
  int nm3 = nm * 3;

  double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
  double *b = (double *) malloc (sizeof (double) * n3);
  double *x = (double *) malloc (sizeof (double) * n3);
  CHECK_MALLOC (uf0, "solve_mix_lub_3f_old");
  CHECK_MALLOC (b, "solve_mix_lub_3f_old");
  CHECK_MALLOC (x, "solve_mix_lub_3f_old");

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mix_lub_3f_old (sys, f, uf0, b);
  free (uf0);

  solve_iter (n3, b, x,
	      atimes_mix_lub_3f_old, (void *)sys,
	      sys->it);

  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
}
