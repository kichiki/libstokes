/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3ft.c,v 4.10 2006/10/22 22:22:45 kichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 3D configuration
 * periodic boundary condition in 3 direction
 * FT version
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
#include "ft.h"
#include "f.h"
#include "ewald.h" // atimes_ewald_3all()
#include "lub.h" // calc_lub_ewald_3ft()

#include "ewald-3ft.h"


/** natural resistance problem **/
/* solve natural resistance problem in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_ewald_3ft (struct stokes * sys,
		     const double *u, const double *o,
		     double *f, double *t)
{
  int np;
  int i;
  int n6;

  double *b;
  double *x;


  sys->version = 1; // FT version
  np = sys->np;
  n6 = np * 6;

  double *u0;
  double *o0;
  u0 = (double *) malloc (sizeof (double) * np * 3);
  o0 = (double *) malloc (sizeof (double) * np * 3);
  b  = (double *) malloc (sizeof (double) * n6);
  x  = (double *) malloc (sizeof (double) * n6);
  if (u0 == NULL ||
      o0 == NULL ||
      b == NULL ||
      x == NULL)
    {
      fprintf (stderr, "libstokes: allocation error"
	       " at solve_res_ewald_3ft()\n");
      exit (1);
    }

  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  set_ft_by_FT (np, b, u0, o0);
  free (u0);
  free (o0);

  /* first guess */
  for (i = 0; i < n6; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n6, b, x,
	      atimes_ewald_3all, (void *) sys,
	      sys->it);
  // for atimes_ewald_3all(), sys->version is 1 (FT)

  set_FT_by_ft (np, f, t, x);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}

/** natural mobility problem **/
/* solve natural mobility problem in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [np * 3] :
 *   t [np * 3] :
 * OUTPUT
 *   u [np * 3] :
 *   o [np * 3] :
 */
void
solve_mob_ewald_3ft (struct stokes * sys,
		     const double *f, const double *t,
		     double *u, double *o)
{
  int np;
  int n6;

  double *b;
  double *x;


  sys->version = 1; // FT version
  np = sys->np;
  n6 = np * 6;

  b = (double *) malloc (sizeof (double) * n6);
  x = (double *) malloc (sizeof (double) * n6);

  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  set_ft_by_FT (np, x, f, t);
  atimes_ewald_3all (n6, x, b, (void *) sys); // sys->version is 1 (FT)
  set_FT_by_ft (np, u, o, b);

  free (b);
  free (x);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, sys->np, u);
  shift_rest_to_labo_O (sys, sys->np, o);
}


/** natural mobility problem with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem under Ewald sum
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
calc_b_mix_ewald_3ft (struct stokes * sys,
		      const double *f, const double *t,
		      const double *uf, const double *of,
		      double *b)
{
  int np, nm;
  int i;
  int i3, i6;
  int j;
  int nf;
  int n3, n6;
  int nm6;

  double *x;
  double *v3_0;


  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FT.\n");
      sys->version = 1;
    }
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  n6 = np * 6;
  nm6 = nm * 6;

  x = (double *) malloc (sizeof (double) * n6);
  v3_0 = (double *) malloc (sizeof (double) * n3);

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set x := [(F,T)_m,(0,0)_f] */
  set_ft_by_FT (nm, x, f, t);
  set_ft_by_FT (nf, x + nm6, v3_0, v3_0);
  atimes_ewald_3all (n6, x, b, (void *) sys); // sys->version is 1 (FT)

  /* set b := M.x - [(0,0)_m,(U,O)_f] */
  for (i = 0; i < nf; ++i)
    {
      i3 = i * 3;
      i6 = (i + nm) * 6;
      for (j = 0; j < 3; ++j)
	{
	  b [i6 + j] -= uf [i3 + j];
	  b [i6 + 3 + j] -= of [i3 + j];
	}
    }

  free (x);
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
atimes_mix_ewald_3ft (int n, const double *x, double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int np3;
  int nm, nf;
  int nf3;
  int nm3, nm6;

  double *z;
  double *v3_0;
  double *u;
  double *o;
  double *ff;
  double *tf;


  sys = (struct stokes *) user_data;
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FT.\n");
      sys->version = 1;
    }
  np = sys->np;
  nm = sys->nm;

  np3 = np * 3;
  nf = np - nm;
  nf3 = nf * 3;
  nm3 = nm * 3;
  nm6 = nm * 6;

  z = (double *) malloc (sizeof (double) * n);
  v3_0 = (double *) malloc (sizeof (double) * np3);
  u = (double *) malloc (sizeof (double) * nm3);
  o = (double *) malloc (sizeof (double) * nm3);
  ff = (double *) malloc (sizeof (double) * nf3);
  tf = (double *) malloc (sizeof (double) * nf3);

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
  atimes_ewald_3all (n, y, z, (void *) sys); // sys->version is 1 (FT)

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
 * under Ewald sum
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
solve_mix_ewald_3ft (struct stokes * sys,
		     const double *f, const double *t,
		     const double *uf, const double *of,
		     double *u, double *o,
		     double *ff, double *tf)
{
  int np, nm;
  int i;
  int n6;
  int nf;
  int nm6;

  double *b;
  double *x;


  sys->version = 1; // FT version
  np = sys->np;
  nm = sys->nm;
  if (np == nm)
    {
      solve_mob_ewald_3ft (sys, f, t, uf, of,
			   u, o, ff, tf);
      return;
    }

  nf = np - nm;
  n6 = np * 6;
  nm6 = nm * 6;

  double *uf0;
  double *of0;
  uf0 = (double *) malloc (sizeof (double) * nf * 3);
  of0 = (double *) malloc (sizeof (double) * nf * 3);
  b = (double *) malloc (sizeof (double) * n6);
  x = (double *) malloc (sizeof (double) * n6);
  if (uf0 == NULL ||
      of0 == NULL ||
      b == NULL ||
      x == NULL)
    {
      fprintf (stderr, "libstokes: allocation error"
	       " at solve_mix_ewald_3ft()\n");
      exit (1);
    }

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mix_ewald_3ft (sys, f, t, uf0, of0, b);
  free (uf0);
  free (of0);

  /* first guess */
  for (i = 0; i < n6; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n6, b, x,
	      atimes_mix_ewald_3ft, (void *) sys,
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
/* solve natural resistance problem with lubrication
 * in FT version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
void
solve_res_lub_ewald_3ft (struct stokes * sys,
			 const double *u, const double *o,
			 double *f, double *t)
{
  int np;
  int i;
  int n6;

  double *b;
  double *x;
  double *lub;


  sys->version = 1; // FT version
  np = sys->np;
  n6 = np * 6;

  double *u0;
  double *o0;
  u0  = (double *) malloc (sizeof (double) * np * 3);
  o0  = (double *) malloc (sizeof (double) * np * 3);
  b = (double *) malloc (sizeof (double) * n6);
  x = (double *) malloc (sizeof (double) * n6);
  lub = (double *) malloc (sizeof (double) * n6);
  if (u0 == NULL ||
      o0 == NULL ||
      b == NULL ||
      x == NULL ||
      lub == NULL)
    {
      fprintf (stderr, "libstokes: allocation error"
	       " at solve_res_lub_ewald_3ft()\n");
      exit (1);
    }

  shift_labo_to_rest_U (sys, np, u, u0);
  shift_labo_to_rest_O (sys, np, o, o0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  set_ft_by_FT (np, b, u0, o0);
  free (u0);
  free (o0);

  calc_lub_ewald_3ft (sys, b, lub);
  atimes_ewald_3all (n6, lub, x, (void *) sys);
  // sys->version is 1 (FT)
  // x[] is used temporarily
  for (i = 0; i < n6; ++i)
    {
      b [i] += x [i];
    }

  /* first guess */
  for (i = 0; i < n6; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n6, b, x,
	      atimes_ewald_3all, (void *) sys,
	      sys->it);
  // for atimes_ewald_3all(), sys->version is 1 (FT)

  set_FT_by_ft (np, f, t, x);

  free (b);
  free (x);
  free (lub);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  // here, no velocity in output, we do nothing
}


/** natural mobility problem with lubrication with fixed particles **/
/* calc b-term (constant term) of (natural) mobility problem
 * with lubrication under Ewald sum
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
calc_b_mix_lub_ewald_3ft (struct stokes * sys,
			  const double *f, const double *t,
			  const double *uf, const double *of,
			  double *b)
{
  int np, nm;
  int i;
  int nf;
  int n3, n6;
  int nm6;

  double *x;
  double *y;
  double *v3_0;


  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FT.\n");
      sys->version = 1;
    }
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n3 = np * 3;
  n6 = np * 6;
  nm6 = nm * 6;

  x = (double *) malloc (sizeof (double) * n6);
  y = (double *) malloc (sizeof (double) * n6);
  v3_0 = (double *) malloc (sizeof (double) * n3);

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set b := [(0,0)_m,(U,O)_f] */
  set_ft_by_FT (nm, b, v3_0, v3_0);
  set_ft_by_FT (nf, b + nm6, uf, of);

  /* set y := L.[(0,0)_m,(U,O)_f] */
  calc_lub_ewald_3ft (sys, b, y);

  /* set x := [(F,T)_m,(0,0)_f] */
  set_ft_by_FT (nm, x, f, t);
  set_ft_by_FT (nf, x + nm6, v3_0, v3_0);

  /* set x := x - y */
  for (i = 0; i < n6; ++i)
    {
      x [i] -= y [i];
    }

  atimes_ewald_3all (n6, x, y, (void *) sys); // sys->version is 1 (FT)

  /* set b := - (I + M.L).[(0,0)_m,(U,O)_f] + M.[(F,T)_m,(0,0)_f] */
  for (i = 0; i < n6; ++i)
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
atimes_mix_lub_ewald_3ft (int n, const double *x,
			  double *y, void * user_data)
{
  struct stokes * sys;

  int i;
  int np;
  int np3;
  int nm, nf;
  int nf3;
  int nm3, nm6;

  double *w;
  double *z;
  double *v3_0;
  double *u;
  double *o;
  double *ff;
  double *tf;


  sys = (struct stokes *) user_data;
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FT.\n");
      sys->version = 1;
    }
  np = sys->np;
  nm = sys->nm;

  np3 = np * 3;
  nf = np - nm;
  nf3 = nf * 3;
  nm3 = nm * 3;
  nm6 = nm * 6;

  w = (double *) malloc (sizeof (double) * n);
  z = (double *) malloc (sizeof (double) * n);
  v3_0 = (double *) malloc (sizeof (double) * np3);
  u = (double *) malloc (sizeof (double) * nm3);
  o = (double *) malloc (sizeof (double) * nm3);
  ff = (double *) malloc (sizeof (double) * nf3);
  tf = (double *) malloc (sizeof (double) * nf3);

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
  calc_lub_ewald_3ft (sys, y, w);

  /* set z := [(0,0)_mobile,(F,T)_fixed] */
  set_ft_by_FT (nm, z, v3_0, v3_0);
  set_ft_by_FT (nf, z + nm6, ff, tf);

  for (i = 0; i < n; ++i)
    {
      w [i] -= z [i];
    }

  atimes_ewald_3all (n, w, z, (void *) sys); // sys->version is 1 (FT)

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
 * with fixed particles in FT version under Ewald sum
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
solve_mix_lub_ewald_3ft (struct stokes * sys,
			 const double *f, const double *t,
			 const double *uf, const double *of,
			 double *u, double *o,
			 double *ff, double *tf)
{
  int np, nm;

  int i;
  int n6;
  int nf;
  int nm6;

  double *b;
  double *x;


  sys->version = 1; // FT version
  np = sys->np;
  nm = sys->nm;
  if (np == nm)
    {
      solve_mob_lub_ewald_3ft (sys, f, t, uf, of,
			       u, o, ff, tf);
      return;
    }

  nf = np - nm;
  n6 = np * 6;
  nm6 = nm * 6;

  double *uf0;
  double *of0;
  uf0 = (double *) malloc (sizeof (double) * nf * 3);
  of0 = (double *) malloc (sizeof (double) * nf * 3);
  b = (double *) malloc (sizeof (double) * n6);
  x = (double *) malloc (sizeof (double) * n6);
  if (uf0 == NULL ||
      of0 == NULL ||
      b == NULL ||
      x == NULL)
    {
      fprintf (stderr, "libstokes: allocation error"
	       " at solve_mix_lub_ewald_3ft()\n");
      exit (1);
    }

  shift_labo_to_rest_U (sys, nf, uf, uf0);
  shift_labo_to_rest_O (sys, nf, of, of0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */

  calc_b_mix_lub_ewald_3ft (sys, f, t, uf0, of0, b);
  free (uf0);
  free (of0);

  /* first guess */
  for (i = 0; i < n6; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n6, b, x,
	       atimes_mix_lub_ewald_3ft, (void *) sys,
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
