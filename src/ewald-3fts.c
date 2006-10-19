/* Beenakker's formulation of Ewald summation technique for RP tensor in 3D
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-3fts.c,v 5.8 2006/10/19 04:14:40 ichiki Exp $
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
#include "ewald.h" // atimes_ewald_3all()
#include "lub.h" // calc_lub_ewald_3fts()

#include "ewald-3fts.h"


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
solve_res_ewald_3fts (struct stokes * sys,
		      const double *u, const double *o, const double *e,
		      double *f, double *t, double *s)
{
  int np;
  int i;
  int n11;

  double *b;
  double *x;


  sys->version = 2; // FTS version
  np = sys->np;

  n11 = np * 11;
  b = (double *) malloc (sizeof (double) * n11);
  x = (double *) malloc (sizeof (double) * n11);

  set_fts_by_FTS (np, b, u, o, e);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_ewald_3all, (void *) sys,
	      sys->it);
  // for atimes_ewald_3all(), sys->version is 2 (FTS)

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


  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }
  np = sys->np;

  n3 = np * 3;
  n5 = np * 5;
  n11 = np * 11;

  x = (double *) malloc (sizeof (double) * n11);
  v5_0 = (double *) malloc (sizeof (double) * n5);

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  set_fts_by_FTS (np, x, f, t, v5_0);
  atimes_ewald_3all (n11, x, b, (void *) sys); // sys->version is 2 (FTS)

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
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes: version is wrong. reset to FTS.\n");
      sys->version = 2;
    }
  np = sys->np;

  n3 = np * 3;
  n5 = np * 5;

  z = (double *) malloc (sizeof (double) * n);
  v5_0 = (double *) malloc (sizeof (double) * n5);
  u = (double *) malloc (sizeof (double) * n3);
  o = (double *) malloc (sizeof (double) * n3);
  s = (double *) malloc (sizeof (double) * n5);

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  set_FTS_by_fts (np, u, o, s, x);

  set_fts_by_FTS (np, y, v5_0, v5_0, s);
  atimes_ewald_3all (n, y, z, (void *) sys); // sys->version is 2 (FTS)

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
solve_mob_ewald_3fts (struct stokes * sys,
		      const double *f, const double *t, const double *e,
		      double *u, double *o, double *s)
{
  int np;
  int i;
  int n11;

  double *b;
  double *x;


  sys->version = 2; // FTS version
  np = sys->np;
  n11 = np * 11;
  b = (double *) malloc (sizeof (double) * n11);
  x = (double *) malloc (sizeof (double) * n11);

  calc_b_mob_ewald_3fts (sys, f, t, e, b);

  /* first guess */
  for (i = 0; i < n11; ++i)
    {
      x [i] = 0.0;
    }

  solve_iter (n11, b, x,
	      atimes_mob_ewald_3fts, (void *) sys,
	      sys->it);

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
calc_b_mix_ewald_3fts (struct stokes * sys,
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

  x = (double *) malloc (sizeof (double) * n11);
  v5_0 = (double *) malloc (sizeof (double) * n5);

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  /* set x := [(F,T,0)_m,(0,0,0)_f] */
  set_fts_by_FTS (nm, x, f, t, v5_0);
  set_fts_by_FTS (nf, x + nm11, v5_0, v5_0, v5_0);
  atimes_ewald_3all (n11, x, b, (void *) sys); // sys->version is 2 (FTS)

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
atimes_mix_ewald_3fts (int n, const double *x, double *y, void * user_data)
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

  z = (double *) malloc (sizeof (double) * n);
  v5_0 = (double *) malloc (sizeof (double) * n5);
  u = (double *) malloc (sizeof (double) * nm3);
  o = (double *) malloc (sizeof (double) * nm3);
  s = (double *) malloc (sizeof (double) * nm5);
  ff = (double *) malloc (sizeof (double) * nf3);
  tf = (double *) malloc (sizeof (double) * nf3);
  sf = (double *) malloc (sizeof (double) * nf5);

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
  atimes_ewald_3all (n, y, z, (void *) sys); // sys->version is 2 (FTS)

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
solve_mix_ewald_3fts (struct stokes * sys,
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


  sys->version = 2; // FTS version
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n11 = np * 11;
  nm11 = nm * 11;

  double *uf0;
  double *of0;
  double *ef0;
  double *e0;
  uf0 = (double *) malloc (sizeof (double) * nf * 3);
  of0 = (double *) malloc (sizeof (double) * nf * 3);
  ef0 = (double *) malloc (sizeof (double) * nf * 5);
  e0  = (double *) malloc (sizeof (double) * nm * 5);
  b = (double *) malloc (sizeof (double) * n11);
  x = (double *) malloc (sizeof (double) * n11);
  if (uf0 == NULL ||
      of0 == NULL ||
      ef0 == NULL ||
      e0  == NULL ||
      b == NULL ||
      x == NULL)
    {
      fprintf (stderr, "libstokes: allocation error"
	       " at solve_mix_ewald_3fts()\n");
      exit (1);
    }

  int i5, ix, iy, iz;
  for (i = 0; i < nf; i ++)
    {
      i5 = i*5;
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;
      uf0 [ix] = uf[ix] - sys->Ui[0];
      uf0 [iy] = uf[iy] - sys->Ui[1];
      uf0 [iz] = uf[iz] - sys->Ui[2];

      of0 [ix] = of[ix] - sys->Oi[0];
      of0 [iy] = of[iy] - sys->Oi[1];
      of0 [iz] = of[iz] - sys->Oi[2];

      ef0 [i5+0] = ef[i5+0] - sys->Ei[0];
      ef0 [i5+1] = ef[i5+1] - sys->Ei[1];
      ef0 [i5+2] = ef[i5+2] - sys->Ei[2];
      ef0 [i5+3] = ef[i5+3] - sys->Ei[3];
      ef0 [i5+4] = ef[i5+4] - sys->Ei[4];
    }
  for (i = 0; i < nm; i ++)
    {
      i5 = i*5;
      e0 [i5+0] = e[i5+0] - sys->Ei[0];
      e0 [i5+1] = e[i5+1] - sys->Ei[1];
      e0 [i5+2] = e[i5+2] - sys->Ei[2];
      e0 [i5+3] = e[i5+3] - sys->Ei[3];
      e0 [i5+4] = e[i5+4] - sys->Ei[4];
    }

  //calc_b_mix_ewald_3fts (sys, f, t, e, uf, of, ef, b);
  calc_b_mix_ewald_3fts (sys, f, t, e0, uf0, of0, ef0, b);
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
	      atimes_mix_ewald_3fts, (void *) sys,
	      sys->it);

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  for (i = 0; i < nm; i ++)
    {
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;
      u [ix] += sys->Ui[0];
      u [iy] += sys->Ui[1];
      u [iz] += sys->Ui[2];

      o [ix] += sys->Oi[0];
      o [iy] += sys->Oi[1];
      o [iz] += sys->Oi[2];
    }

  free (b);
  free (x);
}

/** natural resistance problem with lubrication **/
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
solve_res_lub_ewald_3fts (struct stokes * sys,
			  const double *u, const double *o, const double *e,
			  double *f, double *t, double *s)
{
  int np;
  int i;
  int n11;

  double *b;
  double *x;
  double *lub;


  sys->version = 2; // FTS version
  np = sys->np;
  n11 = np * 11;
  b = (double *) malloc (sizeof (double) * n11);
  x = (double *) malloc (sizeof (double) * n11);
  lub = (double *) malloc (sizeof (double) * n11);

  set_fts_by_FTS (np, b, u, o, e);
  calc_lub_ewald_3fts (sys, b, lub);
  atimes_ewald_3all (n11, lub, x, (void *) sys);
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
	      atimes_ewald_3all, (void *) sys,
	      sys->it);
  // for atimes_ewald_3all(), sys->version is 2 (FTS)

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
calc_b_mix_lub_ewald_3fts (struct stokes * sys,
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

  x = (double *) malloc (sizeof (double) * n11);
  y = (double *) malloc (sizeof (double) * n11);
  v5_0 = (double *) malloc (sizeof (double) * n5);

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

  atimes_ewald_3all (n11, x, y, (void *) sys); // sys->version is 2 (FTS)

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
atimes_mix_lub_ewald_3fts (int n, const double *x,
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

  w = (double *) malloc (sizeof (double) * n);
  z = (double *) malloc (sizeof (double) * n);
  v5_0 = (double *) malloc (sizeof (double) * n5);
  u = (double *) malloc (sizeof (double) * nm3);
  o = (double *) malloc (sizeof (double) * nm3);
  s = (double *) malloc (sizeof (double) * nm5);
  ff = (double *) malloc (sizeof (double) * nf3);
  tf = (double *) malloc (sizeof (double) * nf3);
  sf = (double *) malloc (sizeof (double) * nf5);

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

  atimes_ewald_3all (n, w, z, (void *) sys); // sys->version is 2 (FTS)

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
solve_mix_lub_ewald_3fts (struct stokes * sys,
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


  sys->version = 2; // FTS version
  np = sys->np;
  nm = sys->nm;

  nf = np - nm;
  n11 = np * 11;
  nm11 = nm * 11;

  double *uf0;
  double *of0;
  double *ef0;
  double *e0;
  uf0 = (double *) malloc (sizeof (double) * nf * 3);
  of0 = (double *) malloc (sizeof (double) * nf * 3);
  ef0 = (double *) malloc (sizeof (double) * nf * 5);
  e0  = (double *) malloc (sizeof (double) * nm * 5);
  b = (double *) malloc (sizeof (double) * n11);
  x = (double *) malloc (sizeof (double) * n11);
  if (uf0 == NULL ||
      of0 == NULL ||
      ef0 == NULL ||
      e0  == NULL ||
      b == NULL ||
      x == NULL)
    {
      fprintf (stderr, "libstokes: allocation error"
	       " at solve_mix_lub_ewald_3fts()\n");
      exit (1);
    }

  int i5, ix, iy, iz;
  for (i = 0; i < nf; i ++)
    {
      i5 = i*5;
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;
      uf0 [ix] = uf[ix] - sys->Ui[0];
      uf0 [iy] = uf[iy] - sys->Ui[1];
      uf0 [iz] = uf[iz] - sys->Ui[2];

      of0 [ix] = of[ix] - sys->Oi[0];
      of0 [iy] = of[iy] - sys->Oi[1];
      of0 [iz] = of[iz] - sys->Oi[2];

      ef0 [i5+0] = ef[i5+0] - sys->Ei[0];
      ef0 [i5+1] = ef[i5+1] - sys->Ei[1];
      ef0 [i5+2] = ef[i5+2] - sys->Ei[2];
      ef0 [i5+3] = ef[i5+3] - sys->Ei[3];
      ef0 [i5+4] = ef[i5+4] - sys->Ei[4];
    }
  for (i = 0; i < nm; i ++)
    {
      i5 = i*5;
      e0 [i5+0] = e[i5+0] - sys->Ei[0];
      e0 [i5+1] = e[i5+1] - sys->Ei[1];
      e0 [i5+2] = e[i5+2] - sys->Ei[2];
      e0 [i5+3] = e[i5+3] - sys->Ei[3];
      e0 [i5+4] = e[i5+4] - sys->Ei[4];
    }

  //calc_b_mix_lub_ewald_3fts (sys, f, t, e, uf, of, ef, b);
  calc_b_mix_lub_ewald_3fts (sys, f, t, e0, uf0, of0, ef0, b);
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
	      atimes_mix_lub_ewald_3fts, (void *) sys,
	      sys->it);

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  for (i = 0; i < nm; i ++)
    {
      ix = i*3;
      iy = ix + 1;
      iz = ix + 2;
      u [ix] += sys->Ui[0];
      u [iy] += sys->Ui[1];
      u [iz] += sys->Ui[2];

      o [ix] += sys->Oi[0];
      o [iy] += sys->Oi[1];
      o [iz] += sys->Oi[2];
    }

  free (b);
  free (x);
}
