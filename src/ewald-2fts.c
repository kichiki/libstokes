/* Ewald summation technique under 2D
 * this is a wrapper package for ewald-3fts.c
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-2fts.c,v 1.1 2001/02/02 08:01:53 ichiki Exp $
 *
 * 3 dimensional hydrodynamics, 2D configuration
 * periodic boundary condition in 3 direction,
 * FTS version
 * non-dimension formulation
 */
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include "ewald-3fts.h"

#include <libiter.h> /* gpb(), solve_iter_stab() */
#include "lub.h"
#include "fts.h"

#include "ewald-2fts.h"


static void
calc_res_lub_ewald_3fts_2d (int np,
			    double *u, double *o, double *e,
			    double *f, double *t, double *s);

static void
calc_mob_lub_fix_ewald_3fts_2d (int np, int nm,
				double *f, double *t, double *e,
				double *uf, double *of, double *ef,
				double *u, double *o, double *s,
				double *ff, double *tf, double *sf);
static void
calc_b_mob_lub_fix_ewald_2fts (int np, int nm,
			       double *f, double *t, double *e,
			       double *uf, double *of, double *ef,
			       double *b);
static void
atimes_mob_lub_fix_ewald_2fts (int n, double *x, double *y);
static void
calc_lub_ewald_2fts (int np, double * uoe, double * fts);
static int
cond_lub_2d (double * x1, double * x2);

/** natural resistance problem **/
/* solve natural resistance problem in FTS version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_ewald_2fts (int np,
		     double *u, double *o, double *e,
		     double *f, double *t, double *s)
{
  int i;
  int np3, np5;
  int i2, i3, i5;

  double * u3, * o3, * e3;


  np3 = np * 3;
  np5 = np * 5;
  u3 = malloc (sizeof (double) * np3);
  o3 = malloc (sizeof (double) * np3);
  e3 = malloc (sizeof (double) * np5);
  if (u3 == NULL
      || o3 == NULL
      || e3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_2fts ().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;
      i5 = i * 5;

      u3 [i3 + 0] = u [i2 + 0];
      u3 [i3 + 1] = u [i2 + 1];
      u3 [i3 + 2] = 0.0;

      o3 [i3 + 0] = 0.0;
      o3 [i3 + 1] = 0.0;
      o3 [i3 + 2] = o [i];

      e3 [i5 + 0] = e [i2 + 0]; /* xx */
      e3 [i5 + 1] = e [i2 + 1]; /* xy */
      e3 [i5 + 2] = 0.0; /* xz */
      e3 [i5 + 3] = 0.0; /* yz */
      e3 [i5 + 4] = - e [i2 + 0]; /* yy */
    }
  calc_res_ewald_3fts (np, u3, o3, e3,
		       f, t, s);

  free (u3);
  free (o3);
  free (e3);
}

/** natural mobility problem **/
/* solve natural mobility problem in FTS version under Ewald sum
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [np * 3] : OK, this is 3D form
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
calc_mob_ewald_2fts (int np,
		     double *f, double *t3, double *e,
		     double *u, double *o, double *s)
{
  int i;
  int np3, np5;
  int i2, i3, i5;

  double * f3, * e3;


  np3 = np * 3;
  np5 = np * 5;
  f3 = malloc (sizeof (double) * np3);
  e3 = malloc (sizeof (double) * np5);
  if (f3 == NULL
      || e3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_2fts ().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;
      i5 = i * 5;

      f3 [i3 + 0] = f [i2 + 0];
      f3 [i3 + 1] = f [i2 + 1];
      f3 [i3 + 2] = 0.0;

      e3 [i5 + 0] = e [i2 + 0]; /* xx */
      e3 [i5 + 1] = e [i2 + 1]; /* xy */
      e3 [i5 + 2] = 0.0; /* xz */
      e3 [i5 + 3] = 0.0; /* yz */
      e3 [i5 + 4] = - e [i2 + 0]; /* yy */
    }
  calc_mob_ewald_2fts (np, f3, t3, e3,
		       u, o, s);

  free (f3);
  free (e3);
}

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FTS version
 * under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   e [nm * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   ef [nf * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   s [nm * 5] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 *   sf [nf * 5] :
 */
void
calc_mob_fix_ewald_2fts (int np, int nm,
			 double *f, double *t3, double *e,
			 double *uf, double *of, double *ef,
			 double *u, double *o, double *s,
			 double *ff, double *tf, double *sf)
{
  int i;
  int nm3, nm5;
  int nf, nf3, nf5;
  int i2, i3, i5;

  double * f3, * e3;
  double * uf3, * of3, * ef3;


  nm3 = nm * 3;
  nm5 = nm * 5;
  f3 = malloc (sizeof (double) * nm3);
  e3 = malloc (sizeof (double) * nm5);

  nf = np - nm;
  nf3 = nf * 3;
  nf5 = nf * 5;
  uf3 = malloc (sizeof (double) * nf3);
  of3 = malloc (sizeof (double) * nf3);
  ef3 = malloc (sizeof (double) * nf5);

  if (f3 == NULL
      || e3 == NULL
      || uf3 == NULL
      || of3 == NULL
      || ef3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_fix_ewald_2fts ().\n");
      exit (1);
    }

  for (i = 0; i < nm; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;
      i5 = i * 5;

      f3 [i3 + 0] = f [i2 + 0];
      f3 [i3 + 1] = f [i2 + 1];
      f3 [i3 + 2] = 0.0;

      e3 [i5 + 0] = e [i2 + 0]; /* xx */
      e3 [i5 + 1] = e [i2 + 1]; /* xy */
      e3 [i5 + 2] = 0.0; /* xz */
      e3 [i5 + 3] = 0.0; /* yz */
      e3 [i5 + 4] = - e [i2 + 0]; /* yy */
    }

  for (i = 0; i < nf; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;
      i5 = i * 5;

      uf3 [i3 + 0] = uf [i2 + 0];
      uf3 [i3 + 1] = uf [i2 + 1];
      uf3 [i3 + 2] = 0.0;

      of3 [i3 + 0] = 0.0;
      of3 [i3 + 1] = 0.0;
      of3 [i3 + 2] = of [i];

      ef3 [i5 + 0] = ef [i2 + 0]; /* xx */
      ef3 [i5 + 1] = ef [i2 + 1]; /* xy */
      ef3 [i5 + 2] = 0.0; /* xz */
      ef3 [i5 + 3] = 0.0; /* yz */
      ef3 [i5 + 4] = - ef [i2 + 0]; /* yy */
    }

  calc_mob_fix_ewald_3fts (np, nm,
			   f3, t3, e3, uf3, of3, ef3,
			   u, o, s, ff, tf, sf);

  free (f3);
  free (e3);
  free (uf3);
  free (of3);
  free (ef3);
}

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in FTS version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_lub_ewald_2fts (int np,
			 double *u, double *o, double *e,
			 double *f, double *t, double *s)
{
  int i;
  int np3, np5;
  int i2, i3, i5;

  double * u3, * o3, * e3;


  np3 = np * 3;
  np5 = np * 5;
  u3 = malloc (sizeof (double) * np3);
  o3 = malloc (sizeof (double) * np3);
  e3 = malloc (sizeof (double) * np5);
  if (u3 == NULL
      || o3 == NULL
      || e3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_2fts ().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;
      i5 = i * 5;

      u3 [i3 + 0] = u [i2 + 0];
      u3 [i3 + 1] = u [i2 + 1];
      u3 [i3 + 2] = 0.0;

      o3 [i3 + 0] = 0.0;
      o3 [i3 + 1] = 0.0;
      o3 [i3 + 2] = o [i];

      e3 [i5 + 0] = e [i2 + 0]; /* xx */
      e3 [i5 + 1] = e [i2 + 1]; /* xy */
      e3 [i5 + 2] = 0.0; /* xz */
      e3 [i5 + 3] = 0.0; /* yz */
      e3 [i5 + 4] = - e [i2 + 0]; /* yy */
    }
  calc_res_lub_ewald_3fts_2d (np, u3, o3, e3,
			      f, t, s);

  free (u3);
  free (o3);
  free (e3);
}

/* solve natural resistance problem with lubrication (2D!)
 * in FTS version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 3] :
 *   o [np * 3] :
 *   e [np * 5] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 *   s [np * 5] :
 */
static void
calc_res_lub_ewald_3fts_2d (int np,
			    double *u, double *o, double *e,
			    double *f, double *t, double *s)
{
  int i;
  int n11;

  double *b;
  double *x;
  double *lub;


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
  calc_lub_ewald_2fts (np, b, lub);
  atimes_ewald_3fts (n11, lub, x); // x[] is used temporaly
  for (i = 0; i < n11; ++i)
    b [i] += x [i];

  /* first guess */
  for (i = 0; i < n11; ++i)
    x [i] = 0.0;

  solve_iter_stab (n11, b, x, atimes_ewald_3fts,
		   gpb);

  set_FTS_by_fts (np, f, t, s, x);

  free (b);
  free (x);
  free (lub);
}


/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
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
calc_mob_lub_fix_ewald_2fts (int np, int nm,
			     double *f, double *t3, double *e,
			     double *uf, double *of, double *ef,
			     double *u, double *o, double *s,
			     double *ff, double *tf, double *sf)
{
  int i;
  int nm3, nm5;
  int nf, nf3, nf5;
  int i2, i3, i5;

  double * f3, * e3;
  double * uf3, * of3, * ef3;


  nm3 = nm * 3;
  nm5 = nm * 5;
  f3 = malloc (sizeof (double) * nm3);
  e3 = malloc (sizeof (double) * nm5);

  nf = np - nm;
  nf3 = nf * 3;
  nf5 = nf * 5;
  uf3 = malloc (sizeof (double) * nf3);
  of3 = malloc (sizeof (double) * nf3);
  ef3 = malloc (sizeof (double) * nf5);

  if (f3 == NULL
      || e3 == NULL
      || uf3 == NULL
      || of3 == NULL
      || ef3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_fix_ewald_2fts ().\n");
      exit (1);
    }

  for (i = 0; i < nm; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;
      i5 = i * 5;

      f3 [i3 + 0] = f [i2 + 0];
      f3 [i3 + 1] = f [i2 + 1];
      f3 [i3 + 2] = 0.0;

      e3 [i5 + 0] = e [i2 + 0]; /* xx */
      e3 [i5 + 1] = e [i2 + 1]; /* xy */
      e3 [i5 + 2] = 0.0; /* xz */
      e3 [i5 + 3] = 0.0; /* yz */
      e3 [i5 + 4] = - e [i2 + 0]; /* yy */
    }

  for (i = 0; i < nf; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;
      i5 = i * 5;

      uf3 [i3 + 0] = uf [i2 + 0];
      uf3 [i3 + 1] = uf [i2 + 1];
      uf3 [i3 + 2] = 0.0;

      of3 [i3 + 0] = 0.0;
      of3 [i3 + 1] = 0.0;
      of3 [i3 + 2] = of [i];

      ef3 [i5 + 0] = ef [i2 + 0]; /* xx */
      ef3 [i5 + 1] = ef [i2 + 1]; /* xy */
      ef3 [i5 + 2] = 0.0; /* xz */
      ef3 [i5 + 3] = 0.0; /* yz */
      ef3 [i5 + 4] = - ef [i2 + 0]; /* yy */
    }

  calc_mob_lub_fix_ewald_3fts_2d (np, nm,
				  f3, t3, e3, uf3, of3, ef3,
				  u, o, s, ff, tf, sf);

  free (f3);
  free (e3);
  free (uf3);
  free (of3);
  free (ef3);
}


/* solve natural mobility problem with lubrication (2D!)
 * with fixed particles in FTS version under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
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
static void
calc_mob_lub_fix_ewald_3fts_2d (int np, int nm,
				double *f, double *t, double *e,
				double *uf, double *of, double *ef,
				double *u, double *o, double *s,
				double *ff, double *tf, double *sf)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

  int i;
  int n11;
  int nf;
  int nm11;

  double *b;
  double *x;


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

  calc_b_mob_lub_fix_ewald_2fts (np, nm, f, t, e, uf, of, ef, b);

  /* first guess */
  for (i = 0; i < n11; ++i)
    x [i] = 0.0;

  NUMB_mobile_particles = nm;
  solve_iter_stab (n11, b, x, atimes_mob_lub_fix_ewald_2fts,
		   gpb);

  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  free (b);
  free (x);
}

/* calc b-term (constant term) of (natural) mobility problem
 * with lubrication under Ewald sum
 * where b := -(0,0,e) + M.(f,t,0).
 * INPUT
 *  np : # all particles (not # elements in b[]!)
 *  nm : # mobile particles (not # elements in b[]!)
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
calc_b_mob_lub_fix_ewald_2fts (int np, int nm,
			       double *f, double *t, double *e,
			       double *uf, double *of, double *ef,
			       double *b)
{
  int i;
  int nf;
  int n3, n5, n11;
  int nm11;

  double *x;
  double *y;
  double *v5_0;


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
      fprintf (stderr, "allocation error in calc_b_mob_lub_ewald_2fts ().\n");
      exit (1);
    }

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  /* set b :=  - (I + M.L).[(0,0,E)_m,(U,O,E)_f] */
  set_fts_by_FTS (nm, b, v5_0, v5_0, e);
  set_fts_by_FTS (nf, b + nm11, uf, of, ef);
  calc_lub_ewald_2fts (np, b, y);
  atimes_ewald_3fts (n11, y, x); // x[] is used temporaly
  for (i = 0; i < n11; ++i)
    {
      b [i] = - (b [i] + x [i]);
    }

  /* set x := [(F,T,0)_m,(0,0,0)_f] */
  set_fts_by_FTS (nm, x, f, t, v5_0);
  set_fts_by_FTS (nf, x + nm11, v5_0, v5_0, v5_0);
  atimes_ewald_3fts (n11, x, y);

  /* set b := - (I + M.L).[(0,0,E)_m,(U,O,E)_f] + [(F,T,0)_m,(0,0,0)_f] */
  for (i = 0; i < n11; ++i)
    {
      b [i] += y [i];
    }

  free (x);
  free (y);
  free (v5_0);
}

/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := [(u,o,0)_m,(0,0,0)_f] - M.[(0,0,s)_m,(f,t,s)_f].
 * INPUT
 *  (global) : NUMB_mobile_particles -- this is dirty, though ...
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_lub_fix_ewald_2fts (int n, double *x, double *y)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

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


  np = n / 11;
  nm = NUMB_mobile_particles;
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
      fprintf (stderr, "allocation error in atimes_mob_lub_ewald_2fts ().\n");
      exit (1);
    }

  for (i = 0; i < n5; ++i)
    {
      v5_0 [i] = 0.0;
    }

  /* set (U,O,S)_mobile,(F,T,S)_fixed by x[] */
  set_FTS_by_fts (nm, u, o, s, x);
  set_FTS_by_fts (nf, ff, tf, sf, x + nm11);

  /* set y := (I + M.L).[(U,O,0)_mobile,(0,0,0)_fixed] */
  set_fts_by_FTS (nm, y, u, o, v5_0);
  set_fts_by_FTS (nf, y + nm11, v5_0, v5_0, v5_0);
  calc_lub_ewald_2fts (np, y, w);
  atimes_ewald_3fts (n, w, z);
  for (i = 0; i < n; ++i)
    {
      y [i] += z [i];
    }

  /* set z := [(0,0,S)_mobile,(F,T,S)_fixed] */
  set_fts_by_FTS (nm, w, v5_0, v5_0, s);
  set_fts_by_FTS (nf, w + nm11, ff, tf, sf);
  atimes_ewald_3fts (n, w, z);

  /* set y := (I + M.L).[(U,O,0)_m,(0,0,0)_f] - M.[(0,0,S)_m,(F,T,S)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
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

/* calculate lubrication fts by uoe for all particles
 * under the periodic boundary condition
 * INPUT
 *   (global) pos [np * 3] : position of particles
 *   np : # particles
 *   uoe [np * 11] : velocity, angular velocity, strain
 * OUTPUT
 *   fts [np * 11] : force, torque, stresslet
 */
static void
calc_lub_ewald_2fts (int np, double * uoe, double * fts)
{
  extern double * pos;
  extern double llx [27], llz [27];

  int i, j, k;
  int i3, i11;
  int j3, j11;

  double * tmp_pos;


  tmp_pos = malloc (sizeof (double) * 3);
  if (tmp_pos == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_ewald_3fts().\n");
      exit (1);
    }

  /* clear fts [np * 11] */
  for (i = 0; i < np * 11; ++i)
    fts [i] = 0.0;

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i11 = i * 11;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  j11 = j * 11;
	  /* all image cells */
	  for (k = 0; k < 9; ++k)
	    {
	      tmp_pos [0] = pos [j3 + 0] + llx [k];
	      tmp_pos [1] = pos [j3 + 1] + llz [k];
	      tmp_pos [2] = pos [j3 + 2];
	      if (cond_lub_2d (pos + i3, tmp_pos) == 0)
		{
		  calc_lub_2b (uoe + i11, uoe + j11,
			       pos + i3, tmp_pos,
			       fts + i11, fts + j11);
		}
	    }
	}
    }

  free (tmp_pos);
}

/* condition for lubrication
 * INPUT
 *  x1 [3], x2 [3] : position
 * OUTPUT (return value)
 *  0 : r != 0 and r < 3.0
 *  1 : otherwise
 */
static int
cond_lub_2d (double * x1, double * x2)
{
  double x, y;
  double r2;


  x = x1 [0] - x2 [0];
  y = x1 [1] - x2 [1];

  r2 = x * x
    + y * y;

  if (r2 != 0.0
      && r2 < 9.0) // r = 3.0 is the critical separation for lubrication now.
    return 0;
  else
    return 1;
}
