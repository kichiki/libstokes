/* Ewald summation technique under 2D
 * this is a wrapper package for ewald-3fts.c
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-2ft.c,v 1.1 2001/02/03 14:48:36 ichiki Exp $
 *
 * 3 dimensional hydrodynamics, 2D configuration
 * periodic boundary condition in 3 direction,
 * FT version
 * non-dimension formulation
 */
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include "ewald-3ft.h"

#include <libiter.h> /* gpb(), solve_iter_stab() */
#include "ft.h"

#include "ewald-2ft.h"


static void
calc_res_lub_ewald_3ft_2d (int np,
			   double *u, double *o,
			   double *f, double *t);

static void
calc_mob_lub_fix_ewald_3ft_2d (int np, int nm,
			       double *f, double *t,
			       double *uf, double *of,
			       double *u, double *o,
			       double *ff, double *tf);
static void
calc_b_mob_lub_fix_ewald_2ft (int np, int nm,
			      double *f, double *t,
			      double *uf, double *of,
			      double *b);
static void
atimes_mob_lub_fix_ewald_2ft (int n, double *x, double *y);
static void
calc_lub_ewald_2ft (int np, double * uo, double * ft);
static int
cond_lub_2d (double * x1, double * x2);


/** natural resistance problem **/
/* solve natural resistance problem in FT version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
calc_res_ewald_2ft (int np,
		    double *u, double *o,
		    double *f, double *t)
{
  int i;
  int np3;
  int i2, i3;

  double * u3, * o3;


  np3 = np * 3;
  u3 = malloc (sizeof (double) * np3);
  o3 = malloc (sizeof (double) * np3);
  if (u3 == NULL
      || o3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_2ft ().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;

      u3 [i3 + 0] = u [i2 + 0];
      u3 [i3 + 1] = u [i2 + 1];
      u3 [i3 + 2] = 0.0;

      o3 [i3 + 0] = 0.0;
      o3 [i3 + 1] = 0.0;
      o3 [i3 + 2] = o [i];
    }
  calc_res_ewald_3ft (np, u3, o3,
		      f, t);

  free (u3);
  free (o3);
}

/** natural mobility problem **/
/* solve natural mobility problem in FT version under Ewald sum
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [np * 3] : OK, this is 3D form
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 *   o [np * 3] :
 */
void
calc_mob_ewald_2ft (int np,
		    double *f, double *t3,
		    double *u, double *o)
{
  int i;
  int np3;
  int i2, i3;

  double * f3;


  np3 = np * 3;
  f3 = malloc (sizeof (double) * np3);
  if (f3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_2ft ().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;

      f3 [i3 + 0] = f [i2 + 0];
      f3 [i3 + 1] = f [i2 + 1];
      f3 [i3 + 2] = 0.0;
    }
  calc_mob_ewald_2ft (np, f3, t3,
		      u, o);

  free (f3);
}

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FT version
 * under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
calc_mob_fix_ewald_2ft (int np, int nm,
			double *f, double *t3,
			double *uf, double *of,
			double *u, double *o,
			double *ff, double *tf)
{
  int i;
  int nm3;
  int nf, nf3;
  int i2, i3;

  double * f3;
  double * uf3, * of3;


  nm3 = nm * 3;
  nf = np - nm;
  nf3 = nf * 3;

  f3 = malloc (sizeof (double) * nm3);
  uf3 = malloc (sizeof (double) * nf3);
  of3 = malloc (sizeof (double) * nf3);
  if (f3 == NULL
      || uf3 == NULL
      || of3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_fix_ewald_2ft ().\n");
      exit (1);
    }

  for (i = 0; i < nm; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;

      f3 [i3 + 0] = f [i2 + 0];
      f3 [i3 + 1] = f [i2 + 1];
      f3 [i3 + 2] = 0.0;
    }

  for (i = 0; i < nf; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;

      uf3 [i3 + 0] = uf [i2 + 0];
      uf3 [i3 + 1] = uf [i2 + 1];
      uf3 [i3 + 2] = 0.0;

      of3 [i3 + 0] = 0.0;
      of3 [i3 + 1] = 0.0;
      of3 [i3 + 2] = of [i];
    }

  calc_mob_fix_ewald_3ft (np, nm,
			  f3, t3, uf3, of3,
			  u, o, ff, tf);

  free (f3);
  free (uf3);
  free (of3);
}

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in FT version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
calc_res_lub_ewald_2ft (int np,
			double *u, double *o,
			double *f, double *t)
{
  int i;
  int np3;
  int i2, i3;

  double * u3, * o3;


  np3 = np * 3;
  u3 = malloc (sizeof (double) * np3);
  o3 = malloc (sizeof (double) * np3);
  if (u3 == NULL
      || o3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_2ft ().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;

      u3 [i3 + 0] = u [i2 + 0];
      u3 [i3 + 1] = u [i2 + 1];
      u3 [i3 + 2] = 0.0;

      o3 [i3 + 0] = 0.0;
      o3 [i3 + 1] = 0.0;
      o3 [i3 + 2] = o [i];
    }
  calc_res_lub_ewald_3ft_2d (np, u3, o3,
			     f, t);

  free (u3);
  free (o3);
}

/* solve natural resistance problem with lubrication (2D!)
 * in FT version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 3] :
 *   o [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 *   t [np * 3] :
 */
static void
calc_res_lub_ewald_3ft_2d (int np,
			   double *u, double *o,
			   double *f, double *t)
{
  int i;
  int n6;

  double *b;
  double *x;
  double *lub;


  n6 = np * 6;
  b = malloc (sizeof (double) * n6);
  x = malloc (sizeof (double) * n6);
  lub = malloc (sizeof (double) * n6);
  if (b == NULL
      || x == NULL
      || lub == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_lub_ewald_3ft ().\n");
      exit (1);
    }

  set_ft_by_FT (np, b, u, o);
  calc_lub_ewald_2ft (np, b, lub);
  atimes_ewald_3ft (n6, lub, x); // x[] is used temporaly
  for (i = 0; i < n6; ++i)
    b [i] += x [i];

  /* first guess */
  for (i = 0; i < n6; ++i)
    x [i] = 0.0;

  solve_iter_stab (n6, b, x, atimes_ewald_3ft,
		   gpb);

  set_FT_by_ft (np, f, t, x);

  free (b);
  free (x);
  free (lub);
}


/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [nm * 3] : OK, this is 3D form
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   of [nf * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   o [nm * 3] :
 *   ff [nf * 3] :
 *   tf [nf * 3] :
 */
void
calc_mob_lub_fix_ewald_2ft (int np, int nm,
			    double *f, double *t3,
			    double *uf, double *of,
			    double *u, double *o,
			    double *ff, double *tf)
{
  int i;
  int nm3;
  int nf, nf3;
  int i2, i3;

  double * f3;
  double * uf3, * of3;


  nm3 = nm * 3;
  nf = np - nm;
  nf3 = nf * 3;

  f3 = malloc (sizeof (double) * nm3);
  uf3 = malloc (sizeof (double) * nf3);
  of3 = malloc (sizeof (double) * nf3);
  if (f3 == NULL
      || uf3 == NULL
      || of3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_fix_ewald_2ft ().\n");
      exit (1);
    }

  for (i = 0; i < nm; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;

      f3 [i3 + 0] = f [i2 + 0];
      f3 [i3 + 1] = f [i2 + 1];
      f3 [i3 + 2] = 0.0;
    }

  for (i = 0; i < nf; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;

      uf3 [i3 + 0] = uf [i2 + 0];
      uf3 [i3 + 1] = uf [i2 + 1];
      uf3 [i3 + 2] = 0.0;

      of3 [i3 + 0] = 0.0;
      of3 [i3 + 1] = 0.0;
      of3 [i3 + 2] = of [i];
    }

  calc_mob_lub_fix_ewald_3ft_2d (np, nm,
				 f3, t3, uf3, of3,
				 u, o, ff, tf);

  free (f3);
  free (uf3);
  free (of3);
}


/* solve natural mobility problem with lubrication (2D!)
 * with fixed particles in FT version under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
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
static void
calc_mob_lub_fix_ewald_3ft_2d (int np, int nm,
			       double *f, double *t,
			       double *uf, double *of,
			       double *u, double *o,
			       double *ff, double *tf)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

  int i;
  int n6;
  int nf;
  int nm6;

  double *b;
  double *x;


  nf = np - nm;
  n6 = np * 6;
  nm6 = nm * 6;

  b = malloc (sizeof (double) * n6);
  x = malloc (sizeof (double) * n6);
  if (b == NULL
      || x == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_3ft ().\n");
      exit (1);
    }

  calc_b_mob_lub_fix_ewald_2ft (np, nm, f, t, uf, of, b);

  /* first guess */
  for (i = 0; i < n6; ++i)
    x [i] = 0.0;

  NUMB_mobile_particles = nm;
  solve_iter_stab (n6, b, x, atimes_mob_lub_fix_ewald_2ft,
		   gpb);

  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  free (b);
  free (x);
}

/* calc b-term (constant term) of (natural) mobility problem
 * with lubrication under Ewald sum
 * where b := -(0,0) + M.(f,t).
 * INPUT
 *  np : # all particles (not # elements in b[]!)
 *  nm : # mobile particles (not # elements in b[]!)
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  uf [nf * 3] :
 *  of [nf * 3] :
 * OUTPUT
 *  b [np * 6] : constant vector
 */
static void
calc_b_mob_lub_fix_ewald_2ft (int np, int nm,
			      double *f, double *t,
			      double *uf, double *of,
			      double *b)
{
  int i;
  int nf;
  int n3, n6;
  int nm6;

  double *x;
  double *y;
  double *v3_0;


  nf = np - nm;
  n3 = np * 3;
  n6 = np * 6;
  nm6 = nm * 6;

  x = malloc (sizeof (double) * n6);
  y = malloc (sizeof (double) * n6);
  v3_0 = malloc (sizeof (double) * n3);
  if (x == NULL
      || y == NULL
      || v3_0 == NULL)
    {
      fprintf (stderr, "allocation error in calc_b_mob_lub_ewald_2ft ().\n");
      exit (1);
    }

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set b :=  - (I + M.L).[(0,0)_m,(U,O)_f] */
  set_ft_by_FT (nm, b, v3_0, v3_0);
  set_ft_by_FT (nf, b + nm6, uf, of);
  calc_lub_ewald_2ft (np, b, y);
  atimes_ewald_3ft (n6, y, x); // x[] is used temporaly
  for (i = 0; i < n6; ++i)
    {
      b [i] = - (b [i] + x [i]);
    }

  /* set x := [(F,T)_m,(0,0)_f] */
  set_ft_by_FT (nm, x, f, t);
  set_ft_by_FT (nf, x + nm6, v3_0, v3_0);
  atimes_ewald_3ft (n6, x, y);

  /* set b := - (I + M.L).[(0,0)_m,(U,O)_f] + [(F,T)_m,(0,0)_f] */
  for (i = 0; i < n6; ++i)
    {
      b [i] += y [i];
    }

  free (x);
  free (y);
  free (v3_0);
}

/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := [(u,o)_m,(0,0)_f] - M.[(0,0)_m,(f,t)_f].
 * INPUT
 *  (global) : NUMB_mobile_particles -- this is dirty, though ...
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_lub_fix_ewald_2ft (int n, double *x, double *y)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

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


  np = n / 6;
  np3 = np * 3;
  nm = NUMB_mobile_particles;
  nf = np - nm;
  nf3 = nf * 3;
  nm3 = nm * 3;
  nm6 = nm * 6;

  w = malloc (sizeof (double) * n);
  z = malloc (sizeof (double) * n);
  v3_0 = malloc (sizeof (double) * np3);
  u = malloc (sizeof (double) * nm3);
  o = malloc (sizeof (double) * nm3);
  ff = malloc (sizeof (double) * nf3);
  tf = malloc (sizeof (double) * nf3);
  if (z == NULL
      || v3_0 == NULL
      || u == NULL
      || o == NULL
      || ff == NULL
      || tf == NULL)
    {
      fprintf (stderr, "allocation error in atimes_mob_lub_ewald_2ft ().\n");
      exit (1);
    }

  for (i = 0; i < np3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set (U,O)_mobile,(F,T)_fixed by x[] */
  set_FT_by_ft (nm, u, o, x);
  set_FT_by_ft (nf, ff, tf, x + nm6);

  /* set y := (I + M.L).[(U,O)_mobile,(0,0)_fixed] */
  set_ft_by_FT (nm, y, u, o);
  set_ft_by_FT (nf, y + nm6, v3_0, v3_0);
  calc_lub_ewald_2ft (np, y, w);
  atimes_ewald_3ft (n, w, z);
  for (i = 0; i < n; ++i)
    {
      y [i] += z [i];
    }

  /* set z := [(0,0)_mobile,(F,T)_fixed] */
  set_ft_by_FT (nm, w, v3_0, v3_0);
  set_ft_by_FT (nf, w + nm6, ff, tf);
  atimes_ewald_3ft (n, w, z);

  /* set y := (I + M.L).[(U,O)_m,(0,0)_f] - M.[(0,0)_m,(F,T)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
    }

  free (w);
  free (z);
  free (v3_0);
  free (u);
  free (o);
  free (ff);
  free (tf);
}

/* calculate lubrication ft by uo for all particles
 * under the periodic boundary condition
 * INPUT
 *   (global) pos [np * 3] : position of particles
 *   np : # particles
 *   uo [np * 6] : velocity, angular velocity, strain
 * OUTPUT
 *   ft [np * 6] : force, torque, stresslet
 */
static void
calc_lub_ewald_2ft (int np, double * uo, double * ft)
{
  extern double * pos;
  extern double llx [27], lly [27];

  int i, j, k;
  int i3, i6;
  int j3, j6;

  double * tmp_pos;


  tmp_pos = malloc (sizeof (double) * 3);
  if (tmp_pos == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_ewald_3ft().\n");
      exit (1);
    }

  /* clear ft [np * 6] */
  for (i = 0; i < np * 6; ++i)
    ft [i] = 0.0;

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i6 = i * 6;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  j6 = j * 6;
	  /* all image cells */
	  for (k = 0; k < 9; ++k)
	    {
	      tmp_pos [0] = pos [j3 + 0] + llx [k];
	      tmp_pos [1] = pos [j3 + 1] + lly [k];
	      tmp_pos [2] = pos [j3 + 2];
	      if (cond_lub_2d (pos + i3, tmp_pos) == 0)
		{
		  calc_lub_ft_2b (uo + i6, uo + j6,
				  pos + i3, tmp_pos,
				  ft + i6, ft + j6);
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
