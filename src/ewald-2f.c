/* Ewald summation technique under 2D
 * this is a wrapper package for ewald-3f.c
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-2f.c,v 1.1 2001/02/03 15:05:28 ichiki Exp $
 *
 * 3 dimensional hydrodynamics, 2D configuration
 * periodic boundary condition in 3 direction,
 * F version
 * non-dimension formulation
 */
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include "ewald-3f.h"

#include <libiter.h> /* gpb(), solve_iter_stab() */
#include "f.h"

#include "ewald-2f.h"


static void
calc_res_lub_ewald_3f_2d (int np,
			  double *u,
			  double *f);
static void
calc_mob_lub_fix_ewald_3f_2d (int np, int nm,
			      double *f,
			      double *uf,
			      double *u,
			      double *ff);
static void
calc_b_mob_lub_fix_ewald_2f (int np, int nm,
			     double *f,
			     double *uf,
			     double *b);
static void
atimes_mob_lub_fix_ewald_2f (int n, double *x, double *y);
static void
calc_lub_ewald_2f (int np, double * u, double * f);
static int
cond_lub_2d (double * x1, double * x2);


/** natural resistance problem **/
/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
calc_res_ewald_2f (int np,
		   double *u,
		   double *f)
{
  int i;
  int np3;
  int i2, i3;

  double * u3;


  np3 = np * 3;
  u3 = malloc (sizeof (double) * np3);
  if (u3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_2f ().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;

      u3 [i3 + 0] = u [i2 + 0];
      u3 [i3 + 1] = u [i2 + 1];
      u3 [i3 + 2] = 0.0;
    }
  calc_res_ewald_3f (np, u3,
		     f);

  free (u3);
}

/** natural mobility problem **/
/* solve natural mobility problem in F version under Ewald sum
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 */
void
calc_mob_ewald_2f (int np,
		   double *f,
		   double *u)
{
  int i;
  int np3;
  int i2, i3;

  double * f3;


  np3 = np * 3;
  f3 = malloc (sizeof (double) * np3);
  if (f3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_ewald_2f ().\n");
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
  calc_mob_ewald_2f (np, f3,
		     u);

  free (f3);
}

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in F version
 * under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
calc_mob_fix_ewald_2f (int np, int nm,
		       double *f,
		       double *uf,
		       double *u,
		       double *ff)
{
  int i;
  int nm3;
  int nf, nf3;
  int i2, i3;

  double * f3;
  double * uf3;


  nm3 = nm * 3;
  nf = np - nm;
  nf3 = nf * 3;

  f3 = malloc (sizeof (double) * nm3);
  uf3 = malloc (sizeof (double) * nf3);
  if (f3 == NULL
      || uf3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_fix_ewald_2f ().\n");
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
    }

  calc_mob_fix_ewald_3f (np, nm,
			 f3, uf3,
			 u, ff);

  free (f3);
  free (uf3);
}

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in F version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
calc_res_lub_ewald_2f (int np,
		       double *u,
		       double *f)
{
  int i;
  int np3;
  int i2, i3;

  double * u3;


  np3 = np * 3;
  u3 = malloc (sizeof (double) * np3);
  if (u3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_ewald_2f ().\n");
      exit (1);
    }

  for (i = 0; i < np; ++i)
    {
      i2 = i * 2;
      i3 = i * 3;

      u3 [i3 + 0] = u [i2 + 0];
      u3 [i3 + 1] = u [i2 + 1];
      u3 [i3 + 2] = 0.0;
    }
  calc_res_lub_ewald_3f_2d (np, u3,
			    f);

  free (u3);
}

/* solve natural resistance problem with lubrication (2D!)
 * in F version under Ewald sum
 * INPUT
 *  np : # particles
 *   u [np * 3] :
 * OUTPUT
 *   f [np * 3] :
 */
static void
calc_res_lub_ewald_3f_2d (int np,
			  double *u,
			  double *f)
{
  int i;
  int n3;

  double *b;
  double *x;
  double *lub;


  n3 = np * 3;
  b = malloc (sizeof (double) * n3);
  x = malloc (sizeof (double) * n3);
  lub = malloc (sizeof (double) * n3);
  if (b == NULL
      || x == NULL
      || lub == NULL)
    {
      fprintf (stderr, "allocation error in calc_res_lub_ewald_3f ().\n");
      exit (1);
    }

  set_F_by_f (np, b, u);
  calc_lub_ewald_2f (np, b, lub);
  atimes_ewald_3f (n3, lub, x); // x[] is used temporaly
  for (i = 0; i < n3; ++i)
    b [i] += x [i];

  /* first guess */
  for (i = 0; i < n3; ++i)
    x [i] = 0.0;

  solve_iter_stab (n3, b, x, atimes_ewald_3f,
		   gpb);

  set_F_by_f (np, f, x);

  free (b);
  free (x);
  free (lub);
}


/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in F version under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
calc_mob_lub_fix_ewald_2f (int np, int nm,
			   double *f,
			   double *uf,
			   double *u,
			   double *ff)
{
  int i;
  int nm3;
  int nf, nf3;
  int i2, i3;

  double * f3;
  double * uf3;


  nm3 = nm * 3;
  nf = np - nm;
  nf3 = nf * 3;

  f3 = malloc (sizeof (double) * nm3);
  uf3 = malloc (sizeof (double) * nf3);
  if (f3 == NULL
      || uf3 == NULL)
    {
      fprintf (stderr, "allocation error in calc_mob_fix_ewald_2f ().\n");
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
    }

  calc_mob_lub_fix_ewald_3f_2d (np, nm,
				f3, uf3,
				u, ff);

  free (f3);
  free (uf3);
}


/* solve natural mobility problem with lubrication (2D!)
 * with fixed particles in F version under Ewald sum
 * INPUT
 *  np : # all particles
 *  nm : # mobile particles, so that (np - nm) is # fixed particles
 *   f [nm * 3] :
 *   uf [nf * 3] :
 * OUTPUT
 *   u [nm * 3] :
 *   ff [nf * 3] :
 */
static void
calc_mob_lub_fix_ewald_3f_2d (int np, int nm,
			      double *f,
			      double *uf,
			      double *u,
			      double *ff)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

  int i;
  int n3;
  int nf;
  int nm3;

  double *b;
  double *x;


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

  calc_b_mob_lub_fix_ewald_2f (np, nm, f, uf, b);

  /* first guess */
  for (i = 0; i < n3; ++i)
    x [i] = 0.0;

  NUMB_mobile_particles = nm;
  solve_iter_stab (n3, b, x, atimes_mob_lub_fix_ewald_2f,
		   gpb);

  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

  free (b);
  free (x);
}

/* calc b-term (constant term) of (natural) mobility problem
 * with lubrication under Ewald sum
 * where b := -(0) + M.(f).
 * INPUT
 *  np : # all particles (not # elements in b[]!)
 *  nm : # mobile particles (not # elements in b[]!)
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *  b [np * 3] : constant vector
 */
static void
calc_b_mob_lub_fix_ewald_2f (int np, int nm,
			     double *f,
			     double *uf,
			     double *b)
{
  int i;
  int nf;
  int n3;
  int nm3;

  double *x;
  double *y;
  double *v3_0;


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
      fprintf (stderr, "allocation error in calc_b_mob_lub_ewald_2f ().\n");
      exit (1);
    }

  for (i = 0; i < n3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set b :=  - (I + M.L).[(0,0)_m,(U,O)_f] */
  set_F_by_f (nm, b, v3_0);
  set_F_by_f (nf, b + nm3, uf);
  calc_lub_ewald_2f (np, b, y);
  atimes_ewald_3f (n3, y, x); // x[] is used temporaly
  for (i = 0; i < n3; ++i)
    {
      b [i] = - (b [i] + x [i]);
    }

  /* set x := [(F)_m,(0)_f] */
  set_F_by_f (nm, x, f);
  set_F_by_f (nf, x + nm3, v3_0);
  atimes_ewald_3f (n3, x, y);

  /* set b := - (I + M.L).[(0)_m,(U)_f] + [(F)_m,(0)_f] */
  for (i = 0; i < n3; ++i)
    {
      b [i] += y [i];
    }

  free (x);
  free (y);
  free (v3_0);
}

/* calc atimes of (natural) mobility problem under Ewald sum
 * where A.x := [(u)_m,(0)_f] - M.[(0)_m,(f)_f].
 * INPUT
 *  (global) : NUMB_mobile_particles -- this is dirty, though ...
 *  n : # elements in x[] and b[] (not # particles!)
 *  x [n] :
 * OUTPUT
 *  y [n] :
 */
static void
atimes_mob_lub_fix_ewald_2f (int n, double *x, double *y)
{
  extern int NUMB_mobile_particles; /* this is dirty, though ... */

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


  np = n / 3;
  np3 = np * 3;
  nm = NUMB_mobile_particles;
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
      fprintf (stderr, "allocation error in atimes_mob_lub_ewald_2f ().\n");
      exit (1);
    }

  for (i = 0; i < np3; ++i)
    {
      v3_0 [i] = 0.0;
    }

  /* set (U,O)_mobile,(F,T)_fixed by x[] */
  set_F_by_f (nm, u, x);
  set_F_by_f (nf, ff, x + nm3);

  /* set y := (I + M.L).[(U,O)_mobile,(0,0)_fixed] */
  set_F_by_f (nm, y, u);
  set_F_by_f (nf, y + nm3, v3_0);
  calc_lub_ewald_2f (np, y, w);
  atimes_ewald_3f (n, w, z);
  for (i = 0; i < n; ++i)
    {
      y [i] += z [i];
    }

  /* set z := [(0)_mobile,(F)_fixed] */
  set_F_by_f (nm, w, v3_0);
  set_F_by_f (nf, w + nm3, ff);
  atimes_ewald_3f (n, w, z);

  /* set y := (I + M.L).[(U)_m,(0)_f] - M.[(0)_m,(F)_f] */
  for (i = 0; i < n; ++i)
    {
      y [i] -= z [i];
    }

  free (w);
  free (z);
  free (v3_0);
  free (u);
  free (ff);
}

/* calculate lubrication f by uo for all particles
 * under the periodic boundary condition
 * INPUT
 *   (global) pos [np * 3] : position of particles
 *   np : # particles
 *   u [np * 3] : velocity
 * OUTPUT
 *   f [np * 3] : force
 */
static void
calc_lub_ewald_2f (int np, double * u, double * f)
{
  extern double * pos;
  extern double llx [27], lly [27];

  int i, j, k;
  int i3;
  int j3;

  double * tmp_pos;


  tmp_pos = malloc (sizeof (double) * 3);
  if (tmp_pos == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_ewald_3f().\n");
      exit (1);
    }

  /* clear f [np * 3] */
  for (i = 0; i < np * 3; ++i)
    f [i] = 0.0;

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  j3 = j * 3;
	  /* all image cells */
	  for (k = 0; k < 9; ++k)
	    {
	      tmp_pos [0] = pos [j3 + 0] + llx [k];
	      tmp_pos [1] = pos [j3 + 1] + lly [k];
	      tmp_pos [2] = pos [j3 + 2];
	      if (cond_lub_2d (pos + i3, tmp_pos) == 0)
		{
		  calc_lub_f_2b (u + i3, u + j3,
				 pos + i3, tmp_pos,
				 f + i3, f + j3);
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
