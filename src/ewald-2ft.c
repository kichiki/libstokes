/* Ewald summation technique under 2D
 * this is a wrapper package for ewald-3fts.c
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-2ft.c,v 1.2 2001/02/04 12:19:15 ichiki Exp $
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
  calc_mob_ewald_3ft (np, f3, t3,
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
  calc_res_lub_ewald_3ft (np, u3, o3,
			  f, t);

  free (u3);
  free (o3);
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

  calc_mob_lub_fix_ewald_3ft (np, nm,
			      f3, t3, uf3, of3,
			      u, o, ff, tf);

  free (f3);
  free (uf3);
  free (of3);
}
