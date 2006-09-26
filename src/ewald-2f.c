/* Ewald summation technique under 2D
 * this is a wrapper package for ewald-3f.c
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-2f.c,v 1.4 2006/09/26 05:38:58 ichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 2D configuration
 * periodic boundary condition in 3 direction
 * F version
 * non-dimension formulation
 */
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include "stokes.h" /* struct stokeks */
#include "ewald-3f.h"

#include "ewald-2f.h"


/** natural resistance problem **/
/* solve natural resistance problem in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
calc_res_ewald_2f (struct stokes * sys,
		   const double *u,
		   double *f)
{
  int np;
  int i;
  int np3;
  int i2, i3;

  double * u3;


  np = sys->np;
  np3 = np * 3;
  u3 = (double *) malloc (sizeof (double) * np3);
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
  calc_res_ewald_3f (sys, u3, f);

  free (u3);
}

/** natural mobility problem **/
/* solve natural mobility problem in F version under Ewald sum
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 */
void
calc_mob_ewald_2f (struct stokes * sys,
		   const double *f,
		   double *u)
{
  int np;
  int i;
  int np3;
  int i2, i3;

  double * f3;


  np = sys->np;
  np3 = np * 3;
  f3 = (double *) malloc (sizeof (double) * np3);
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
  calc_mob_ewald_3f (sys, f3, u);

  free (f3);
}

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in F version
 * under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
calc_mob_fix_ewald_2f (struct stokes * sys,
		       const double *f,
		       const double *uf,
		       double *u,
		       double *ff)
{
  int np, nm;
  int i;
  int nm3;
  int nf, nf3;
  int i2, i3;

  double * f3;
  double * uf3;


  np = sys->np;
  nm = sys->nm;

  nm3 = nm * 3;
  nf = np - nm;
  nf3 = nf * 3;

  f3  = (double *) malloc (sizeof (double) * nm3);
  uf3 = (double *) malloc (sizeof (double) * nf3);
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

  calc_mob_fix_ewald_3f (sys, f3, uf3, u, ff);

  free (f3);
  free (uf3);
}

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication
 * in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 */
void
calc_res_lub_ewald_2f (struct stokes * sys,
		       const double *u,
		       double *f)
{
  int np;
  int i;
  int np3;
  int i2, i3;

  double * u3;


  np = sys->np;
  np3 = np * 3;
  u3 = (double *) malloc (sizeof (double) * np3);
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
  calc_res_lub_ewald_3f (sys, u3, f);

  free (u3);
}

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in F version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   f [nm * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   uf [nf * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 * OUTPUT
 *   u [nm * 3] : results are given in 3D form
 *   ff [nf * 3] :
 */
void
calc_mob_lub_fix_ewald_2f (struct stokes * sys,
			   const double *f,
			   const double *uf,
			   double *u,
			   double *ff)
{
  int np, nm;
  int i;
  int nm3;
  int nf, nf3;
  int i2, i3;

  double * f3;
  double * uf3;


  np = sys->np;
  nm = sys->nm;

  nm3 = nm * 3;
  nf = np - nm;
  nf3 = nf * 3;

  f3  = (double *) malloc (sizeof (double) * nm3);
  uf3 = (double *) malloc (sizeof (double) * nf3);
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

  calc_mob_lub_fix_ewald_3f (sys, f3, uf3, u, ff);

  free (f3);
  free (uf3);
}
