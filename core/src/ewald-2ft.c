/* Solvers for 2 dimensional (monolayer) FT version problems
 * this is a wrapper package for ewald-3fts.c
 * Copyright (C) 2001-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-2ft.c,v 1.6 2007/03/07 22:25:23 kichiki Exp $
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
#include "stokes.h" /* struct stokeks */
#include "ewald-3ft.h"

#include "ewald-2ft.h"


/** natural resistance problem **/
/* solve natural resistance problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
solve_res_2ft (struct stokes * sys,
	       const double *u, const double *o,
	       double *f, double *t)
{
  int np;
  int i;
  int np3;
  int i2, i3;

  double * u3, * o3;


  np = sys->np;
  np3 = np * 3;
  u3 = (double *) malloc (sizeof (double) * np3);
  o3 = (double *) malloc (sizeof (double) * np3);
  if (u3 == NULL
      || o3 == NULL)
    {
      fprintf (stderr, "allocation error in solve_res_2ft ().\n");
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
  solve_res_3ft (sys, u3, o3, f, t);

  free (u3);
  free (o3);
}

/** natural mobility problem **/
/* solve natural mobility problem in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [np * 3] : OK, this is 3D form
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 *   o [np * 3] :
 */
void
solve_mob_2ft (struct stokes * sys,
	       const double *f, const double *t3,
	       double *u, double *o)
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
      fprintf (stderr, "allocation error in solve_mob_2ft ().\n");
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
  solve_mob_3ft (sys, f3, t3, u, o);

  free (f3);
}

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
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
solve_mix_2ft (struct stokes * sys,
	       const double *f, const double *t3,
	       const double *uf, const double *of,
	       double *u, double *o,
	       double *ff, double *tf)
{
  int np, nm;
  int i;
  int nm3;
  int nf, nf3;
  int i2, i3;

  double * f3;
  double * uf3, * of3;


  np = sys->np;
  nm = sys->nm;

  nm3 = nm * 3;
  nf = np - nm;
  nf3 = nf * 3;

  f3  = (double *) malloc (sizeof (double) * nm3);
  uf3 = (double *) malloc (sizeof (double) * nf3);
  of3 = (double *) malloc (sizeof (double) * nf3);
  if (f3 == NULL
      || uf3 == NULL
      || of3 == NULL)
    {
      fprintf (stderr, "allocation error in solve_mix_2ft ().\n");
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

  solve_mix_3ft (sys,
		 f3, t3, uf3, of3,
		 u, o, ff, tf);

  free (f3);
  free (uf3);
  free (of3);
}

/** natural resistance problem with lubrication **/
/* solve natural resistance problem with lubrication in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 */
void
solve_res_lub_2ft (struct stokes * sys,
		   const double *u, const double *o,
		   double *f, double *t)
{
  int np;
  int i;
  int np3;
  int i2, i3;

  double * u3, * o3;


  np = sys->np;
  np3 = np * 3;
  u3 = (double *) malloc (sizeof (double) * np3);
  o3 = (double *) malloc (sizeof (double) * np3);
  if (u3 == NULL
      || o3 == NULL)
    {
      fprintf (stderr, "allocation error in solve_res_2ft ().\n");
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
  solve_res_lub_3ft (sys, u3, o3, f, t);

  free (u3);
  free (o3);
}

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FT version
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
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
solve_mix_lub_2ft (struct stokes * sys,
		   const double *f, const double *t3,
		   const double *uf, const double *of,
		   double *u, double *o,
		   double *ff, double *tf)
{
  int np, nm;
  int i;
  int nm3;
  int nf, nf3;
  int i2, i3;

  double * f3;
  double * uf3, * of3;


  np = sys->np;
  nm = sys->nm;

  nm3 = nm * 3;
  nf = np - nm;
  nf3 = nf * 3;

  f3  = (double *) malloc (sizeof (double) * nm3);
  uf3 = (double *) malloc (sizeof (double) * nf3);
  of3 = (double *) malloc (sizeof (double) * nf3);
  if (f3 == NULL
      || uf3 == NULL
      || of3 == NULL)
    {
      fprintf (stderr, "allocation error in solve_mix_2ft ().\n");
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

  solve_mix_lub_3ft (sys,
		     f3, t3, uf3, of3,
		     u, o, ff, tf);

  free (f3);
  free (uf3);
  free (of3);
}
