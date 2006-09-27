/* Ewald summation technique under 2D
 * this is a wrapper package for ewald-3fts.c
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald-2fts.c,v 1.5 2006/09/27 00:03:29 ichiki Exp $
 *
 * 3 dimensional hydrodynamics
 * 2D configuration
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
#include "stokes.h" /* struct stokeks */
#include "ewald-3fts.h"

#include "ewald-2fts.h"


/** natural resistance problem **/
/* solve natural resistance problem in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_ewald_2fts (struct stokes * sys,
		     const double *u, const double *o, const double *e,
		     double *f, double *t, double *s)
{
  int np;
  int i;
  int np3, np5;
  int i2, i3, i5;

  double * u3, * o3, * e3;


  np = sys->np;
  np3 = np * 3;
  np5 = np * 5;
  u3 = (double *) malloc (sizeof (double) * np3);
  o3 = (double *) malloc (sizeof (double) * np3);
  e3 = (double *) malloc (sizeof (double) * np5);
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
  calc_res_ewald_3fts (sys, u3, o3, e3,
		       f, t, s);

  free (u3);
  free (o3);
  free (e3);
}

/** natural mobility problem **/
/* solve natural mobility problem in FTS version under Ewald sum
 *  sys : system parameters
 *   f [np * 2] : f_x, f_y are given (and f_z = 0 is assumed)
 *   t3 [np * 3] : OK, this is 3D form
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   u [np * 3] : results are given in 3D form
 *   o [np * 3] :
 *   s [np * 5] :
 */
void
calc_mob_ewald_2fts (struct stokes * sys,
		     const double *f, const double *t3, const double *e,
		     double *u, double *o, double *s)
{
  int np;
  int i;
  int np3, np5;
  int i2, i3, i5;

  double * f3, * e3;


  np = sys->np;
  np3 = np * 3;
  np5 = np * 5;
  f3 = (double *) malloc (sizeof (double) * np3);
  e3 = (double *) malloc (sizeof (double) * np5);
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
  calc_mob_ewald_3fts (sys, f3, t3, e3,
		       u, o, s);

  free (f3);
  free (e3);
}

/** natural mobility problem with fixed particles **/
/* solve natural mobility problem with fixed particles in FTS version
 * under Ewald sum
 * INPUT
 *  sys : system parameters
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
calc_mob_fix_ewald_2fts (struct stokes * sys,
			 const double *f, const double *t3, const double *e,
			 const double *uf, const double *of, const double *ef,
			 double *u, double *o, double *s,
			 double *ff, double *tf, double *sf)
{
  int np, nm;
  int i;
  int nm3, nm5;
  int nf, nf3, nf5;
  int i2, i3, i5;

  double * f3, * e3;
  double * uf3, * of3, * ef3;


  np = sys->np;
  nm = sys->nm;

  nm3 = nm * 3;
  nm5 = nm * 5;
  f3 = (double *) malloc (sizeof (double) * nm3);
  e3 = (double *) malloc (sizeof (double) * nm5);

  nf = np - nm;
  nf3 = nf * 3;
  nf5 = nf * 5;
  uf3 = (double *) malloc (sizeof (double) * nf3);
  of3 = (double *) malloc (sizeof (double) * nf3);
  ef3 = (double *) malloc (sizeof (double) * nf5);

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

  calc_mob_fix_ewald_3fts (sys,
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
 *  sys : system parameters
 *   u [np * 2] : u_x, u_y are given (and u_z = 0 is assumed)
 *   o [np * 1] : o_z is given (and o_x = o_y = 0 is assumed)
 *   e [np * 2] : e_xx, e_xy are given (e_?z = e_z? = 0 is assumed)
 * OUTPUT
 *   f [np * 3] : results are given in 3D form
 *   t [np * 3] :
 *   s [np * 5] :
 */
void
calc_res_lub_ewald_2fts (struct stokes * sys,
			 const double *u, const double *o, const double *e,
			 double *f, double *t, double *s)
{
  int np;
  int i;
  int np3, np5;
  int i2, i3, i5;

  double * u3, * o3, * e3;


  np = sys->np;
  np3 = np * 3;
  np5 = np * 5;
  u3 = (double *) malloc (sizeof (double) * np3);
  o3 = (double *) malloc (sizeof (double) * np3);
  e3 = (double *) malloc (sizeof (double) * np5);
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
  calc_res_lub_ewald_3fts (sys, u3, o3, e3,
			   f, t, s);

  free (u3);
  free (o3);
  free (e3);
}

/** natural mobility problem with lubrication with fixed particles **/
/* solve natural mobility problem with lubrication
 * with fixed particles in FTS version under Ewald sum
 * INPUT
 *  sys : system parameters
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
calc_mob_lub_fix_ewald_2fts (struct stokes * sys,
			     const double *f, const double *t3,
			     const double *e,
			     const double *uf, const double *of,
			     const double *ef,
			     double *u, double *o, double *s,
			     double *ff, double *tf, double *sf)
{
  int np, nm;
  int i;
  int nm3, nm5;
  int nf, nf3, nf5;
  int i2, i3, i5;

  double * f3, * e3;
  double * uf3, * of3, * ef3;


  np = sys->np;
  nm = sys->nm;

  nm3 = nm * 3;
  nm5 = nm * 5;
  f3 = (double *) malloc (sizeof (double) * nm3);
  e3 = (double *) malloc (sizeof (double) * nm5);

  nf = np - nm;
  nf3 = nf * 3;
  nf5 = nf * 5;
  uf3 = (double *) malloc (sizeof (double) * nf3);
  of3 = (double *) malloc (sizeof (double) * nf3);
  ef3 = (double *) malloc (sizeof (double) * nf5);

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

  calc_mob_lub_fix_ewald_3fts (sys,
			       f3, t3, e3, uf3, of3, ef3,
			       u, o, s, ff, tf, sf);

  free (f3);
  free (e3);
  free (uf3);
  free (of3);
  free (ef3);
}
