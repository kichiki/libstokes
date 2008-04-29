/* Solvers for no-hydrodynamic interaction problems
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: noHI.c,v 1.1 2008/04/29 03:21:03 kichiki Exp $
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
#include "memory-check.h" // macro CHECK_MALLOC

#include "stokes.h" /* struct stokeks */
#include "f.h"   // shift_labo_to_rest_U()
#include "ft.h"  // shift_labo_to_rest_O()
#include "fts.h" // shift_labo_to_rest_E()

#include "noHI.h"


/* solve natural mobility problem with fixed particles in F version
 * without HI
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  uf [nf * 3] :
 * OUTPUT
 *  u [nm * 3] :
 *  ff [nf * 3] :
 */
void
solve_mix_3f_noHI (struct stokes * sys,
		   const double *f,
		   const double *uf,
		   double *u,
		   double *ff)
{
  if (sys->version != 0)
    {
      fprintf (stderr, "libstokes solve_mix_3f_noHI :"
	       " the version is wrong. reset to F\n");
      sys->version = 0;
    }

  int np = sys->np;
  int nm = sys->nm;
  int i;

  // for fixed particles
  int nf = np - nm;
  if (nf > 0)
    {
      double *uf0 = (double *)malloc (sizeof (double) * nf * 3);
      CHECK_MALLOC (uf0, "solve_mix_3f_noHI");

      shift_labo_to_rest_U (sys, nf, uf, uf0);
      /* the main calculation is done in the the fluid-rest frame;
       * u(x)=0 as |x|-> infty */
      for (i = 0; i < nf * 3; i ++)
	{
	  ff[i] = uf0[i];
	}
      free (uf0);
    }

  // for mobile particles
  for (i = 0; i < nm * 3; i ++)
    {
      u[i] = f[i];
    }
  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
}
/* solve natural mobility problem with fixed particles in FT version
 * without HI
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  uf [nf * 3] :
 *  of [nf * 3] :
 * OUTPUT
 *  u [nm * 3] :
 *  o [nm * 3] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 */
void
solve_mix_3ft_noHI (struct stokes * sys,
		    const double *f, const double *t,
		    const double *uf, const double *of,
		    double *u, double *o,
		    double *ff, double *tf)
{
  if (sys->version != 1)
    {
      fprintf (stderr, "libstokes solve_mix_3ft_noHI :"
	       " the version is wrong. reset to FT\n");
      sys->version = 1;
    }

  int np = sys->np;
  int nm = sys->nm;
  int i;

  // for fixed particles
  int nf = np - nm;
  if (nf > 0)
    {
      double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
      double *of0 = (double *) malloc (sizeof (double) * nf * 3);
      CHECK_MALLOC (uf0, "solve_mix_3ft_noHI");
      CHECK_MALLOC (of0, "solve_mix_3ft_noHI");

      shift_labo_to_rest_U (sys, nf, uf, uf0);
      shift_labo_to_rest_O (sys, nf, of, of0);
      /* the main calculation is done in the the fluid-rest frame;
       * u(x)=0 as |x|-> infty */
      for (i = 0; i < nf * 3; i ++)
	{
	  ff[i] = uf0[i];
	  tf[i] = of0[i] / 0.75;
	}
      free (uf0);
      free (of0);
    }

  // for mobile particles
  for (i = 0; i < nm * 3; i ++)
    {
      u[i] = f[i];
      o[i] = t[i] * 0.75;
    }
  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}
/* solve natural mobility problem with fixed particles in FTS version
 * without HI
 * for both periodic and non-periodic boundary conditions
 * INPUT
 *  sys : system parameters
 *  f [nm * 3] :
 *  t [nm * 3] :
 *  e [nm * 5] : in the labo frame.
 *  uf [nf * 3] : in the labo frame.
 *  of [nf * 3] : in the labo frame.
 *  ef [nf * 5] : in the labo frame.
 * OUTPUT
 *  u [nm * 3] : in the labo frame.
 *  o [nm * 3] : in the labo frame.
 *  s [nm * 5] :
 *  ff [nf * 3] :
 *  tf [nf * 3] :
 *  sf [nf * 5] :
 */
void
solve_mix_3fts_noHI (struct stokes * sys,
		     const double *f, const double *t, const double *e,
		     const double *uf, const double *of, const double *ef,
		     double *u, double *o, double *s,
		     double *ff, double *tf, double *sf)
{
  if (sys->version != 2)
    {
      fprintf (stderr, "libstokes solve_mix_3fts_noHI :"
	       " the version is wrong. reset to FTS\n");
      sys->version = 2;
    }

  int np = sys->np;
  int nm = sys->nm;
  int i;

  // for fixed particles
  int nf = np - nm;
  if (nf > 0)
    {
      double *uf0 = (double *) malloc (sizeof (double) * nf * 3);
      double *of0 = (double *) malloc (sizeof (double) * nf * 3);
      double *ef0 = (double *) malloc (sizeof (double) * nf * 5);
      CHECK_MALLOC (uf0, "solve_mix_3fts");
      CHECK_MALLOC (of0, "solve_mix_3fts");
      CHECK_MALLOC (ef0, "solve_mix_3fts");

      shift_labo_to_rest_U (sys, nf, uf, uf0);
      shift_labo_to_rest_O (sys, nf, of, of0);
      shift_labo_to_rest_E (sys, nf, ef, ef0);
      /* the main calculation is done in the the fluid-rest frame;
       * u(x)=0 as |x|-> infty */
      for (i = 0; i < nf * 3; i ++)
	{
	  ff[i] = uf0[i];
	  tf[i] = of0[i] / 0.75;
	}
      for (i = 0; i < nf * 5; i ++)
	{
	  sf[i] = ef0[i] / 0.9;
	}
      free (uf0);
      free (of0);
      free (ef0);
    }

  // for mobile particles
  double *e0 = (double *) malloc (sizeof (double) * nm * 5);
  CHECK_MALLOC (e0, "solve_mix_3fts");

  shift_labo_to_rest_E (sys, nm, e, e0);
  /* the main calculation is done in the the fluid-rest frame;
   * u(x)=0 as |x|-> infty */
  for (i = 0; i < nm * 3; i ++)
    {
      u[i] = f[i];
      o[i] = t[i] * 0.75;
    }
  for (i = 0; i < nm * 5; i ++)
    {
      s[i] = e0[i] / 0.9;
    }
  free (e0);

  /* for the interface, we are in the labo frame, that is
   * u(x) is given by the imposed flow field as |x|-> infty */
  shift_rest_to_labo_U (sys, nm, u);
  shift_rest_to_labo_O (sys, nm, o);
}
