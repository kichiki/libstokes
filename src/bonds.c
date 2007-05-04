/* bond interaction between particles
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds.c,v 1.3 2007/05/04 01:20:47 kichiki Exp $
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
#include <stdlib.h>
#include <math.h> // sqrt

#include "stokes.h" // struct stokes
#include "memory-check.h" // macro CHECK_MALLOC

#include "bonds.h" // struct bonds


struct bond_pairs *
bond_pairs_init (void)
{
  struct bond_pairs *pairs = NULL;
  pairs = (struct bond_pairs *)malloc (sizeof (struct bond_pairs));

  pairs->n = 0;
  pairs->ia = NULL;
  pairs->ib = NULL;

  return (pairs);
}

void
bond_pairs_free (struct bond_pairs *pairs)
{
  if (pairs == NULL) return;
  if (pairs->ia != NULL) free (pairs->ia);
  if (pairs->ib != NULL) free (pairs->ib);
  free (pairs);
}

void
bond_pairs_add (struct bond_pairs *pairs,
		int ia, int ib)
{
  pairs->n++;

  int n;
  n = pairs->n;
  pairs->ia = (int *)realloc (pairs->ia, sizeof (int) * n);
  pairs->ib = (int *)realloc (pairs->ib, sizeof (int) * n);

  // set n as the newly added element
  n--;
  pairs->ia [n] = ia;
  pairs->ib [n] = ib;
}


struct bonds *
bonds_init (void)
{
  struct bonds *bonds = NULL;
  bonds = (struct bonds *)malloc (sizeof (struct bonds));

  bonds->ntypes = 0;
  bonds->k  = NULL;
  bonds->r0 = NULL;
  bonds->pairs = NULL;

  return (bonds);
}

void
bonds_free (struct bonds *bonds)
{
  if (bonds == NULL) return;
  if (bonds->k  != NULL) free (bonds->k);
  if (bonds->r0 != NULL) free (bonds->r0);
  if (bonds->pairs != NULL)
    {
      int i;
      for (i = 0; i < bonds->ntypes; i ++)
	{
	  bond_pairs_free (bonds->pairs [i]);
	}
      free (bonds->pairs);
    }
  free (bonds);
}

void
bonds_add_type (struct bonds *bonds,
		double k, double r0)
{
  bonds->ntypes ++;

  int n;
  n = bonds->ntypes;
  bonds->k  = (double *)realloc (bonds->k,  sizeof (double) * n);
  bonds->r0 = (double *)realloc (bonds->r0, sizeof (double) * n);
  bonds->pairs
    = (struct bond_pairs **)realloc (bonds->pairs,
				     sizeof (struct bond_pairs *) * n);

  // set n as the newly added element
  n--;
  bonds->k  [n] = k;
  bonds->r0 [n] = r0;

  bonds->pairs [n] = bond_pairs_init ();
}

/*
 * INPUT
 *  b          : struct bond
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
bonds_calc_force (struct bonds *bonds,
		  struct stokes *sys,
		  double *f,
		  int flag_add)
{
  int i;

  if (flag_add == 0)
    {
      // clear the force
      for (i = 0; i < sys->nm * 3; i ++)
	{
	  f [i] = 0.0;
	}
    }

  for (i = 0; i < bonds->ntypes; i ++)
    {
      double k  = bonds->k  [i];
      double r0 = bonds->r0 [i];
      struct bond_pairs *pairs = bonds->pairs [i];

      int j;
      for (j = 0; j < pairs->n; j ++)
	{
	  int ia = pairs->ia [j];
	  int ib = pairs->ib [j];
	  // skip if both particles are fixed
	  if (ia >= sys->nm && ib >= sys->nm) continue;

	  int ia3 = ia * 3;
	  int ib3 = ib * 3;
	  double x = sys->pos [ia3+0] - sys->pos [ib3+0];
	  double y = sys->pos [ia3+1] - sys->pos [ib3+1];
	  double z = sys->pos [ia3+2] - sys->pos [ib3+2];
	  double r = sqrt (x*x + y*y + z*z);
	  double ex = x / r;
	  double ey = y / r;
	  double ez = z / r;

	  double fac = k * (r - r0);

	  // F_a = - k (r-r0) (R_a - R_b)/|R_a - R_b|
	  if (ia < sys->nm)
	    {
	      f[ia3+0] += - fac * ex;
	      f[ia3+1] += - fac * ey;
	      f[ia3+2] += - fac * ez;
	    }

	  if (ib < sys->nm)
	    {
	      f[ib3+0] += fac * ex;
	      f[ib3+1] += fac * ey;
	      f[ib3+2] += fac * ez;
	    }
	}
    }
}


void
fprint_bonds (FILE *out, struct bonds *bonds)
{
  int i;
  for (i = 0; i < bonds->ntypes; i ++)
    {
      fprintf (out, "bond-type %d: k = %f, r0 = %f, number of pairs = %d\n",
	       i, bonds->k[i], bonds->r0[i], bonds->pairs[i]->n);
      int j;
      for (j = 0; j < bonds->pairs[i]->n; j ++)
	{
	  fprintf (out, "\t pair %d : %d, %d\n",
		   j,
		   bonds->pairs[i]->ia[j],
		   bonds->pairs[i]->ib[j]);
	}
    }
}
