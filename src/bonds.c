/* bond interaction between particles
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds.c,v 1.7 2007/11/30 06:39:06 kichiki Exp $
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


/* initialize struct bonds
 * INPUT
 * OUTPUT
 *  returned value : struct bonds
 */
struct bonds *
bonds_init (void)
{
  struct bonds *bonds = NULL;
  bonds = (struct bonds *)malloc (sizeof (struct bonds));

  bonds->n = 0;
  bonds->type  = NULL;
  bonds->fene  = NULL;
  bonds->p1    = NULL;
  bonds->p2    = NULL;
  bonds->pairs = NULL;

  return (bonds);
}

void
bonds_free (struct bonds *bonds)
{
  if (bonds == NULL) return;
  if (bonds->type  != NULL) free (bonds->type);
  if (bonds->fene  != NULL) free (bonds->fene);
  if (bonds->p1    != NULL) free (bonds->p1);
  if (bonds->p2    != NULL) free (bonds->p2);
  if (bonds->pairs != NULL)
    {
      int i;
      for (i = 0; i < bonds->n; i ++)
	{
	  bond_pairs_free (bonds->pairs [i]);
	}
      free (bonds->pairs);
    }
  free (bonds);
}

/* add a spring into bonds
 * INPUT
 *  bonds  : struct bonds
 *  type   : type of the spring
 *  fene   : 0 == (p1,p2) are (A^{sp}, L_{s})
 *           1 == (p1, p2) = (N_{K,s}, b_{K})
 *  p1, p2 : spring parameters
 * OUTPUT
 *  bonds  :
 */
void
bonds_add_type (struct bonds *bonds,
		int type, int fene, double p1, double p2)
{
  bonds->n ++;

  int n;
  n = bonds->n;
  bonds->type = (int *)realloc (bonds->type,  sizeof (int) * n);
  bonds->fene = (int *)realloc (bonds->fene,  sizeof (int) * n);
  bonds->p1   = (double *)realloc (bonds->p1, sizeof (double) * n);
  bonds->p2   = (double *)realloc (bonds->p2, sizeof (double) * n);
  bonds->pairs
    = (struct bond_pairs **)realloc (bonds->pairs,
				     sizeof (struct bond_pairs *) * n);

  // set n as the newly added element
  n--;
  bonds->type [n] = type;
  bonds->fene [n] = fene;
  bonds->p1   [n] = p1;
  bonds->p2   [n] = p2;

  bonds->pairs [n] = bond_pairs_init ();
}

/* set FENE spring parameters for run
 * INPUT
 *  bonds : p1 (N_{K,s}) and p2 (b_{K}) are used.
 *  a     : length scale in the simulation
 *  pe    : peclet number
 * OUTPUT
 *  bonds->p1[] := A^{sp} = 3a / pe b_{K}
 *  bonds->p2[] := Ls / a = N_{K,s} b_{K} / a 
 */
void
bonds_set_FENE (struct bonds *bonds,
		double a, double pe)
{
  int i;
  for (i = 0; i < bonds->n; i ++)
    {
      if (bonds->fene[i] == 1)
	{
	  // bond i is FENE spring
	  double N_Ks = bonds->p1[i];
	  double b_K  = bonds->p2[i];
	  bonds->p1[i] = 3.0 * a / (pe * b_K);
	  bonds->p2[i] = N_Ks * b_K / a;

	  // now, (p1, p2) = (A^{sp}, L_{s}).
	  bonds->fene[i] = 0;
	}
    }
}


static double 
bonds_ILC (double x)
{
  // coth (x) = cosh (x) / sinh (x)
  return (cosh (x) / sinh (x) - 1.0 / x);
}

static void
bonds_set_force_ij (struct bonds *bonds,
		    struct stokes *sys,
		    int bond_index,
		    int ia, int ib, 
		    double x, double y, double z,
		    double *f)
{
  double r2 = x*x + y*y + z*z;
  double r = sqrt (r2);
  double ex = x / r;
  double ey = y / r;
  double ez = z / r;

  double Asp = bonds->p1[bond_index];
  double Ls  = bonds->p2[bond_index];

  double fr;
  double rLs = r / Ls;
  double rLs2;

  int bond_type = bonds->type[bond_index];
  if ((bond_type != 0 && bond_type != 5) // FENE chain
      && rLs >= 1.0)
    {
      fprintf (stderr, "bonds: extension %e exceeds the max %e for (%d,%d)\n",
	       r, Ls, ia, ib);
      fprintf (stderr, "bond_type = %d\n", bond_type);
      exit (1);
    }

  switch (bond_type)
    {
    case 1: // Wormlike chain (WLC)
      rLs2 = (1.0 - rLs) * (1.0 - rLs);
      fr = (2.0/3.0) * Asp * (0.25 / rLs2 - 0.25 + rLs);
      break;

    case 2: // inverse Langevin chain (ILC)
      fr = Asp / 3.0 * bonds_ILC (rLs);
      break;

    case 3: // Cohen's Pade approximation for ILC
      rLs2 = rLs * rLs;
      fr = Asp * rLs * (1.0 - rLs2 / 3.0) / (1.0 - rLs2);
      break;

    case 4: // Warner spring
      rLs2 = rLs * rLs;
      fr = Asp * rLs / (1.0 - rLs2);
      break;

    case 5: // Hookean
      fr = Asp * rLs;
      break;

    case 0: // Hookean
    default:
      fr = Asp * (r - Ls);
      break;
    }

  /* F_a = - fr (R_a - R_b)/|R_a - R_b|
   * where fr > 0 corresponds to the attractive
   *   and fr < 0 corresponds to the repulsive
   */
  if (ia < sys->nm)
    {
      int ia3 = ia * 3;
      f[ia3+0] += - fr * ex;
      f[ia3+1] += - fr * ey;
      f[ia3+2] += - fr * ez;
    }
  if (ib < sys->nm)
    {
      int ib3 = ib * 3;
      f[ib3+0] += fr * ex;
      f[ib3+1] += fr * ey;
      f[ib3+2] += fr * ez;
    }
}
 
/* search the closest image in 27 periodic images
 * INPUT
 *  sys : struct stokes
 *  x, y, z : relative distance for the pair
 * OUTPUT
 *  k : image index in [0, 27)
 */
static int
search_close_image (struct stokes *sys,
		    double x, double y, double z)
{
  int k0 = 0;
  double r2 = x*x + y*y + z*z;

  int k;
  for (k = 1; k < 27; k ++)
    {
      double xx = x - sys->llx[k];
      double yy = y - sys->lly[k];
      double zz = z - sys->llz[k];
      double rr2 = xx*xx + yy*yy + zz*zz;
      if (rr2 < r2) k0 = k;
    }

  return (k0);
}


/*
 * INPUT
 *  bonds      : struct bond
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

  for (i = 0; i < bonds->n; i ++)
    {
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

	  if (sys->periodic == 0)
	    {
	      bonds_set_force_ij (bonds, sys,
				  i,
				  ia, ib, x, y, z,
				  f);
	    }
	  else
	    {
	      int k = search_close_image (sys, x, y, z);
	      double xx = x - sys->llx[k];
	      double yy = y - sys->lly[k];
	      double zz = z - sys->llz[k];

	      bonds_set_force_ij (bonds, sys,
				  i,
				  ia, ib, xx, yy, zz,
				  f);
	    }
	}
    }
}


void
fprint_bonds (FILE *out, struct bonds *bonds)
{
  int i;
  for (i = 0; i < bonds->n; i ++)
    {
      fprintf (out, "bond-index %d: type = %d, k = %f, r0 = %f,"
	       " number of pairs = %d\n",
	       i, bonds->type[i], bonds->p1[i], bonds->p2[i], 
	       bonds->pairs[i]->n);
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


/**
 * SWIG utility routine
 * For examplean, expected usage in python by SWIG:
 *   n = stokes.bonds_get_pairs_n(bonds, i)
 *   ia = stokes.iarray(n)
 *   ib = stokes.iarray(n)
 *   stokes.bonds_get_pairs(bonds, i, ia, ib)
 * then, you have arrays ia[n] and ib[n].
 */

/* to get the number of pairs for the bond "i"
 */
int
bonds_get_pairs_n (struct bonds *b, int i)
{
  return (b->pairs[i]->n);
}
/* 
 * INPUT
 *  b : struct bonds
 *  i : index of the bond
 *  ia[n] : "a" particle index of j-th pair for the bond "i",
 *  ib[n] : "b" particle index of j-th pair for the bond "i",
 *          where j runs from 0 to (n-1) and
 *          n = b->pairs[i]->n is the number of pairs for the bond "i".
 *          before calling, allocate ia and ib with (sizeof(int) * n).
 * OUTPUT
 *  ia[n], ib[n] : 
 */
void
bonds_get_pairs (struct bonds *b, int i,
		 int *ia, int *ib)
{
  int j;
  for (j = 0; j < b->pairs[i]->n; j ++)
    {
      ia[j] = b->pairs[i]->ia[j];
      ib[j] = b->pairs[i]->ib[j];
    }
}
