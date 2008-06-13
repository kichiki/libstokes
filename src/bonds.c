/* bond interaction between particles
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds.c,v 1.13 2008/06/13 03:04:52 kichiki Exp $
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
#include "memory-check.h" // macro CHECK_MALLOC

#include "stokes.h" // struct stokes

#include "bonds.h" // struct bonds


struct bond_pairs *
bond_pairs_init (void)
{
  struct bond_pairs *pairs
    = (struct bond_pairs *)malloc (sizeof (struct bond_pairs));
  CHECK_MALLOC (pairs, "bond_pairs_init");

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
  bonds->nex   = NULL;

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
  if (bonds->nex   != NULL) free (bonds->nex);
  free (bonds);
}

/* add a spring into bonds
 * INPUT
 *  bonds  : struct bonds
 *  type   : type of the spring
 *  fene   : 0 == (p1,p2) are (A^{sp}, L_{s})
 *           1 == (p1, p2) = (N_{K,s}, b_{K}) or
 *                (p1, p2) = (k, r0) for dWLC (type == 6).
 *                in the latter case, potential is given by
 *                (k/2) * (kT / r0^2) * (r-r0)^2
 *  p1, p2 : spring parameters
 *  nex    : number of excluded particles in the chain
 * OUTPUT
 *  bonds  :
 */
void
bonds_add_type (struct bonds *bonds,
		int type, int fene, double p1, double p2,
		int nex)
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
  bonds->nex  = (int *)realloc (bonds->nex, sizeof (int) * n);


  // set n as the newly added element
  n--;
  bonds->type [n] = type;
  bonds->fene [n] = fene;
  bonds->p1   [n] = p1;
  bonds->p2   [n] = p2;

  bonds->pairs [n] = bond_pairs_init ();

  bonds->nex  [n] = nex;
}

/* set FENE spring parameters for run
 * INPUT
 *  bonds : p1 (N_{K,s}) and p2 (b_{K}) are used
 *          (bonds->fene[i] == 1 is expected), or 
 *          p1 (k) and p2 (r0) are used for dWLC spring (type == 6).
 *            in this case, potential is given by
 *            (k/2) * (kT / r0^2) * (r-r0)^2
 *  length : length scale in the simulation
 *  peclet : peclet number
 * OUTPUT
 *  bonds->p1[] := A^{sp} = 3 length / (peclet b_{K})
 *  bonds->p2[] := Ls / length = N_{K,s} b_{K} / length
 *    or, for dWLC spring, the conversions are given by 
 *  bonds->p1[] := A^{sp} = k / (pe * (r0/length)^2)
 *  bonds->p2[] := Ls / length = r0 / length
 */
void
bonds_set_FENE (struct bonds *bonds,
		double length, double peclet)
{
  int i;
  for (i = 0; i < bonds->n; i ++)
    {
      if (bonds->fene[i] == 1)
	{
	  if (bonds->type[i] == 6)
	    {
	      // dWLC spring
	      double k  = bonds->p1[i];
	      double r0 = bonds->p2[i];

	      // first scale r0
	      r0 /= length;
	      bonds->p1[i] = k / (peclet * r0 * r0);
	      bonds->p2[i] = r0;
	    }
	  else
	    {
	      // bond i is FENE spring
	      double N_Ks = bonds->p1[i];
	      double b_K  = bonds->p2[i];
	      bonds->p1[i] = 3.0 * length / (peclet * b_K);
	      bonds->p2[i] = N_Ks * b_K / length;
	    }
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
  if ((bond_type != 0 && bond_type != 5 && bond_type != 6) // FENE chain
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
    case 6: // dWLC
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
 *  x, y, z : relative distance for the pair in the primary cell
 * OUTPUT
 *  k : image index in [0, 27)
 *  x, y, z : relative distance for the closest pair
 */
static int
search_close_image (struct stokes *sys,
		    double *x, double *y, double *z)
{
  int k0 = 0;
  double x0 = (*x);
  double y0 = (*y);
  double z0 = (*z);
  double r2 = x0*x0 + y0*y0 + z0*z0;

  int k;
  for (k = 1; k < 27; k ++)
    {
      double xx = x0 + (double)sys->ilx[k] * sys->lx;
      double yy = y0 + (double)sys->ily[k] * sys->ly;
      double zz = z0 + (double)sys->ilz[k] * sys->lz;

      // shift for shear
      if (sys->shear_mode == 1)
	{
	  xx += (double)sys->ily[k] * sys->shear_shift;
	}
      else if (sys->shear_mode == 2)
	{
	  xx += (double)sys->ilz[k] * sys->shear_shift;
	}

      double rr2 = xx*xx + yy*yy + zz*zz;
      if (rr2 < r2)
	{
	  k0 = k;
	  (*x) = xx;
	  (*y) = yy;
	  (*z) = zz;
	  r2 = rr2;
	}
    }

  return (k0);
}


/*
 * INPUT
 *  bonds      : struct bonds
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

	  if (sys->periodic != 0)
	    {
	      search_close_image (sys, &x, &y, &z);
	    }
	  bonds_set_force_ij (bonds, sys,
			      i,
			      ia, ib, x, y, z,
			      f);
	}
    }
}


void
bonds_print (FILE *out, struct bonds *bonds)
{
  int i;
  for (i = 0; i < bonds->n; i ++)
    {
      fprintf (out, "bond-index %d:\n", i);
      fprintf (out, " type = %d\n", bonds->type[i]);
      fprintf (out, " k = %f\n", bonds->p1[i]);
      fprintf (out, " r0 = %f\n", bonds->p2[i]);
      fprintf (out, " nex = %d\n", bonds->nex[i]);
      fprintf (out, " number of pairs = %d\n", bonds->pairs[i]->n);
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


/**
 * exclusion list for lubrication due to the bonding
 */
struct list_ex *
list_ex_init (int np)
{
  struct list_ex *ex = (struct list_ex *)malloc (sizeof (struct list_ex));
  CHECK_MALLOC (ex, "list_ex_init");

  ex->np = np;
  ex->n = (int *)malloc (sizeof (int) * np);
  CHECK_MALLOC (ex->n, "list_ex_init");
  ex->i = (int **)malloc (sizeof (int *) * np);
  CHECK_MALLOC (ex->i, "list_ex_init");
  
  int i;
  for (i = 0; i < np; i ++)
    {
      ex->n[i] = 0;
      ex->i[i] = NULL;
    }

  return (ex);
}

void
list_ex_add (struct list_ex *ex, int j, int k)
{
  // self is not in the excluded list
  if (j == k) return;

  // check the duplication
  int i;
  for (i = 0; i < ex->n[j]; i ++)
    {
      if (ex->i[j][i] == k)
	{
	  // k is already in the list
	  return;
	}
    }

  // add k in the excluded list for particle j.
  ex->n[j] ++;
  ex->i[j] = (int *)realloc (ex->i[j], sizeof (int) * ex->n[j]);
  ex->i[j][ex->n[j]-1] = k;
}

void
list_ex_free (struct list_ex *ex)
{
  if (ex == NULL) return;
  int i;
  for (i = 0; i < ex->np; i ++)
    {
      if (ex->i[i] != NULL) free (ex->i[i]);
    }
  if (ex->i != NULL) free (ex->i);
  if (ex->n != NULL) free (ex->n);
}

struct list_ex *
list_ex_copy (struct list_ex *ex0)
{
  if (ex0 == NULL) return (NULL);

  struct list_ex *ex = list_ex_init (ex0->np);
  CHECK_MALLOC (ex, "list_ex_copy");

  int i;
  for (i = 0; i < ex0->np; i ++)
    {
      int j;
      for (j = 0; j < ex0->n[i]; j ++)
	{
	  list_ex_add (ex, i, ex0->i[i][j]);
	}
    }

  return (ex);
}

/* 
 * OUTPUT
 *  returned value : 0 == ia is NOT found in the list nn[n]
 *                   1 == ia IS found in the list nn[n]
 */
static int
check_nn_list (int n, const int *nn, int ia)
{
  int flag = 0;
  int k;
  for (k = 0; k < n; k ++)
    {
      if (nn[k] == ia)
	{
	  flag = 1;
	  break;
	}
    }
  return (flag);
}

/* construct the excluded list by struct bonds
 */
void
list_ex_set_by_bonds (struct list_ex *ex, const struct bonds *b)
{
  if (b == NULL) return;

  // i is the bond type
  int i;
  for (i = 0; i < b->n; i ++)
    {
      if (b->nex[i] == 0)
	{
	  // no exclusion
	  continue;
	}

      // first, count the number of particles in the chain
      struct bond_pairs *pairs = b->pairs [i];
      int n = 0;
      int *nn = NULL;
      int j;
      for (j = 0; j < pairs->n; j ++)
	{
	  int ia = pairs->ia [j];
	  int ib = pairs->ib [j];
	  // check the particle is in the list
	  if (check_nn_list (n, nn, ia) == 0)
	    {
	      // ia is new
	      n ++;
	      nn = (int *)realloc (nn, sizeof (int) * n);
	      CHECK_MALLOC (nn, "list_ex_set_by_bonds");
	      nn[n-1] = ia;
	    }
	  if (check_nn_list (n, nn, ib) == 0)
	    {
	      // ib is new
	      n ++;
	      nn = (int *)realloc (nn, sizeof (int) * n);
	      CHECK_MALLOC (nn, "list_ex_set_by_bonds");
      	      nn[n-1] = ib;
	    }
	}
      /* now n is the number of particles
       * and nn[n] is the particle indices.
       */

      int nex = b->nex[i];
      if (nex < 0) nex = n-1;
      // this is equivalent to exclude all particles in the chain


      // make the nearest-neighbor list in "ex"
      for (j = 0; j < pairs->n; j ++)
	{
	  int ia = pairs->ia [j];
	  int ib = pairs->ib [j];

	  list_ex_add (ex, ia, ib);
	  list_ex_add (ex, ib, ia);
	}
      /* now ex is the list of the nearest neighbors,
       * which is equivalent to the excluded list with nex = 1.
       */

      // extend the list to the level of nex.
      int m = 1;
      while (m < nex)
	{
	  // reference should be the old one
	  struct list_ex *ex0 = list_ex_copy (ex);
	  CHECK_MALLOC (ex0, "list_ex_set_by_bonds");

	  // loop for the particles in the chain
	  int k;
	  for (k = 0; k < n; k ++)
	    {
	      int ia = nn[k]; // the particle now considering
	      for (j = 0; j < ex0->n[ia]; j ++)
		{
		  int ib = ex0->i[ia][j]; // ia's j-th excluded particle
		  int l;
		  for (l = 0; l < ex0->n[ib]; l ++)
		    {
		      int ic = ex0->i[ib][l]; // l-th excluded particle for ib
		      if (ic != ia)
			{
			  // add ic in the excluded list for ia.
			  list_ex_add (ex, ia, ic);
			  /* note that the duplication is checked
			   * in list_ex_add()
			   */
			}
		    }
		}
	    }
	  list_ex_free (ex0);
	  m ++;
	}

      // done for the i-th chain.
      free (nn);
    }
}

/* check whether j is excluded for i
 * INPUT
 *  ex : struct list_ex
 *  i  : particle now we are considering
 *  j  : particle whether it is in the list or not.
 * OUTPUT
 *  returned value : 0 (false); j is NOT in the excluded list for i.
 *                   1 (true);  j IS in the excluded list for i.
 *                   NOTE, the self for i in some chain is EXCLUDED,
 *                   while the self for particles is NOT excluded.
 */
int
list_ex_check (struct list_ex *ex, int i, int j)
{
  if (ex->n[i] == 0) return 0; /* i is not in some chain, so 
				* false; j is NOT in the excluded list.
				*/
  if (i == j) return 1; // if i is in some chain, the self is excluded.

  int k;
  for (k = 0; k < ex->n[i]; k ++)
    {
      if (ex->i[i][k] == j) return 1; // true; j IS in the excluded list.
    }

  return 0;
}
