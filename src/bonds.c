/* bond interaction between particles
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds.c,v 1.16 2008/07/27 00:50:02 kichiki Exp $
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

#include "bonds.h"


struct BONDS *
BONDS_init (void)
{
  struct BONDS *b
    = (struct BONDS *)malloc (sizeof (struct BONDS));
  CHECK_MALLOC (b, "BONDS_init");

  b->n = 0;

  b->type = NULL;
  b->fene = NULL;

  b->p1 = NULL;
  b->p2 = NULL;
  b->p3 = NULL;

  b->ia = NULL;
  b->ib = NULL;

  return (b);
}

void
BONDS_free (struct BONDS *b)
{
  if (b == NULL) return;

  if (b->type != NULL) free (b->type);
  if (b->fene != NULL) free (b->fene);
  if (b->p1 != NULL) free (b->p1);
  if (b->p2 != NULL) free (b->p2);
  if (b->p3 != NULL) free (b->p3);
  if (b->ia != NULL) free (b->ia);
  if (b->ib != NULL) free (b->ib);

  free (b);
}

/*
 * INPUT
 *  ia, ib : (global) particle index, that is, they are in the range [0, NP), 
 *           where NP is the total number of particles.
 */
void
BONDS_append (struct BONDS *b,
	      int type,
	      int fene,
	      double p1,
	      double p2,
	      double p3,
	      int ia,
	      int ib)
{
  if (ia == ib)
    {
      fprintf (stderr, "# BONDS_append() : invalid ia,ib = %d,%d\n",
	       ia, ib);
      exit (1);
    }
  // sort as ia < ib
  if (ia > ib)
    {
      int ic = ia;
      ia = ib;
      ib = ic;
    }

  b->n ++;

  b->type = (int *)realloc (b->type, sizeof (int) * b->n);
  b->fene = (int *)realloc (b->fene, sizeof (int) * b->n);
  b->p1 = (double *)realloc (b->p1, sizeof (double) * b->n);
  b->p2 = (double *)realloc (b->p2, sizeof (double) * b->n);
  b->p3 = (double *)realloc (b->p3, sizeof (double) * b->n);
  b->ia = (int *)realloc (b->ia, sizeof (int) * b->n);
  b->ib = (int *)realloc (b->ib, sizeof (int) * b->n);
  CHECK_MALLOC (b->type, "BONDS_append");
  CHECK_MALLOC (b->fene, "BONDS_append");
  CHECK_MALLOC (b->p1, "BONDS_append");
  CHECK_MALLOC (b->p2, "BONDS_append");
  CHECK_MALLOC (b->p3, "BONDS_append");
  CHECK_MALLOC (b->ia, "BONDS_append");
  CHECK_MALLOC (b->ib, "BONDS_append");

  int n = b->n - 1;
  b->type[n] = type;
  b->fene[n] = fene;
  b->p1[n] = p1;
  b->p2[n] = p2;
  b->p3[n] = p3;
  b->ia[n] = ia;
  b->ib[n] = ib;
}


static void
BONDS_exchange_bonds (struct BONDS *b,
		      int i, int j)
{
  int type = b->type[i];
  int fene = b->fene[i];
  double p1 = b->p1[i];
  double p2 = b->p2[i];
  double p3 = b->p3[i];
  int ia = b->ia[i];
  int ib = b->ib[i];

  b->type[i] = b->type[j];
  b->fene[i] = b->fene[j];
  b->p1[i] = b->p1[j];
  b->p2[i] = b->p2[j];
  b->p3[i] = b->p3[j];
  b->ia[i] = b->ia[j];
  b->ib[i] = b->ib[j];

  b->type[j] = type;
  b->fene[j] = fene;
  b->p1[j] = p1;
  b->p2[j] = p2;
  b->p3[j] = p3;
  b->ia[j] = ia;
  b->ib[j] = ib;
}

void
BONDS_sort_by_ia (struct BONDS *b)
{
  int i;
  for (i = 0; i < b->n - 1; i ++)
    {
      if (b->ia[i] > b->ia[i+1])
	{
	  BONDS_exchange_bonds (b, i, i+1);
	  int k;
	  for (k = i - 1;
	       k >= 0 && b->ia[k] > b->ia[k+1];
	       k --)
	    {
	      BONDS_exchange_bonds (b, k, k+1);
	    }
	}
    }
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
BONDS_set_FENE (struct BONDS *b,
		double length, double peclet)
{
  int i;
  for (i = 0; i < b->n; i ++)
    {
      if (b->fene[i] == 1)
	{
	  if (b->type[i] == 6)
	    {
	      // dWLC spring
	      double k  = b->p1[i];
	      double r0 = b->p2[i];

	      // first scale r0
	      r0 /= length;
	      b->p1[i] = k / (peclet * r0 * r0);
	      b->p2[i] = r0;
	    }
	  else
	    {
	      // bond i is FENE spring
	      double N_Ks = b->p1[i];
	      double b_K  = b->p2[i];
	      b->p1[i] = 3.0 * length / (peclet * b_K);
	      b->p2[i] = N_Ks * b_K / length;
	    }
	  // now, (p1, p2) = (A^{sp}, L_{s}).
	  b->fene[i] = 0;
	}
    }
}

static double 
BONDS_ILC (double x)
{
  // coth (x) = cosh (x) / sinh (x)
  return (cosh (x) / sinh (x) - 1.0 / x);
}


/* return force function (scalar part) of the spring
 */
double
BONDS_fr_i (struct BONDS *b,
	    int bond_index,
	    double Q)
{
  double fr;

  double A  = b->p1[bond_index];
  double Q0 = b->p2[bond_index];
  double s  = b->p3[bond_index]; // only for FENE-Fraenkel

  double hatQ = Q / Q0;
  double x2;

  int bond_type = b->type[bond_index];
  if ((bond_type != 0 &&
       bond_type != 5 &&
       bond_type != 6 &&
       bond_type != 7) // FENE chain
      && hatQ >= 1.0)
    {
      fprintf (stderr, "BONDS_fr_i:"
	       " extension %e exceeds the max %e for the bond %d\n",
	       Q, Q0, bond_index);
      fprintf (stderr, "bond_type = %d\n", bond_type);
      exit (1);
    }

  switch (bond_type)
    {
    case 1: // Wormlike chain (WLC)
      x2 = (1.0 - hatQ) * (1.0 - hatQ);
      fr = (2.0/3.0) * A * (0.25 / x2 - 0.25 + hatQ);
      break;

    case 2: // inverse Langevin chain (ILC)
      fr = A / 3.0 * BONDS_ILC (hatQ);
      break;

    case 3: // Cohen's Pade approximation for ILC
      x2 = hatQ * hatQ;
      fr = A * hatQ * (1.0 - x2 / 3.0) / (1.0 - x2);
      break;

    case 4: // Warner spring
      x2 = hatQ * hatQ;
      fr = A * hatQ / (1.0 - x2);
      break;

    case 5: // Hookean
      fr = A * hatQ;
      break;

    case 7: // FENE-Fraenkel
      // note that (p1,p2,p3) = (H,r0,s)
      x2 = (1.0 - hatQ) / s; // s is the tolerance
      x2 = x2 * x2;
      fr = A * (hatQ - 1.0) / (1.0 - x2);
      break;

    case 0: // Hookean with natural length == Fraenkel
    case 6: // dWLC
      // note that here (p1,p2) = (k, r0), not (N_{K,s}, b_{K})
      fr = A * (Q - Q0);
      break;

    default:
      fprintf (stderr, "bonds: invalid spring type %d\n", bond_type);
      exit (1);
      break;
    }

  return (fr);
}

/* calc spring force for the bond "i"
 * INPUT
 *  bonds : struct BONDS
 *  ib    : bond index
 *  q[3]  : connector vector for the bond "i"
 * OUTPUT
 *  f[3]  : force due to the bond "i"
 */
void
BONDS_calc_force_spring_i (struct BONDS *bonds,
			   int ib,
			   const double *q,
			   double *f)
{
  double Q = sqrt(q[0] * q[0]
		  + q[1] * q[1]
		  + q[2] * q[2]);
  double ex = q[0] / Q;
  double ey = q[1] / Q;
  double ez = q[2] / Q;

  double fr = BONDS_fr_i (bonds, ib, Q);

  f[0] = fr * ex;
  f[1] = fr * ey;
  f[2] = fr * ez;
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
 *  b          : struct BONDS
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
BONDS_calc_force (struct BONDS *b,
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

  for (i = 0; i < b->n; i ++)
    {
      int ia = b->ia [i];
      int ib = b->ib [i];
      // skip if both particles are fixed
      if (ia >= sys->nm && ib >= sys->nm) continue;

      int ia3 = ia * 3;
      int ib3 = ib * 3;
      double x = sys->pos [ia3  ] - sys->pos [ib3  ];
      double y = sys->pos [ia3+1] - sys->pos [ib3+1];
      double z = sys->pos [ia3+2] - sys->pos [ib3+2];

      if (sys->periodic != 0)
	{
	  search_close_image (sys, &x, &y, &z);
	}

      double r = sqrt (x*x + y*y + z*z);
      x /= r;
      y /= r;
      z /= r;

      double fr = BONDS_fr_i (b, i, r);

      /* F_a = - fr (R_a - R_b)/|R_a - R_b|
       * where fr > 0 corresponds to the attractive
       *   and fr < 0 corresponds to the repulsive
       */
      if (ia < sys->nm)
	{
	  int ia3 = ia * 3;
	  f[ia3+0] += - fr * x;
	  f[ia3+1] += - fr * y;
	  f[ia3+2] += - fr * z;
	}
      if (ib < sys->nm)
	{
	  int ib3 = ib * 3;
	  f[ib3+0] += fr * x;
	  f[ib3+1] += fr * y;
	  f[ib3+2] += fr * z;
	}
    }
}
