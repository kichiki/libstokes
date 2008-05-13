/* excluded-volume interactions
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: excluded-volume.c,v 1.5 2008/05/13 01:10:18 kichiki Exp $
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

#include "excluded-volume.h"


/* initialize struct EV
 * INPUT
 *  bonds  : struct bonds (either fene=0 or fene=1 is fine).
 *  length : unit of length given by "length" in SCM (dimensional number)
 *  peclet : peclet number (with respect to "length")
 *  r2     : square of the max distance for F^{EV}
 *  v[n]   : EV parameters for each spring.
 *           the index should correspond to that in bonds.
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV, where l and A are defined by 
 *      ev->l[i] = sqrt((2/3) hat(ls)^2)
 *      ev->A[i] = (1/pi)^{3/2} hat(v) N_Ks^2 / (peclet ev->l[i]^5)
 *    where
 *      ls^2     = N_Ks * b_K^2 / 3 (dimensional, so as b_K)
 *      hat(ls)  = ls / length      (dimensionless)
 *      hat(v)   = v / length^3     (dimensionless)
 *    with which the force F^EV is given by 
 *      F^{EV}_{i} = A * r_{ij} * exp (- r_{ij}^2 / l^2)
 *    (note: r_{ij} is also dimensionless scaled by length)
 */
struct EV *
EV_init (const struct bonds *bonds,
	 double length, double peclet,
	 double r2, const double *v,
	 int np)
{
  struct EV *ev = (struct EV *)malloc (sizeof (struct EV));
  CHECK_MALLOC (ev, "EV_init");

  ev->r2 = r2;
  ev->n = bonds->n;
  ev->l = (double *)malloc (sizeof (double) * ev->n);
  ev->A = (double *)malloc (sizeof (double) * ev->n);
  ev->ch = (int *)malloc (sizeof (int) * np);
  CHECK_MALLOC (ev->l,  "EV_init");
  CHECK_MALLOC (ev->A,  "EV_init");
  CHECK_MALLOC (ev->ch, "EV_init");

  /* make the table for each chain */
  double length3 = length * length * length;
  double pi_m32 = pow (M_PI, -1.5); // (1 / pi)^{3/2}

  int i;
  for (i = 0; i < ev->n; i ++)
    {
      if (v[i] == 0.0)
	{
	  // no excluded-volume interaction for the spring "i"
	  ev->l[i] = 0.0;
	  ev->A[i] = 0.0;
	}
      else if (bonds->type[i] == 6)
	{
	  // EV cannot be defined for dWLC chain
	  fprintf (stderr, "EV_init:"
		   " EV for spring %d of dWLC chain couldn't define.\n",
		   i);
	  ev->l[i] = 0.0;
	  ev->A[i] = 0.0;
	}
      else
	{
	  double N_Ks, b_K;
	  if (bonds->fene[i] == 1)
	    {
	      // given parameters
	      N_Ks = bonds->p1[i];
	      b_K  = bonds->p2[i];

	      b_K /= length; // dimensionless
	    }
	  else
	    {
	      // given parameters
	      double Asp = bonds->p1[i];
	      double Ls  = bonds->p2[i];

	      b_K = 3.0 / (peclet * Asp); // dimensionless
	      N_Ks = Ls / b_K;
	    }

	  double ls = sqrt (N_Ks * b_K * b_K / 3.0); // dimensionless
	  ev->l[i] = ls / sqrt (1.5);
	  ev->A[i] = pi_m32 / peclet
	    * (v[i] / length3)
	    * N_Ks * N_Ks
	    * pow (ev->l[i], -5.0);
	  // with these, F^EV(r) = A * r * exp (- r^2 / l^2)
	}
    }
  
  /* make the particle table */
  for (i = 0; i < np; i ++)
    {
      ev->ch[i] = -1; // initialized by the unassigment
    }
  for (i = 0; i < bonds->n; i ++)
    {
      // now i is the chain (bond) type
      struct bond_pairs *pairs = bonds->pairs [i];

      int j;
      for (j = 0; j < pairs->n; j ++)
	{
	  int ia = pairs->ia [j];
	  int ib = pairs->ib [j];
	  ev->ch[ia] = i; // i is the chain (bond) type.
	  ev->ch[ib] = i; // i is the chain (bond) type.
	}
    }


  return (ev);
}

void
EV_free (struct EV *ev)
{
  if (ev == NULL) return;
  if (ev->l  != NULL) free (ev->l);
  if (ev->A  != NULL) free (ev->A);
  if (ev->ch != NULL) free (ev->ch);
  free (ev);
}

/* retrieve the EV parameters for particles i and j
 * with the combination rules 
 *   l(12) = (l1 + l2) / 2
 *   A(12) = (A1 * A2)^{1/2}
 * INPUT
 *  ev   : struct EV
 *  i, j : particle index (0 ~ np-1)
 * OUTPUT
 *  *l : 
 *  *A : 
 */
void
EV_get_coefficients (struct EV *ev,
		     int i, int j,
		     double *l, double *A)
{
  int ic = ev->ch[i];
  int jc = ev->ch[j];

  if (ic == jc)
    {
      *l = ev->l[ic];
      *A = ev->A[ic];
    }
  else
    {
      *l = 0.5 * (ev->l[ic] + ev->l[jc]);
      *A = sqrt (ev->A[ic] * ev->A[jc]);
    }
}


static void
EV_set_force_ij (struct stokes *sys,
		 struct EV *ev,
		 int i, int j, double r2, double x, double y, double z,
		 double *f)
{
  // define fr
  double l, A;
  EV_get_coefficients (ev, i, j, &l, &A);
  double l2 = l * l; // where 1 / l2 = (3 / 2 ls^2)
  double fr = A * exp (- r2 / l2);
  /* NOTE: fr has "r" and the unit vector e has "1/r", so that
   * here we implement fr' = A * exp (- rl2)
   * and the vector f' = fr' * (x, y, z).
   */

  /* F_a = fr * (R_a - R_b)/|R_a - R_b|
   * where fr > 0 corresponds to the repulsive
   *   and fr < 0 corresponds to the attractive
   */
  if (i < sys->nm)
    {
      int i3 = i * 3;
      f[i3+0] += fr * x;
      f[i3+1] += fr * y;
      f[i3+2] += fr * z;
    }
  if (j < sys->nm)
    {
      int j3 = j * 3;
      f[j3+0] += - fr * x;
      f[j3+1] += - fr * y;
      f[j3+2] += - fr * z;
    }
}

/*
 * INPUT
 *  ev         : struct EV
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_calc_force (struct EV *ev,
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

  double x, y, z;
  double r2;
  for (i = 0; i < sys->nm; i ++)
    {
      int ia3 = i * 3;
      int j;
      // loop for each pair (i <= j)
      for (j = i; j < sys->np; j ++)
	{
	  int ib3 = j * 3;
	  x = sys->pos [ia3+0] - sys->pos [ib3+0];
	  y = sys->pos [ia3+1] - sys->pos [ib3+1];
	  z = sys->pos [ia3+2] - sys->pos [ib3+2];

	  if (i != j)
	    {
	      r2 = x*x + y*y + z*z;
	      if (r2 < ev->r2)
		{
		  EV_set_force_ij (sys, ev, i, j, r2, x, y, z,
				   f);
		}
	    }
	  if (sys->periodic != 0)
	    {
	      int k;
	      for (k = 1; k < 27; k ++) // excluding the primary cell (k=0)
		{
		  double xx = x + (double)sys->ilx[k] * sys->lx;
		  double yy = y + (double)sys->ily[k] * sys->ly;
		  double zz = z + (double)sys->ilz[k] * sys->lz;

		  // shift for shear
		  if (sys->shear_mode == 1)
		    {
		      xx += (double)sys->ily[k] * sys->shear_shift;
		    }
		  else if (sys->shear_mode == 2)
		    {
		      xx += (double)sys->ilz[k] * sys->shear_shift;
		    }

		  r2 = xx*xx + yy*yy + zz*zz;
		  if (r2 < ev->r2)
		    {
		      EV_set_force_ij (sys, ev, i, j, r2, xx, yy, zz,
				       f);
		    }
		}
	    }
	}
    }
}
