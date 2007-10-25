/* excluded-volume interactions
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: excluded-volume.c,v 1.1 2007/10/25 05:57:05 kichiki Exp $
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
 *  a, pe  : parameters for bonds parameters
 *  r2     : square of the max distance for F^{EV}
 *  v[n]   : EV parameters for each spring.
 *           the index should correspond to that in bonds.
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV, where l and A are defined by 
 *      ev->l[i] characteristic distance = sqrt{(1/3) N_{K,s} b_{K}^2},
 *      ev->A[i] prefactor = (9/2) A^{sp} z',
 *    and
 *      A^{sp} = 3 a / Pe b_{K},
 *      z' = (N_{K,s}/2 pi)^{3/2} (v/l_s^3)
 *         = (3 / 2 pi b_{K}^2)^{3/2} v.
 */
struct EV *
EV_init (const struct bonds *bonds, double a, double pe,
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
  CHECK_MALLOC (ev->l, "EV_init");
  CHECK_MALLOC (ev->A, "EV_init");
  CHECK_MALLOC (ev->ch, "EV_init");

  /* make the table for each chain */
  double N_Ks, b_K;
  double Asp, Ls;
  double zz;
  double l, A;
  double c = pow (1.5 / M_PI, 1.5); // (3 / 2 pi)^{3/2}
  int i;
  for (i = 0; i < ev->n; i ++)
    {
      if (bonds->fene[i] == 0)
	{
	  // given parameters
	  N_Ks = bonds->p1[i];
	  b_K  = bonds->p2[i];

	  l = sqrt (N_Ks * b_K * b_K / 3.0);
	  Asp = 3.0 * a / (pe * b_K);
	  zz = c * v[i] / (b_K * b_K * b_K);
	  A = 4.5 * Asp * zz;
	}
      else
	{
	  // given parameters
	  Asp = bonds->p1[i];
	  Ls  = bonds->p2[i];

	  b_K = Asp * pe / (3.0 * a);
	  N_Ks = Ls * a / b_K;
	  l = sqrt (N_Ks * b_K * b_K / 3.0);
	  zz = c * v[i] / (b_K * b_K * b_K);
	  A = 4.5 * Asp * zz;
	}
      ev->l[i] = l;
      ev->A[i] = A;
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
  if (ev->l != NULL) free (ev->l);
  if (ev->A != NULL) free (ev->A);
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
  double r = sqrt (r2);
  double ex = x / r;
  double ey = y / r;
  double ez = z / r;

  // define fr
  double fr, l, A;
  EV_get_coefficients (ev, i, j, &l, &A);
  r /= l;
  fr = A * r * exp (-1.5 * r * r);

  /* F_a = fr * (R_a - R_b)/|R_a - R_b|
   * where fr > 0 corresponds to the repulsive
   *   and fr < 0 corresponds to the attractive
   */
  if (i < sys->nm)
    {
      int i3 = i * 3;
      f[i3+0] += fr * ex;
      f[i3+1] += fr * ey;
      f[i3+2] += fr * ez;
    }
  if (j < sys->nm)
    {
      int j3 = j * 3;
      f[j3+0] += - fr * ex;
      f[j3+1] += - fr * ey;
      f[j3+2] += - fr * ez;
    }
}

/*
 * INPUT
 *  sys        : struct stokes (only nm and pos are used)
 *  ev         : struct EV
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_calc_force (struct stokes *sys,
	       struct EV *ev,
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
      for (j = 0; j < sys->np; j ++)
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
		  double xx = x - sys->llx[k];
		  double yy = y - sys->lly[k];
		  double zz = z - sys->llz[k];

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
