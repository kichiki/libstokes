/* angle interaction between particles
 * Copyright (C) 2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include <math.h> // M_PI, sqrt()

#include "memory-check.h" // macro CHECK_MALLOC

#include "angles.h" // struct angle


/**
 * struct angles
 */
struct angles *
angles_init (void)
{
  struct angles *ang = (struct angles *)malloc (sizeof (struct angles));
  CHECK_MALLOC (ang, "angles_init");
  
  ang->n = 0;
  ang->a = NULL;

  return (ang);
}

void
angles_free (struct angles *ang)
{
  if (ang == NULL) return;
  if (ang->a != NULL) free (ang->a);
  free (ang);
}

/*
 * INPUT
 *  ia, ib, ic : particle indices (ib is the center particle)
 *  k  : potential factor
 *  t0 : natural angle (in radian)
 *  scale : flag for scale (0 == k is scaled,
 *                          1 == k is not scaled yet.)
 */
void
angles_add (struct angles *ang,
	    int ia, int ib, int ic,
	    double k, double t0, int scale)
{
  ang->n ++;

  int n = ang->n;
  ang->a = (struct angle *)realloc (ang->a, sizeof (struct angle) * n);
  CHECK_MALLOC (ang->a, "angles_add");

  // set n as the newly added element
  n--;
  struct angle *a = ang->a + n;

  // potential factor
  a->k  = k;
  // natural angle (theta_0)
  a->t0 = t0;
  // scale flag
  a->scale = scale;

  // particle indices
  a->ia = ia;
  a->ib = ib;
  a->ic = ic;
}

/* scale parameter by the Peclet number
 * INPUT
 *  ang : struct angles
 *        (ang->a [n])->k is scaled as (k / pe)
 *  a     : length scale in the simulation
 *  pe    : peclet number
 * OUTPUT
 *  (ang->a [i])->k : scaled as (k / pe)
 */
void
angles_scale_k (struct angles *ang,
		double a, double pe)
{
  int i;
  struct angle *ai;
  for (i = 0, ai = ang->a;
       i < ang->n;
       i++, ai++)
    {
      if (ai->scale == 1)
	{
	  ai->k /= pe;
	  ai->scale = 0; // scaled
	}
    }
}

/** copy from bonds.c **/
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
 *  ang        : struct angles
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
angles_calc_force (struct angles *ang,
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

  // calc angles
  double Ax[3];
  double Az[3];

  double X[3]; /* X := x_a - x_b */
  double Z[3]; /* Z := x_c - x_b */
  struct angle *a;
  for (i = 0, a = ang->a;
       i < ang->n;
       i++, a++)
    {
      /* note ib is the center particle */
      int ix3 = 3 * a->ia;
      int iy3 = 3 * a->ib;
      int iz3 = 3 * a->ic;

      X[0] = sys->pos [ix3+0] - sys->pos [iy3+0];
      X[1] = sys->pos [ix3+1] - sys->pos [iy3+1];
      X[2] = sys->pos [ix3+2] - sys->pos [iy3+2];

      Z[0] = sys->pos [iz3+0] - sys->pos [iy3+0];
      Z[1] = sys->pos [iz3+1] - sys->pos [iy3+1];
      Z[2] = sys->pos [iz3+2] - sys->pos [iy3+2];

      if (sys->periodic != 0)
	{
	  search_close_image (sys, &X[0], &X[1], &X[2]);
	  search_close_image (sys, &Z[0], &Z[1], &Z[2]);
	}

      double XX;
      double ZZ;

      XX
	= X[0] * X[0]
	+ X[1] * X[1]
	+ X[2] * X[2];
      XX = sqrt (XX);
      X[0] /= XX;
      X[1] /= XX;
      X[2] /= XX;

    REDO_angles_calc_force:
      ZZ
	= Z[0] * Z[0]
	+ Z[1] * Z[1]
	+ Z[2] * Z[2];
      ZZ = sqrt (ZZ);
      Z[0] /= ZZ;
      Z[1] /= ZZ;
      Z[2] /= ZZ;

      // cos theta = - (X . Z)/(|X| |Z|)
      double ct
	= - (X[0] * Z[0]
	     + X[1] * Z[1]
	     + X[2] * Z[2]);

      
      if (ct == 1.0 && a->t0 == 0.0)
	{
	  // force is zero
	  continue;
	}

      double st;// sin theta
      double t; // theta
      if (ct >= +1.0 || ct <= -1.0)
	{
	  // vectors AB and BC are parallel...
	  fprintf (stderr, "# angles_calc_force() "
		   ": cos(theta) = %e, parallel bonds for angle"
		   " for the particles %d, %d, %d\n",
		   ct,
		   a->ia, a->ib, a->ic);

	  // recover to the real vector
	  Z[0] *= ZZ;
	  Z[1] *= ZZ;
	  Z[2] *= ZZ;

	  // give infinitesimal artificial displacement for the vector BC
	  double tiny = 1.0e-7;
	  Z[0] += tiny * (drand48 () - 0.5);
	  Z[1] += tiny * (drand48 () - 0.5);
	  Z[2] += tiny * (drand48 () - 0.5);

	  goto REDO_angles_calc_force;
	}
      else
	{
	  t = acos (ct);
	  st = sqrt (1.0 - ct * ct);
	  // st = sin (t);
	}

      double factor = - a->k * (t - a->t0) / st;

      /* Ax = A(-Z, x) */
      Ax[0] = factor / XX* (Z[0]
			    - X[0] * X[0] * Z[0]
			    - X[0] * X[1] * Z[1]
			    - X[0] * X[2] * Z[2]);
      Ax[1] = factor / XX* (Z[1]
			    - X[1] * X[0] * Z[0]
			    - X[1] * X[1] * Z[1]
			    - X[1] * X[2] * Z[2]);
      Ax[2] = factor / XX* (Z[2]
			    - X[2] * X[0] * Z[0]
			    - X[2] * X[1] * Z[1]
			    - X[2] * X[2] * Z[2]);

      /* Az = A(-X, Z) */
      Az[0] = factor / ZZ* (X[0]
			    - Z[0] * Z[0] * X[0]
			    - Z[0] * Z[1] * X[1]
			    - Z[0] * Z[2] * X[2]);
      Az[1] = factor / ZZ* (X[1]
			    - Z[1] * Z[0] * X[0]
			    - Z[1] * Z[1] * X[1]
			    - Z[1] * Z[2] * X[2]);
      Az[2] = factor / ZZ* (X[2]
			    - Z[2] * Z[0] * X[0]
			    - Z[2] * Z[1] * X[1]
			    - Z[2] * Z[2] * X[2]);

      // X (ia)
      if (a->ia < sys->nm)
	{
	  f[ix3+0] += Ax[0];
	  f[ix3+1] += Ax[1];
	  f[ix3+2] += Ax[2];
	}
      // Y (ib) == the center particle
      if (a->ib < sys->nm)
	{
	  f[iy3+0] += - Ax[0] - Az[0];
	  f[iy3+1] += - Ax[1] - Az[1];
	  f[iy3+2] += - Ax[2] - Az[2];
	}
      // Z (ic)
      if (a->ic < sys->nm)
	{
	  // top particle
	  f[iz3+0] += Az[0];
	  f[iz3+1] += Az[1];
	  f[iz3+2] += Az[2];
	}
    }
}
