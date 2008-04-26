/* excluded-volume interactions by Debye-Huckel
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-dh.c,v 1.2 2008/04/26 05:09:39 kichiki Exp $
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

#include "ev-dh.h"


/* initialize struct EV
 * INPUT
 *  a  : characteristic length (in the same dimension for rd below, usually nm)
 *  pe : peclet number
 *  rd : Debye length (in the same dimension for a above, usually nm)
 *  T  : temperature in Kelvin.
 *  e  : dielectric constant of the solution
 *  r2 : square of the max distance for F^{EV_DH}
 *       (in the same dimension squared for a above, usually nm^2)
 *  np : number of particles
 * OUTPUT
 *  returned value : struct EV_DH,
 *      where only the system parameters (a_sys, rd) are defined.
 *      nu[i] is zero cleared.
 */
struct EV_DH *
EV_DH_init (double a, double pe, double rd, double T, double e,
	    double r2,
	    int np)
{
  struct EV_DH *ev_dh = (struct EV_DH *)malloc (sizeof (struct EV_DH));
  CHECK_MALLOC (ev_dh, "EV_DH_init");

  ev_dh->r2 = r2 / (a * a);

  // Boltzmann constant
  double kB = 1.3806504e-23; // [N m / K]
  // elementary charge
  double ec = 1.602176487e-19; // [C]
  // permittivity of vacuum
  double e0 = 8.8541878176e-12; // [F / m] = [C^2 / N m^2]

  double kT = kB * T; // [N m]
  double e2by4pie0 = ec * ec / (4.0 * M_PI * e * e0); // [N m^2]

  ev_dh->a_sys = e2by4pie0 / (pe * kT * a); // dimensionless
  ev_dh->rd    = rd / a; // dimensionless

  ev_dh->n     = np;
  ev_dh->nu = (double *)malloc (sizeof (double) * np);
  CHECK_MALLOC (ev_dh->nu, "EV_DH_init");

  // zero clear
  int i;
  for (i = 0; i < np; i ++)
    {
      ev_dh->nu [i] = 0.0;
    }

  return (ev_dh);
}

void
EV_DH_free (struct EV_DH *ev_dh)
{
  if (ev_dh == NULL) return;
  if (ev_dh->nu != NULL) free (ev_dh->nu);
  free (ev_dh);
}

static void
EV_DH_set_force_ij (struct stokes *sys,
		    struct EV_DH *ev_dh,
		    int i, int j, double r2, double x, double y, double z,
		    double *f)
{
  double r = sqrt (r2);
  double ex = x / r;
  double ey = y / r;
  double ez = z / r;

  double kr = r / ev_dh->rd;

  // define fr
  double fr
    = ev_dh->a_sys
    * ev_dh->nu[i] * ev_dh->nu[j]
    / r2 * (1.0 + kr) * exp (-kr);

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
 *  ev_dh      : struct EV_DH
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_DH_calc_force (struct EV_DH *ev_dh,
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
      for (j = 0; j < sys->np; j ++)
	{
	  int ib3 = j * 3;
	  x = sys->pos [ia3+0] - sys->pos [ib3+0];
	  y = sys->pos [ia3+1] - sys->pos [ib3+1];
	  z = sys->pos [ia3+2] - sys->pos [ib3+2];

	  if (i != j)
	    {
	      r2 = x*x + y*y + z*z;
	      if (r2 < ev_dh->r2)
		{
		  EV_DH_set_force_ij (sys, ev_dh, i, j, r2, x, y, z,
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
		  if (r2 < ev_dh->r2)
		    {
		      EV_DH_set_force_ij (sys, ev_dh, i, j, r2, xx, yy, zz,
					  f);
		    }
		}
	    }
	}
    }
}
