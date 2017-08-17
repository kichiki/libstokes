/* excluded-volume interactions in Lennard-Jones type
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-LJ.c,v 1.2 2008/06/13 03:06:04 kichiki Exp $
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

#include "ev-LJ.h"


/* initialize struct EV_LJ
 * INPUT
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV_LJ,
 *      where LJ parameters are set by zero.
 */
struct EV_LJ *
EV_LJ_init (int np)
{
  struct EV_LJ *ev_LJ = (struct EV_LJ *)malloc (sizeof (struct EV_LJ));
  CHECK_MALLOC (ev_LJ, "EV_LJ_init");

  ev_LJ->n    = np;
  ev_LJ->flag = (int *)malloc (sizeof (int) * np);
  ev_LJ->e    = (double *)malloc (sizeof (double) * np);
  ev_LJ->r0   = (double *)malloc (sizeof (double) * np);
  CHECK_MALLOC (ev_LJ->flag, "EV_LJ_init");
  CHECK_MALLOC (ev_LJ->e,    "EV_LJ_init");
  CHECK_MALLOC (ev_LJ->r0,   "EV_LJ_init");

  // zero clear
  int i;
  for (i = 0; i < np; i ++)
    {
      ev_LJ->flag[i] = 0;
      ev_LJ->e   [i] = 0.0;
      ev_LJ->r0  [i] = 0.0;
    }

  return (ev_LJ);
}

void
EV_LJ_free (struct EV_LJ *ev_LJ)
{
  if (ev_LJ == NULL) return;
  if (ev_LJ->flag != NULL) free (ev_LJ->flag);
  if (ev_LJ->e    != NULL) free (ev_LJ->e);
  if (ev_LJ->r0   != NULL) free (ev_LJ->r0);
  free (ev_LJ);
}

/* scale LJ parameter e for runs
 * INPUT
 *  ev_LJ  : struct EV_LJ
 *  peclet : peclet number
 */
void
EV_LJ_scale (struct EV_LJ *ev_LJ,
	     double peclet)
{
  if (ev_LJ == NULL) return;

  int i;
  for (i = 0; i < ev_LJ->n; i ++)
    {
      if (ev_LJ->flag[i] == 0)
	{
	  ev_LJ->flag[i] = 1;
	  ev_LJ->e [i] /= peclet * ev_LJ->r0[i];
	}
    }
}


/*
 * INPUT
 *  x : depth of overlap
 *      therefore, the LJ distance is given by 
 *        r = r0 - x.
 *      x have to be smaller than r0.
 *      zero is returned if x < 0 (no overlap).
 * OUTPUT
 *  (returned value)
 */
static double
EV_LJ_fr (double e, double r0, double x)
{
  if (x < 0.0) return (0.0);

  // NOTE: r0 is converted to hat(r0) = r0 / length
  double r = r0 - x;
  if (r < 0.0)
    {
      fprintf (stderr, "EV_LJ_fr: too much overlap"
	       " (x=%e > r0=%e)"
	       " r=r0-x=%e is adjusted to 1.0e-12.\n",
	       x, r0, r);
      r = 1.0e-12; // give some small value
    }

  double ri = r0 / r;
  double ri6 = pow (ri, 6.0);
  double ri7 = ri6 * ri;

  // NOTE: e is converted to hat(e)/(peclet * hat(r0))
  double fr = 12.0 * e * ri7 * (ri6 - 1.0);

  return (fr);
}

static void
EV_LJ_coef_ij (struct EV_LJ *ev_LJ,
	       int i, int j,
	       double *e, double *r0)
{
  *e  = sqrt (ev_LJ->e[i] * ev_LJ->e[j]);
  *r0 = 0.5 * (ev_LJ->r0[i] + ev_LJ->r0[j]);
}

static void
EV_LJ_set_force_ij (struct stokes *sys,
		    struct EV_LJ *ev_LJ,
		    int i, int j,
		    double overlap,
		    double r, double x, double y, double z,
		    double *f)
{
  // get the LJ coef. for (i-j) pair
  double e;
  double r0;
  EV_LJ_coef_ij (ev_LJ, i, j, &e, &r0);

  // define fr
  double fr = EV_LJ_fr (e, r0, overlap);

  double ex = x / r;
  double ey = y / r;
  double ez = z / r;


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
 *  ev_LJ      : struct EV_LJ
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_LJ_calc_force (struct EV_LJ *ev_LJ,
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
  double r;
  double overlap;
  for (i = 0; i < sys->nm; i ++)
    {
      double ai;
      if (sys->a == NULL) ai = 1.0;
      else                ai = sys->a[i];
      int i3 = i * 3;
      int j;
      // loop for each pair (i <= j)
      for (j = i; j < sys->np; j ++)
	{
	  double aj;
	  if (sys->a == NULL) aj = 1.0;
	  else                aj = sys->a[j];
	  int j3 = j * 3;
	  x = sys->pos [i3+0] - sys->pos [j3+0];
	  y = sys->pos [i3+1] - sys->pos [j3+1];
	  z = sys->pos [i3+2] - sys->pos [j3+2];

	  if (i != j)
	    {
	      r = sqrt (x*x + y*y + z*z);
	      overlap = ai + aj - r;
	      if (overlap > 0.0)
		{
		  EV_LJ_set_force_ij (sys, ev_LJ, i, j,
				      overlap, r, x, y, z,
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

		  r = sqrt (xx*xx + yy*yy + zz*zz);
		  overlap = ai + aj - r;
		  if (overlap > 0.0)
		    {
		      EV_LJ_set_force_ij (sys, ev_LJ, i, j,
					  overlap, r, x, y, z,
					  f);
		    }
		}
	    }
	}
    }
}
