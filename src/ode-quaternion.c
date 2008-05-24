/* ODE utility routines for angle by quaternion
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ode-quaternion.c,v 1.8 2008/05/24 05:58:56 kichiki Exp $
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h> // sin(), cos()
#include <gsl/gsl_errno.h> // GSL_SUCCESS
#include "memory-check.h" // macro CHECK_MALLOC

#include "dgetri_c.h" /* lapack_inv_() */
#include "bonds.h" // bonds_calc_force()

#include "ode.h" // solve_mix_3all()
#include "ode-quaternion.h"


void
quaternion_W (const double *Q, double *W)
{
  W[ 0] = +Q[3];
  W[ 1] = +Q[2];
  W[ 2] = -Q[1];
  W[ 3] = -Q[0];

  W[ 4] = -Q[2];
  W[ 5] = +Q[3];
  W[ 6] = +Q[0];
  W[ 7] = -Q[1];

  W[ 8] = +Q[1];
  W[ 9] = -Q[0];
  W[10] = +Q[3];
  W[11] = -Q[2];

  W[12] = +Q[0];
  W[13] = +Q[1];
  W[14] = +Q[2];
  W[15] = +Q[3];
}

void
quaternion_Winv_lapack (const double *Q, double *W)
{
  quaternion_W (Q, W);
  lapack_inv_ (4, W);
}

void
quaternion_Winv (const double *Q, double *W)
{
  double den = Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];

  W[ 0] = +Q[3] / den;
  W[ 1] = -Q[2] / den;
  W[ 2] = +Q[1] / den;
  W[ 3] = +Q[0] / den;

  W[ 4] = +Q[2] / den;
  W[ 5] = +Q[3] / den;
  W[ 6] = -Q[0] / den;
  W[ 7] = +Q[1] / den;

  W[ 8] = -Q[1] / den;
  W[ 9] = +Q[0] / den;
  W[10] = +Q[3] / den;
  W[11] = +Q[2] / den;

  W[12] = -Q[0] / den;
  W[13] = -Q[1] / den;
  W[14] = -Q[2] / den;
  W[15] = +Q[3] / den;
}

void
quaternion_dQdt (const double *Q, const double *O,
		 double *dQdt)
{
  double Winv[16];
  quaternion_Winv (Q, Winv);

  int i, j;
  for (i = 0; i < 4; i++)
    {
      dQdt[i] = 0.0;
      for (j = 0; j < 3; j ++)
	{
	  dQdt[i] += 0.5 * Winv[i*4+j] * O[j];
	}
    }
}

void
quaternion_by_Euler (double phi, double theta, double psi,
		     double *Q)
{
  double ct = cos (theta);
  double st = sin (theta);
  double cm = cos (0.5*(phi - psi));
  double cp = cos (0.5*(phi + psi));
  double sm = sin (0.5*(phi - psi));
  double sp = sin (0.5*(phi + psi));

  Q[0] = st * cm;
  Q[1] = st * sm;
  Q[2] = ct * sp;
  Q[3] = ct * cp;
}


/* calc dydt for center and angle by gsl_odeiv with hydrodynamics
 * without inertia (zero Stokes number)
 * INPUT
 *  t   : current time
 *  y[nm*3 + nm*4] : center position and quaternion of particles at time t
 *  params : (struct ode_params*)ode.
 *           the following parameters are used here;
 *           ode->sys       : (struct stokes *)
 *           ode->F [np*3]
 *           ode->T [np*3]
 *           ode->E [np*5]
 *           ode->uf [np*3]
 *           ode->of [np*3]
 *           ode->ef [np*5]
 *           ode->flag_mat
 *           ode->flag_lub
 *           ode->bonds     : (struct bonds *)
 * OUTPUT
 *  dydt[] := (d/dt) x(t) = U(t)
 *         where U(t) := R^-1 . F, the terminal velocity
 */
int
dydt_Q_hydro (double t, const double *y, double *dydt,
	      void *params)
{
  // get the parameters
  struct ode_params *ode = (struct ode_params *)params;


  // prepare variables for mobile particles
  int nm = ode->sys->nm;
  int nm3 = nm * 3;
  int nm5 = nm * 5;
  double *O = (double *)malloc (sizeof (double) * nm3);
  double *S = (double *)malloc (sizeof (double) * nm5);
  CHECK_MALLOC (O, "dydt_hydro");
  CHECK_MALLOC (S, "dydt_hydro");

  // prepare variables for fixed particles
  int nf = ode->sys->np - nm;
  int nf3 = nf * 3;
  int nf5 = nf * 5;
  double *ff = NULL;
  double *tf = NULL;
  double *sf = NULL;
  if (nf > 0)
    {
      ff = (double *)malloc (sizeof (double) * nf3);
      tf = (double *)malloc (sizeof (double) * nf3);
      sf = (double *)malloc (sizeof (double) * nf5);
      CHECK_MALLOC (ff, "dydt_hydro");
      CHECK_MALLOC (tf, "dydt_hydro");
      CHECK_MALLOC (sf, "dydt_hydro");
    }


  double *fm = (double *)malloc (sizeof (double) * nm3);
  CHECK_MALLOC (fm, "dydt_hydro_st");

  int i;
  for (i = 0; i < nm3; i ++)
    {
      fm[i] = ode->F[i];
    }


  // set the current configuration into struct stokes "sys"
  stokes_set_pos_mobile (ode->sys, y);
  stokes_set_shear_shift (ode->sys, t, ode->t0, ode->s0);

  if (ode->bonds->n > 0)
    {
      // calc force on the mobile particles
      bonds_calc_force (ode->bonds, ode->sys,
			fm, 1/* add */);
    }

  // set dydt = U = R^-1 . F
  solve_mix_3all (ode->sys,
		  ode->flag_noHI,
		  ode->flag_lub, ode->flag_mat,
		  fm, ode->T, ode->E, ode->uf, ode->of, ode->ef,
		  dydt, O, S, ff, tf, sf);

  // calc dQdt
  double       *dQdt = dydt + nm3;
  const double *Q    = y    + nm3;
  for (i = 0; i < nm; i ++)
    {
      int i3 = i * 3;
      int i4 = i * 4;
      quaternion_dQdt (Q + i4, O + i3, dQdt + i4);
    }

  // house-keeping
  free (fm);
  free (O);
  free (S);
  if (nf > 0)
    {
      free (ff);
      free (tf);
      free (sf);
    }

  return GSL_SUCCESS;
}

/* calc dydt for center and angle by gsl_odeiv with hydrodynamics
 * with inertia (scaled Stokes number)
 * INPUT
 *  t   : current time
 *  y[nm*3 + nm*3 + nm*4] :
 *           position, (real) velocity, and quaternion of particles at time t
 *  params : (struct ode_params*)ode.
 *           the following parameters are used here;
 *           ode->sys       : (struct stokes *)
 *           ode->F [np*3]
 *           ode->T [np*3]
 *           ode->E [np*5]
 *           ode->uf [np*3]
 *           ode->of [np*3]
 *           ode->ef [np*5]
 *           ode->flag_mat
 *           ode->flag_lub
 *           ode->stokes
 *           ode->bonds     : (struct bonds *)
 * OUTPUT
 *  dydt[] := (d/dt) y(t), where y(t) = (x(t), U(t)).
 *         (d/dt) x(t) = U(t)
 *         (d/dt) U(t) = (1/stokes)(-U(t) + V(t)),
 *         where V(t) := R^-1 . F, the terminal velocity
 */
int
dydt_Q_hydro_st (double t, const double *y, double *dydt,
		 void *params)
{
  // get the parameters
  struct ode_params *ode = (struct ode_params *)params;


  // prepare variables for mobile particles
  int nm = ode->sys->nm;
  int nm3 = nm * 3;
  int nm5 = nm * 5;
  double *V = (double *)malloc (sizeof (double) * nm3);
  double *O = (double *)malloc (sizeof (double) * nm3);
  double *S = (double *)malloc (sizeof (double) * nm5);
  CHECK_MALLOC (V, "dydt_hydro_st");
  CHECK_MALLOC (O, "dydt_hydro_st");
  CHECK_MALLOC (S, "dydt_hydro_st");

  // prepare variables for fixed particles
  int nf = ode->sys->np - nm;
  int nf3 = nf * 3;
  int nf5 = nf * 5;
  double *ff = NULL;
  double *tf = NULL;
  double *sf = NULL;
  if (nf > 0)
    {
      ff = (double *)malloc (sizeof (double) * nf3);
      tf = (double *)malloc (sizeof (double) * nf3);
      sf = (double *)malloc (sizeof (double) * nf5);
      CHECK_MALLOC (ff, "dydt_hydro_st");
      CHECK_MALLOC (tf, "dydt_hydro_st");
      CHECK_MALLOC (sf, "dydt_hydro_st");
    }

  double *fm = (double *)malloc (sizeof (double) * nm3);
  CHECK_MALLOC (fm, "dydt_hydro_st");

  int i;
  for (i = 0; i < nm3; i ++)
    {
      fm[i] = ode->F[i];
    }

  // set the current configuration into struct stokes "sys"
  stokes_set_pos_mobile (ode->sys, y);
  stokes_set_shear_shift (ode->sys, t, ode->t0, ode->s0);

  if (ode->bonds->n > 0)
    {
      // calc force on the mobile particles
      bonds_calc_force (ode->bonds, ode->sys,
			fm, 1/* add */);
    }

  // calculate the terminal velocity V[]
  solve_mix_3all (ode->sys,
		  ode->flag_noHI,
		  ode->flag_lub, ode->flag_mat,
		  fm, ode->T, ode->E, ode->uf, ode->of, ode->ef,
		  V, O, S, ff, tf, sf);

  for (i = 0; i < nm3; i ++)
    {
      dydt [i]       =   y[nm3 + i];
      dydt [nm3 + i] = (-y[nm3 + i] + V[i]) / ode->st;
    }

  // calc dQdt
  int nm6 = nm * 6;
  double       *dQdt = dydt + nm6;
  const double *Q    = y    + nm6;
  for (i = 0; i < nm; i ++)
    {
      int i3 = i * 3;
      int i4 = i * 4;
      quaternion_dQdt (Q + i4, O + i3, dQdt + i4);
    }

  // house-keeping
  free (fm);
  free (V);
  free (O);
  free (S);
  if (nf > 0)
    {
      free (ff);
      free (tf);
      free (sf);
    }

  return GSL_SUCCESS;
}

