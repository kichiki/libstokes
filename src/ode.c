/* ODE utility routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ode.c,v 1.6 2007/12/05 03:47:25 kichiki Exp $
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
#include <gsl/gsl_errno.h> // GSL_SUCCESS
#include "memory-check.h"

#include "stokes.h" // struct stokes
#include "ewald-3f.h"
#include "ewald-3ft.h"
#include "ewald-3fts.h"
#include "ewald-3f-matrix.h"
#include "ewald-3ft-matrix.h"
#include "ewald-3fts-matrix.h"
#include "bonds.h" // struct bonds
#include "coll.h"
#include "periodicity.h" // check_periodic()


#include "ode.h"


/* wrapper for solving mob and mix problems in ODE routines
 */
void
solve_mix_3all (struct stokes * sys,
		int flag_lub,
		int flag_mat,
		const double *fm, const double *tm, const double *em,
		const double *uf, const double *of, const double *ef,
		double *um, double *om, double *sm,
		double *ff, double *tf, double *sf)
{
  if (sys->version == 0) // F version
    {
      if (flag_lub == 0)
	{
	  if (flag_mat == 0)
	    {
	      solve_mix_3f
		(sys,
		 fm, uf,
		 um, ff);
	    }
	  else
	    {
	      solve_mix_3f_matrix
		(sys,
		 fm, uf,
		 um, ff);
	    }
	}
      else
	{
	  if (flag_mat == 0)
	    {
	      solve_mix_lub_3f
		(sys,
		 fm, uf,
		 um, ff);
	    }
	  else
	    {
	      solve_mix_lub_3f_matrix
		(sys,
		 fm, uf,
		 um, ff);
	    }
	}
    }
  else if (sys->version == 1) // FT version
    {
      if (flag_lub == 0)
	{
	  if (flag_mat == 0)
	    {
	      solve_mix_3ft
		(sys,
		 fm, tm, uf, of,
		 um, om, ff, tf);
	    }
	  else
	    {
	      solve_mix_3ft_matrix
		(sys,
		 fm, tm, uf, of,
		 um, om, ff, tf);
	    }
	}
      else
	{
	  if (flag_mat == 0)
	    {
	      solve_mix_lub_3ft
		(sys,
		 fm, tm, uf, of,
		 um, om, ff, tf);
	    }
	  else
	    {
	      solve_mix_lub_3ft_matrix
		(sys,
		 fm, tm, uf, of,
		 um, om, ff, tf);
	    }
	}
    }
  else // FTS version
    {
      if (flag_lub == 0)
	{
	  if (flag_mat == 0)
	    {
	      solve_mix_3fts
		(sys,
		 fm, tm, em, uf, of, ef,
		 um, om, sm, ff, tf, sf);
	    }
	  else
	    {
	      solve_mix_3fts_matrix
		(sys,
		 fm, tm, em, uf, of, ef,
		 um, om, sm, ff, tf, sf);
	    }
	}
      else
	{
	  if (flag_mat == 0)
	    {
	      solve_mix_lub_3fts
		(sys,
		 fm, tm, em, uf, of, ef,
		 um, om, sm, ff, tf, sf);
	    }
	  else
	    {
	      solve_mix_lub_3fts_matrix
		(sys,
		 fm, tm, em, uf, of, ef,
		 um, om, sm, ff, tf, sf);
	    }
	}
    }
}

/* set the parameters to struct ode_params
 * INPUT
 *  (struct stokes *)sys
 *  F [np*3]
 *  T [np*3]
 *  E [np*5]
 *  uf [np*3]
 *  of [np*3]
 *  ef [np*5]
 *  (int) flag_mat
 *  (int) flag_lub
 *  (double) stokes
 *  (struct bonds *)bonds
 *  (double) gamma (for the bond relaxation scheme)
 * OUTPUT :
 *  (struct ode_params) params
 */
struct ode_params *
ode_params_init (struct stokes *sys,
		 double *F,
		 double *T,
		 double *E,
		 double *uf,
		 double *of,
		 double *ef,
		 int flag_lub,
		 int flag_mat,
		 double st,
		 struct bonds *bonds,
		 double gamma)
{
  struct ode_params *params
    = (struct ode_params *)malloc (sizeof (struct ode_params));
  CHECK_MALLOC (params, "ode_params_init");

  params->sys       = sys;
  params->F         = F;
  params->T         = T;
  params->E         = E;
  params->uf        = uf;
  params->of        = of;
  params->ef        = ef;
  params->flag_lub  = flag_lub;
  params->flag_mat  = flag_mat;
  params->st        = st;
  params->bonds     = bonds;
  params->gamma     = gamma;

  return (params);
}


void 
ode_params_free (struct ode_params *params)
{
  if (params != NULL) free (params);
}


/* calc dydt for gsl_odeiv ONLY with bond interactions for relaxation
 * this is equivalent to dydt_relax_bond with a constant friction, where
 *  dx/dt = U = f_bond / gamma
 * INPUT
 *  t   : current time
 *  y[] : position of particles at time t
 *  params : (struct ode_params*)ode.
 *           the following parameters are used here;
 *           ode->sys       : (struct stokes *)
 *           ode->bonds     : (struct bonds *)
 *           ode->gamma     : the friction coef
 * OUTPUT
 *  f[] := (d/dt) y (t), that is, the velocity of particles at time t
 */
int
dydt_relax_bond (double t, const double *y, double *f,
		 void *params)
{
  // get the parameters
  struct ode_params *ode = (struct ode_params *)params;


  stokes_set_pos_mobile (ode->sys, y);

  bonds_calc_force (ode->bonds, ode->sys, f, 0/* zero-clear */);

  // scale by the friction coefficient
  int i;
  for (i = 0; i < ode->sys->nm * 3; i ++)
    {
      f [i] /= ode->gamma;
    }

  return GSL_SUCCESS;
}


/* calc dydt for gsl_odeiv with hydrodynamics
 * with inertia (scaled Stokes number)
 * INPUT
 *  t   : current time
 *  y[] : position and (real) velocity of particles at time t
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
dydt_hydro_st (double t, const double *y, double *dydt,
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

  if (ode->bonds->n > 0)
    {
      // calc force on the mobile particles
      bonds_calc_force (ode->bonds, ode->sys,
			fm, 1/* add */);
    }

  // calculate the terminal velocity V[]
  solve_mix_3all (ode->sys,
		  ode->flag_lub, ode->flag_mat,
		  fm, ode->T, ode->E, ode->uf, ode->of, ode->ef,
		  V, O, S, ff, tf, sf);

  for (i = 0; i < nm3; i ++)
    {
      dydt [i]       =   y[nm3 + i];
      dydt [nm3 + i] = (-y[nm3 + i] + V[i]) / ode->st;
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

/* calc dydt for gsl_odeiv with hydrodynamics
 * without inertia (zero Stokes number)
 * INPUT
 *  t   : current time
 *  y[] : position of particles at time t
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
dydt_hydro (double t, const double *y, double *dydt,
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

  if (ode->bonds->n > 0)
    {
      // calc force on the mobile particles
      bonds_calc_force (ode->bonds, ode->sys,
			fm, 1/* add */);
    }

  // set dydt = U = R^-1 . F
  solve_mix_3all (ode->sys,
		  ode->flag_lub, ode->flag_mat,
		  fm, ode->T, ode->E, ode->uf, ode->of, ode->ef,
		  dydt, O, S, ff, tf, sf);

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
