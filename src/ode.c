/* ODE utility routines
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ode.c,v 1.13 2008/05/24 05:59:27 kichiki Exp $
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
#include "memory-check.h" // macro CHECK_MALLOC

#include "stokes.h" // struct stokes
#include "ewald-3f.h"
#include "ewald-3ft.h"
#include "ewald-3fts.h"
#include "ewald-3f-matrix.h"
#include "ewald-3ft-matrix.h"
#include "ewald-3fts-matrix.h"
#include "bonds.h" // struct bonds
#include "noHI.h" // solve_mix_3[f|ft|fts]_noHI()

#include "ode.h"


/* old wrapper for solving mob and mix problems in ODE routines
 * without flag_noHI in the argument
 */
static void
solve_mix_3all_HI (struct stokes *sys,
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

/* wrapper for solving mob and mix problems in ODE routines
 * for noHI version
 */
static void
solve_mix_3all_noHI (struct stokes *sys,
		     const double *fm, const double *tm, const double *em,
		     const double *uf, const double *of, const double *ef,
		     double *um, double *om, double *sm,
		     double *ff, double *tf, double *sf)
{
  if (sys->version == 0)
    {
      // F version
      solve_mix_3f_noHI
	(sys,
	 fm, uf,
	 um, ff);
    }
  else if (sys->version == 1)
    {
      // FT version
      solve_mix_3ft_noHI
	(sys,
	 fm, tm, uf, of,
	 um, om, ff, tf);
    }
  else
    {
      // FTS version
      solve_mix_3fts_noHI
	(sys,
	 fm, tm, em, uf, of, ef,
	 um, om, sm, ff, tf, sf);
    }
}


/**
 ** conversion between Global scaling and Particle-wise scaling
 **/
static void
scaling_G_to_P_fm (struct stokes *sys,
		   const double *fm_G,
		   double *fm_P)
{
  int i;
  for (i = 0; i < sys->nm; i ++)
    {
      int ix = i*3;
      int iy = ix + 1;
      int iz = ix + 2;

      fm_P[ix] = fm_G[ix];
      fm_P[iy] = fm_G[iy];
      fm_P[iz] = fm_G[iz];
    }
}
static void
scaling_G_to_P_tm (struct stokes *sys,
		   const double *tm_G,
		   double *tm_P)
{
  if (sys->version <= 0)
    {
      // F version
      return;
    }

  int i;
  for (i = 0; i < sys->nm; i ++)
    {
      int ix = i*3;
      int iy = ix + 1;
      int iz = ix + 2;

      tm_P[ix] = tm_G[ix];
      tm_P[iy] = tm_G[iy];
      tm_P[iz] = tm_G[iz];
    }
}
static void
scaling_G_to_P_em (struct stokes *sys,
		   const double *em_G,
		   double *em_P)
{
  if (sys->version <= 1)
    {
      // F or FT versions
      return;
    }

  int i;
  for (i = 0; i < sys->nm; i ++)
    {
      double a = sys->a[i];
      double a3 = a * a * a;

      int i1 = i*5;
      int i2 = i1 + 1;
      int i3 = i1 + 2;
      int i4 = i1 + 3;
      int i5 = i1 + 4;

      em_P[i1] = em_G[i1] * a3;
      em_P[i2] = em_G[i2] * a3;
      em_P[i3] = em_G[i3] * a3;
      em_P[i4] = em_G[i4] * a3;
      em_P[i5] = em_G[i5] * a3;
    }
}
static void
scaling_G_to_P_uf (struct stokes *sys,
		   const double *uf_G,
		   double *uf_P)
{
  int i;
  for (i = sys->nm; i < sys->np; i ++)
    {
      double a = sys->a[i];

      int ix = i*3;
      int iy = ix + 1;
      int iz = ix + 2;

      uf_P[ix] = uf_G[ix] * a;
      uf_P[iy] = uf_G[iy] * a;
      uf_P[iz] = uf_G[iz] * a;
    }
}
static void
scaling_G_to_P_of (struct stokes *sys,
		   const double *of_G,
		   double *of_P)
{
  if (sys->version <= 0)
    {
      // F version
      return;
    }

  int i;
  for (i = sys->nm; i < sys->np; i ++)
    {
      double a = sys->a[i];
      double a3 = a * a * a;

      int ix = i*3;
      int iy = ix + 1;
      int iz = ix + 2;

      of_P[ix] = of_G[ix] * a3;
      of_P[iy] = of_G[iy] * a3;
      of_P[iz] = of_G[iz] * a3;
    }
}
static void
scaling_G_to_P_ef (struct stokes *sys,
		   const double *ef_G,
		   double *ef_P)
{
  if (sys->version <= 1)
    {
      // F or FT versions
      return;
    }

  int i;
  for (i = sys->nm; i < sys->np; i ++)
    {
      double a = sys->a[i];
      double a3 = a * a * a;

      int i1 = i*5;
      int i2 = i1 + 1;
      int i3 = i1 + 2;
      int i4 = i1 + 3;
      int i5 = i1 + 4;

      ef_P[i1] = ef_G[i1] * a3;
      ef_P[i2] = ef_G[i2] * a3;
      ef_P[i3] = ef_G[i3] * a3;
      ef_P[i4] = ef_G[i4] * a3;
      ef_P[i5] = ef_G[i5] * a3;
    }
}

/* convert Global (G) scaling to Particle-wise (P) scaling
 */
static void
scaling_G_to_P (struct stokes *sys,
		const double *fm_G,
		const double *tm_G,
		const double *em_G,
		const double *uf_G,
		const double *of_G,
		const double *ef_G,
		double *fm_P,
		double *tm_P,
		double *em_P,
		double *uf_P,
		double *of_P,
		double *ef_P)
{
  scaling_G_to_P_fm (sys, fm_G, fm_P);
  scaling_G_to_P_tm (sys, tm_G, tm_P);
  scaling_G_to_P_em (sys, em_G, em_P);

  scaling_G_to_P_uf (sys, uf_G, uf_P);
  scaling_G_to_P_of (sys, of_G, of_P);
  scaling_G_to_P_ef (sys, ef_G, ef_P);
}

static void
scaling_P_to_G_um (struct stokes *sys,
		   const double *um_P,
		   double *um_G)
{
  int i;
  for (i = 0; i < sys->nm; i ++)
    {
      double a = sys->a[i];

      int ix = i*3;
      int iy = ix + 1;
      int iz = ix + 2;

      um_G[ix] = um_P[ix] / a;
      um_G[iy] = um_P[iy] / a;
      um_G[iz] = um_P[iz] / a;
    }
}
static void
scaling_P_to_G_om (struct stokes *sys,
		   const double *om_P,
		   double *om_G)
{
  if (sys->version <= 0)
    {
      // F version
      return;
    }

  int i;
  for (i = 0; i < sys->nm; i ++)
    {
      double a = sys->a[i];
      double a3 = a * a * a;

      int ix = i*3;
      int iy = ix + 1;
      int iz = ix + 2;

      om_G[ix] = om_P[ix] / a3;
      om_G[iy] = om_P[iy] / a3;
      om_G[iz] = om_P[iz] / a3;
    }
}
static void
scaling_P_to_G_sm (struct stokes *sys,
		   const double *sm_P,
		   double *sm_G)
{
  if (sys->version <= 1)
    {
      // F or FT versions
      return;
    }

  int i;
  for (i = 0; i < sys->nm; i ++)
    {
      int i1 = i*5;
      int i2 = i1 + 1;
      int i3 = i1 + 2;
      int i4 = i1 + 3;
      int i5 = i1 + 4;

      sm_G[i1] = sm_P[i1];
      sm_G[i2] = sm_P[i2];
      sm_G[i3] = sm_P[i3];
      sm_G[i4] = sm_P[i4];
      sm_G[i5] = sm_P[i5];
    }
}
static void
scaling_P_to_G_ff (struct stokes *sys,
		   const double *ff_P,
		   double *ff_G)
{
  int i;
  for (i = sys->nm; i < sys->np; i ++)
    {
      int ix = i*3;
      int iy = ix + 1;
      int iz = ix + 2;

      ff_G[ix] = ff_P[ix];
      ff_G[iy] = ff_P[iy];
      ff_G[iz] = ff_P[iz];
    }
}
static void
scaling_P_to_G_tf (struct stokes *sys,
		   const double *tf_P,
		   double *tf_G)
{
  if (sys->version <= 0)
    {
      // F version
      return;
    }

  int i;
  for (i = sys->nm; i < sys->np; i ++)
    {
      int ix = i*3;
      int iy = ix + 1;
      int iz = ix + 2;

      tf_G[ix] = tf_P[ix];
      tf_G[iy] = tf_P[iy];
      tf_G[iz] = tf_P[iz];
    }
}
static void
scaling_P_to_G_sf (struct stokes *sys,
		   const double *sf_P,
		   double *sf_G)
{
  if (sys->version <= 1)
    {
      // F or FT versions
      return;
    }

  int i;
  for (i = sys->nm; i < sys->np; i ++)
    {
      int i1 = i*5;
      int i2 = i1 + 1;
      int i3 = i1 + 2;
      int i4 = i1 + 3;
      int i5 = i1 + 4;

      sf_G[i1] = sf_P[i1];
      sf_G[i2] = sf_P[i2];
      sf_G[i3] = sf_P[i3];
      sf_G[i4] = sf_P[i4];
      sf_G[i5] = sf_P[i5];
    }
}
/* convert Global (G) scaling to Particle-wise (P) scaling
 */
static void
scaling_P_to_G (struct stokes *sys,
		const double *um_P,
		const double *om_P,
		const double *sm_P,
		const double *ff_P,
		const double *tf_P,
		const double *sf_P,
		double *um_G,
		double *om_G,
		double *sm_G,
		double *ff_G,
		double *tf_G,
		double *sf_G)
{
  scaling_P_to_G_um (sys, um_P, um_G);
  scaling_P_to_G_om (sys, om_P, om_G);
  scaling_P_to_G_sm (sys, sm_P, sm_G);

  scaling_P_to_G_ff (sys, ff_P, ff_G);
  scaling_P_to_G_tf (sys, tf_P, tf_G);
  scaling_P_to_G_sf (sys, sf_P, sf_G);
}

/* wrapper for solving mob and mix problems in ODE routines
 * INPUT
 *  sys : struct stokes
 *  flag_noHI : 0 == with hydrodynamic interaction
 *              1 == without hydrodynamic interaction
 *  flag_lub  : 0 == without lubrication
 *              1 == with lubrication
 *  flag_mat  : 0 == O(N^3) approach using matrix explicitly
 *              1 == O(N^2) approach using matrix implicitly
 *  fm, tm, em : F, T, and E for mobile particles (given parameters)
 *  uf, of, ef : U, O, and E for fixed particles (given parameters)
 * OUTPUT
 *  um, om, sm : U, O, and S for mobile particles
 *  ff, tf, sf : F, T, and S for fixed particles
 */
void
solve_mix_3all (struct stokes *sys,
		int flag_noHI,
		int flag_lub,
		int flag_mat,
		const double *fm, const double *tm, const double *em,
		const double *uf, const double *of, const double *ef,
		double *um, double *om, double *sm,
		double *ff, double *tf, double *sf)
{
  if (sys->a == NULL)
    {
      // monodisperse
      if (flag_noHI == 0)
	{
	  // with HI
	  solve_mix_3all_HI (sys, flag_lub, flag_mat,
			     fm, tm, em, uf, of, ef,
			     um, om, sm, ff, tf, sf);
	}
      else
	{
	  // without HI
	  solve_mix_3all_noHI (sys,
			       fm, tm, em, uf, of, ef,
			       um, om, sm, ff, tf, sf);
	}
    }
  else 
    {
      // polydisperse
      int nm = sys->nm;
      int nm3 = nm * 3;
      int nm5 = nm * 5;
      int nf = sys->np - nm;
      int nf3 = nf * 3;
      int nf5 = nf * 5;

      // allocate memory for the Particle-wise scaling
      double *fm_P = NULL;
      double *tm_P = NULL;
      double *em_P = NULL;
      double *uf_P = NULL;
      double *of_P = NULL;
      double *ef_P = NULL;

      double *um_P = NULL;
      double *om_P = NULL;
      double *sm_P = NULL;
      double *ff_P = NULL;
      double *tf_P = NULL;
      double *sf_P = NULL;

      // at least, F version
      fm_P = (double *)malloc (sizeof (double) * nm3);
      um_P = (double *)malloc (sizeof (double) * nm3);
      CHECK_MALLOC (fm_P, "solve_mix_3all");
      CHECK_MALLOC (um_P, "solve_mix_3all");
      if (nf > 0)
	{
	  uf_P = (double *)malloc (sizeof (double) * nf3);
	  ff_P = (double *)malloc (sizeof (double) * nf3);
	  CHECK_MALLOC (uf_P, "solve_mix_3all");
	  CHECK_MALLOC (ff_P, "solve_mix_3all");
	}
      if (sys->version > 0) // at least, FT version
	{
	  tm_P = (double *)malloc (sizeof (double) * nm3);
	  om_P = (double *)malloc (sizeof (double) * nm3);
	  CHECK_MALLOC (tm_P, "solve_mix_3all");
	  CHECK_MALLOC (om_P, "solve_mix_3all");
	  if (nf > 0)
	    {
	      of_P = (double *)malloc (sizeof (double) * nf3);
	      tf_P = (double *)malloc (sizeof (double) * nf3);
	      CHECK_MALLOC (of_P, "solve_mix_3all");
	      CHECK_MALLOC (tf_P, "solve_mix_3all");
	    }
	  if (sys->version > 1) // FTS version
	    {
	      em_P = (double *)malloc (sizeof (double) * nm5);
	      sm_P = (double *)malloc (sizeof (double) * nm5);
	      CHECK_MALLOC (em_P, "solve_mix_3all");
	      CHECK_MALLOC (sm_P, "solve_mix_3all");
	      if (nf > 0)
		{
		  ef_P = (double *)malloc (sizeof (double) * nf5);
		  sf_P = (double *)malloc (sizeof (double) * nf5);
		  CHECK_MALLOC (ef_P, "solve_mix_3all");
		  CHECK_MALLOC (sf_P, "solve_mix_3all");
		}
	    }
	}

      // convert the scaling from Global to Particle-wise
      scaling_G_to_P (sys,
		      fm, tm, em, uf, of, ef,
		      fm_P, tm_P, em_P, uf_P, of_P, ef_P);

      if (flag_noHI == 0)
	{
	  // with HI
	  solve_mix_3all_HI (sys, flag_lub, flag_mat,
			     fm_P, tm_P, em_P, uf_P, of_P, ef_P,
			     um_P, om_P, sm_P, ff_P, tf_P, sf_P);
	}
      else
	{
	  solve_mix_3all_noHI (sys,
			       fm_P, tm_P, em_P, uf_P, of_P, ef_P,
			       um_P, om_P, sm_P, ff_P, tf_P, sf_P);
	}

      // convert the scaling from Particle-wise to Global
      scaling_P_to_G (sys,
		      um_P, om_P, sm_P, ff_P, tf_P, sf_P,
		      um, om, sm, ff, tf, sf);

      if (fm_P != NULL) free (fm_P);
      if (tm_P != NULL) free (tm_P);
      if (em_P != NULL) free (em_P);
      if (uf_P != NULL) free (uf_P);
      if (of_P != NULL) free (of_P);
      if (ef_P != NULL) free (ef_P);

      if (um_P != NULL) free (um_P);
      if (om_P != NULL) free (om_P);
      if (sm_P != NULL) free (sm_P);
      if (ff_P != NULL) free (ff_P);
      if (tf_P != NULL) free (tf_P);
      if (sf_P != NULL) free (sf_P);
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
 *  (int) flag_noHI
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
		 int flag_noHI,
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
  params->flag_noHI = flag_noHI;
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

/* set the reference for cell-shift (shear_mode = 1 and 2)
 * INTPUT
 *  t0 : reference time for s0
 *  s0 : cell shift at time t0 (for shear_mode = 1 or 2)
 * OUTPUT
 *  ode->t0, ode->s0 :
 */
void
ode_set_shear_shift_ref (struct ode_params *ode,
			 double t0, double s0)
{
  ode->t0 = t0;
  ode->s0 = s0;
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
  stokes_set_shear_shift (ode->sys, t, ode->t0, ode->s0);

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
