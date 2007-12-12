/* implicit Brownian dynamics algorithms
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bd-imp.c,v 1.3 2007/12/12 06:27:19 kichiki Exp $
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
#include <math.h>
#include "memory-check.h" // CHECK_MALLOC
#include <gsl/gsl_multiroots.h>

#include <stokes.h> // stokes_set_pos()
#include <ode.h> // solve_mix_3all()
#include <ode-quaternion.h> // quaternion_dQdt()

#include <bonds.h> // bonds_calc_force ()
#include <excluded-volume.h> // EV_calc_force ()

#include <brownian.h> // struct BD_params

#include "bd-imp.h"


/* initialize struct BD_imp
 * INPUT
 *  BD : struct BD_params
 *       note that BDimp->BD is just a pointer to BD in the argument.
 *  itmax : max of iteration for the root-finding
 *  eps   : tolerance for the root-finding
 */
struct BD_imp *
BD_imp_init (struct BD_params *BD,
	     int itmax, double eps)
{
  struct BD_imp *BDimp = (struct BD_imp *)malloc (sizeof (struct BD_imp));
  CHECK_MALLOC (BDimp, "BD_imp_init");

  BDimp->BD = BD;
  BDimp->itmax = itmax;
  BDimp->eps = eps;

  // dt is set by 0 here. use BD_imp_set_dt() to set later.
  BDimp->dt = 0.0;

  int nm = BD->sys->nm;

  // x0[] and q0[] are set by 0 here. use BD_imp_set_xq() to set later.
  BDimp->x0 = (double *)malloc (sizeof (double) * nm * 3);
  CHECK_MALLOC (BDimp->x0, "BD_imp_init");
  int i;
  for (i = 0; i < nm * 3; i ++)
    {
      BDimp->x0[i] = 0.0;
    }

  int nz; // dimension for BDimp->z[]
  if (BD->sys->version == 0)
    {
      nz = nm * 3;
      BDimp->q0 = NULL;
    }
  else
    {
      nz = nm * 6; // z[] is for F and T (not the quaternion)
      BDimp->q0 = (double *)malloc (sizeof (double) * nm * 4);
      CHECK_MALLOC (BDimp->q0, "BD_imp_init");
      for (i = 0; i < nm * 4; i ++)
	{
	  BDimp->q0[i] = 0.0;
	}
    }

  // z[] is set by 0 here. it is defined by BD_imp_set_xq().
  BDimp->z = (double *)malloc (sizeof (double) * nz);
  CHECK_MALLOC (BDimp->z, "BD_imp_init");

  // initialize GSL stuff
  BDimp->T = gsl_multiroot_fsolver_hybrid;

  BDimp->F = (gsl_multiroot_function *)malloc (sizeof (gsl_multiroot_function));
  CHECK_MALLOC (BDimp->F, "BD_imp_init");

  BDimp->F->f = BD_imp_JGdP00_func;
  if (BDimp->BD->sys->version == 0)
    {
      BDimp->F->n = nm * 3;
    }
  else
    {
      BDimp->F->n = nm * 7; // quaternion is used here
    }
  BDimp->F->params = BDimp;

  BDimp->S = gsl_multiroot_fsolver_alloc (BDimp->T, BDimp->F->n);
  CHECK_MALLOC (BDimp->S, "BD_imp_init");

  BDimp->guess = gsl_vector_alloc (BDimp->F->n);
  CHECK_MALLOC (BDimp->guess, "BD_imp_init");

  // working area
  BDimp->FTS = FTS_init (BDimp->BD->sys);
  CHECK_MALLOC (BDimp->FTS, "BD_imp_init");

  BDimp->pos = (double *)malloc (sizeof (double) * nm * 3);
  CHECK_MALLOC (BDimp->pos, "BD_imp_init");
  BDimp->q   = (double *)malloc (sizeof (double) * nm * 4);
  CHECK_MALLOC (BDimp->q,   "BD_imp_init");

  return (BDimp);
}

void
BD_imp_free (struct BD_imp *BDimp)
{
  if (BDimp == NULL) return;

  if (BDimp->x0 != NULL) free (BDimp->x0);
  if (BDimp->q0 != NULL) free (BDimp->q0);

  if (BDimp->z  != NULL) free (BDimp->z);

  if (BDimp->F  != NULL) free (BDimp->F);
  if (BDimp->S  != NULL) gsl_multiroot_fsolver_free (BDimp->S);
  if (BDimp->guess != NULL) gsl_vector_free (BDimp->guess);

  if (BDimp->FTS != NULL) FTS_free (BDimp->FTS);
  if (BDimp->pos != NULL) free (BDimp->pos);
  if (BDimp->q   != NULL) free (BDimp->q);

  free (BDimp);
}

/* set configuration x[] and q[],
 * and calculate the Brownian force vector BDimp->z[] for it.
 * note that even though q0==NULL, BDimp->q[] is set for FT and FTS versions.
 */
void
BD_imp_set_xq (struct BD_imp *BDimp,
	       const double *x, const double *q)
{
  int i;
  for (i = 0; i < BDimp->BD->sys->nm * 3; i ++)
    {
      BDimp->x0[i] = x[i];
    }

  if (BDimp->BD->sys->version > 0)
    {
      if (q == NULL)
	{
	  for (i = 0; i < BDimp->BD->sys->nm; i ++)
	    {
	      BDimp->q0[i*4+0] = 0.0;
	      BDimp->q0[i*4+1] = 0.0;
	      BDimp->q0[i*4+2] = 0.0;
	      BDimp->q0[i*4+3] = 1.0;
	    }
	}
      else
	{
	  for (i = 0; i < BDimp->BD->sys->nm * 4; i ++)
	    {
	      BDimp->q0[i] = q[i];
	    }
	}
    }

  /**
   * update the brownian force by the new config x[]
   */
  stokes_set_pos_mobile (BDimp->BD->sys, x);
  calc_brownian_force (BDimp->BD, BDimp->z);
}

/* set the time step,
 * and pre-factor BDimp->fact of the Brownian force to BDimp->z[].
 */
void
BD_imp_set_dt (struct BD_imp *BDimp,
	       double dt)
{
  BDimp->dt = dt;
  if (BDimp->BD->peclet > 0.0)
    {
      BDimp->fact = sqrt(2.0 / (BDimp->BD->peclet * dt));
    }
  else
    {
      BDimp->fact = 0.0;
    }
}


/*
 * INTPUT
 *  x[n] : = (x[nm*3])          for F version
 *         = (x[nm*3], q[nm*4]) for FT and FTS versions
 *  p    : (struct BD_imp *)
 * OUTPUT
 *  f[n] := -x + x0 + dt * (uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0)))
 *  f[n] := x - x0 - dt * (uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0)))
 */
int
BD_imp_JGdP00_func (const gsl_vector *x, void *p,
		    gsl_vector *f)
{
  struct BD_imp *BDimp = (struct BD_imp *)p;
  struct BD_params *BD = BDimp->BD;
  struct stokes *sys = BD->sys;

  int i;

  int nm = sys->nm;
  int nm3 = nm * 3;
  int nm4 = nm3 + nm;
  int nm5 = nm4 + nm;

  /**
   * set the forces for x0[] (and independent of the config)
   */
  // set force (and torque)
  if (sys->version == 0) // F version
    {
      for (i = 0; i < nm3; i ++)
	{
	  BDimp->FTS->f[i] = BD->F[i] + BDimp->fact * BDimp->z[i];
	}
    }
  else // FT and FTS version
    {
      for (i = 0; i < nm3; i ++)
	{
	  BDimp->FTS->f[i] = BD->F[i] + BDimp->fact * BDimp->z[i];
	}
      for (i = 0; i < nm3; i ++)
	{
	  BDimp->FTS->t[i] = BD->T[i] + BDimp->fact * BDimp->z[nm3 + i];
	}
      if (sys->version == 2)
	{
	  // FTS version
	  for (i = 0; i < nm5; i ++)
	    {
	      BDimp->FTS->e[i] = BD->E[i];
	    }
	}
    }

  /* extract pos[nm3] and q[nm4] from gsl_vector *x */
  for (i = 0; i < nm3; i ++)
    {
      BDimp->pos[i] = gsl_vector_get (x, i);
    }
  if (sys->version > 0)
    {
      for (i = 0; i < nm4; i ++)
	{
	  BDimp->q[i] = gsl_vector_get (x, nm3 + i);
	}
    }

  /**
   * set the force for x[]
   * F^P should be evaluated for the new config x[]
   */
  stokes_set_pos_mobile (sys, BDimp->pos); // first nm3 is the position

  if (BD->bonds->n > 0)
    {
      // calc bond (spring) force
      bonds_calc_force (BD->bonds, sys,
			BDimp->FTS->f, 1/* add */);
    }
  if (BD->ev != NULL)
    {
      // calc EV force
      EV_calc_force (BD->ev, sys,
		     BDimp->FTS->f, 1/* add */);
    }

  /**
   * solve (u, o) with F^P(x) for the configuration to x0[]
   */
  // reset the configuration to x0[]
  stokes_set_pos_mobile (sys, BDimp->x0);
  solve_mix_3all (sys,
		  BD->flag_lub, BD->flag_mat,
		  BDimp->FTS->f, BDimp->FTS->t, BDimp->FTS->e,
		  BD->uf, BD->of, BD->ef,
		  BDimp->FTS->u, BDimp->FTS->o, BDimp->FTS->s,
		  BDimp->FTS->ff, BDimp->FTS->tf, BDimp->FTS->sf);
  /* note that FTS->[u,o] are in the labo frame */

  /**
   * adjust the imposed flow with the new config x[].
   *
   * now FTS->u = U (in the labo frame)
   *            = u(x0) + R^{-1}.(F^ext(x0) + F^B(x0) + F^P(x))
   * so that FTS->u += (u(x) - u(x0))
   *                => u(x) + R^{-1}.(F^ext(x0) + F^B(x0) + F^P(x))
   *                   ^^^^
   * note that only u depends on the position (not o (nor e))
   * thereofre, there is no correction
   * (that is, the components in (u(x) - u(x0)) are zero)
   */
  for (i = 0; i < nm; i ++)
    {
      int ix = i*3;
      int iy = ix + 1;
      int iz = ix + 2;

      double dx = BDimp->pos[ix] - BDimp->x0[ix];
      double dy = BDimp->pos[iy] - BDimp->x0[iy];
      double dz = BDimp->pos[iz] - BDimp->x0[iz];

      // O\times (x - x0)
      double uOx = sys->Oi[1] * dz - sys->Oi[2] * dy;
      double uOy = sys->Oi[2] * dx - sys->Oi[0] * dz;
      double uOz = sys->Oi[0] * dy - sys->Oi[1] * dx;

      // E.(x - x0)
      double Ezz =
	- sys->Ei[0]
	- sys->Ei[4];
      double uEx
	= sys->Ei[0] * dx // xx . x
	+ sys->Ei[1] * dy // xy . y
	+ sys->Ei[2] * dz;// xz . z
      double uEy
	= sys->Ei[1] * dx // yx . x
	+ sys->Ei[4] * dy // yy . y
	+ sys->Ei[3] * dz;// yz . z
      double uEz
	= sys->Ei[2] * dx // zx . x
	+ sys->Ei[3] * dy // zy . y
	+ Ezz        * dz;// zz . z

      BDimp->FTS->u [ix] += uOx + uEx;
      BDimp->FTS->u [iy] += uOy + uEy;
      BDimp->FTS->u [iz] += uOz + uEz;
    }
  /* now FTS->u = u(x) + U(x0, F^ext(x0), F^B(x0), F^P(x)) */

  /**
   * form the vector f[] for the root-finding routines
   */
  for (i = 0; i < nm3; i ++)
    {
      double fi = 
	/*
	- BDimp->pos[i]
	+ BDimp->x0[i]
	+ BDimp->dt * BDimp->FTS->u[i];
	*/
	BDimp->pos[i] - BDimp->x0[i]
	- BDimp->dt * BDimp->FTS->u[i];
      gsl_vector_set (f, i, fi);
    }
  if (sys->version > 0)
    {
      // FT or FTS version
      for (i = 0; i < nm; i ++)
	{
	  int ix = i * 3;
	  int iq = i * 4;
	  double dQdt[4];
	  quaternion_dQdt (BDimp->q + iq, BDimp->FTS->o + ix, dQdt);
	  int ii;
	  for (ii = 0; ii < 4; ii ++)
	    {
	      double fi =
		/*
		- BDimp->q[iq + ii]
		+ BDimp->q0[iq + ii]
		+ BDimp->dt * dQdt[ii];
		*/
		BDimp->q[iq + ii] - BDimp->q0[iq + ii]
		- BDimp->dt * dQdt[ii];
	      gsl_vector_set (f, nm3 + iq + ii, fi);
	    }
	}
    }

  return (GSL_SUCCESS);
}

/* set gsl_vector
 * if q == NULL, q = (0,0,0,1) is set.
 */
void
BD_imp_set_guess (const struct BD_imp *BDimp,
		  const double *x,
		  const double *q,
		  gsl_vector *guess)
{
  int nm3 = BDimp->BD->sys->nm * 3;
  int i;
  for (i = 0; i < nm3; i ++)
    {
      gsl_vector_set (guess, i, x[i]);
    }
  if (BDimp->BD->sys->version > 0)
    {
      if (q == NULL)
	{
	  for (i = 0; i < BDimp->BD->sys->nm; i ++)
	    {
	      gsl_vector_set (guess, nm3 + i*4 + 0, 0.0);
	      gsl_vector_set (guess, nm3 + i*4 + 1, 0.0);
	      gsl_vector_set (guess, nm3 + i*4 + 2, 0.0);
	      gsl_vector_set (guess, nm3 + i*4 + 3, 1.0);
	    }
	}
      else
	{
	  for (i = 0; i < BDimp->BD->sys->nm * 4; i ++)
	    {
	      gsl_vector_set (guess, nm3 + i, q[i]);
	    }
	}
    }
}

void
BD_imp_get_root (const struct BD_imp *BDimp,
		 gsl_vector *root,
		 double *x,
		 double *q)
{
  int nm3 = BDimp->BD->sys->nm * 3;
  int i;
  for (i = 0; i < nm3; i ++)
    {
      x[i] = gsl_vector_get (root, i);
    }
  if (q != NULL)
    {
      for (i = 0; i < BDimp->BD->sys->nm * 4; i ++)
	{
	  q[i] = gsl_vector_get (root, nm3 + i);
	}
    }
}


/* evolve position of particles by semi-implicit scheme
 * ref: Jendrejack et al (2000) J. Chem. Phys. vol.113 p.2894.
 * INPUT
 *  BDimp   : struct BD_imp
 *  x[nm*3] : positions of particles   at t = t0
 *  q[nm*4] : quaternions of particles at t = t0 (only for FT and FTS)
 *            if NULL is given, just ignored.
 *  dt      : time step (scaled by a/U)
 * OUTPUT
 *  x[nm*3] : updated positions of particles at t = t0 + dt
 *  q[nm*4] : quaternions of particles       at t = t0 + dt
 *            (only if q[] is given for FT and FTS)
 *  returned value : the integrated time duration
 */
double
BD_evolve_JGdP00 (struct BD_imp *BDimp,
		  double *x, double *q,
		  double dt)
{
  /**
   * update BD_imp data
   */
  BD_imp_set_xq (BDimp, x, q);
  if (dt != BDimp->dt)
    {
      BD_imp_set_dt (BDimp, dt);
    }

  /**
   * set the initial guess
   */
  BD_imp_set_guess (BDimp, x, q, BDimp->guess);
  gsl_multiroot_fsolver_set (BDimp->S, BDimp->F, BDimp->guess);

  int status;
  int iter = 0;
  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (BDimp->S);

      if (status)   /* check if solver is stuck */
	break;

      status = gsl_multiroot_test_residual (BDimp->S->f, BDimp->eps);
    }
  while (status == GSL_CONTINUE && iter < BDimp->itmax);

  /**
   * retreive the solution
   */
  gsl_vector *root = gsl_multiroot_fsolver_root (BDimp->S);
  BD_imp_get_root (BDimp, root, x, q);

  /** check
  int i;
  for (i = 0; i < F->n; i ++)
    {
      fprintf (stdout, "#JGdP00 solution x[%d] = %e %e %e\n",
	       i,
	       x[i], 
	       gsl_vector_get (guess, i),
	       fabs (x[i] - gsl_vector_get (guess, i)));
    }
  **/

  //gsl_vector_free (root); // don't we need that?

  return (dt);
}


/* wrapper for BD_imp_evolve()
 * INPUT
 *  *t    : (input) current time
 *  t_out : output time
 *  *dt   : (input) current inner time-step
 *  y[n]  : (input) current configuration at (*t),
 *          where n = nm*3 for F version
 *                    (y[] in the first (nm*3) elements = velocity),
 *                n = nm*3 + nm*4 with quaternion
 *                    (y[] in the first (nm*3) elements = velocity,
 *                     y[] in the next (nm*4) elements = quaternion).
 * OUTPUT
 *  *t    : (output) output time (= t_out)
 *  *dt   : (output) current (updated, if necessary) inner time-step
 *  y[n]  : (output) updated configuration at t_out
 */
void
BD_imp_ode_evolve (struct BD_imp *BDimp,
		   double *t, double t_out, double *dt,
		   double *y)
{
  struct BD_params *BD = BDimp->BD;

  int nm3 = BD->sys->nm * 3;
  // asign local pointers x[] and q[] by y[].
  double *x = y;
  double *q = NULL;
  if (BD->flag_Q != 0)
    {
      q = y + nm3;
    }

  // local time step
  double dt_local = *dt;

  // loop until *t == t_out
  do
    {
      if ((*t) + dt_local > t_out)
	{
	  dt_local = t_out - (*t);
	}

      double dt_done = BD_evolve_JGdP00 (BDimp, x, q, dt_local);

      //(*t) += dt_local;
      (*t) += dt_done;
    }
  while ((*t) < t_out);
}
