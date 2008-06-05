/* implicit Brownian dynamics algorithms by NITSOL
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bd-imp-nitsol.c,v 1.1 2008/06/05 03:20:48 kichiki Exp $
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

#include <stokes.h> // stokes_set_pos()
#include <ode.h> // solve_mix_3all()
#include <ode-quaternion.h> // quaternion_dQdt()

#include <brownian.h> // struct BD_params
#include <bd-imp.h> // struct BD_imp

#include "nitsol.h"
#include "bd-imp-nitsol.h"

// BLAS functions
double
ddot_(int* N, 
      double* X, int* incX, 
      double* Y, int* incY);
double
dnrm2_(int* N, 
       double* X, int* incX);


/* jacobian and preconditioner
 * at this time, nothing will be done for both
 */
void BD_imp_NITSOL_jacv (int *n, double *xcur, double *fcur,
			 int *ijob, double *v, double *z,
			 double *rpar, int *ipar, int *itrmjv)
{
  if (*ijob == 0)
    {
      // calc Jacobian-vector product
      return;
    }
  else if (*ijob == 1)
    {
      // apply the preconditioner
      return;
    }
  else
    {
      fprintf (stderr, "# NITSOL_jacv: invalid ijob %d\n", *ijob);
      return;
    }
}

/* a wrapper of the nonlinear equations for the semi-implicit algorithms
 * (both JGdP and siPC) for NITSOL
 * INTPUT
 *  xcur[n] : = (x[nm*3])          for F version
 *            = (x[nm*3], q[nm*4]) for FT and FTS versions
 *  rpar    : (struct BD_imp *)BDimp
 *  ipar    : not used.
 * OUTPUT
 *  itrmf   : 0 means success
 *  fcur[n] := x - x0
 *           - dt * (uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for JGdP, or
 *  fcur[n] := x - x0
 *           - (dt/2) * (U^pr
 *                       + uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for siPC, where U^pr is the predictor obtained by Euler scheme as 
 *  U^pr = uinf(x0) + M(x0).(F^E + F^P(x0) + F^B(x0)).
 */
void
BD_imp_NITSOL_f (int *n, double *xcur, double *fcur,
		 double *rpar, int *ipar, int *itrmf)
{
  // let's say that (double *)rpar corresponds to (struct BD_imp *)BDimp
  struct BD_imp *BDimp = (struct BD_imp *)rpar;

  BD_imp_f (BDimp, xcur, fcur);

  *itrmf = 0;
}


/* set gsl_vector
 * if q == NULL, q = (0,0,0,1) is set.
 */
static void
BD_imp_NITSOL_set_guess (const struct BD_imp *BDimp,
			 const double *x,
			 const double *q,
			 double *x_nitsol)
{
  int nm3 = BDimp->BD->sys->nm * 3;
  int i;
  for (i = 0; i < nm3; i ++)
    {
      x_nitsol[i] = x[i];
    }
  if (BDimp->BD->sys->version > 0)
    {
      if (q == NULL)
	{
	  for (i = 0; i < BDimp->BD->sys->nm; i ++)
	    {
	      x_nitsol[nm3 + i*4 + 0] = 0.0;
	      x_nitsol[nm3 + i*4 + 1] = 0.0;
	      x_nitsol[nm3 + i*4 + 2] = 0.0;
	      x_nitsol[nm3 + i*4 + 3] = 1.0;
	    }
	}
      else
	{
	  for (i = 0; i < BDimp->BD->sys->nm * 4; i ++)
	    {
	      x_nitsol[nm3 + i] = q[i];
	    }
	}
    }
}

static void
BD_imp_NITSOL_get_root (const struct BD_imp *BDimp,
			const double *x_nitsol,
			double *x,
			double *q)
{
  int nm3 = BDimp->BD->sys->nm * 3;
  int i;
  for (i = 0; i < nm3; i ++)
    {
      x[i] = x_nitsol[i];
    }
  if (q != NULL)
    {
      for (i = 0; i < BDimp->BD->sys->nm * 4; i ++)
	{
	  q[i] = x_nitsol[nm3 + i];
	}
    }
}

void
BD_imp_NITSOL_wrap (struct BD_imp *BDimp,
		    double *x, double *q)
{
  // maximum Krylov subspace dimension
  int maxkd = 50;

  int input[10];
  // Initialize all inputs to zero (=> default options). 
  int i;
  for (i = 0; i < 10; i ++)
    {
      input[i] = 0;
    }

  // Jacobian : input[1]
  input[1] = 0; // approximation for J*v
  //input[1] = 1; // analytical J*v

  // order of finite-difference formula for the Jacobian
  //input[7] = 1; // 1st order
  //input[7] = 2; // 2nd order
  input[7] = 4; // 4th order

  // precondition : intput[4]
  input[4] = 0; // no preconditioning?
  //input[4] = 1; // right preconditioning

  // method : input[2]
  input[2] = 0; // GMRES
  //input[2] = 1; // BiCGSTAB
  //input[2] = 2; // TFQMR

  // max krylov subspace dimension : input[3]
  input[3] = maxkd;

  // forcing term : input[9]
  input[9] = 0;
  //input[9] = 1;
  //input[9] = 2; // require choice2_exp (alpha), choice2_coef (gamma)
  //input[9] = 3; // require etafixed

  nitprint_.iplvl  = 1; // level 1
  nitprint_.ipunit = 6; // stdout


  int n;
  if (BDimp->BD->sys->version == 0)
    {
      n = BDimp->BD->sys->nm * 3;
    }
  else
    {
      n = BDimp->BD->sys->nm * 7;
    }

  int lrwork = n * (maxkd + 10) + maxkd * (maxkd + 3);
  double *rwork = (double *)malloc (sizeof (double) * lrwork);
  CHECK_MALLOC (rwork, "BD_imp_NITSOL_wrap");


  double *x_nitsol = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x_nitsol, "BD_imp_NITSOL_wrap");

  BD_imp_NITSOL_set_guess (BDimp, x, q, x_nitsol);


  double ftol   = 1.0e-6;
  double stptol = 1.0e-6;
  int iterm;
  int info[6];
  nitsol_ (&n,
	   x_nitsol,
	   BD_imp_NITSOL_f,
	   BD_imp_NITSOL_jacv,
	   &ftol,
	   &stptol,
	   input,
	   info,
	   rwork,
	   (double *)BDimp, // rpar
	   NULL, // ipar
	   &iterm,
	   ddot_, dnrm2_);

  // check
  fprintf (stdout, "# No. function evaluations:     %d\n", info[0]);
  fprintf (stdout, "# No. J*v evaluations:          %d\n", info[1]);
  fprintf (stdout, "# No. P(inverse)*v evaluations: %d\n", info[2]);
  fprintf (stdout, "# No. linear iterations:        %d\n", info[3]);
  fprintf (stdout, "# No. nonlinear iterations:     %d\n", info[4]);
  fprintf (stdout, "# No. backtracks:               %d\n", info[5]);

  BD_imp_NITSOL_get_root (BDimp, x_nitsol, x, q);

  free (x_nitsol);
  free (rwork);
}


/**
 * semi-implicit algorithm by
 * Jendrejack et al (2000) J.Chem.Phys. vol 113 p.2894.
 */
/* evolve position of particles by semi-implicit scheme
 * with NITSOL
 * ref: Jendrejack et al (2000) J. Chem. Phys. vol.113 p.2894.
 * INPUT
 *  t       : current time
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
BD_evolve_NITSOL_JGdP00 (double t,
			 struct BD_imp *BDimp,
			 double *x, double *q,
			 double dt)
{
  double dt_local = dt;
  struct overlap ol;

  // set sys->shear_shift if necessary
  // JGdP is 1-step scheme so that we only need to set "shift" once
  stokes_set_shear_shift (BDimp->BD->sys, t,
			  BDimp->BD->t0, BDimp->BD->s0);

  /**
   * update BD_imp data
   */
 BD_evolve_NITSOL_JGdP_REDO:
  BD_imp_set_xq (BDimp, x, q);
 BD_evolve_NITSOL_JGdP_REDO_scale:
  if (dt_local != BDimp->dt)
    {
      BD_imp_set_dt (BDimp, dt_local);
    }

  /**
   * solve the nonlinear equations by NITSOL
   */
  BD_imp_NITSOL_wrap (BDimp, x, q);


  /* dt-ajustment process -- overlap check */
  if (BDimp->BD->sys->rmin == 0.0 // if rmin is defined, skip overlap check
      && check_overlap (BDimp->BD->sys, x, BDimp->BD->rmin, &ol) > 0)
    {
      //dt_local = reset_dt_by_ol (BDimp->BD->sys, dt_local, x, &ol);
      dt_local *= 0.5;
      if (dt_local < BDimp->BD->dt_lim)
	{
	  // just reject the random force z[]
	  /* BD->pos is still set by the initial config x[] */
	  fprintf (stderr, "# NITSOL_JGdP: overlap. "
		   "dt = %e > %e => reconstruct FB\n",
		   dt_local, BDimp->BD->dt_lim);
	  // reset dt to the original step
	  dt_local = dt;
	  goto BD_evolve_NITSOL_JGdP_REDO;
	}
      else
	{
	  //fact = BD_params_get_fact (BD, dt_local);
	  /* adjustment of BDimp->fact is done by BD_imp_set_dt() above
	   * just after the BD_evolve_NITSOL_JGdP_REDO_scale
	   */
	  fprintf (stderr, "# NITSOL_JGdP: overlap. "
		   "dt = %e > %e => reconstruct FB\n",
		   dt_local, BDimp->BD->dt_lim);
	  goto BD_evolve_NITSOL_JGdP_REDO_scale;
	}
    }

  return (dt_local);
}


/**
 * semi-implicit predictor-corrector algorithm
 */
/* evolve position of particles by semi-implicit predictor-corrector
 * with NITSOL
 * INPUT
 *  t       : current time
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
BD_evolve_NITSOL_imp_PC (double t,
			 struct BD_imp *BDimp,
			 double *x, double *q,
			 double dt)
{
  double dt_local = dt;
  struct overlap ol;

  /**
   * update BD_imp data
   */
 BD_evolve_NITSOL_imp_PC_REDO:
  // set sys->shear_shift if necessary
  // at the predictor step
  stokes_set_shear_shift (BDimp->BD->sys, t,
			  BDimp->BD->t0, BDimp->BD->s0);
  BD_imp_set_xq (BDimp, x, q);
 BD_evolve_NITSOL_imp_PC_REDO_scale:
  if (dt_local != BDimp->dt)
    {
      BD_imp_set_dt (BDimp, dt_local);
    }
  // set the predictor in BDimp->uP[] and BDimp->qP[]
  BD_imp_set_P (BDimp);

  /* dt-ajustment process -- overlap check */
  if (BDimp->BD->sys->rmin == 0.0 // if rmin is defined, skip overlap check
      && check_overlap (BDimp->BD->sys, BDimp->xP, BDimp->BD->rmin, &ol) > 0)
    {
      //dt_local = reset_dt_by_ol (BDimp->BD->sys, dt_local, x, &ol);
      dt_local *= 0.5;
      if (dt_local < BDimp->BD->dt_lim)
	{
	  // just reject the random force z[]
	  /* BD->pos is still set by the initial config x[] */
	  fprintf (stderr, "# NITSOL_imp_PC: overlap in the predictor step. "
		   "dt = %e > %e => reconstruct FB\n",
		   dt_local, BDimp->BD->dt_lim);
	  // reset dt to the original step
	  dt_local = dt;
	  goto BD_evolve_NITSOL_imp_PC_REDO;
	}
      else
	{
	  //fact = BD_params_get_fact (BD, dt_local);
	  /* adjustment of BDimp->fact is done by BD_imp_set_dt() above
	   * just after the BD_evolve_NITSOL_imp_PC_REDO_scale
	   */
	  fprintf (stderr, "# NITSOL_imp_PC: overlap in the predictor step. "
		   "dt = %e > %e => reconstruct FB\n",
		   dt_local, BDimp->BD->dt_lim);
	  goto BD_evolve_NITSOL_imp_PC_REDO_scale;
	}
    }


  // set sys->shear_shift if necessary
  // at the corrector step
  stokes_set_shear_shift (BDimp->BD->sys, t + dt_local,
			  BDimp->BD->t0, BDimp->BD->s0);
  /**
   * solve the nonlinear equations by NITSOL
   */
  BD_imp_NITSOL_wrap (BDimp, x, q);


  /* dt-ajustment process -- overlap check */
  if (BDimp->BD->sys->rmin == 0.0 // if rmin is defined, skip overlap check
      && check_overlap (BDimp->BD->sys, x, BDimp->BD->rmin, &ol) > 0)
    {
      //dt_local = reset_dt_by_ol (BDimp->BD->sys, dt_local, x, &ol);
      dt_local *= 0.5;
      if (dt_local < BDimp->BD->dt_lim)
	{
	  // just reject the random force z[]
	  /* BD->pos is still set by the initial config x[] */
	  fprintf (stderr, "# NITSOL_imp_PC: overlap in the corrector step. "
		   "dt = %e > %e => reconstruct FB\n",
		   dt_local, BDimp->BD->dt_lim);
	  // reset dt to the original step
	  dt_local = dt;
	  goto BD_evolve_NITSOL_imp_PC_REDO;
	}
      else
	{
	  //fact = BD_params_get_fact (BD, dt_local);
	  /* adjustment of BDimp->fact is done by BD_imp_set_dt() above
	   * just after the BD_evolve_NITSOL_imp_PC_REDO_scale
	   */
	  fprintf (stderr, "# NITSOL_imp_PC: overlap in the corrector step. "
		   "dt = %e > %e => reconstruct FB\n",
		   dt_local, BDimp->BD->dt_lim);
	  goto BD_evolve_NITSOL_imp_PC_REDO_scale;
	}
    }

  return (dt_local);
}
