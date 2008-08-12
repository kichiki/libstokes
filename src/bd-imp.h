/* header file for bd-imp.c --
 * implicit Brownian dynamics algorithms
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bd-imp.h,v 1.6 2008/08/12 05:31:10 kichiki Exp $
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
#ifndef	_BD_IMP_H_
#define	_BD_IMP_H_


#include <gsl/gsl_multiroots.h>
#include <nitsol_c.h> // struct NITSOL
#include <brownian.h> // struct BD_params

struct BD_imp {
  struct BD_params *BD;
  double dt;
  double *x0;
  double *q0;
  double fact;
  double *z;

  int solver; /* 0 == GSL
	       * 1 == NITSOL
	       * 2 == fastSI
	       */

  // GSL stuff
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_function *F;
  gsl_multiroot_fsolver  *S;
  gsl_vector *guess;

  // NITSOL stuff
  struct NITSOL *nit;

  int itmax;
  double eps;


  // working area used in BD_imp_JGdP00_func()
  struct FTS *FTS;
  double *pos;
  double *q;

  // for predictor-corrector algorithm
  int flag_PC;
  double *xP;
  double *qP;
  double *uP;
  double *dQdtP;
};


/* initialize struct BD_imp
 * INPUT
 *  BD : struct BD_params
 *       note that BDimp->BD is just a pointer to BD in the argument.
 *  solver : 0 == GSL solver
 *           1 == NITSOL
 *           2 == fastSI
 *  itmax : max of iteration for the root-finding
 *  eps   : tolerance for the root-finding
 */
struct BD_imp *
BD_imp_init (struct BD_params *BD,
	     int solver,
	     int itmax, double eps);

void
BD_imp_free (struct BD_imp *BDimp);

/* set configuration x[] and q[],
 * and calculate the Brownian force vector BDimp->z[] for it.
 * note that even though q0==NULL, BDimp->q[] is set for FT and FTS versions.
 */
void
BD_imp_set_xq (struct BD_imp *BDimp,
	       const double *x, const double *q);

/* set the time step,
 * and pre-factor BDimp->fact of the Brownian force to BDimp->z[].
 */
void
BD_imp_set_dt (struct BD_imp *BDimp,
	       double dt);

/* calculate forces on each particle including
 *   constant force in BDimp->BD->F[],T[]
 *   Brownian force BDimp->fact * BDimp->z[],
 *   bond force F^bond(pos[]) by ...
 *   EV force F^EV(pos[]) by ...
 * INPUT
 *  pos[] : position for F^bond and F^EV
 *  BDimp->BD->FTS->F[],T[],E[] for constant force
 *  BDimp->fact, BDimp->z[] for F^B
 * OUTPUT
 *  BDimp->FTS->f[],t[],e[]
 */
void
BD_imp_calc_forces (struct BD_imp *BDimp,
		    const double *pos);

/* adjust the imposed flow with the new config x[].
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
void
BD_imp_adj_uinf (struct BD_imp *BDimp,
		 double *x0, double *x);


/* form the nonlinear equations for the semi-implicit algorithms
 * (both JGdP and siPC)
 * INTPUT
 *  x[n] : = (x[nm*3])          for F version
 *         = (x[nm*3], q[nm*4]) for FT and FTS versions
 *  p    : (struct BD_imp *)
 * OUTPUT
 *  f[n] := x - x0
 *        - dt * (uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for JGdP, or
 *  f[n] := x - x0
 *        - (dt/2) * (U^pr
 *                    + uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for siPC, where U^pr is the predictor obtained by Euler scheme as 
 *  U^pr = uinf(x0) + M(x0).(F^E + F^P(x0) + F^B(x0)).
 */
void
BD_imp_f (struct BD_imp *BDimp, const double *x,
	  double *f);


/* wrapper of BD_imp_f() for GSL-MULTROOT routine
 * this is applicable for both JGdP and siPC
 */
int
BD_imp_GSL_MULTIROOT_func (const gsl_vector *x, void *p,
			   gsl_vector *f);


/**
 * semi-implicit algorithm by
 * Jendrejack et al (2000) J.Chem.Phys. vol 113 p.2894.
 */
/* evolve position of particles by semi-implicit scheme
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
BD_evolve_JGdP00 (double t,
		  struct BD_imp *BDimp,
		  double *x, double *q,
		  double dt);


/**
 * semi-implicit predictor-corrector algorithm
 */
/* evolve position of particles by semi-implicit predictor-corrector
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
BD_evolve_imp_PC (double t,
		  struct BD_imp *BDimp,
		  double *x, double *q,
		  double dt);


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
		   double *y);


#endif /* !_BD_IMP_H_ */
