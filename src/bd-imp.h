/* header file for bd-imp.c --
 * implicit Brownian dynamics algorithms
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bd-imp.h,v 1.3 2007/12/22 04:30:55 kichiki Exp $
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
#include <brownian.h> // struct BD_params

struct BD_imp {
  struct BD_params *BD;
  double dt;
  double *x0;
  double *q0;
  double fact;
  double *z;

  // GSL stuff
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_function *F;
  gsl_multiroot_fsolver  *S;
  gsl_vector *guess;
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
 *  itmax : max of iteration for the root-finding
 *  eps   : tolerance for the root-finding
 */
struct BD_imp *
BD_imp_init (struct BD_params *BD,
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

/**
 * semi-implicit algorithm by
 * Jendrejack et al (2000) J.Chem.Phys. vol 113 p.2894.
 */
/* form the nonlinear equations for the algorithm by Jendrejack et al (2000)
 * INTPUT
 *  x[n] : = (x[nm*3])          for F version
 *         = (x[nm*3], q[nm*4]) for FT and FTS versions
 *  p    : (struct BD_imp *)
 * OUTPUT
 *  f[n] := x - x0 - dt * (uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0)))
 */
int
BD_imp_JGdP00_func (const gsl_vector *x, void *p,
		    gsl_vector *f);

/* set gsl_vector
 * if q == NULL, q = (0,0,0,1) is set.
 */
void
BD_imp_set_guess (const struct BD_imp *BDimp,
		  const double *x,
		  const double *q,
		  gsl_vector *guess);
void
BD_imp_get_root (const struct BD_imp *BDimp,
		 gsl_vector *root,
		 double *x,
		 double *q);

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
/* set the predictor for the configuration (BDimp->x0, BDimp->q0).
 * this must be called after setting both BD_imp_set_xq() and BD_imp_set_dt().
 */
void
BD_imp_set_P (struct BD_imp *BDimp);

/* form the nonlinear equations for semi-implicit predictor-corrector
 * INTPUT
 *  x[n] : = (x[nm*3])          for F version
 *         = (x[nm*3], q[nm*4]) for FT and FTS versions
 *  p    : (struct BD_imp *)
 * OUTPUT
 *  f[n] := x - x0
 *        - (dt/2) * (U^pr
 *                    + uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  where U^pr is the predictor obtained by Euler scheme as 
 *   U^pr = uinf(x0) + M(x0).(F^E + F^P(x0) + F^B(x0)).
 */
int
BD_imp_PC_func (const gsl_vector *x, void *p,
		gsl_vector *f);

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
