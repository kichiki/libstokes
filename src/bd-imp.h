/* header file for bd-imp.c --
 * implicit Brownian dynamics algorithms
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bd-imp.h,v 1.1 2007/11/07 04:41:29 kichiki Exp $
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
};


struct BD_imp *
BD_imp_init (struct BD_params *BD,
	     double dt,
	     const double *x0,
	     const double *q0);

void
BD_imp_free (struct BD_imp *BDimp);

void
BD_imp_set_xq (struct BD_imp *BDimp,
	       const double *x, const double *q);

void
BD_imp_set_dt (struct BD_imp *BDimp,
	       double dt);


/*
 * INTPUT
 *  x[n] : = (x[nm*3])          for F version
 *         = (x[nm*3], q[nm*4]) for FT and FTS versions
 *  p    : (struct BD_imp *)
 * OUTPUT
 *  f[n] := -x + x0 + dt * (uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0)))
 */
int
BD_imp_JGdP00_func (const gsl_vector *x, void *p,
		    gsl_vector *f);

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
 *  BD      : struct BD_params (sys, rng, flag_lub, flag_mat,
 *                              flag_Q, F, T, E, peclet are used.)
 *  x[nm*3] : positions of particles   at t = t0
 *  q[nm*4] : quaternions of particles at t = t0 (only for FT and FTS)
 *            if NULL is given, just ignored.
 *  dt      : time step (scaled by a/U)
 * OUTPUT
 *  x[nm*3] : updated positions of particles at t = t0 + dt
 *  q[nm*4] : quaternions of particles       at t = t0 + dt
 *            (only if q[] is given for FT and FTS)
 */
void
BD_evolve_JGdP00 (struct BD_params *BD,
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
BD_imp_ode_evolve (struct BD_params *BD,
		   double *t, double t_out, double *dt,
		   double *y);


#endif /* !_BD_IMP_H_ */
