/* header file for bd-imp-fast.c --
 * fast semi-implicit Brownian dynamics algorithms
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bd-imp-fast.h,v 1.1 2008/07/27 00:52:14 kichiki Exp $
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
#ifndef	_BD_IMP_FAST_H_
#define	_BD_IMP_FAST_H_

#include <gsl/gsl_multiroots.h>


/* calc imposed flow difference for a spring at the mid-point
 * basically a copy from shift_rest_to_labo_U() in f.c.
 * INPUT
 *  b : struct BD_imp
 *  q[3] : connector vector at which the imposed flow difference is evaluated
 * OUTPUT
 *  u[3] : overwritten
 */
void
fastSI_calc_imposed (struct stokes *sys,
		     const double *q,
		     double *u);

/* calc imposed flow difference for a spring at the mid-point
 * basically a copy from shift_rest_to_labo_U() in f.c.
 * INPUT
 *  b : struct BD_imp
 *  q1[3] q2[3] : two connector vectors
 *                the imposed flow difference is evaluated at the mid-point
 *                (1/2)(q1 + q2).
 * OUTPUT
 *  u[3] : overwritten
 */
void
fastSI_calc_imposed_midpoint (struct stokes *sys,
			      const double *q1,
			      const double *q2,
			      double *u);

double
calc_B (struct BD_imp *b, int ibond);

/* calculate P = V_a - V_b, where V_a = sum_{j} M_{a,j} F_{j}.
 * INPUT
 *  sys : struct stokes
 *  flag_noHI : 1 == no HI
 *              0 == with HI
 *  bonds : struct BONDS
 *  f[np] : forces on the particles
 *  ib : bond index considering
 * OUTPUT
 *  P[3]  : zero cleared and set.
 */
void
calc_Pstar (struct stokes *sys, int flag_noHI,
	    struct BONDS *bonds,
	    const double *f,
	    int ib,
	    double *P);

/* calc right-hand side vector for the spring "i" with HI case
 * for general configuration
 * INPUT
 *  b        : struct BD_imp
 *             **NOTE** b->z[] should be in the mobility mode
 *                      and stored in the connector index.
 *  q0[ng*3] : initial connector vectors for imposed velocity
 *  q1[ng*3] : updated connector vectors for imposed velocity
 *  q2[ng*3] : the last-step connector vectors
 *  q [ng*3] : the connector vectors updated up to (i-1)-th spring
 *             where ng is the number of particles belongs to "ig"
 *  ib       : the bond index for native bond (not COM), that is,
 *             ib = 0, 1, ..., ng-2. (COM is at ib=(ng-1).)
 *  iq       : index for connector divided by 3
 * OUTPUT
 *  r[3]     : the right-hand side vector for the cubic equation
 */
void
fastSI_rhs (struct BD_imp *b,
	    const double *q0,
	    const double *q1,
	    const double *q2,
	    const double *q,
	    int ib,
	    int iq,
	    double *r);

/* calculate P^COM = (1/N)sum_i V_i, where V_i = sum_{j} M_{i,j} F_{j}.
 * INPUT
 *  sys : struct stokes
 *  flag_noHI : 1 == no HI
 *              0 == with HI
 *  group : struct BONDS_GROUP for the considering group
 *  f[np] : forces on the particles
 * OUTPUT
 *  P[3]  : zero cleared and set.
 */
void
calc_Pstar_COM (struct stokes *sys, int flag_noHI,
		struct BONDS_GROUP *group,
		const double *f,
		double *P);
/* calc right-hand side vector for the spring "i" with HI case
 * for general configuration
 * INPUT
 *  b        : struct BD_imp
 *             **NOTE** b->z[] should be in the mobility mode
 *                      and stored in the connector index.
 *  q0[ng*3] : initial connector vectors for imposed velocity
 *  q1[ng*3] : updated connector vectors for imposed velocity
 *  q2[ng*3] : the last-step connector vectors
 *  q [ng*3] : the connector vectors updated up to (i-1)-th spring
 *             where ng is the number of particles belongs to "ig"
 *  ig       : the group index
 *  iq       : index for connector divided by 3
 *             this should be the COM component
 * OUTPUT
 *  qCOM [3] : the updated COM
 */
void
fastSI_update_COM (struct BD_imp *b,
		   const double *q0,
		   const double *q1,
		   const double *q2,
		   const double *q,
		   int ig,
		   int iq,
		   double *qCOM);

/* calculate V_i = sum_{j} M_{i,j} F_{j}.
 * INPUT
 *  sys : struct stokes
 *  flag_noHI : 1 == no HI
 *              0 == with HI
 *  ip       : particle index of the isolated particle
 *  f[np] : forces on the particles
 * OUTPUT
 *  P[3]  : zero cleared and set.
 */
void
calc_Pstar_isolated (struct stokes *sys, int flag_noHI,
		     int ip,
		     const double *f,
		     double *P);
/* calc right-hand side vector for the spring "i" with HI case
 * for general configuration
 * INPUT
 *  b        : struct BD_imp
 *             **NOTE** b->z[] should be in the mobility mode
 *                      and stored in the connector index.
 *  q0[ng*3] : initial connector vectors for imposed velocity
 *  q1[ng*3] : updated connector vectors for imposed velocity
 *  q2[ng*3] : the last-step connector vectors
 *  q [ng*3] : the connector vectors updated up to (i-1)-th spring
 *             where ng is the number of particles belongs to "ig"
 *  ip       : particle index of the isolated particle
 *  iq       : index for connector divided by 3
 *             this should be the component of isolated particle
 * OUTPUT
 *  qp [3]   : the updated connector of the isolated particle
 */
void
fastSI_update_isolated (struct BD_imp *b,
			const double *q0,
			const double *q1,
			const double *q2,
			const double *q,
			int ip,
			int iq,
			double *qp);

double
fastSI_solve_cubic_in_range (double a, double b, double c, double d,
			     double xmin, double xmax);

/* solve the cubic equation
 * INPUT
 *  b  : struct BD_imp
 *  ib : bond index for (struct BONDS) b->BD->bonds
 *  R  : magnitude of the right-hand side vector
 */
double
fastSI_solve_cubic (struct BD_imp *b,
		    int ib,
		    double R);

/* 
 * INPUT
 *  b    : struct BD_imp
 *  q0[] : initial connector
 *  qC[] : the predictor for the initial corrector step (flag_first == 1)
 *         the previous corrector for update steps (flag_first == 0)
 *  flag_first : 1 for the initial corrector
 *               0 for the update steps
 * OUTPUT
 *  returned value : residual (only for flag_first == 0)
 *  q[]  : updated corrector
 */
double
fastSI_solve_corrector (struct BD_imp *b,
			const double *q0,
			const double *qC,
			int flag_first,
			double *q);

/* update connector q[] by Euler scheme
 * this is used to build predictor in fastSI algorithm
 * INPUT
 *  b : struct BD_imp
 *      **NOTE** b->z[] should be in the mobility mode
 *               and stored in the connector index.
 */
void
fastSI_solve_Euler (struct BD_imp *b,
		    const double *q0,
		    double *q);

void
fastSI_solve (struct BD_imp *b,
	      const double *q0,
	      double *q);


/**
 * for nonlinear solvers
 */
/* set the initial connector vector q0[np*3] in b->x0[]
 */
void
fastSI_set_Q0 (struct BD_imp *b,
	       const double *q0);

/* calc nonlinear equation f[] = 0, where
 * f[i] = q0[i]
 *       + dt * [uimp((q0[i]+q[i])/2)
 *               + (1/a[i+1]) * F(q[i+1])
 *               - (1/a[i+1] + 1/a[i]) * F(q[i])
 *               + (1/a[i]) * F(q[i-1])]
 *       + Y[i]
 *       - q[i]
 * INTPUT
 *  b : struct BD_imp
 *      NOTE that b->x0[np*3] should be set by q0[] before calling
 */
void
fastSI_f (struct BD_imp *b,
	  const double *q,
	  double *f);


/**
 * NITSOL stuff
 */
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
fastSI_NITSOL_f (int *n, double *xcur, double *fcur,
		 double *rpar, int *ipar, int *itrmf);
void
fastSI_NITSOL_wrap (struct BD_imp *b,
		    const double *q0,
		    double *q);


/**
 * GSL related stuff
 */
/* wrapper of fastSI_f() for GSL-MULTROOT routine
 */
int
fastSI_GSL_MULTIROOT_func (const gsl_vector *x, void *p,
			   gsl_vector *f);

void
fastSI_GSL_MULTIROOT_wrap (struct BD_imp *b,
			   const double *q0,
			   double *q);


void
fastSI_calc_FB (struct BD_imp *b,
		double *X);

/* evolve position of particles by semi-implicit predictor-corrector
 * INPUT
 *  t       : current time
 *  b       : struct BD_imp
 *            b->flag_solver == 0 : GSL-multiroot
 *                           == 1 : NITSOL
 *                           == 2 : fastSI
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
fastSI_evolve (double t,
	       struct BD_imp *b,
	       double *x,
	       double dt);


#endif /* !_BD_IMP_FAST_H_ */
