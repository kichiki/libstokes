/* ODE utility routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ode.h,v 1.4 2007/12/05 03:47:49 kichiki Exp $
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
#ifndef	_ODE_H_
#define	_ODE_H_

#include "stokes.h" // struct stokes

/* wrapper for solving mob and mix problems in ODE routines
 */
void
solve_mix_3all (struct stokes * sys,
		int flag_lub,
		int flag_mat,
		const double *fm, const double *tm, const double *em,
		const double *uf, const double *of, const double *ef,
		double *um, double *om, double *sm,
		double *ff, double *tf, double *sf);


struct ode_params
{
  struct stokes *sys;
  double *F;
  double *T;
  double *E;
  double *uf;
  double *of;
  double *ef;
  int flag_mat;
  int flag_lub;
  double st;
  struct bonds *bonds;
  double gamma;
};


/* set the parameters to struct ode_params
 * INPUT
 *  (struct stokes *)sys -- initialize before calling!
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
		 double gamma);

void 
ode_params_free (struct ode_params *params);


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
		 void *params);

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
	       void *params);

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
	    void *params);


/** test routines **/

void
ODE_st_euler (double t, int n, double *y,
	      double t_out,
	      int (*f_dydt)(double, const double *, double *, void *),
	      void *params);

void
ODE_plain_euler (double t, int n, double *y,
		 double t_out,
		 void *params);


#endif /* !_ODE_H_ */
