/* header file for ode-quaternion.c --
 * ODE utility routines for angle by quaternion
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ode-quaternion.h,v 1.3 2008/07/17 02:18:53 kichiki Exp $
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
#ifndef	_ODE_QUATERNION_H_
#define	_ODE_QUATERNION_H_


void
quaternion_W (const double *Q, double *W);

void
quaternion_Winv_lapack (const double *Q, double *W);

void
quaternion_Winv (const double *Q, double *W);

void
quaternion_by_Euler (double phi, double theta, double psi,
		     double *Q);

void
quaternion_dQdt (const double *Q, const double *O,
		 double *dQdt);

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
 *           ode->bonds     : (struct BONDS *)
 * OUTPUT
 *  dydt[] := (d/dt) x(t) = U(t)
 *         where U(t) := R^-1 . F, the terminal velocity
 */
int
dydt_Q_hydro (double t, const double *y, double *dydt,
	      void *params);

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
 *           ode->bonds     : (struct BONDS *)
 * OUTPUT
 *  dydt[] := (d/dt) y(t), where y(t) = (x(t), U(t)).
 *         (d/dt) x(t) = U(t)
 *         (d/dt) U(t) = (1/stokes)(-U(t) + V(t)),
 *         where V(t) := R^-1 . F, the terminal velocity
 */
int
dydt_Q_hydro_st (double t, const double *y, double *dydt,
		 void *params);


#endif /* !_ODE_QUATERNION_H_ */
