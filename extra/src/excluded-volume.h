/* header file for excluded-volume.c --
 * excluded-volume interactions
 * Copyright (C) 2007-2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_EXCLUDED_VOLUME_H_
#define	_EXCLUDED_VOLUME_H_


#include <libstokes-core.h> // struct stokes

struct EV {
  double r2; // square of the max distance for F^{EV}

  /* table for each particle */
  int n;     // number of particles
  double *l; /* := sqrt((2/3) hat(ls)^2) */
  double *A; /* := (1/pi)^{3/2} hat(v) N_Ks^2 / (peclet ev->l[i]^5)
	      * where
	      *   ls^2     = N_Ks * b_K^2 / 3 (dimensional, so as b_K)
	      *   hat(ls)  = ls / length      (dimensionless)
	      *   hat(v)   = v / length^3     (dimensionless)
	      * with which the force F^EV is given by 
	      *   F^{EV}_{i} = A * r_{ij} * exp (- r_{ij}^2 / l^2)
	      * (note: r_{ij} is also dimensionless scaled by length)
	      */
};


/* initialize struct EV
 * INPUT
 *  np     : number of particles
 *  length : unit of length given by "length" in SCM (dimensional number)
 *  peclet : peclet number (with respect to "length")
 *  r2     : square of the max distance for F^{EV}
 *  v[n]   : EV parameters for each spring in the dimension of "length" cubed
 *  NKs[n] : Kuhn steps for a spring belongs to the particle
 *  bK[n]  : (dimensional) Kuhn length in the dimension of "length"
 * OUTPUT
 *  returned value : struct EV, where l and A are defined by 
 *      ev->l[i] = sqrt((2/3) hat(ls)^2)
 *      ev->A[i] = (1/pi)^{3/2} hat(v) N_Ks^2 / (peclet ev->l[i]^5)
 *    where
 *      ls^2     = N_Ks * b_K^2 / 3 (dimensional, so as b_K)
 *      hat(ls)  = ls / length      (dimensionless)
 *      hat(v)   = v / length^3     (dimensionless)
 *    with which the force F^EV is given by 
 *      F^{EV}_{i} = A * r_{ij} * exp (- r_{ij}^2 / l^2)
 *    (note: r_{ij} is also dimensionless scaled by length)
 */
struct EV *
EV_init (int np,
	 double length, double peclet,
	 double r2,
	 const double *v,
	 const double *NKs,
	 const double *bK);

void
EV_free (struct EV *ev);


/* retrieve the EV parameters for particles i and j
 * with the combination rules 
 *   l(12) = (l1 + l2) / 2
 *   A(12) = (A1 * A2)^{1/2}
 * INPUT
 *  ev   : struct EV
 *  i, j : particle index (0 ~ np-1)
 * OUTPUT
 *  *l : 
 *  *A : 
 */
void
EV_get_coefficients (struct EV *ev,
		     int i, int j,
		     double *l, double *A);


/*
 * INPUT
 *  ev         : struct EV
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_calc_force (struct EV *ev,
	       struct stokes *sys,
	       double *f,
	       int flag_add);


#endif /* !_EXCLUDED_VOLUME_H_ */
