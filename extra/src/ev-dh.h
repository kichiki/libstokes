/* header file for ev-dh.c --
 * excluded-volume interactions
 * Copyright (C) 2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_EV_DH_H_
#define	_EV_DH_H_

#include <libstokes-core.h> // struct stokes

#include "grid.h"   // struct RYUON_grid


struct EV_DH {
  /* system parameters */
  double r_cutoff; /* the max distance for F^{EV} in length */
  double r2;       /* square of r_cutoff */
  double a_sys;    /* := (1/pe)(1/kT)(e^2/4pi e0)(1/a) dimensionless number */
  double rd;       /* debye length scaled by the characteristic length */

  /* parameters for each chain */
  /* currently this is implemented particle-wise for simplicity */
  int n;
  double *nu; /* nu := (nu) l0 / e, the dimensionless charge density, where
	       *       (nu) is the line density of charge [e/nm]
	       *       l0   is the bond length [nm]
	       *       e    is the elementary charge [C] */

  /* RYUON_grid */
  int flag_grid; /* == 0, plain particle-particle loop, O(N^2)
		  * == 1, near-particle loop, O(N)
		  */
  struct RYUON_grid *grid;
};


/* initialize struct EV
 * INPUT
 *  length : characteristic length
 *           (in the same dimension for rd below, usually nm)
 *  peclet : peclet number
 *  rd     : Debye length (in the same dimension for a above, usually nm)
 *  T      : temperature in Kelvin.
 *  e      : dielectric constant of the solution (dimensionless number)
 *  eps    : to determine cut-off distance through eps = exp (-r_cutoff / rd) 
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV_DH,
 *      where only the system parameters (a_sys, rd) are defined.
 *      nu[i] is zero cleared.
 */
struct EV_DH *
EV_DH_init (double length, double peclet,
	    double rd, double T, double e,
	    double eps,
	    int np);

void
EV_DH_free (struct EV_DH *ev_dh);


void
EV_DH_set_force_ij (struct stokes *sys,
		    struct EV_DH *ev_dh,
		    int i, int j, double r2, double x, double y, double z,
		    double *f);


/*
 * INPUT
 *  ev_dh      : struct EV_DH
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_DH_calc_force (struct EV_DH *ev_dh,
		  struct stokes *sys,
		  double *f,
		  int flag_add);


#endif /* !_EV_DH_H_ */
