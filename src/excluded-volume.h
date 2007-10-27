/* header file for excluded-volume.c --
 * excluded-volume interactions
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: excluded-volume.h,v 1.2 2007/10/27 03:20:15 kichiki Exp $
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


struct EV {
  double r2; // square of the max distance for F^{EV}

  /* table for chain type */
  int n;     // number of chain types
  double *l; // characteristic distance = (1/3) N_{K,s} b_{K}^2
  double *A; /* prefactor = (9/2) A^{sp} z', where
	      *   A^{sp} = 3 a / Pe b_{K},
	      *   z' = (N_{K,s}/2 pi)^{3/2} (v/l_s^3)
	      *      = (3 / 2 pi b_{K}^2)^{3/2} v
	      */

  /* table for particles */
  int *ch;   /* chain type for each particle
	      * negative value == no assignement to the chain
	      */
};


#include "bonds.h" // struct bonds

/* initialize struct EV
 * INPUT
 *  bonds  : struct bonds (either fene=0 or fene=1 is fine).
 *  a, pe  : parameters for bonds parameters
 *  r2     : square of the max distance for F^{EV}
 *  v[n]   : EV parameters for each spring.
 *           the index should correspond to that in bonds.
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV, where l and A are defined by 
 *      ev->l[i] characteristic distance = (1/3) N_{K,s} b_{K}^2,
 *      ev->A[i] prefactor = (9/2) A^{sp} z',
 *    and
 *      A^{sp} = 3 a / Pe b_{K},
 *      z' = (N_{K,s}/2 pi)^{3/2} (v/l_s^3)
 *         = (3 / 2 pi b_{K}^2)^{3/2} v.
 */
struct EV *
EV_init (const struct bonds *bonds, double a, double pe,
	 double r2, const double *v,
	 int np);

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
