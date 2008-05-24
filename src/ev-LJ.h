/* header file for ev-LJ.c --
 * excluded-volume interactions in Lennard-Jones type
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-LJ.h,v 1.1 2008/05/24 05:47:12 kichiki Exp $
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
#ifndef	_EV_LJ_H_
#define	_EV_LJ_H_

#include "stokes.h" // struct stokes


struct EV_LJ {
  /* LJ parameters, where the potential is given in DIMENSIONAL FORM as 
   *  U(r) = e ((r0/r)^12 - 2(r0/r)^6)
   * Note that r0 = 2^{1/6} sigma is the distance of potential minimum.
   * The repulsive force is evaluated by taking r = r0 at the contact point.
   * The parameters e and r0 should be dimensionless numbers defined by 
   *         e = hat(e) for flag == 0, or
   *         e = hat(e) / (peclet * hat(r0)) for flag == 1,
   *         r0 = hat(r0) = "r0" / length
   * where
   *         hat(e)  = "e" / kT
   *         hat(r0) = "r0" / length
   * and "e" and "r0" are dimensional numbers.
   * The dimensionless force (scalar part) is then given by
   *  hat(F) = 12 e (r0/r)^7 ((r0/r)^6 - 1)
   * where e is in the form of flag == 1.
   */
  /* parameters for each particle */
  int n;      // number of particles
  int    *flag;
  double *e;  /* always dimensionless
	       * either e/kT              for flag==0
	       * or     e/(kT Pe hat(r0)) for flag==1
	       */
  double *r0; // always dimensionless (scaled by "length")
};


/* initialize struct EV_LJ
 * INPUT
 *  np     : number of particles
 * OUTPUT
 *  returned value : struct EV_LJ,
 *      where LJ parameters are set by zero.
 */
struct EV_LJ *
EV_LJ_init (int np);

void
EV_LJ_free (struct EV_LJ *ev_LJ);

/* scale LJ parameter e for runs
 * INPUT
 *  ev_LJ  : struct EV_LJ
 *  peclet : peclet number
 */
void
EV_LJ_scale (struct EV_LJ *ev_LJ,
	     double peclet);

/*
 * INPUT
 *  ev_LJ      : struct EV_LJ
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_LJ_calc_force (struct EV_LJ *ev_LJ,
		  struct stokes *sys,
		  double *f,
		  int flag_add);


#endif /* !_EV_DH_H_ */
