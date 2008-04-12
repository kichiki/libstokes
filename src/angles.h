/* header file for angles.c --
 * angle interaction between particles
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: angles.h,v 1.1 2008/04/12 18:17:30 kichiki Exp $
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
#ifndef	_ANGLES_H_
#define	_ANGLES_H_


#include "stokes.h" // struct stokes

struct angle {
  double k;  // potential factor
  double t0; // natural angle (theta_0) in radian

  // particle indices
  int ia;
  int ib;
  int ic;
};

struct angles {
  /* table for angle type */
  int n;           // number of angles
  struct angle *a; // angles
};


/**
 * struct angles
 */
struct angles *
angles_init (void);

void
angles_free (struct angles *ang);

/*
 * INPUT
 *  ia, ib, ic : particle indices (ib is the center particle)
 *  k  : potential factor
 *  t0 : natural angle (in radian)
 */
void
angles_add (struct angles *ang,
	    int ia, int ib, int ic,
	    double k, double t0);

/*
 * INPUT
 *  ang        : struct angles
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
angles_calc_force (struct angles *ang,
		   struct stokes *sys,
		   double *f,
		   int flag_add);


#endif /* !_ANGLES_H_ */
