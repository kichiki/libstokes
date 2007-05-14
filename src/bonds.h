/* header file for bonds.c --
 * bond interaction between particles
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds.h,v 1.2 2007/05/14 00:19:59 kichiki Exp $
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
#ifndef	_BONDS_H_
#define	_BONDS_H_


#include "stokes.h" // struct stokes


struct bond_pairs {
  int n;   // number of pairs
  int *ia; // particle a for the pair
  int *ib; // particle b for the pair
};

struct bonds {
  int ntypes;                // number of bond-types
  double *k;                 // spring constant for the bond-type
  double *r0;                // equilibrium distance for the bond-type
  struct bond_pairs **pairs; // pairs for the bond-type
};


struct bond_pairs *
bond_pairs_init (void);

void
bond_pairs_free (struct bond_pairs *pairs);

void
bond_pairs_add (struct bond_pairs *pairs,
		int ia, int ib);


struct bonds *
bonds_init (void);

void
bonds_free (struct bonds *bonds);

void
bonds_add_type (struct bonds *bonds,
		double k, double r0);

/*
 * INPUT
 *  b          : struct bond
 *  sys        : struct stokes
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
bonds_calc_force (struct bonds *bonds,
		  struct stokes *sys,
		  double *f,
		  int flag_add);

void
fprint_bonds (FILE *out, struct bonds *bonds);


#endif /* !_BONDS_H_ */
