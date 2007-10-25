/* header file for bonds.c --
 * bond interaction between particles
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds.h,v 1.3 2007/10/25 05:54:23 kichiki Exp $
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
  /* table for bond type */
  int n;      // number of bond types
  int *type;  /* type of the spring
	       * 0 : Hookean (p1 = k, spring constant,
	       *              p2 = r0, natural length)
	       * 1 : Wormlike chain (WLC)
	       * 2 : inverse Langevin chain (ILC)
	       * 3 : Cohen's Pade approx for ILC
	       * 4 : Werner spring (approx for ILC)
	       * 5 : another Hookean
	       *     where for these FENE chains,
	       *     p1 = N_{K,s} the Kuhn steps for a spring
	       *     p2 = b_{K}   the Kuhn length [nm]
	       */
  int *fene;   /* flag for the parameters p1 and p2
		* 0 : (p1, p2) = (A^{sp}, L_{s})
		* 1 : (p1, p2) = (N_{K,s}, b_{K})
		*/
  double *p1;  // the first parameter (k or N_{K,s})
  double *p2;  // the second parameter (r0 or b_{K})

  struct bond_pairs **pairs; // pairs for the bond
};


void
bond_pairs_free (struct bond_pairs *pairs);

void
bond_pairs_add (struct bond_pairs *pairs,
		int ia, int ib);


/* initialize struct bonds
 * INPUT
 * OUTPUT
 *  returned value : struct bonds
 */
struct bonds *
bonds_init (void);

void
bonds_free (struct bonds *bonds);

/* add a spring into bonds
 * INPUT
 *  bonds  : struct bonds
 *  type   : type of the spring
 *  fene   : 0 == (p1,p2) are (A^{sp}, L_{s})
 *           1 == (p1, p2) = (N_{K,s}, b_{K})
 *  p1, p2 : spring parameters
 * OUTPUT
 *  bonds  :
 */
void
bonds_add_type (struct bonds *bonds,
		int type, int fene, double p1, double p2);

/* set FENE spring parameters for run
 * INPUT
 *  bonds : p1 (N_{K,s}) and p2 (b_{K}) are used.
 *  a     : length scale in the simulation
 *  pe    : peclet number
 * OUTPUT
 *  bonds->p1[] := A^{sp} = 3a / pe b_{K}
 *  bonds->p2[] := Ls / a = N_{K,s} b_{K} / a 
 */
void
bonds_set_FENE (struct bonds *bonds,
		double a, double pe);

/*
 * INPUT
 *  bonds      : struct bond
 *  sys        : struct stokes (only nm and pos are used)
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
