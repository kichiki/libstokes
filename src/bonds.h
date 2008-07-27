/* header file for bonds.c --
 * bond interaction between particles
 * Copyright (C) 2007-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds.h,v 1.12 2008/07/27 00:51:03 kichiki Exp $
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


// for each group
struct BONDS {
  int n; // number of bonds
  int *type; /* type[n] : type of the spring for each bond
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
	      * 6 : discrete Wormlike chain (dWLC), where
	      *     p1 = k, dimensionless spring constant, 
	      *     p2 = r0, the natural length [nm],
	      *     the potential is given by
	      *     U(r) = (k/2) * (kT / r0^2) * (r-r0)^2
	      * 7 : FENE-Fraenkel
	      *     p1, p2 are the same for FENE chains.
	      *     p3 = s, the tolerance
	      */
  int *fene; /* fene[n] : flag for the parameters p1 and p2
	      * 0 : (p1, p2) = (A^{sp}, L_{s})
	      * 1 : (p1, p2) = (N_{K,s}, b_{K})
	      */
  double *p1; // p1[n] : the first parameter (k or N_{K,s})
  double *p2; // p2[n] : the second parameter (r0 or b_{K})
  double *p3; // p3[n] : the third parameter (tol for FENE-Fraenkel)

  int *ia;    // ia[n] : particle index of one end of the bond
  int *ib;    // ib[n] : particle index of the other end of the bond
};


struct BONDS *
BONDS_init (void);

void
BONDS_free (struct BONDS *b);

/*
 * INPUT
 *  ia, ib : (global) particle index, that is, they are in the range [0, NP), 
 *           where NP is the total number of particles.
 */
void
BONDS_append (struct BONDS *b,
	      int type,
	      int fene,
	      double p1,
	      double p2,
	      double p3,
	      int ia,
	      int ib);


void
BONDS_sort_by_ia (struct BONDS *b);


/* set FENE spring parameters for run
 * INPUT
 *  bonds : p1 (N_{K,s}) and p2 (b_{K}) are used
 *          (bonds->fene[i] == 1 is expected), or 
 *          p1 (k) and p2 (r0) are used for dWLC spring (type == 6).
 *            in this case, potential is given by
 *            (k/2) * (kT / r0^2) * (r-r0)^2
 *  length : length scale in the simulation
 *  peclet : peclet number
 * OUTPUT
 *  bonds->p1[] := A^{sp} = 3 length / (peclet b_{K})
 *  bonds->p2[] := Ls / length = N_{K,s} b_{K} / length
 *    or, for dWLC spring, the conversions are given by 
 *  bonds->p1[] := A^{sp} = k / (pe * (r0/length)^2)
 *  bonds->p2[] := Ls / length = r0 / length
 */
void
BONDS_set_FENE (struct BONDS *b,
		double length, double peclet);

/* return force function (scalar part) of the spring
 */
double
BONDS_fr_i (struct BONDS *b,
	    int bond_index,
	    double Q);

/* calc spring force for the bond "i"
 * INPUT
 *  bonds : struct BONDS
 *  ib    : bond index
 *  q[3]  : connector vector for the bond "i"
 * OUTPUT
 *  f[3]  : force due to the bond "i"
 */
void
BONDS_calc_force_spring_i (struct BONDS *bonds,
			   int ib,
			   const double *q,
			   double *f);

/*
 * INPUT
 *  b          : struct BONDS
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
BONDS_calc_force (struct BONDS *b,
		  struct stokes *sys,
		  double *f,
		  int flag_add);


#endif /* !_BONDS_H_ */
