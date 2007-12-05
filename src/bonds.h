/* header file for bonds.c --
 * bond interaction between particles
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds.h,v 1.5 2007/12/05 03:43:13 kichiki Exp $
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

  int *nex;    // number of excluded particles in the chain
};

struct list_ex {
  int np;  // total number of particles (must be equal to sys->np)
  int *n;  // n[np] : number of excluded particles for each particles
  int **i; // i[j][k] : k-th particle index to exclude for particle j.
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
 *  nex    : number of excluded particles in the chain
 * OUTPUT
 *  bonds  :
 */
void
bonds_add_type (struct bonds *bonds,
		int type, int fene, double p1, double p2,
		int nex);

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


/**
 * SWIG utility routine
 * For examplean, expected usage in python by SWIG:
 *   n = stokes.bonds_get_pairs_n(bonds, i)
 *   ia = stokes.iarray(n)
 *   ib = stokes.iarray(n)
 *   stokes.bonds_get_pairs(bonds, i, ia, ib)
 * then, you have arrays ia[n] and ib[n].
 */

/* to get the number of pairs for the bond "i"
 */
int
bonds_get_pairs_n (struct bonds *b, int i);
/* 
 * INPUT
 *  b : struct bonds
 *  i : index of the bond
 *  ia[n] : "a" particle index of j-th pair for the bond "i",
 *  ib[n] : "b" particle index of j-th pair for the bond "i",
 *          where j runs from 0 to (n-1) and
 *          n = b->pairs[i]->n is the number of pairs for the bond "i".
 *          before calling, allocate ia and ib with (sizeof(int) * n).
 * OUTPUT
 *  ia[n], ib[n] : 
 */
void
bonds_get_pairs (struct bonds *b, int i,
		 int *ia, int *ib);


/**
 * exclusion list for lubrication due to the bonding
 */
struct list_ex *
list_ex_init (int np);

void
list_ex_add (struct list_ex *ex, int j, int k);

void
list_ex_free (struct list_ex *ex);

struct list_ex *
list_ex_copy (struct list_ex *ex0);

/* construct the excluded list by struct bonds
 */
void
list_ex_set_by_bonds (struct list_ex *ex, const struct bonds *b);

/* check whether j is excluded for i
 * INPUT
 *  ex : struct list_ex
 *  i  : particle now we are considering
 *  j  : particle whether it is in the list or not.
 * OUTPUT
 *  returned value : 0 (false); j is NOT in the excluded list for i.
 *                   1 (true);  j IS in the excluded list for i.
 *                   NOTE, the self for i in some chain is EXCLUDED,
 *                   while the self for particles is NOT excluded.
 */
int
list_ex_check (struct list_ex *ex, int i, int j);


#endif /* !_BONDS_H_ */
