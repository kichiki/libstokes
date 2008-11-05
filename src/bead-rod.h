/* header file for bead-rod.c --
 * bead-rod model
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bead-rod.h,v 1.4 2008/11/05 02:14:07 kichiki Exp $
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
#ifndef	_BEAD_ROD_H_
#define	_BEAD_ROD_H_


#include "stokes.h"

struct BeadRod {
  struct stokes *sys;

  int nc; // number of constraints

  // particle indices consisting of the rods
  double *a2;  // a2[nc] : square rod distance
  int *ia; // ia[nc] : particle at one end of the rod
  int *ib; // ib[nc] : particle at the other end of the rod

  int verbose; // 0 == quiet, 1 == verbose

  int scheme; /* 0 == iterative (Liu 1989)
	       * 1 == NITSOL (Newton-GMRES)
	       */
  double eps; // tolerance

  struct NITSOL *nit; // for scheme == 1
  struct iter *it;    // for scheme == 0

  // dynamic properties
  double dt; // time difference between u[] and uu[] below
  double d1; // = (2 dt / zeta)
  double d2; // = (2 dt / zeta)^2
  // note : the linear-term coefficient is (2 * d1)

  double *u;  // u [nc * 3], initial connector at t
  double *uu; // uu[nc * 3], connector without constraint at t+dt
};


/* initialize struct BeadRod
 * INPUT
 *  sys   : struct stokes
 *  nc    : number of constraints
 *  a[nc] : distances for each constraint
 *  ia[nc], ib[nc] : particle indices for each constraint
 */
struct BeadRod *
BeadRod_init (struct stokes *sys,
	      int nc,
	      const double *a,
	      const int *ia,
	      const int *ib);

void
BeadRod_free (struct BeadRod *br);

/* set scheme for BeadRod solver
 * INPUT
 *  br : struct BeadRod
 *  scheme : 0 == iterative (Liu 1989)
 *           1 == NITSOL (Newton-GMRES)
 *  eps    : tolerance
 */
void
BeadRod_set_scheme (struct BeadRod *br,
		    int scheme,
		    double eps);


/* set coefficients for linear and quadratic terms
 * INPUT
 *  br : struct BeadRod
 *  dt : 
 *  zeta :
 * OUTPUT
 *  br->d1 = 2 dt / zeta
 *  br->d2 = (2 dt / zeta)^2
 */
void
BeadRod_set_coefs (struct BeadRod *br,
		   double dt,
		   double zeta);

/* convert bead vectors r[] to connector vectors u[]
 * INPUT
 *  n  : number of beads
 *  nc : number of constraints (rods)
 *  r[3*n]    : position of each beads
 * OUTPUT
 *  u[3*nc] : connector vectors (NOT scaled by the distance)
 */
void
BeadRod_bead_to_connector (int nc,
			   const int *ia,
			   const int *ib,
			   const double *r,
			   double *u);

void
BeadRod_set_u_by_r (struct BeadRod *br,
		    const double *r);

void
BeadRod_set_uu_by_r (struct BeadRod *br,
		     const double *r);

/* return (i,j) element of the Rouse matrix
 * = 2 \delta_{i,j} - \delta_{ia[i], ib[j]} - \delta_{ib[i], ia[j]}
 * INPUT
 *  br : struct BeadRod
 *  i, j : constraint indices
 */
double
BeadRod_Rouse_matrix (struct BeadRod *br,
		      int i, int j);

/* nonlinear equation
 * INPUT
 *  nc : number of contraints (number of rods)
 *  x[nc] : Lagrange multipliers gamma[n]
 * OUTPUT
 *  f[nc] := (2 dt/zeta)^2(gamma_j A_{ij} u_j)^2
 *          -(4 dt/zeta) gamma_j A_{ij} u'_i u_j
 *          +(|u'_i|^2 - 1)
 */
void
BeadRod_NLEQ_for_gamma (int nc, const double *gamma, double *f, void *params);


/* update gamma iteratively
 * INPUT
 *  br : struct BeadRod
 *  x[nc] : Lagrange multipliers gamma[n]
 * OUTPUT
 *  (4 dt/zeta) gamma_j A_{ij} u'_i u_j
 *        = (2 dt/zeta)^2(gamma_j A_{ij} u_j)^2
 *          +(|u'_i|^2 - 1)
 */
void
BeadRod_update_gamma (struct BeadRod *br,
		      const double *gamma0,
		      double *gamma);


/* solve gamma iteratively : Liu (1989)
 * INPUT
 *  br : struct BeadRod
 * OUTPUT
 *  gamma[nc] :
 */
void
BeadRod_solve_gamma_iter (struct BeadRod *br,
			  double *gamma);


/* a wrapper of the nonlinear equations for NITSOL
 * INTPUT
 *  xcur[n] : gamma
 *  rpar    : (struct BeadRod *)br
 *  ipar    : not used.
 * OUTPUT
 *  itrmf   : 0 means success
 *  fcur[n] :
 */
void
BeadRod_NITSOL_f (int *n, double *xcur, double *fcur,
		  double *rpar, int *ipar, int *itrmf);


/* solve gamma for the constraints
 * either by linear approximation iteratively (Liu 1989)
 * or by Newton-GMRES method in NITSOL
 * INPUT
 *  br : struct BeadRod
 *  gamma[nc] : all elements should be initialized.
 * OUTPUT
 *  gamma[nc] :
 */
void
BeadRod_solve_gamma (struct BeadRod *br,
		     double *gamma);


/*
 * INPUT
 *  br : struct BeadRod
 *  gamma[nc] : Lagrange multiplier
 *  n  : number of beads
 *  a  : prefactor for dr[n*3], that is, a * dr[] + (correction) is returned.
 *       so that if you give the unconstraint position in dr[]
 *       and give a = 1, then the resultant dr[] is the constraint position.
 *       if a = 0, dr[] given is ignored and
 *       dr[] in return is just the correction
 * OUTPUT
 *  dr[n*3] := a * dr(input) + correction.
 */
void
BeadRod_constraint_displacement (struct BeadRod *br,
				 const double *gamma,
				 int n,
				 double a,
				 double *dr);

/* adjust the displacement according to the constraints
 * INPUT
 *  br      : struct BeadRod
 *  n       : number of beads
 *  r [n*3] : initial positions at t
 *  rr[n*3] : unconstraint positions at t + dt (THIS IS OVERWRITTEN)
 * OUTPUT
 *  rr[n*3] : corrected positions is stored.
 */
void
BeadRod_adjust_by_constraints (struct BeadRod *br,
			       int n,
			       const double *r,
			       double *rr);


/**
 * HI version
 */
/* calc vector g[nc * 3] with Lagrange multiplier gamma[nc] 
 * and the correction displacement dr[n * 3]
 * INPUT
 *  br     : struct BeadRod
 *  gamma[nc] : Lagrange multiplier
 * OUTPUT
 *  g[nc*3]   :
 */
void
BeadRod_calc_g (struct BeadRod *br,
		const double *gamma,
		double *g);

/*
 * INPUT
 *  br : struct BeadRod
 *  gamma[nc] : Lagrange multiplier
 *  a  : prefactor for dr[n*3], that is, a * dr[] + (correction) is returned.
 *       so that if you give the unconstraint position in dr[]
 *       and give a = 1, then the resultant dr[] is the constraint position.
 *       if a = 0, dr[] given is ignored and
 *       dr[] in return is just the correction
 * OUTPUT
 *  dr[nm*3] := a * dr(input) + correction.
 */
void
BeadRod_calc_dr (struct BeadRod *br,
		 const double *gamma,
		 double a,
		 double *dr);

/* nonlinear equation
 * INPUT
 *  nc : number of contraints (number of rods)
 *  x[nc] : Lagrange multipliers gamma[n]
 * OUTPUT
 *  f[nc] := (2 dt/zeta)^2(gamma_j A_{ij} u_j)^2
 *          -(4 dt/zeta) gamma_j A_{ij} u'_i u_j
 *          +(|u'_i|^2 - 1)
 */
void
BeadRod_NLEQ_for_gamma_ (struct BeadRod *br,
			 const double *gamma,
			 double *f);

/* update gamma iteratively
 * INPUT
 *  br : struct BeadRod
 *  x[nc] : Lagrange multipliers gamma[n]
 * OUTPUT
 *  (4 dt/zeta) gamma_j A_{ij} u'_i u_j
 *        = (2 dt/zeta)^2(gamma_j A_{ij} u_j)^2
 *          +(|u'_i|^2 - 1)
 */
void
BeadRod_update_gamma_ (struct BeadRod *br,
		       const double *gamma0,
		       double *gamma);

/* solve gamma iteratively : Liu (1989)
 * INPUT
 *  br : struct BeadRod
 * OUTPUT
 *  gamma[nc] :
 */
void
BeadRod_solve_gamma_iter_ (struct BeadRod *br,
			   double *gamma);


#endif /* !_BEAD_ROD_H_ */
