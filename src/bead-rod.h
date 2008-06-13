/* header file for bead-rod.c --
 * bead-rod model
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bead-rod.h,v 1.2 2008/06/13 05:09:15 kichiki Exp $
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


struct BeadRod {
  int nc; // number of constraints

  // particle indices consisting of the rods
  double *a;  // a[nc] : rod distance
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
  double c1; // coefficient of linear term (4 dt/zeta)
  double c2; // coefficient of quadratic term (2 dt/zeta)^2

  double *u;  // u(t)
  double *uu; // u'(t+dt)
};

/* initialize struct BeadRod
 * INPUT
 *  nc    : number of constraints
 *  a[nc] : distances for each constraint (NULL for the uniform case a=1)
 *  ia[nc], ib[nc] : particle indices for each constraint
 */
struct BeadRod *
BeadRod_init (int nc,
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
 *  br->c1 = 4 dt / zeta
 *  br->c2 = (2 dt / zeta)^2
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
 *  a[nc]  : rod distances (NULL for a = 1.0)
 * OUTPUT
 *  u[3*nc] : connector vectors
 */
void
BeadRod_bead_to_connector (int nc,
			   const int *ia,
			   const int *ib,
			   const double *a,
			   const double *r,
			   double *u);

void
BeadRod_set_u_by_r (struct BeadRod *br,
		    const double *r);

void
BeadRod_set_uu_by_r (struct BeadRod *br,
		     const double *r);


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
BeadRod_solve_iter_gamma (struct BeadRod *br,
			  double *gamma);

/*
 * INPUT
 *  br : struct BeadRod
 *  gamma[nc] : Lagrange multiplier
 *  n  : number of beads
 * OUTPUT
 *  dr[n*3] : constraint correction for displacement
 */
void
BeadRod_constraint_displacement (struct BeadRod *br,
				 const double *gamma,
				 int n,
				 double *dr);


/**
 * NITSOL stuff
 */
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

/* solve gamma by NITSOL
 * INPUT
 *  br  : struct BeadRod
 * OUTPUT
 *  gamma[nc] :
 */
void
BeadRod_solve_gamma_by_NITSOL (struct BeadRod *br,
			       double *gamma);


/**
 * top level non-linear equasiont solver for gamma
 */
void
BeadRod_solve_gamma (struct BeadRod *br,
		     double *gamma);


/* adjust the displacement according to the constraints
 * INPUT
 *  br      : struct BeadRod
 *  n       : number of beads
 *  rr[n*3] : unconstraint positions at t + dt
 *  r [n*3] : initial positions at t
 * OUTPUT
 *  r [n*3] : corrected positions is stored.
 */
void
BeadRod_adjust_by_constraints (struct BeadRod *br,
			       int n,
			       const double *rr,
			       double *r);


#endif /* !_BEAD_ROD_H_ */
