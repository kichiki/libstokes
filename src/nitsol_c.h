/* header file for nitsol_c.c --
 * C wrappers for NITSOL
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: nitsol_c.h,v 1.2 2008/06/07 02:58:21 kichiki Exp $
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
#ifndef	_NITSOL_C_H_
#define	_NITSOL_C_H_


struct NITSOL {
  // size of the problem
  int n;

  // tolerance values
  double ftol;
  double stptol;

  // input parameters
  int input[10];

  // work area
  int lrwork;
  double *rwork;

  // parameters to the functions
  double *rpar;
  int *ipar;

  // functions
  void (*f)(int *n, double *xcur, double *fcur,
	    double *rpar, int *ipar, int *itrmf);
  void (*jacv)(int *n, double *xcur, double *fcur,
	       int *ijob, double *v, double *z,
	       double *rpar, int *ipar, int *itrmjv);
  double (*dinpr)(int* N, 
		  double* X, int* incX, 
		  double* Y, int* incY);
  double (*dnorm)(int* N, double* X, int* incX);

  // other parameters
  int iplvl;
  int ipunit;

  double choice1_exp;
  double choice2_exp;
  double choice2_coef;
  double eta_cutoff;
  double etamax;
  double thmin;
  double thmax;
  double etafixed;

  // output parameters
  int iterm;
  int info[6];
};


struct NITSOL *
NITSOL_init (void);

void
NITSOL_free (struct NITSOL *nit);

void
NITSOL_set_n (struct NITSOL *nit,
	      int n);

void
NITSOL_set_GMRES (struct NITSOL *nit,
		  int restart);

void
NITSOL_set_BiCGSTAB (struct NITSOL *nit);

void
NITSOL_set_TFQMR (struct NITSOL *nit);

/*
 * INPUT
 *  p_flag  : 0 == no preconditioning
 *            1 == right preconditioning
 *  j_flag  : 0 == approximated for J*v
 *            1 == analytical J*v
 *  j_order : 1, 2, or 4 for the approximation of J*v
 *            (ignored for j_flag == 1)
 *  jacv  : jacobian and preconditioning function
 *          (ignored for p_flag == 0 and j_flag == 0)
 */
void
NITSOL_set_jacv (struct NITSOL *nit,
		 int p_flag,
		 int j_flag,
		 int j_order,
		 void (*jacv)(int *n, double *xcur, double *fcur,
			      int *ijob, double *v, double *z,
			      double *rpar, int *ipar, int *itrmjv));

void
NITSOL_set_tol (struct NITSOL *nit,
		double ftol, double stptol);

/*
 * i do not know the way of specifying parameters is fine (for choices 1, 2, 3)
 * INPUT
 *  flag : choice 0
 *         choice 1, p1 is for exp (alpha)
 *         choice 2, p1 and p2 are for exp (alpha) and coef (gamma)
 *         choice 3, p1 is for etafixed
 *  p1, p2 : parameters
 */
void
NITSOL_set_forcing (struct NITSOL *nit,
		    int flag,
		    double p1, double p2);

/* set iplvl and ipunit
 * INPUT
 *  iplvl = 0 => no printout
 *        = 1 => iteration numbers and F-norms
 *        = 2 => ... + some stats, step norms, and linear model norms
 *        = 3 => ... + some Krylov solver and backtrack information
 *        = 4 => ... + more Krylov solver and backtrack information
 *  ipunit = printout unit number, e.g., ipunit = 6 => standard output. 
 *           NOTE: If ipunit = 0 on input, then it is set to 6 below.
 */
void
NITSOL_set_iplvl (struct NITSOL *nit,
		  int iplvl, int ipunit);


/**
 * parsing functions for the information about NITSOL
 */
void
NITSOL_parse_input (struct NITSOL *nit,
		    FILE *out);

void
NITSOL_parse_info (struct NITSOL *nit,
		   FILE *out);

void
NITSOL_parse_nitinfo (FILE *out);

void
NITSOL_parse_iterm (struct NITSOL *nit,
		    FILE *out);


/**
 * solver routine via struct NITSOL
 */
void
NITSOL_solve (struct NITSOL *nit,
	      double *x);


#endif /* !_NITSOL_C_H_ */
