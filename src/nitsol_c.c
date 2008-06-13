/* C wrappers for NITSOL
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: nitsol_c.c,v 1.3 2008/06/13 02:56:08 kichiki Exp $
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
#include <stdio.h>
#include <stdlib.h>
#include "memory-check.h" // macro CHECK_MALLOC

#include <nitsol.h>
#include "nitsol_c.h" // struct NITSOL


// BLAS functions
double
ddot_(int* N, 
      double* X, int* incX, 
      double* Y, int* incY);
double
dnrm2_(int* N, 
       double* X, int* incX);


struct NITSOL *
NITSOL_init (void)
{
  struct NITSOL *nit = (struct NITSOL *)malloc (sizeof (struct NITSOL));
  CHECK_MALLOC (nit, "NITSOL_init");

  nit->n = 0;

  nit->ftol = 0.0;
  nit->stptol = 0.0;

  int i;
  for (i = 0; i < 10; i ++)
    {
      nit->input[i] = 0;
    }

  nit->lrwork = 0;
  nit->rwork = NULL;

  nit->rpar = NULL;
  nit->ipar = NULL;

  nit->f = NULL;
  nit->jacv = NULL;
  nit->dinpr = NULL;
  nit->dnorm = NULL;

  /*
c If diagnostic information is desired, include this common block in the 
c main program and set iplvl and ipunit according to the following: 
c
c     iplvl = 0 => no printout
c           = 1 => iteration numbers and F-norms
c           = 2 => ... + some stats, step norms, and linear model norms
c           = 3 => ... + some Krylov solver and backtrack information
c           = 4 => ... + more Krylov solver and backtrack information
  */
  nit->iplvl = 0;
  nit->ipunit = 6; // stdout

  nit->choice1_exp  = 0.0;
  nit->choice2_exp  = 0.0;
  nit->choice2_coef = 0.0;
  nit->eta_cutoff   = 0.0;
  nit->etamax       = 0.0;
  nit->thmin        = 0.0;
  nit->thmax        = 0.0;
  nit->etafixed     = 0.0;

  nit->iterm = 0;
  for (i = 0; i < 6; i ++)
    {
      nit->info[i] = 0;
    }

  return (nit);
}

void
NITSOL_free (struct NITSOL *nit)
{
  if (nit == NULL) return;
  if (nit->rwork != NULL) free (nit->rwork);
  free (nit);
}

/* allocate work area
 * based on the following explanation in bratu.f
c>>> Alternative parameter statements for different circumstances -- HFW. 
c The following is always safe but may require a little unnecessary storage.
      parameter ( LRWORK=MAXN*(MAXKD+10)+MAXKD*(MAXKD+3))
c The following can be used if the compiler allows the "max". 
c      parameter ( LRWORK=max(11*MAXN,MAXN*(MAXKD+5)+MAXKD*(MAXKD+3)) )
c The following can be used if MAXKD > 5.
c      parameter ( LRWORK=MAXN*(MAXKD+5)+MAXKD*(MAXKD+3))
 */
static void
NITSOL_alloc_rwork (struct NITSOL *nit)
{
  if (nit->rwork != NULL) free(nit->rwork);

  if (nit->n == 0) return;

  int maxkd = nit->input[3];
  if (maxkd == 0)
    {
      nit->lrwork = nit->n * 11;
    }
  else
    {
      int n1 = nit->n * 11;
      int n2 = nit->n * (maxkd + 5) + maxkd * (maxkd + 3);
      if (n1 > n2) nit->lrwork = n1;
      else         nit->lrwork = n2;
    }

  nit->rwork = (double *)malloc (sizeof (double) * nit->lrwork);
  CHECK_MALLOC (nit->rwork, "NITSOL_set_GMRES");
}


void
NITSOL_set_n (struct NITSOL *nit,
	      int n)
{
  nit->n = n;

  NITSOL_alloc_rwork (nit);
}


void
NITSOL_set_GMRES (struct NITSOL *nit,
		  int restart)
{
  nit->input[2] = 0; // GMRES
  nit->input[3] = restart;

  NITSOL_alloc_rwork (nit);
}

void
NITSOL_set_BiCGSTAB (struct NITSOL *nit)
{
  nit->input[2] = 1; // BiCGSTAB
  nit->input[3] = 0;

  NITSOL_alloc_rwork (nit);
}

void
NITSOL_set_TFQMR (struct NITSOL *nit)
{
  nit->input[2] = 2; // TFQMR
  nit->input[3] = 0;

  NITSOL_alloc_rwork (nit);
}


/* do nothing for jacobian and preconditioner
 */
static void
NITSOL_no_jacv (int *n, double *xcur, double *fcur,
		int *ijob, double *v, double *z,
		double *rpar, int *ipar, int *itrmjv)
{
  if (*ijob == 0)
    {
      // calc Jacobian-vector product
      return;
    }
  else if (*ijob == 1)
    {
      // apply the preconditioner
      return;
    }
  else
    {
      fprintf (stderr, "# NITSOL_no_jacv: invalid ijob %d\n", *ijob);
      return;
    }
}

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
			      double *rpar, int *ipar, int *itrmjv))
{
  nit->input[4] = p_flag; /* 0 == no preconditioning
			  * 1 == right preconditioning
			  */

  nit->input[1] = j_flag; /* 0 == approximated for J*v
			  * 1 == analytical J*v
			  */
  if (nit->input[1] == 0)
    {
      if (j_order <= 0 || j_order == 3 || j_order > 4)
	{
	  fprintf (stderr, "# NITSOL_set_Jacobian: invalid order %d."
		   " force to set 4.\n", j_order);
	  j_order = 4;
	}
      nit->input[7] = j_order;
    }
  else
    {
      nit->input[7] = 0;
    }

  // set function jacv()
  if (p_flag != 0 || j_flag != 0)
    {
      nit->jacv = jacv;
    }
  else
    {
      //nit->jacv = NULL;
      nit->jacv = NITSOL_no_jacv;
    }
}

void
NITSOL_set_tol (struct NITSOL *nit,
		double ftol, double stptol)
{
  nit->ftol = ftol;
  nit->stptol = stptol;
}

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
		    double p1, double p2)
{
  // zero clear for the forcing parameters
  nit->choice1_exp  = 0.0;
  nit->choice2_exp  = 0.0;
  nit->choice2_coef = 0.0;
  nit->eta_cutoff   = 0.0;
  nit->etamax       = 0.0;
  nit->thmin        = 0.0;
  nit->thmax        = 0.0;
  nit->etafixed     = 0.0;

  if (flag == 1)
    {
      nit->choice1_exp  = p1;
    }
  else if (flag == 2)
    {
      nit->choice2_exp  = p1; // alpha
      nit->choice2_coef = p2; // gamma
    }
  else if (flag == 3)
    {
      nit->etafixed     = p1;
    }
  else if (flag != 0)
    {
      fprintf (stderr, "# NITSOL_set_forcing: invalid forcing type %d."
	       " force to set 0\n", flag);
      flag = 0;
    }
  nit->input[9] = flag;
}

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
		  int iplvl, int ipunit)
{
  nit->iplvl = iplvl;
  nit->ipunit = ipunit;
}


/* set routines for calculating norm by BLAS routines (ddot_ and dnrm2_)
 */
void
NITSOL_set_norm_by_BLAS (struct NITSOL *nit)
{
  nit->dinpr = ddot_;
  nit->dnorm = dnrm2_;
}


/**
 * parsing functions for the information about NITSOL
 */
void
NITSOL_parse_input (struct NITSOL *nit,
		    FILE *out)
{
  // Optional every-user input:
  // input(1) = nnimax = maximum number of nonlinear iterations (default 200).
  fprintf (out, "# NITSOL: nnimax = %d\n", nit->input[0]);


  /*
     input(2) = ijacv = flag for determining the method of J*v evaluation:
                  0 => finite-difference evaluation (default) 
                  1 => analytic evaluation
  */
  if (nit->input[1] == 0)
    {
      fprintf (out, "# NITSOL: J*v is calculated by finite-difference.\n");
    }
  else if (nit->input[1] == 1)
    {
      fprintf (out, "# NITSOL: J*v is calculated analytically.\n");
    }
  else
    {
      fprintf (out, "# NITSOL: invalid ijacv parameter %d\n", nit->input[1]);
    }

  /*
     input(3) = ikrysl = flag for determining the Krylov solver: 
                  0 => GMRES (default)
                  1 => BiCGSTAB
                  2 => TFQMR

                For brief descriptions of the solvers plus references, 
                see the subroutines nitgm, nitstb, and nittfq. 
  */
  if (nit->input[2] == 0)
    {
      fprintf (out, "# NITSOL: GMRES method for the Krylov solver\n");
    }
  else if (nit->input[2] == 1)
    {
      fprintf (out, "# NITSOL: BiCGSTAB method for the Krylov solver\n");
    }
  else if (nit->input[2] == 2)
    {
      fprintf (out, "# NITSOL: TFQMR method for the Krylov solver\n");
    }
  else
    {
      fprintf (out, "# NITSOL: invalid ikrysl parameter %d\n", nit->input[2]);
    }

  /*
     input(4) = kdmax = maximum Krylov subspace dimension when GMRES is used 
                (default 20). 
  */
  fprintf (out, "# NITSOL: kdmax = %d\n", nit->input[3]);

  /*
     input(5) = irpre = flag for right preconditioning: 
                  0 => no right preconditioning
                  1 => right preconditioning
  */
  if (nit->input[4] == 0)
    {
      fprintf (out, "# NITSOL: no right preconditioning\n");
    }
  else if (nit->input[4] == 1)
    {
      fprintf (out, "# NITSOL: right preconditioning\n");
    }
  else
    {
      fprintf (out, "# NITSOL: invalid irpre parameter %d\n", nit->input[4]);
    }


  // Optional experienced user input:
  /*
     input(6) = iksmax = maximum allowable number of iterations per call 
                to the Krylov solver routine (default 1000). 
  */
  fprintf (out, "# NITSOL: iksmax = %d\n", nit->input[5]);

  /*
     input(7) = iresup = residual update flag when GMRES is used; on 
                restarts, the residual is updated as follows: 
                  0 => linear combination (default) 
                  1 => direct evaluation
                The first is cheap (one n-vector saxpy) but may lose 
                accuracy with extreme residual reduction; the second 
                retains accuracy better but costs one J*v product per 
                restart. 
  */
  if (nit->input[6] == 0)
    {
      fprintf (out, "# NITSOL: linear conbination for residual update for GMRES.\n");
    }
  else if (nit->input[6] == 1)
    {
      fprintf (out, "# NITSOL: direct evaluation for residual update for GMRES.\n");
    }
  else
    {
      fprintf (out, "# NITSOL: invalid iresup parameter %d\n", nit->input[6]);
    }

  /*
     input(8) = ifdord = order of the finite-difference formula (sometimes) 
                used when input(2) = ijacv = 0. When input(2) = ijacv = 0, 
                this must be 0, 1, 2, or 4 on input; otherwise, it is 
                irrelevant. With input(2) = ijacv = 0, the precise 
                meaning is as follows: 
 
                If GMRES is used, then ifdord matters only if input(7) = 
                iresup = 1, in which case it determines the order of 
                the finite-difference formula used in evaluating the 
                initial residual at each GMRES restart (default 2); if 
                ifdord = 0 on input, then it is set to 2 below. NOTE: This 
                only affects initial residuals at restarts; first-order 
                differences are always used within each GMRES cycle. Using 
                higher-order differences at restarts only should give 
                the same accuracy as if higher-order differences were 
                used throughout; see K. Turner and H. F. Walker, "Efficient 
                high accuracy solutions with GMRES(m)," SIAM J. Sci. 
                Stat. Comput., 13 (1992), pp. 815--825. 
                
                If BiCGSTAB or TFQMR is used, then ifdord determines the 
                order of the finite-difference formula used at each 
                iteration (default 1); if ifdord = 0 on input, then it 
                is set to 1 below. 
  */
  fprintf (out, "# NITSOL: ifdord = %d\n", nit->input[7]);

  /*
     input(9) = ibtmax = maximum allowable number of backtracks (step 
                reductions) per call to nitbt (default 10). 
 
                USAGE NOTE: Backtracking can be turned off by setting 
        		ibtmax = -1. Other negative values of ibtmax are not 
                valid. 
  */
  fprintf (out, "# NITSOL: ibtmax = %d\n", nit->input[8]);

  /*
     input(10) = ieta = flag determining the forcing term eta as follows: 
                  0 => abs( ||fcur|| - ||fprev+Jprev*sprev|| )/||fprev||
                       (default) 
                  1 => (||fcur||/||fprev||)**2 
                  2 => gamma*(||fcur||/||fprev||)**alpha 
                       for user-supplied gamma in (0,1] and alpha in (1,2] 
                  3 => fixed (constant) eta in (0,1), either 0.1 (default) 
        		       or specified by the user (see USAGE NOTE below) 
                Here, fcur = current f, fprev = previous f, etc. The Krylov 
                iterations are terminated when an iterate s satisfies 
                an inexact Newton condition ||F + J*s|| .le. eta*||F||.
 
                USAGE NOTE: If input(10) = ieta = 2, then alpha and gamma 
                must be set in common block nitparam.h as described below. 
        		If input(10) = ieta = 3, then the desired constant eta may 
        		be similarly set in nitparam.h if a value other than the 
        		default of 0.1 is desired. 
                
                The first three expressions above are from S. C. Eisenstat 
                and H. F. Walker, "Choosing the forcing terms in an inexact 
                Newton method", SIAM J. Scientific Computing, 17 (1996), 
                pp. 16--32. (They may be modified according to certain 
                safeguards in subroutine nitdrv.) The first gives convergence 
                that is q-superlinear and of r-order (1+sqrt(5))/2. The 
                second gives convergence that is r-quadratic and of q-order 
                p for every p in [1,2). The third gives convergence that is 
                of q-order alpha when gamma < 1 and, when gamma = 1, of 
                r-order alpha and q-order p for every p in [1,alpha). The 
                fourth gives q-linear convergence with asymptotic rate 
                constant eta in a certain norm; see R. S. Dembo, S. C. 
        		Eisenstat, and T. Steihaug, "Inexact Newton methods", 
                SIAM J. Numer. Anal., 18 (1982), pp. 400-408. 
 
                Of these four choices, the 1st is usually satisfactory, 
                the 2nd or 3rd is sometimes preferred, and the 4th may be 
                useful in some situations, e.g., it may be desirable to 
                choose a fairly large fixed eta in (0,1), such as eta = .1, 
                when numerical inaccuracy prevents the Krylov solver 
                from obtaining much residual reduction. 
  */
  if (nit->input[9] == 0)
    {
      fprintf (out, "# NITSOL: type 0 forcing"
	       " abs( ||fcur|| - ||fprev+Jprev*sprev|| )/||fprev||\n");
    }
  else if (nit->input[9] == 1)
    {
      fprintf (out, "# NITSOL: type 1 forcing"
	       " (||fcur||/||fprev||)**2\n");
    }
  else if (nit->input[9] == 2)
    {
      fprintf (out, "# NITSOL: type 2 forcing"
	       " gamma*(||fcur||/||fprev||)**alpha \n");
    }
  else if (nit->input[9] == 3)
    {
      fprintf (out, "# NITSOL: type 3 forcing by fixed (constant) eta.\n");
    }
  else
    {
      fprintf (out, "# NITSOL: invalid ieta parameter %d\n", nit->input[9]);
    }
}

void
NITSOL_parse_info (struct NITSOL *nit,
		   FILE *out)
{
  /*
c Further explanation of info: On output, the components of info are 
c as follows: 
c
c     info(1)   = nfe (number of function evaluations)
c     info(2)   = njve (number of J*v evaluations)
c     info(3)   = nrpre (number of P(inverse)*v evaluations)
c     info(4)   = nli (number of linear iterations)
c     info(5)   = nni (number of nonlinear iterations)
c     info(6)   = nbt (number of backtracks)
   */
  fprintf (out, "# NITSOL: No. function evaluations     %d\n", nit->info[0]);
  fprintf (out, "# NITSOL: No. J*v evaluations          %d\n", nit->info[1]);
  fprintf (out, "# NITSOL: No. P(inverse)*v evaluations %d\n", nit->info[2]);
  fprintf (out, "# NITSOL: No. linear iterations        %d\n", nit->info[3]);
  fprintf (out, "# NITSOL: No. nonlinear iterations     %d\n", nit->info[4]);
  fprintf (out, "# NITSOL: No. backtracks               %d\n", nit->info[5]);
}

void
NITSOL_parse_nitinfo (FILE *out)
{
  /*
  For passing information about the nonlinear iterations to user-supplied 
  subroutines: 
 
      include 'nitinfo.h'
 
  If information on the current state of the nonlinear iteration is
  desired in a user-supplied subroutine (for example, deciding 
  whether to update a preconditioner), include this common block
  in the subroutine. The variables are as follows: 
  */
  /*
     instep - inexact Newton step number. 
  */
  fprintf (out, "# NITSOL: inexact Newton step number %d\n",
	   nitinfo_.instep);
  /*
     newstep - set to 0 at the beginning of an inexact Newton step.
               This may be checked in a user-supplied jacv to decide
               whether to update the preconditioner.  If you test on
               newstep .eq. 0 to determine whether to take some 
               special action at the beginning of a nonlinear iteration, 
               you must also set newstep to some nonzero value to
               subsequently avoid taking that action unnecessarily. 
  */
  if (nitinfo_.newstep == 0)
    {
      fprintf (out, "# NITSOL: at the beginning of an inexact Newton step\n");
    }
  else
    {
      fprintf (out, "# NITSOL: not at the beginning of an inexact Newton step\n");
    }

  /*
     krystat - status of the Krylov iteration; same as itrmks (see 
               the nitsol documentation). 
  */
  fprintf (out, "# NITSOL: krystat = %d\n", nitinfo_.krystat);

  /*
     avrate  - average rate of convergence of the Krylov solver during
               the previous inexact Newton step.  This may be checked
               in a user-supplied jacv to decide when to update the
               preconditioner.
  */
  fprintf (out, "# NITSOL: avrate = %e\n", nitinfo_.avrate);

  /*
     fcurnrm - ||f(xcur)||. 
   */
  fprintf (out, "# NITSOL: fcurnrm = %e\n", nitinfo_.fcurnrm);
}

void
NITSOL_parse_iterm (struct NITSOL *nit,
		    FILE *out)
{
  /*
c  iterm   = termination flag; values have the following meanings: 
c             -k => illegal value in input(k). 
c              0 => normal termination: ||F||.le.ftol or ||step||.le.stptol.
c              1 => nnimax nonlinear iterations reached without success. 
c              2 => failure to evaluate F.
c              3 => in nitjv, J*v failure. 
c              4 => in nitjv, P(inverse)*v failure. 
c              5 => in nitdrv, insufficient initial model norm reduction 
c                   for adequate progress. NOTE: This can occur for several 
c                   reasons; examine itrmks on return from the Krylov 
c                   solver for further information. (This will be printed out 
c                   if iplvl .ge. 3, see the discussion of optional 
c                   common blocks below). 
c              6 => in nitbt, failure to reach an acceptable step through 
c                   backtracking. 
  */
  if (nit->iterm < 0)
    {
      fprintf (out, "# NITSOL: illeagal value in input[%d]\n",
	       -(nit->iterm)-1);
    }
  else
    {
      switch (nit->iterm)
	{
	case 0:
	  fprintf (out, "# NITSOL: normally terminated.\n");
	  break;
	case 1:
	  fprintf (out, "# NITSOL: nnimax nonlinear iterations (%d)"
		   " reached without success.\n",
		   nit->input[0]);
	  break;
	case 2:
	  fprintf (out, "# NITSOL: failure to evaluate F.\n");
	  break;
	case 3:
	  fprintf (out, "# NITSOL: in nitjv, J*v failure.\n");
	  break;
	case 4:
	  fprintf (out, "# NITSOL: in nitjv, P(inverse)*v failure.\n");
	  break;
	case 5:
	  fprintf (out, "# NITSOL: in nitdrv, insufficient initial model norm reduction.\n");
	  break;
	case 6:
	  fprintf (out, "# NITSOL: in nitbt, failure to reach an acceptable step through backtracking.\n");
	  break;
	default:
	  fprintf (out, "# NITSOL: invalid iterm %d\n", nit->iterm);
	  break;
	}
    }
}


/**
 * solver routine via struct NITSOL
 */
void
NITSOL_solve (struct NITSOL *nit,
	      double *x)
{
  nitprint_.iplvl  = nit->iplvl;
  nitprint_.ipunit = nit->ipunit;

  nitparam_.choice1_exp  = nit->choice1_exp;
  nitparam_.choice2_exp  = nit->choice2_exp;
  nitparam_.choice2_coef = nit->choice2_coef;
  nitparam_.eta_cutoff   = nit->eta_cutoff;
  nitparam_.etamax       = nit->etamax;
  nitparam_.thmin        = nit->thmin;
  nitparam_.thmax        = nit->thmax;
  nitparam_.etafixed     = nit->etafixed;

  //NITSOL_parse_input (nit, stdout);
  nitsol_ (&(nit->n),
	   x,
	   nit->f,
	   nit->jacv,
	   &(nit->ftol),
	   &(nit->stptol),
	   nit->input,
	   nit->info,
	   nit->rwork,
	   nit->rpar,
	   nit->ipar,
	   &(nit->iterm),
	   nit->dinpr,
	   nit->dnorm);

  if (nit->iterm != 0)
    {
      NITSOL_parse_iterm (nit, stdout);
      exit (1);
    }

}
