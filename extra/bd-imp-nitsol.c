/* implicit Brownian dynamics algorithms by NITSOL
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bd-imp-nitsol.c,v 1.2 2008/06/06 03:45:50 kichiki Exp $
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
#include <math.h>
#include "memory-check.h" // CHECK_MALLOC

#include <stokes.h> // stokes_set_pos()
#include <ode.h> // solve_mix_3all()
#include <ode-quaternion.h> // quaternion_dQdt()

#include <brownian.h> // struct BD_params
#include <bd-imp.h> // struct BD_imp

#include <nitsol.h>
#include "nitsol_c.h" // struct NITSOL
#include "bd-imp-nitsol.h"

// BLAS functions
double
ddot_(int* N, 
      double* X, int* incX, 
      double* Y, int* incY);
double
dnrm2_(int* N, 
       double* X, int* incX);


/* jacobian and preconditioner
 * at this time, nothing will be done for both
 */
void BD_imp_NITSOL_jacv (int *n, double *xcur, double *fcur,
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
      fprintf (stderr, "# NITSOL_jacv: invalid ijob %d\n", *ijob);
      return;
    }
}

/* a wrapper of the nonlinear equations for the semi-implicit algorithms
 * (both JGdP and siPC) for NITSOL
 * INTPUT
 *  xcur[n] : = (x[nm*3])          for F version
 *            = (x[nm*3], q[nm*4]) for FT and FTS versions
 *  rpar    : (struct BD_imp *)BDimp
 *  ipar    : not used.
 * OUTPUT
 *  itrmf   : 0 means success
 *  fcur[n] := x - x0
 *           - dt * (uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for JGdP, or
 *  fcur[n] := x - x0
 *           - (dt/2) * (U^pr
 *                       + uinf(x) + M(x0).(F^E + F^P(x) + F^B(x0))),
 *  for siPC, where U^pr is the predictor obtained by Euler scheme as 
 *  U^pr = uinf(x0) + M(x0).(F^E + F^P(x0) + F^B(x0)).
 */
void
BD_imp_NITSOL_f (int *n, double *xcur, double *fcur,
		 double *rpar, int *ipar, int *itrmf)
{
  // let's say that (double *)rpar corresponds to (struct BD_imp *)BDimp
  struct BD_imp *BDimp = (struct BD_imp *)rpar;

  BD_imp_f (BDimp, xcur, fcur);

  *itrmf = 0;
}


/* set gsl_vector
 * if q == NULL, q = (0,0,0,1) is set.
 */
static void
BD_imp_NITSOL_set_guess (const struct BD_imp *BDimp,
			 const double *x,
			 const double *q,
			 double *x_nitsol)
{
  int nm3 = BDimp->BD->sys->nm * 3;
  int i;
  for (i = 0; i < nm3; i ++)
    {
      x_nitsol[i] = x[i];
    }
  if (BDimp->BD->sys->version > 0)
    {
      if (q == NULL)
	{
	  for (i = 0; i < BDimp->BD->sys->nm; i ++)
	    {
	      x_nitsol[nm3 + i*4 + 0] = 0.0;
	      x_nitsol[nm3 + i*4 + 1] = 0.0;
	      x_nitsol[nm3 + i*4 + 2] = 0.0;
	      x_nitsol[nm3 + i*4 + 3] = 1.0;
	    }
	}
      else
	{
	  for (i = 0; i < BDimp->BD->sys->nm * 4; i ++)
	    {
	      x_nitsol[nm3 + i] = q[i];
	    }
	}
    }
}

static void
BD_imp_NITSOL_get_root (const struct BD_imp *BDimp,
			const double *x_nitsol,
			double *x,
			double *q)
{
  int nm3 = BDimp->BD->sys->nm * 3;
  int i;
  for (i = 0; i < nm3; i ++)
    {
      x[i] = x_nitsol[i];
    }
  if (q != NULL)
    {
      for (i = 0; i < BDimp->BD->sys->nm * 4; i ++)
	{
	  q[i] = x_nitsol[nm3 + i];
	}
    }
}

void
BD_imp_NITSOL_wrap (struct BD_imp *BDimp,
		    double *x, double *q)
{
  double *x_nitsol = (double *)malloc (sizeof (double) * BDimp->nit->n);
  CHECK_MALLOC (x_nitsol, "BD_imp_NITSOL_wrap");

  BD_imp_NITSOL_set_guess (BDimp, x, q, x_nitsol);

  NITSOL_solve (BDimp->nit, x_nitsol);

  BD_imp_NITSOL_get_root (BDimp, x_nitsol, x, q);

  free (x_nitsol);
}
