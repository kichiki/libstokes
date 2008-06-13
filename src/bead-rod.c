/* bead-rod model
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bead-rod.c,v 1.2 2008/06/13 05:09:43 kichiki Exp $
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
#include <math.h> // sqrt
#include "memory-check.h" // macro CHECK_MALLOC

#include <libiter.h>
#include <nitsol_c.h>
#include "bead-rod.h"


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
	      const int *ib)
{
  struct BeadRod *br = (struct BeadRod *)malloc (sizeof (struct BeadRod));
  CHECK_MALLOC (br, "BeadRod_init");

  br->nc = nc;

  br->c1 = 0.0;
  br->c2 = 0.0;

  int i;
  if (a == NULL)
    {
      br->a  = NULL; // default : all the rod distances are 1
    }
  else
    {
      br->a = (double *)malloc (sizeof (double) * nc);
      CHECK_MALLOC (br->a, "BeadRod_init");
      for (i = 0; i < nc; i ++)
	{
	  br->a[i] = a[i];
	}
    }

  br->ia = (int *)malloc (sizeof (int) * nc);
  br->ib = (int *)malloc (sizeof (int) * nc);
  CHECK_MALLOC (br->ia, "BeadRod_init");
  CHECK_MALLOC (br->ib, "BeadRod_init");
  for (i = 0; i < nc; i ++)
    {
      br->ia[i] = ia[i];
      br->ib[i] = ib[i];
    }

  br->u  = (double *)malloc (sizeof (double) * nc * 3);
  br->uu = (double *)malloc (sizeof (double) * nc * 3);
  CHECK_MALLOC (br->u,  "BeadRod_init");
  CHECK_MALLOC (br->uu, "BeadRod_init");
  // zero clear
  for (i = 0; i < nc; i ++)
    {
      br->u[i]  = 0.0;
      br->uu[i] = 0.0;
    }


  br->verbose = 0;

  br->scheme = 0; // default : Liu (1989)
  br->nit = NULL;
  br->it  = NULL;

  return (br);
}

void
BeadRod_free (struct BeadRod *br)
{
  if (br == NULL) return;

  if (br->a  != NULL) free (br->a);
  if (br->ia != NULL) free (br->ia);
  if (br->ib != NULL) free (br->ib);

  if (br->u  != NULL) free (br->u);
  if (br->uu != NULL) free (br->uu);

  if (br->nit!= NULL) NITSOL_free (br->nit);
  if (br->it != NULL) iter_free (br->it);
  free (br);
}

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
		    double eps)
{
  br->scheme = scheme;
  br->eps    = eps;

  if (scheme == 0)
    {
      if (br->nit != NULL)
	{
	  NITSOL_free (br->nit);
	  br->nit = NULL;
	}

      br->it = iter_init ("gmres", 2000, 20,
			  eps,
			  br->nc,
			  NULL, // initial guess[]
			  1,    // use guess[]
			  br->verbose,
			  stdout);
      CHECK_MALLOC (br->it, "BeadRod_set_scheme");
    }
  else
    {
      if (br->it != NULL)
	{
	  iter_free (br->it);
	  br->it = NULL;
	}

      br->nit = NITSOL_init();
      CHECK_MALLOC (br->nit, "BeadRod_set_scheme");

      NITSOL_set_n (br->nit, br->nc);
      NITSOL_set_GMRES (br->nit, 50);
      NITSOL_set_jacv (br->nit,
		       0,   // p_flag
		       0,   // j_flag
		       4,   // j_order
		       NULL); // do nothing for jacv
      br->nit->f = BeadRod_NITSOL_f;
      NITSOL_set_forcing (br->nit, 0, 0.0, 0.0);

      NITSOL_set_iplvl (br->nit,
			0,  // iplvl
			6); // stdout

      NITSOL_set_norm_by_BLAS (br->nit);

      NITSOL_set_tol (br->nit,
		      eps, // ftol
		      eps);// stptol

      // parameter for f() and jacv()
      br->nit->rpar = (double *)br;
    }
}


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
		   double zeta)
{
  double d = 2.0 * dt / zeta;
  br->c1 = 2.0 * d;
  br->c2 = d * d;
}

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
			   double *u)
{
  int i;
  for (i = 0; i < nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      int iax = ia[i] * 3;
      int iay = iax + 1;
      int iaz = iax + 2;

      int ibx = ib[i] * 3;
      int iby = ibx + 1;
      int ibz = ibx + 2;

      if (a == NULL)
	{
	  u[ix] = r[iax] - r[ibx];
	  u[iy] = r[iay] - r[iby];
	  u[iz] = r[iaz] - r[ibz];
	}
      else
	{
	  u[ix] = (r[iax] - r[ibx]) / a[i];
	  u[iy] = (r[iay] - r[iby]) / a[i];
	  u[iz] = (r[iaz] - r[ibz]) / a[i];
	}
    }
}

void
BeadRod_set_u_by_r (struct BeadRod *br,
		    const double *r)
{
  BeadRod_bead_to_connector (br->nc,
			     br->ia,
			     br->ib,
			     br->a,
			     r,
			     br->u);
}

void
BeadRod_set_uu_by_r (struct BeadRod *br,
		     const double *r)
{
  BeadRod_bead_to_connector (br->nc,
			     br->ia,
			     br->ib,
			     br->a,
			     r,
			     br->uu);
}

/* return (i,j) element of the Rouse matrix
 * = 2 \delta_{i,j} - \delta_{ia[i], ib[j]} - \delta_{ib[i], ia[j]}
 * INPUT
 *  br : struct BeadRod
 *  i, j : constraint indices
 */
static double
BeadRod_Rouse_matrix (struct BeadRod *br,
		      int i, int j)
{
  double A = 0.0;
  if (i < 0 || i >= br->nc ||
      j < 0 || j >= br->nc) return (A);

  if (i == j) A += 2.0;
  if (br->ia[i] == br->ib[j]) A -= 1.0;
  if (br->ib[i] == br->ia[j]) A -= 1.0;

  return (A);
}

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
BeadRod_NLEQ_for_gamma (int nc, const double *gamma, double *f, void *params)
{
  struct BeadRod *br = (struct BeadRod *)params;
  if (nc != br->nc)
    {
      fprintf (stderr, "BeadRod_NLEQ_for_gamma: nc mismatches %d != %d\n",
	       nc, br->nc);
      exit (1);
    }

  int i;
  for (i = 0; i < nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // vector g[i] := gamma_j A_{ij} u_j
      double gx = 0.0;
      double gy = 0.0;
      double gz = 0.0;

      int j;
      for (j = 0; j < nc; j ++)
	{
	  double A = BeadRod_Rouse_matrix (br, i, j);
	  if (A != 0.0)
	    {
	      int jx = j * 3;
	      gx += gamma[j] * A * br->u[jx+0];
	      gy += gamma[j] * A * br->u[jx+1];
	      gz += gamma[j] * A * br->u[jx+2];
	    }
	}

      // linear term
      double lin
	= br->uu[ix] * gx
	+ br->uu[iy] * gy
	+ br->uu[iz] * gz;

      f[i] = br->c2 * (  gx * gx
		       + gy * gy
		       + gz * gz)
	- br->c1 * lin
	+ (  br->uu[ix] * br->uu[ix]
	   + br->uu[iy] * br->uu[iy]
	   + br->uu[iz] * br->uu[iz])
	- 1.0;
    }
}

static void
BeadRod_atimes (int n, const double *x,
		double *b, void *param)
{
  struct BeadRod *br = (struct BeadRod *)param;
  if (n != br->nc)
    {
      fprintf (stderr, "BeadRod_atimes: nc mismatches %d != %d\n",
	       n, br->nc);
      exit (1);
    }

  /* C_{ij} gamma_j = f[i]
   * where C_{ij} = c1 * A_{ij} u'_i u_j
   * NOTE: A_{ij} is tri-diagonal matrix
   *       A = [ 2 -1  0 ...   ]
   *           [-1  2 -1 ...   ]
   *           [ 0 -1  2 -1 ...]
   *           [...............]
   *           [    0 -1  2  -1]
   *           [       0 -1   2]
   */
  int i;
  for (i = 0; i < br->nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // quatratic term
      double gx = 0.0;
      double gy = 0.0;
      double gz = 0.0;

      int j;
      for (j = 0; j < br->nc; j ++)
	{
	  double A = BeadRod_Rouse_matrix (br, i, j);
	  if (A != 0.0)
	    {
	      int jx = j * 3;
	      gx += x[j] * A * br->u[jx+0];
	      gy += x[j] * A * br->u[jx+1];
	      gz += x[j] * A * br->u[jx+2];
	    }
	}

      // linear term
      b[i] = br->c1 * 
	(  br->uu[ix] * gx
	 + br->uu[iy] * gy
	 + br->uu[iz] * gz);
    }
}


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
		      double *gamma)
{
  int nc = br->nc;

  double *f = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (f, "BeadRod_update_gamma");

  int i;
  for (i = 0; i < nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // quatratic term
      double gx = 0.0;
      double gy = 0.0;
      double gz = 0.0;

      int j;
      for (j = 0; j < nc; j ++)
	{
	  double A = BeadRod_Rouse_matrix (br, i, j);
	  if (A != 0.0)
	    {
	      int jx = j * 3;
	      gx += gamma0[j] * A * br->u[jx+0];
	      gy += gamma0[j] * A * br->u[jx+1];
	      gz += gamma0[j] * A * br->u[jx+2];
	    }
	}

      f[i] = br->c2 * (gx * gx + gy * gy + gz * gz)
	+ (  br->uu[ix] * br->uu[ix]
	   + br->uu[iy] * br->uu[iy]
	   + br->uu[iz] * br->uu[iz])
	- 1.0;
    }

  // solve the linear set of equations for gamma
  solve_iter (nc, f, gamma,
	      BeadRod_atimes,
	      (void *)br,
	      br->it);

  free (f);
}

/* solve gamma iteratively : Liu (1989)
 * INPUT
 *  br : struct BeadRod
 * OUTPUT
 *  gamma[nc] :
 */
void
BeadRod_solve_iter_gamma (struct BeadRod *br,
			  double *gamma)
{
  double eps2 = br->eps * br->eps;

  int nc = br->nc; // number of rods (constraints)
  double *g0 = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (g0, "BeadRod_solve_iter_gamma");

  // initial condition
  int i;
  for (i = 0; i < nc; i ++)
    {
      g0[i] = 0.0;
    }

  double err2;
  int iter = 0;
  do
    {
      BeadRod_update_gamma (br, g0, gamma);

      // estimate error
      err2 = 0.0;
      for (i = 0; i < nc; i ++)
	{
	  double d = gamma[i] - g0[i];
	  err2 += d * d;

	  // for the next update
	  g0[i] = gamma[i];
	}
      err2 /= (double)nc;

      if (br->verbose != 0)
	{
	  fprintf (stdout, " BeadRod_solve_iter_gamma:"
		   " iter = %d, err^2 = %e\n", iter, err2);
	}
      iter ++;
    }
  while (err2 > eps2);

  free (g0);
}


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
				 double *dr)
{
  double c = - 0.5 * br->c1; // = - 2 dt / zeta

  // particle loop (nu)
  int i;
  for (i = 0; i < n; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;
      dr[ix] = 0.0;
      dr[iy] = 0.0;
      dr[iz] = 0.0;

      int j;
      for (j = 0; j < br->nc; j ++)
	{
	  double cag;
	  // \overline{B}_{j,nu} = \delta_{ia[j],nu} - \delta_{ib[j],nu}
	  if (br->ia[j] == i)
	    {
	      if (br->a == NULL)
		{
		  cag = c * gamma[j];
		}
	      else
		{
		  cag = c * br->a[j] * gamma[j];
		}
	      int jx = j * 3;
	      dr[ix] += cag * br->u[jx  ];
	      dr[iy] += cag * br->u[jx+1];
	      dr[iz] += cag * br->u[jx+2];
	    }
	  else if (br->ib[j] == i)
	    {
	      if (br->a == NULL)
		{
		  cag = c * gamma[j];
		}
	      else
		{
		  cag = c * br->a[j] * gamma[j];
		}
	      int jx = j * 3;
	      dr[ix] -= cag * br->u[jx  ];
	      dr[iy] -= cag * br->u[jx+1];
	      dr[iz] -= cag * br->u[jx+2];
	    }
	}
    }
}


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
		  double *rpar, int *ipar, int *itrmf)
{
  // let's say that (double *)rpar corresponds to (struct BeadRod *)params
  void *params = (void *)rpar;

  BeadRod_NLEQ_for_gamma (*n, xcur, fcur, params);

  *itrmf = 0;
}

/* solve gamma by NITSOL
 * INPUT
 *  br  : struct BeadRod
 * OUTPUT
 *  gamma[nc] :
 */
void
BeadRod_solve_gamma_by_NITSOL (struct BeadRod *br,
			       double *gamma)
{
  NITSOL_solve (br->nit, gamma);
}

/**
 * top level non-linear equasiont solver for gamma
 */
void
BeadRod_solve_gamma (struct BeadRod *br,
		     double *gamma)
{
  if (br->scheme == 0)
    {
      BeadRod_solve_iter_gamma (br, gamma);
    }
  else
    {
      NITSOL_solve (br->nit, gamma);
    }
}


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
			       double *r)
{
  //return;
  BeadRod_set_u_by_r  (br, r);
  BeadRod_set_uu_by_r (br, rr);

  double *gamma = (double *)malloc (sizeof (double) * br->nc);
  CHECK_MALLOC (gamma, "BeadRod_adjust_by_constraints");
  BeadRod_solve_gamma (br, gamma);

  BeadRod_constraint_displacement (br, gamma, n, r);
  free (gamma);

  int i;
  for (i = 0; i < n*3; i ++)
    {
      r[i] += rr[i];
    }
}
