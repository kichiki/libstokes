/* bead-rod model
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bead-rod.c,v 1.3 2008/07/16 16:44:54 kichiki Exp $
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
#include "stokes.h"
#include "nitsol_c.h"

#include "f.h"         // matrix_f_atimes()
#include "non-ewald.h" // scalars_nonewald()

#include "bead-rod.h"


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
	      const int *ib)
{
  struct BeadRod *br = (struct BeadRod *)malloc (sizeof (struct BeadRod));
  CHECK_MALLOC (br, "BeadRod_init");

  br->sys = sys;

  br->nc = nc;

  int i;
  br->a2 = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (br->a2, "BeadRod_init");
  for (i = 0; i < nc; i ++)
    {
      br->a2[i] = a[i] * a[i];
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

  br->verbose = 0;

  br->scheme = 0; // default : Liu (1989)
  br->eps = 0.0;
  br->nit = NULL;
  br->it  = NULL;

  br->dt = 0.0;
  br->d1 = 0.0;
  br->d2 = 0.0;

  br->u  = (double *)malloc (sizeof (double) * nc * 3);
  br->uu = (double *)malloc (sizeof (double) * nc * 3);
  CHECK_MALLOC (br->u,  "BeadRod_init");
  CHECK_MALLOC (br->uu, "BeadRod_init");
  // zero clear
  for (i = 0; i < nc*3; i ++)
    {
      br->u [i] = 0.0;
      br->uu[i] = 0.0;
    }

  return (br);
}

void
BeadRod_free (struct BeadRod *br)
{
  if (br == NULL) return;

  if (br->a2 != NULL) free (br->a2);
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
      /*
      br->it = iter_init ("gmres", 2000, 20,
			  eps,
			  br->nc,
			  NULL, // initial guess[]
			  1,    // use guess[]
			  br->verbose,
			  stdout);
      */
      br->it = iter_init ("otmk", 2000, 20,
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
 *  br->d1 = 2 dt / zeta
 *  br->d2 = (2 dt / zeta)^2
 */
void
BeadRod_set_coefs (struct BeadRod *br,
		   double dt,
		   double zeta)
{
  br->dt = dt;
  br->d1 = 2.0 * dt / zeta;
  br->d2 = br->d1 * br->d1;
}

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

      u[ix] = r[iax] - r[ibx];
      u[iy] = r[iay] - r[iby];
      u[iz] = r[iaz] - r[ibz];
    }
}

void
BeadRod_set_u_by_r (struct BeadRod *br,
		    const double *r)
{
  BeadRod_bead_to_connector (br->nc,
			     br->ia,
			     br->ib,
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
			     r,
			     br->uu);
}

/* return (i,j) element of the Rouse matrix
 * = 2 \delta_{i,j} - \delta_{ia[i], ib[j]} - \delta_{ib[i], ia[j]}
 * INPUT
 *  br : struct BeadRod
 *  i, j : constraint indices
 */
//static double
double
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

      f[i] = br->d2
	* (gx * gx
	   + gy * gy
	   + gz * gz)
	- 2.0 * br->d1 * lin
	+ (  br->uu[ix] * br->uu[ix]
	   + br->uu[iy] * br->uu[iy]
	   + br->uu[iz] * br->uu[iz])
	- br->a2[i];
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
      b[i] = 2.0 * br->d1 * 
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

      f[i] = br->d2 * (gx * gx + gy * gy + gz * gz)
	+ (  br->uu[ix] * br->uu[ix]
	   + br->uu[iy] * br->uu[iy]
	   + br->uu[iz] * br->uu[iz])
	- br->a2[i];
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
BeadRod_solve_gamma_iter (struct BeadRod *br,
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
      /*
      for (i = 0; i < nc; i ++)
	{
	  fprintf (stdout, "# iter %d : g0 = %e, g = %e\n",
		   i, g0[i], gamma[i]);
	}
      */

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


/* solve gamma for the constraints
 * either by linear approximation iteratively (Liu 1989)
 * or by Newton-GMRES method in NITSOL
 * INPUT
 *  br : struct BeadRod
 * OUTPUT
 *  gamma[nc] :
 */
void
BeadRod_solve_gamma (struct BeadRod *br,
		     double *gamma)
{
  if (br->scheme == 0)
    {
      BeadRod_solve_gamma_iter (br, gamma);
    }
  else
    {
      NITSOL_solve (br->nit, gamma);
    }
}




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
				 double *dr)
{
  // particle loop (nu)
  int i;
  for (i = 0; i < n; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      dr[ix] = a * dr[ix];
      dr[iy] = a * dr[iy];
      dr[iz] = a * dr[iz];

      // constraint loop
      int j;
      for (j = 0; j < br->nc; j ++)
	{
	  double cag;
	  // \overline{B}_{j,nu} = \delta_{ia[j],nu} - \delta_{ib[j],nu}
	  if (br->ia[j] == i)
	    {
	      cag = - br->d1 * gamma[j];
	      int jx = j * 3;
	      dr[ix] += cag * br->u[jx  ];
	      dr[iy] += cag * br->u[jx+1];
	      dr[iz] += cag * br->u[jx+2];
	    }
	  else if (br->ib[j] == i)
	    {
	      cag = - br->d1 * gamma[j];
	      int jx = j * 3;
	      dr[ix] -= cag * br->u[jx  ];
	      dr[iy] -= cag * br->u[jx+1];
	      dr[iz] -= cag * br->u[jx+2];
	    }
	}
    }
}



/**
 * HI version
 */
// borrowed from atimes_nonewald_3all() in non-ewald.c
/* calc U += M(i,j) . F
 * INPUT
 *  sys : struct stokes
 *  i, j : particle indices
 *  f[3] : force on particle j
 *  a    : prefactors
 * OUTPUT
 *  u[3] := a * M(i,j).f[]
 */
static void
mobility_F_atimes (struct stokes *sys,
		   int i, int j,
		   const double *f,
		   double a,
		   double *u)
{
  if (sys->periodic != 0)
    {
      fprintf (stderr, "# mobility_F_atimes: "
	       "periodic is not implemented yet\n");
      exit (1);
    }

  // used in the slip case
  double ai;
  double aj;
  if (sys->a == NULL)
    {
      ai = 1.0;
      aj = 1.0;
    }
  else
    {
      ai = sys->a [i];
      aj = sys->a [j];
    }


  // zero clear
  u[0] = 0;
  u[1] = 0;
  u[2] = 0;

  if (i == j)
    {
      // self
      double self_a = 1.0;
      if (sys->slip != NULL) // slip
	{
	  self_a = 1.0  * sys->slip_G32 [i]; // Lambda(3,2) = 1/Lambda(2,3)
	}
      matrix_f_atimes (f, u,
		       0.0, 0.0, 0.0,
		       self_a, self_a); // a part
    }
  else
    {
      // no-HI
      /*
      double mob [22];

      double xx = sys->pos [j*3  ] - sys->pos [i*3  ];
      double yy = sys->pos [j*3+1] - sys->pos [i*3+1];
      double zz = sys->pos [j*3+2] - sys->pos [i*3+2];
      double rr = xx * xx + yy * yy + zz * zz;

      double r = sqrt (rr);
      double rmin = (ai + aj) * sys->rmin;
      if (r < rmin) r = rmin;

      double ex = xx / r;
      double ey = yy / r;
      double ez = zz / r;

      if (sys->slip == NULL // no-slip
	  && sys->a == NULL) // monodisperse
	{
	  scalars_nonewald (sys->version, r, mob);
	}
      else if (sys->slip == NULL) // no-slip polydisperse
	{
	  scalars_nonewald_poly (0, // F version
				 r,
				 sys->a[i], sys->a[j],
				 mob);
	  scalars_mob_poly_scale_SD (0, // F version
				     sys->a[i],
				     mob);
	  // now mob is in the SD form
	}
      else // slip (for both mono- and poly-disperse systems)
	{
	  // use the effective radius in slip->a[] defined by
	  // slip->a[i] = sys->a[i] * sqrt (Lambda^(i)(0,2)).
	  // for both monodisperse and polydisperse.
	  // (slip->a[] is defined for both cases properly.)
	  scalars_nonewald_poly (0, // F version
				 r,
				 sys->slip_a[i], sys->slip_a[j],
				 mob);
	  // scale is done by the real radius
	  scalars_mob_poly_scale_SD (0, // F version
				     //sys->a[i],
				     ai,
				     mob);
	  // now scalars are in the SD form
	}

      // note that interaction (i,j) should be for (U[i], F[j])
      matrix_f_atimes (f, u,
		       ex, ey, ez,
		       mob[0], mob[1]); // xa, ya
      */
    }

  /* convert the resultant velocity into the Global scaling
   * and multiply the factor "a"
   */
  double c = a / ai;
  u[0] *= c;
  u[1] *= c;
  u[2] *= c;
}

/* calc \mathcal{M}_{i,j}\cdot\bm{u}_{j}
 * INPUT
 *  br     : struct BeadRod
 *  ci, cj : constraint indices
 *  u[3]   : connector vector for constraint j
 * OUTPUT
 *  mu[3]  : the resultant vector
 */
static void
BeadRod_mobility_atimes (struct BeadRod *br,
			 int ci, int cj,
			 const double *u,
			 double *mu)
{
  struct stokes *sys = br->sys;

  // zero clear
  mu[0] = 0.0;
  mu[1] = 0.0;
  mu[2] = 0.0;

  double mu_[3];
  mobility_F_atimes (sys, br->ia[ci], br->ia[cj], u, +1.0, mu_);
  mu[0] += mu_[0];
  mu[1] += mu_[1];
  mu[2] += mu_[2];

  mobility_F_atimes (sys, br->ia[ci], br->ib[cj], u, -1.0, mu_);
  mu[0] += mu_[0];
  mu[1] += mu_[1];
  mu[2] += mu_[2];

  mobility_F_atimes (sys, br->ib[ci], br->ia[cj], u, -1.0, mu_);
  mu[0] += mu_[0];
  mu[1] += mu_[1];
  mu[2] += mu_[2];

  mobility_F_atimes (sys, br->ib[ci], br->ib[cj], u, +1.0, mu_);
  mu[0] += mu_[0];
  mu[1] += mu_[1];
  mu[2] += mu_[2];
}
			 
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
		double *g)
{
  double c = - 2.0 * br->dt;

  int i;
  for (i = 0; i < br->nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;
      
      g[ix] = 0.0;
      g[iy] = 0.0;
      g[iz] = 0.0;

      int j;
      for (j = 0; j < br->nc; j ++)
	{
	  double mu[3];
	  BeadRod_mobility_atimes (br,
				   i, j,
				   br->u + j * 3,
				   mu);
	  g[ix] += gamma[j] * mu[0];
	  g[iy] += gamma[j] * mu[1];
	  g[iz] += gamma[j] * mu[2];
	} 
      // c = - 2 dt
      g[ix] *= c;
      g[iy] *= c;
      g[iz] *= c;
    }
}

/* calc \widetilde{\mathcal{M}}_{\nu,j}\cdot\bm{u}_{j}
 * INPUT
 *  br     : struct BeadRod
 *  i  : particle index (\nu)
 *  cj : constraint index (j)
 *  u[3]   : connector vector for constraint j
 * OUTPUT
 *  mu[3]  : the resultant vector
 */
static void
BeadRod_tilde_mobility_atimes (struct BeadRod *br,
			       int i, int cj,
			       const double *u,
			       double *mu)
{
  struct stokes *sys = br->sys;

  // zero clear
  mu[0] = 0.0;
  mu[1] = 0.0;
  mu[2] = 0.0;

  double mu_[3];
  mobility_F_atimes (sys, i, br->ia[cj], u, +1.0, mu_);
  mu[0] += mu_[0];
  mu[1] += mu_[1];
  mu[2] += mu_[2];

  mobility_F_atimes (sys, i, br->ib[cj], u, -1.0, mu_);
  mu[0] += mu_[0];
  mu[1] += mu_[1];
  mu[2] += mu_[2];
}

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
		 double *dr)
{
  struct stokes *sys = br->sys;

  double c = - 2.0 * br->dt;

  // particle loop (nu)
  int i;
  for (i = 0; i < sys->nm; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;
      
      dr[ix] = a * dr[ix];
      dr[iy] = a * dr[iy];
      dr[iz] = a * dr[iz];

      double dri[3] = {0.0, 0.0, 0.0};
      // constraint loop
      int j;
      for (j = 0; j < br->nc; j ++)
	{
	  double mu[3];
	  BeadRod_tilde_mobility_atimes (br,
					 i, // particle index
					 j, // constraint index
					 br->u + j * 3,
					 mu);
	  dri[0] += gamma[j] * mu[0];
	  dri[1] += gamma[j] * mu[1];
	  dri[2] += gamma[j] * mu[2];
	} 
      // c = - 2 dt
      dri[0] *= c;
      dri[1] *= c;
      dri[2] *= c;

      dr[ix] += dri[0];
      dr[iy] += dri[1];
      dr[iz] += dri[2];
    }
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
BeadRod_NLEQ_for_gamma_ (struct BeadRod *br,
			 const double *gamma,
			 double *f)
{
  double *g = (double *)malloc (sizeof (double) * br->nc * 3);
  CHECK_MALLOC (g, "BeadRod_NLEQ_for_gamma_");

  BeadRod_calc_g (br, gamma, g);

  int i;
  for (i = 0; i < br->nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      f[i]
	= (  g[ix] * g[ix]
	   + g[iy] * g[iy]
	   + g[iz] * g[iz])
	+ 2.0 * (  br->uu[ix] * g[ix]
		 + br->uu[iy] * g[iy]
		 + br->uu[iz] * g[iz])
	+ (  br->uu[ix] * br->uu[ix]
	   + br->uu[iy] * br->uu[iy]
	   + br->uu[iz] * br->uu[iz])
	- br->a2[i];
    }

  free (g);
}

/*
 * INPUT
 *  x[n] : gamma[nc]
 */
static void
BeadRod_atimes_ (int n, const double *x,
		 double *b, void *param)
{
  struct BeadRod *br = (struct BeadRod *)param;
  if (n != br->nc)
    {
      fprintf (stderr, "BeadRod_atimes: nc mismatches %d != %d\n",
	       n, br->nc);
      exit (1);
    }

  double *g = (double *)malloc (sizeof (double) * br->nc * 3);
  CHECK_MALLOC (g, "BeadRod_atimes_");

  BeadRod_calc_g (br, x, g);

  int i;
  for (i = 0; i < br->nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // linear term
      b[i] =
	- 2.0 * (  br->uu[ix] * g[ix]
		 + br->uu[iy] * g[iy]
		 + br->uu[iz] * g[iz]);
    }

  free (g);
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
BeadRod_update_gamma_ (struct BeadRod *br,
		       const double *gamma0,
		       double *gamma)
{
  int nc = br->nc;

  double *f = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (f, "BeadRod_update_gamma");

  double *g = (double *)malloc (sizeof (double) * br->nc * 3);
  CHECK_MALLOC (g, "BeadRod_update_gamma_");

  BeadRod_calc_g (br, gamma, g);

  int i;
  for (i = 0; i < nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // without the linear term
      f[i]
	= (  g[ix] * g[ix]
	   + g[iy] * g[iy]
	   + g[iz] * g[iz])
	+ (  br->uu[ix] * br->uu[ix]
	   + br->uu[iy] * br->uu[iy]
	   + br->uu[iz] * br->uu[iz])
	- br->a2[i];
    }

  // solve the linear set of equations for gamma
  solve_iter (nc, f, gamma,
	      BeadRod_atimes_,
	      (void *)br,
	      br->it);

  free (f);
  free (g);
}

/* solve gamma iteratively : Liu (1989)
 * INPUT
 *  br : struct BeadRod
 * OUTPUT
 *  gamma[nc] :
 */
void
BeadRod_solve_gamma_iter_ (struct BeadRod *br,
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
      BeadRod_update_gamma_ (br, g0, gamma);
      /*
      for (i = 0; i < nc; i ++)
	{
	  fprintf (stdout, "# iter %d : g0 = %e, g = %e\n",
		   i, g0[i], gamma[i]);
	}
      */

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


/**
 * top-level routine
 */
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
			       double *rr)
{
  BeadRod_set_u_by_r  (br, r);
  BeadRod_set_uu_by_r (br, rr);

  double *gamma = (double *)malloc (sizeof (double) * br->nc);
  CHECK_MALLOC (gamma, "BeadRod_adjust_by_constraints");
  BeadRod_solve_gamma (br, gamma);
  //BeadRod_solve_gamma_iter (br, gamma);
  //BeadRod_solve_gamma_iter_ (br, gamma); // linear scheme with HI version

  free (gamma);

  //BeadRod_constraint_displacement (br, gamma, n, 1.0, rr);
  BeadRod_calc_dr (br, gamma, 1.0, rr); // HI version

  // rr = rr + dr is set.
}


