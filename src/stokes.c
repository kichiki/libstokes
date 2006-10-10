/* structure for system parameters of stokes library.
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes.c,v 2.6 2006/10/10 17:14:55 ichiki Exp $
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
#include <stdlib.h>
#include <math.h> /* log() */
#include <libiter.h>

#include "stokes.h"


/* table for lattice summation
 * Note: now ewald-table is default (always constructed) so that
 * those are called only from the top level routine stokes_init().
 */
static void
stokes_ewald_table_free (struct stokes * sys)
{
  if (sys == NULL) return;

  if (sys->rlx != NULL) free (sys->rlx);
  if (sys->rly != NULL) free (sys->rly);
  if (sys->rlz != NULL) free (sys->rlz);

  if (sys->ex != NULL) free (sys->ex);
  if (sys->ey != NULL) free (sys->ey);
  if (sys->ez != NULL) free (sys->ez);
  if (sys->k  != NULL) free (sys->k);
  if (sys->k1 != NULL) free (sys->k1);
  if (sys->k2 != NULL) free (sys->k2);
  if (sys->k3 != NULL) free (sys->k3);
  if (sys->ya != NULL) free (sys->ya);
  if (sys->yb != NULL) free (sys->yb);
  if (sys->yc != NULL) free (sys->yc);
  if (sys->yg != NULL) free (sys->yg);
  if (sys->yh != NULL) free (sys->yh);
  if (sys->ym != NULL) free (sys->ym);

  sys->rlx = NULL;
  sys->rly = NULL;
  sys->rlz = NULL;

  sys->ex = NULL;
  sys->ey = NULL;
  sys->ez = NULL;
  sys->k  = NULL;
  sys->k1 = NULL;
  sys->k2 = NULL;
  sys->k3 = NULL;
  sys->ya = NULL;
  sys->yb = NULL;
  sys->yc = NULL;
  sys->yg = NULL;
  sys->yh = NULL;
  sys->ym = NULL;

  sys->nr = 0;
  sys->nk = 0;
  sys->flag_table = 0;
}

static void
ewald_r_append (struct stokes * e,
		double rlx, double rly, double rlz)
{
  e->nr ++;

  e->rlx = (double *) realloc (e->rlx, sizeof (double) * e->nr);
  e->rly = (double *) realloc (e->rly, sizeof (double) * e->nr);
  e->rlz = (double *) realloc (e->rlz, sizeof (double) * e->nr);

  int i = e->nr - 1;
  e->rlx [i] = rlx;
  e->rly [i] = rly;
  e->rlz [i] = rlz;
}

static void
ewald_k_append (struct stokes * e,
		double ex, double ey, double ez,
		double k, double k1, double k2, double k3,
		double ya,
		double yb,
		double yc,
		double yg,
		double yh,
		double ym)
{
  e->nk ++;

  e->ex = (double *) realloc (e->ex, sizeof (double) * e->nk);
  e->ey = (double *) realloc (e->ey, sizeof (double) * e->nk);
  e->ez = (double *) realloc (e->ez, sizeof (double) * e->nk);

  e->k  = (double *) realloc (e->k,  sizeof (double) * e->nk);
  e->k1 = (double *) realloc (e->k1, sizeof (double) * e->nk);
  e->k2 = (double *) realloc (e->k2, sizeof (double) * e->nk);
  e->k3 = (double *) realloc (e->k3, sizeof (double) * e->nk);

  e->ya = (double *) realloc (e->ya, sizeof (double) * e->nk);
  e->yb = (double *) realloc (e->yb, sizeof (double) * e->nk);
  e->yc = (double *) realloc (e->yc, sizeof (double) * e->nk);
  e->yg = (double *) realloc (e->yg, sizeof (double) * e->nk);
  e->yh = (double *) realloc (e->yh, sizeof (double) * e->nk);
  e->ym = (double *) realloc (e->ym, sizeof (double) * e->nk);

  int i = e->nk - 1;
  e->ex [i] = ex;
  e->ey [i] = ey;
  e->ez [i] = ez;

  e->k  [i] = k;
  e->k1 [i] = k1;
  e->k2 [i] = k2;
  e->k3 [i] = k3;

  e->ya [i] = ya;
  e->yb [i] = yb;
  e->yc [i] = yc;
  e->yg [i] = yg;
  e->yh [i] = yh;
  e->ym [i] = ym;
}

/* make ewald-summation table
 * INPUT
 *  (struct stokes *) sys :
 *  cutlim   :
 * OUTPUT
 *  (struct stokes *) sys : 
 */
static void
stokes_ewald_table_make (struct stokes * sys, double cutlim)
{
  double ya; 
  double yb;
  double yc;
  double yg;
  double yh;
  double ym;

  double ex, ey, ez;

  double xx, yy, zz, rr;
  double rlx, rly, rlz;

  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double k;
  double kexp;


  if (sys->flag_table == 1)
    {
      stokes_ewald_table_free (sys);
    }

  double rl;
  double rmax;
  rmax = sqrt (sys->rmax2);
  rl = sqrt (sys->lx*sys->lx + sys->ly*sys->ly + sys->lz*sys->lz);
  /* first Ewald part ( real space ) */
  for (m1 = - sys->rmaxx; m1 <= sys->rmaxx; m1++)
    {
      rlx = sys->lx * (double) m1;
      xx = rlx; // this is for minimum distance
      if (m1 > 0) xx -= 0.5 * sys->lx;
      if (m1 < 0) xx += 0.5 * sys->lx;
      for (m2 = - sys->rmaxy; m2 <= sys->rmaxy; m2++)
	{
	  rly = sys->ly * (double) m2;
	  yy = rly; // this is for minimum distance
	  if (m2 > 0) yy -= 0.5 * sys->ly;
	  if (m2 < 0) yy += 0.5 * sys->ly;
	  for (m3 = - sys->rmaxz; m3 <= sys->rmaxz; m3++)
	    {
	      rlz = sys->lz * (double) m3;
	      zz = rlz; // this is for minimum distance
	      if (m3 > 0) zz -= 0.5 * sys->lz;
	      if (m3 < 0) zz += 0.5 * sys->lz;
  
	      rr = sqrt (xx * xx + yy * yy + zz * zz);
	      if (rr <= rmax + rl)
		{
		  ewald_r_append (sys,
				  rlx, rly, rlz);
		}
	    }
	}
    }

  /* Second Ewald part ( reciprocal space ) */
  for (m1 = - sys->kmaxx; m1 <= sys->kmaxx; m1++)
    {
      k1 = 2.0 * M_PI * (double) m1 / sys->lx;
      for (m2 = - sys->kmaxy; m2 <= sys->kmaxy; m2++)
	{
	  k2 = 2.0 * M_PI * (double) m2 / sys->ly;
	  for (m3 = - sys->kmaxz; m3 <= sys->kmaxz; m3++)
	    {
	      k3 = 2.0 * M_PI * (double) m3 / sys->lz;
	      if (m1 != 0 || m2 != 0 || m3 != 0)
		{
		  kk = k1 * k1 + k2 * k2 + k3 * k3;
		  k = sqrt (kk);
		  if (k <= sys->kmax)
		    {
		      k4z = kk / 4.0 / sys->xi2;
		      kexp = sys->pivol
			* (1.0 + k4z * (1.0 + 2.0 * k4z))
			/ kk * exp (- k4z);

		      ex = k1 / k;
		      ey = k2 / k;
		      ez = k3 / k;

		      ya = 6.0 * (1.0 - kk / 3.0) * kexp;
		      yb = 3.0 * k * kexp;
		      yc = 3.0 / 2.0 * kk * kexp;
		      yg = 3.0 * (1.0 - 4.0 / 15.0 * kk) * k * kexp;
		      yh = 3.0 / 2.0 * kk * kexp;
		      ym = 3.0 * (1.0 - kk / 5.0) * kk * kexp;

		      if (sys->version == 0) // F version
			{
			  ewald_k_append (sys,
					  ex, ey, ez,
					  k, k1, k2, k3,
					  ya,
					  0.0,
					  0.0,
					  0.0,
					  0.0,
					  0.0);
			}
		      else if (sys->version == 1) // FT version
			{
			  ewald_k_append (sys,
					  ex, ey, ez,
					  k, k1, k2, k3,
					  ya,
					  yb,
					  yc,
					  0.0,
					  0.0,
					  0.0);
			}
		      else
			{
			  ewald_k_append (sys,
					  ex, ey, ez,
					  k, k1, k2, k3,
					  ya,
					  yb,
					  yc,
					  yg,
					  yh,
					  ym);
			}
		    }
		}
	    }
	}
    }

  sys->flag_table = 1;
}



/* all elements are zero-cleared
 */
struct stokes *
stokes_init (void)
{
  struct stokes * sys = NULL;

  sys = (struct stokes *) malloc (sizeof (struct stokes));

  sys->np = 0;
  sys->nm = 0;
  sys->pos = NULL;

  sys->version = 0;

  /* # of cell in real space */
  sys->rmax2 = 0.0;
  sys->rmaxx = 0;
  sys->rmaxy = 0;
  sys->rmaxz = 0;

  /* # of cell in reciprocal space */
  sys->kmax = 0.0;
  sys->kmaxx = 0;
  sys->kmaxy = 0;
  sys->kmaxz = 0;

  sys->xi = 0.0;
  sys->xi2 = 0.0;
  sys->xiaspi = 0.0;
  sys->xia2 = 0.0;

  sys->pivol = 0.0;

  sys->lx = 0.0;
  sys->ly = 0.0;
  sys->lz = 0.0;

  int i;
  for (i = 0; i < 27; i ++)
    {
      sys->llx [i] = 0.0;
      sys->lly [i] = 0.0;
      sys->llz [i] = 0.0;
    }

  // self part
  sys->self_a = 0.0;
  sys->self_c = 0.0;
  sys->self_m = 0.0;

  // table for lattice summation
  sys->flag_table = 0;
  // real space
  sys->nr = 0;
  sys->rlx = NULL;
  sys->rly = NULL;
  sys->rlz = NULL;
  // reciprocal space
  sys->nk = 0;
  sys->ex = NULL;
  sys->ey = NULL;
  sys->ez = NULL;
  sys->k  = NULL;
  sys->k1 = NULL;
  sys->k2 = NULL;
  sys->k3 = NULL;
  sys->ya = NULL;
  sys->yb = NULL;
  sys->yc = NULL;
  sys->yg = NULL;
  sys->yh = NULL;
  sys->ym = NULL;

  /* for lubrication */
  sys->lubcut = 0.0;

  /* for xi program */
  sys->cpu1 = 0.0;
  sys->cpu2 = 0.0;
  sys->cpu3 = 0.0;

  /* for iterative solvers */
  sys->it = NULL;

  return (sys);
}

void
stokes_free (struct stokes * sys)
{
  if (sys != NULL)
    {
      if (sys->pos != NULL) free (sys->pos);

      if (sys->rlx != NULL) free (sys->rlx);
      if (sys->rly != NULL) free (sys->rly);
      if (sys->rlz != NULL) free (sys->rlz);

      if (sys->ex != NULL) free (sys->ex);
      if (sys->ey != NULL) free (sys->ey);
      if (sys->ez != NULL) free (sys->ez);
      if (sys->k  != NULL) free (sys->k);
      if (sys->k1 != NULL) free (sys->k1);
      if (sys->k2 != NULL) free (sys->k2);
      if (sys->k3 != NULL) free (sys->k3);
      if (sys->ya != NULL) free (sys->ya);
      if (sys->yb != NULL) free (sys->yb);
      if (sys->yc != NULL) free (sys->yc);
      if (sys->yg != NULL) free (sys->yg);
      if (sys->yh != NULL) free (sys->yh);
      if (sys->ym != NULL) free (sys->ym);

      if (sys->it != NULL) iter_free (sys->it);
      free (sys);
    }
}

/* set np and nm and allocate the memory for pos[np*3]
 */
void
stokes_set_np (struct stokes * sys,
	       int np, int nm)
{
  sys->np = np;
  sys->nm = nm;

  if (sys->pos != NULL) free (sys->pos);
  sys->pos = (double *) malloc (sizeof (double) * np * 3);
}


void
stokes_set_ll (struct stokes * sys,
	       double lx, double ly, double lz)
{
  sys->pivol = M_PI / lx / ly / lz;

  sys->lx = lx;
  sys->ly = ly;
  sys->lz = lz;

  /* 0 is always the primary cell */
  sys->llx [ 0] = 0.0; sys->lly [ 0] = 0.0; sys->llz [ 0] = 0.0;
  sys->llx [ 1] = -lx; sys->lly [ 1] = 0.0; sys->llz [ 1] = 0.0;
  sys->llx [ 2] = 0.0; sys->lly [ 2] = -ly; sys->llz [ 2] = 0.0;
  sys->llx [ 3] = +lx; sys->lly [ 3] = 0.0; sys->llz [ 3] = 0.0;
  sys->llx [ 4] = 0.0; sys->lly [ 4] = +ly; sys->llz [ 4] = 0.0;
  sys->llx [ 5] = -lx; sys->lly [ 5] = -ly; sys->llz [ 5] = 0.0;
  sys->llx [ 6] = -lx; sys->lly [ 6] = +ly; sys->llz [ 6] = 0.0;
  sys->llx [ 7] = +lx; sys->lly [ 7] = -ly; sys->llz [ 7] = 0.0;
  sys->llx [ 8] = +lx; sys->lly [ 8] = +ly; sys->llz [ 8] = 0.0;
  /* up to 8, the monolayer (xy plain)
   * note: y is the vertical direction for 2D monolayer case! */

  sys->llx [ 9] = 0.0; sys->lly [ 9] = 0.0; sys->llz [ 9] = - lz;
  sys->llx [10] = -lx; sys->lly [10] = 0.0; sys->llz [10] = - lz;
  sys->llx [11] = 0.0; sys->lly [11] = -ly; sys->llz [11] = - lz;
  sys->llx [12] = +lx; sys->lly [12] = 0.0; sys->llz [12] = - lz;
  sys->llx [13] = 0.0; sys->lly [13] = +ly; sys->llz [13] = - lz;
  sys->llx [14] = -lx; sys->lly [14] = -ly; sys->llz [14] = - lz;
  sys->llx [15] = -lx; sys->lly [15] = +ly; sys->llz [15] = - lz;
  sys->llx [16] = +lx; sys->lly [16] = -ly; sys->llz [16] = - lz;
  sys->llx [17] = +lx; sys->lly [17] = +ly; sys->llz [17] = - lz;

  sys->llx [18] = 0.0; sys->lly [18] = 0.0; sys->llz [18] = + lz;
  sys->llx [19] = -lx; sys->lly [19] = 0.0; sys->llz [19] = + lz;
  sys->llx [20] = 0.0; sys->lly [20] = -ly; sys->llz [20] = + lz;
  sys->llx [21] = +lx; sys->lly [21] = 0.0; sys->llz [21] = + lz;
  sys->llx [22] = 0.0; sys->lly [22] = +ly; sys->llz [22] = + lz;
  sys->llx [23] = -lx; sys->lly [23] = -ly; sys->llz [23] = + lz;
  sys->llx [24] = -lx; sys->lly [24] = +ly; sys->llz [24] = + lz;
  sys->llx [25] = +lx; sys->lly [25] = -ly; sys->llz [25] = + lz;
  sys->llx [26] = +lx; sys->lly [26] = +ly; sys->llz [26] = + lz;
}

/** copied from PBC/ewald-3.c 1.6 2001/02/06 03:47:16 ichiki Exp **/
static double
ewald_real (double r, double xi)
{
  double xir, xir2;
  double r2;

  /* the safety factor is '1.15' */
  xir = xi * r / 1.15;
  xir2 = xir * xir;
  r2 = r * r;
  return (
	  erfc (xir)
	  * (0.75
	     +
	     1.5 / r2) / r
	  +
	  exp (- xir2)
	  * xi / sqrt (M_PI)
	  * (4.5
	     +
	     3.0 * xir2
	     +
	     (3.0
	      + xir2
	      * (14.0 + xir2
		 * (4.0 + xir2)))
	     / r2)
	  );
}

static double
ewald_recip (double k, double xi)
{
  double kx, kx2;
  double k2;

  /* the safety factor is '1.15' */
  kx = k / xi / 1.15;
  kx2 = kx * kx;
  k2 = k * k;
  return (
	  6.0 * M_PI
	  * exp (- kx2 / 4.0)
	  * (1.0 + k2 / 3.0) / k2
	  * (1.0 + kx2 * (1.0 / 4.0 + kx2 / 8.0))
	  );
}

void
stokes_set_xi (struct stokes * sys,
	       double xi, double cutlim)
{
  sys->xi  = xi;
  sys->xi2 = xi * xi;
  sys->xia2   = sys->xi2;
  sys->xiaspi = xi / sqrt (M_PI);

  /* diagonal part ( self part ) */
  sys->self_a = 1.0  - sys->xiaspi * (6.0 - 40.0 / 3.0 * sys->xia2);
  sys->self_c = 0.75 - sys->xiaspi * sys->xia2 * 10.0;
  sys->self_m = 0.9  - sys->xiaspi * sys->xia2 * (12.0 - 30.24 * sys->xia2);


  // set lattice for ewald summation
  double rmax;
  double kmax;
  double red_factor = 0.99;
  double dr, dk;
  double lmin, lmax;
  /* l is minimum among lx, ly, lz */
  lmin = sys->lx;
  if (lmin > sys->ly) lmin = sys->ly;
  if (lmin > sys->lz) lmin = sys->lz;
  dr = lmin;
  lmax = sys->lx;
  if (lmax < sys->ly) lmax = sys->ly;
  if (lmax < sys->lz) lmax = sys->lz;
  dk = 2.0 * M_PI / lmax;
  /* in real space */
  rmax = 10.0 * sqrt (- log (cutlim)) / sys->xi;
  for (;
       (ewald_real (rmax * red_factor, sys->xi)
	+ 6.0 * ewald_real (rmax * red_factor + dr, sys->xi))
	 < cutlim;
       rmax *= red_factor);
  sys->rmax2 = rmax * rmax;

  sys->rmaxx = (int) (rmax / sys->lx) + 1;
  sys->rmaxy = (int) (rmax / sys->ly) + 1;
  sys->rmaxz = (int) (rmax / sys->lz) + 1;

  /* in reciprocal space */
  kmax = 10.0 * sqrt (- log (cutlim)) * 2.0 * sys->xi;
  for (;
       (ewald_recip (kmax * red_factor, sys->xi)
	+ 6.0 * ewald_recip (kmax * red_factor + dk, sys->xi))
	 < cutlim;
       kmax *= red_factor);
  sys->kmax = kmax;

  sys->kmaxx = (int) (kmax * sys->lx / 2.0 / M_PI);
  sys->kmaxy = (int) (kmax * sys->ly / 2.0 / M_PI);
  sys->kmaxz = (int) (kmax * sys->lz / 2.0 / M_PI);

  /* now ewald summation lattice table is default */
  stokes_ewald_table_make (sys, cutlim);
}


double
xi_by_tratio (struct stokes * sys,
	      double tratio)
{
  double xi;

  xi = pow (tratio, 1.0 / 6.0)
    * sqrt (M_PI)
    / pow (sys->lx * sys->ly * sys->lz, 1.0 / 3.0);

  return (xi);
}


/* set iter param
 * INPUT
 *   solver : string indicating the solver
 *            sta, sta2, gpb, otmk, or gmres (default)
 *   eps and log10_eps
 *   max (and restart)
 *   debug = 0 : no debug info
 *         = 1 : iteration numbs and residue
 *   out   : FILE * to output debug info.
 */
void
stokes_set_iter (struct stokes * sys,
		 const char * solver,
		 int max,
		 int restart,
		 double eps,
		 int debug,
		 FILE * out)
{
  sys->iter = iter_init (solver,
			 max, restart, eps,
			 debug, out);
}
