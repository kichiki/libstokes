/* test code for bead-rod.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bead-rod.c,v 1.3 2008/11/05 02:29:00 kichiki Exp $
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
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include "memory-check.h" // macro CHECK_MALLOC

#include "check.h" // compare_max()

#include <nitsol_c.h> // NITSOL_set_iplvl()

#include <bead-rod.h>
#include "check-bead-rod.h"

// in bead-rod.c
double
BeadRod_Rouse_matrix (struct BeadRod *br,
		      int i, int j);


/* set particle list for constraint ia[] and ib[] for straight chain
 * INPUT
 *  br : struct BeadRod
 * OUTPUT
 *  br->ia[]
 *  br->ib[]
 */
void
BeadRod_set_connections_straight (struct BeadRod *br)
{
  if (br->ia != NULL) free (br->ia);
  if (br->ib != NULL) free (br->ib);

  br->ia = (int *)malloc (sizeof (int) * br->nc);
  br->ib = (int *)malloc (sizeof (int) * br->nc);
  CHECK_MALLOC (br->ia, "BeadRod_set_connections_straight");
  CHECK_MALLOC (br->ib, "BeadRod_set_connections_straight");

  int i;
  for (i = 0; i < br->nc; i ++)
    {
      br->ia[i] = i+1;
      br->ib[i] = i;
    }
}

struct BeadRod *
BeadRod_init_for_test (int nc,
		       double dt, double zeta,
		       double dr)
{
  int n = nc + 1; // straight config is assumed.
  double *r0 = (double *)malloc (sizeof (double) * n * 3);
  double *rr = (double *)malloc (sizeof (double) * n * 3);
  CHECK_MALLOC (r0, "BeadRod_init_for_test");
  CHECK_MALLOC (rr, "BeadRod_init_for_test");

  // radius of big circle
  double R = (double)nc / (2.0 * M_PI);
  double dx = 2.0 * asin (0.5 / R);
  int i;
  for (i = 0; i < n; i ++)
    {
      double x = (double)i * dx;
      r0[i*3+0] = R * cos (x);
      r0[i*3+1] = R * sin (x);
      r0[i*3+2] = 0.0;
    }

  // adjust the center of mass
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  for (i = 0; i < n; i ++)
    {
      cx += r0[i*3+0];
      cy += r0[i*3+1];
      cz += r0[i*3+2];
    }
  cx /= (double)n;
  cy /= (double)n;
  cz /= (double)n;
  for (i = 0; i < n; i ++)
    {
      r0[i*3+0] -= cx;
      r0[i*3+1] -= cy;
      r0[i*3+2] -= cz;
    }

  srand48(1);
  for (i = 0; i < (n * 3); i ++)
    {
      rr[i] = r0[i] + dr * 2.0 * (drand48() - 0.5);
      // deviation is in [-dr, +dr)
    }
  // now r0[] and rr[] are ready

  // straight conformation
  double *a = (double *)malloc (sizeof (double) * nc);
  int *ia = (int *)malloc (sizeof (int) * nc);
  int *ib = (int *)malloc (sizeof (int) * nc);
  CHECK_MALLOC (ia, "BeadRod_set_connections_straight");
  CHECK_MALLOC (ib, "BeadRod_set_connections_straight");
  for (i = 0; i < nc; i ++)
    {
      a[i] = 1.0;
      ia[i] = i+1;
      ib[i] = i;
    }

  struct BeadRod *br = BeadRod_init (NULL, // struct stokes
				     nc, a, ia, ib);
  CHECK_MALLOC (br, "BeadRod_init_for_test");
  free (a);
  free (ia);
  free (ib);

  BeadRod_set_connections_straight (br);
  BeadRod_set_coefs (br, dt, zeta);

  BeadRod_set_u_by_r (br, r0);
  BeadRod_set_uu_by_r (br, rr);

  free (r0);
  free (rr);

  return (br);
}


int
check_BeadRod_constraint_displacement (int n,
				       int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BeadRod_constraint_displacement (n=%d)\n", n);
    }

  int check = 0;
  double max = 0.0;

  int nc = n - 1;
  struct BeadRod *br
    = BeadRod_init_for_test (nc,
			     0.1, // dt
			     1.0, // zeta
			     0.1);// dr

  double *gamma = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (gamma, "check_BeadRod_constraint_displacement");
  int i;
  for (i = 0; i < nc; i ++)
    {
      gamma[i] = 1.0 + (double)i;
    }

  double *dr = (double *)calloc (n*3, sizeof (double));
  CHECK_MALLOC (dr, "check_BeadRod_constraint_displacement");
  BeadRod_constraint_displacement (br, gamma, n, 0.0, dr);

  double *du = (double *)malloc (sizeof (double) * nc * 3);
  CHECK_MALLOC (du, "check_BeadRod_constraint_displacement");
  BeadRod_bead_to_connector (nc,
			     br->ia, br->ib,
			     dr, du);


  char label[80];
  double c = - br->d1; // = - 2 dt / zeta
  for (i = 0; i < nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      int j;
      int jx;
      int jy;
      int jz;

      j = i;
      jx = j * 3;
      jy = jx + 1;
      jz = jx + 2;
      double du2x = c * 2.0 * gamma[j] * br->u[jx]; // A_{i,i}
      double du2y = c * 2.0 * gamma[j] * br->u[jy]; // A_{i,i}
      double du2z = c * 2.0 * gamma[j] * br->u[jz]; // A_{i,i}
      if (i > 0)
	{
	  j = i - 1;
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;
	  du2x += c * (-1.0) * gamma[j] * br->u[jx]; // A_{i,i-1}
	  du2y += c * (-1.0) * gamma[j] * br->u[jy]; // A_{i,i-1}
	  du2z += c * (-1.0) * gamma[j] * br->u[jz]; // A_{i,i-1}
	}
      if (i < (nc - 1))
	{
	  j = i + 1;
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;
	  du2x += c * (-1.0) * gamma[j] * br->u[jx]; // A_{i,i+1}
	  du2y += c * (-1.0) * gamma[j] * br->u[jy]; // A_{i,i+1}
	  du2z += c * (-1.0) * gamma[j] * br->u[jz]; // A_{i,i+1}
	}

      sprintf (label, " du_x[%d]", i);
      check += compare_max (du[ix], du2x, label,  verbose, tiny, &max);
      sprintf (label, " du_y[%d]", i);
      check += compare_max (du[iy], du2y, label,  verbose, tiny, &max);
      sprintf (label, " du_z[%d]", i);
      check += compare_max (du[iz], du2z, label,  verbose, tiny, &max);
    }

  BeadRod_free (br);
  free (gamma);
  free (dr);
  free (du);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

int
check_BeadRod_solve_iter_gamma (int n, double eps,
				int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BeadRod_solve_iter_gamma (n=%d, eps=%e)\n",
	       n, eps);
    }

  int check = 0;
  double max = 0.0;


  int nc = n - 1;
  struct BeadRod *br
    = BeadRod_init_for_test (nc,
			     0.1, // dt
			     1.0, // zeta
			     0.1);// dr
  br->verbose = verbose;

  double *gamma = (double *)calloc (sizeof (double), nc);
  CHECK_MALLOC (gamma, "check_BeadRod_solve_iter_gamma");

  BeadRod_set_scheme (br,
		      0, // Liu (1989)
		      eps);
  BeadRod_solve_gamma (br, gamma);

  double *f = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (f, "check_BeadRod_solve_iter_gamma");

  BeadRod_NLEQ_for_gamma (nc, gamma, f, (void *)br);

  char label[80];
  int i;
  for (i = 0; i < nc; i ++)
    {
      sprintf (label, " f[%d]", i);
      check += compare_max (1.0, f[i]+1.0, label,  verbose, tiny, &max);
    }


  BeadRod_free (br);
  free (gamma);
  free (f);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

int
check_BeadRod_solve_gamma_by_NITSOL (int n, double eps,
				     int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BeadRod_solve_gamma_by_NITSOL (n=%d, eps=%e)\n",
	       n, eps);
    }

  int check = 0;
  double max = 0.0;


  int nc = n - 1;
  struct BeadRod *br
    = BeadRod_init_for_test (nc,
			     0.1, // dt
			     1.0, // zeta
			     0.1);// dr
  BeadRod_set_scheme (br,
		      1, // NITSOL
		      eps);
  NITSOL_set_iplvl (br->nit,
		    1,  // iplvl
		    6); // stdout


  double *gamma = (double *)calloc (sizeof (double), nc);
  CHECK_MALLOC (gamma, "check_BeadRod_solve_iter_gamma");

  // NOTE: gamma[] should be initialized!
  BeadRod_solve_gamma (br, gamma);

  double *f = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (f, "check_BeadRod_solve_iter_gamma");

  BeadRod_NLEQ_for_gamma (nc, gamma, f, (void *)br);

  char label[80];
  int i;
  for (i = 0; i < nc; i ++)
    {
      sprintf (label, " f[%d]", i);
      check += compare_max (1.0, f[i]+1.0, label,  verbose, tiny, &max);
    }


  BeadRod_free (br);
  free (gamma);
  free (f);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


static double
norm (const double *r)
{
  return (sqrt (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));
}

static double
diff_norm (const double *r1, const double *r2)
{
  double d[3];
  d[0] = r1[0] - r2[0];
  d[1] = r1[1] - r2[1];
  d[2] = r1[2] - r2[2];
  return (norm (d));
}

int
check_BeadRod_solve_gamma (int n, double eps,
			   int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BeadRod_solve_gamma (n=%d, eps=%e)\n",
	       n, eps);
    }

  int check = 0;
  double max = 0.0;


  // make straight configuration and it's deviation
  int nc = n - 1;
  double *a = (double *)malloc (sizeof (double) * nc);
  int *ia = (int *)malloc (sizeof (int) * nc);
  int *ib = (int *)malloc (sizeof (int) * nc);
  CHECK_MALLOC (a,  "check_BeadRod_solve_gamma");
  CHECK_MALLOC (ia, "check_BeadRod_solve_gamma");
  CHECK_MALLOC (ib, "check_BeadRod_solve_gamma");

  srand48 (1);
  int i;
  for (i = 0; i < nc; i ++)
    {
      a[i] = 1.0 + 0.2 * (drand48() - 0.5);
      ia[i] = i+1;
      ib[i] = i;
    }

  double *r0 = (double *)malloc (sizeof (double) * n*3);
  double *rr = (double *)malloc (sizeof (double) * n*3);
  CHECK_MALLOC (r0,  "check_BeadRod_solve_gamma");
  CHECK_MALLOC (rr,  "check_BeadRod_solve_gamma");

  r0[0] = 0.0;
  r0[1] = 0.0;
  r0[2] = 0.0;

  rr[0] = 0.2 * (drand48() - 0.5);
  rr[1] = 0.2 * (drand48() - 0.5);
  rr[2] = 0.2 * (drand48() - 0.5);

  for (i = 0; i < nc; i ++)
    {
      r0[(i+1)*3+0] = r0[i*3+0] + a[i];
      r0[(i+1)*3+1] = 0.0;
      r0[(i+1)*3+2] = 0.0;

      rr[(i+1)*3+0] = r0[(i+1)*3+0] + 0.2 * (drand48() - 0.5);
      rr[(i+1)*3+1] = r0[(i+1)*3+1] + 0.2 * (drand48() - 0.5);
      rr[(i+1)*3+2] = r0[(i+1)*3+2] + 0.2 * (drand48() - 0.5);
    }

  struct BeadRod *br
    = BeadRod_init (NULL, // struct stokes
		    nc, a, ia, ib);
  BeadRod_set_scheme (br,
		      //1, // nitsol
		      0, // linear
		      eps);
  if (br->scheme == 1)
    {
      NITSOL_set_iplvl (br->nit,
			0,  // iplvl
			6); // stdout
    }

  double dt = 0.1;
  double zeta = 1.0;
  BeadRod_set_coefs (br, dt, zeta);

  BeadRod_set_u_by_r  (br, r0);
  BeadRod_set_uu_by_r (br, rr);

  double *gamma = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (gamma, "check_BeadRod_solve");
  BeadRod_solve_gamma (br, gamma);

  double *f = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (f, "check_BeadRod_solve");
  BeadRod_NLEQ_for_gamma (nc, gamma, f, (void *)br);

  double *dr = (double *)calloc (n*3, sizeof (double));
  double *r1 = (double *)malloc (sizeof (double) * n*3);
  CHECK_MALLOC (dr, "check_BeadRod_solve");
  CHECK_MALLOC (r1, "check_BeadRod_solve");
  BeadRod_constraint_displacement (br, gamma, n, 0.0, dr);
  for (i = 0; i < n*3; i ++)
    {
      r1[i] = rr[i] + dr[i];
    }

  /* now we have positions
   *  r0[n] : initial positions (satisfying the constraints)
   *  rr[n] : unconstraints positions
   *  r1[n] : corrections to rr[n] to satisfy the constraints
   * also we have f[nc], the nonlinear equations which should be 0
   */
  char label[80];
  //fprintf (stdout, "# i, d0, dd, d1, a\n");
  for (i = 0; i < nc; i ++)
    {
      //double d0 = diff_norm (r0 + i*3, r0 + (i+1)*3);
      //double dd = diff_norm (rr + i*3, rr + (i+1)*3);
      double d1 = diff_norm (r1 + i*3, r1 + (i+1)*3);
      /*
      fprintf (stdout, "%d : %e %e %e =? %e : %e\n",
	       i, d0, dd, d1, a[i], fabs (d1-a[i]));
      */
      sprintf (label, " |u[%d]|^2", i);
      check += compare_max (d1, a[i], label,  verbose, tiny, &max);
    }

  /**
   * the following is further detailed check
   * -- now everything is OK, so I comment out

  // let's calculate g1[nc*3] by br->u[nc*3] here
  // borrowed from BeadRod_NLEQ_for_gamma() in bead-rod.c
  double *g1 = (double *)malloc (sizeof (double) * nc*3);
  CHECK_MALLOC (g1, "check_BeadRod_solve");
  for (i = 0; i < nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // vector g[i] := gamma_j A_{ij} u_j
      g1[ix] = 0.0;
      g1[iy] = 0.0;
      g1[iz] = 0.0;

      int j;
      for (j = 0; j < nc; j ++)
	{
	  double A = BeadRod_Rouse_matrix (br, i, j);
	  if (A != 0.0)
	    {
	      int jx = j * 3;
	      g1[ix] += gamma[j] * A * br->u[jx+0];
	      g1[iy] += gamma[j] * A * br->u[jx+1];
	      g1[iz] += gamma[j] * A * br->u[jx+2];
	    }
	}
      g1[ix] *= - br->d1;
      g1[iy] *= - br->d1;
      g1[iz] *= - br->d1;
    }

  // let's recover g2[nc*3] from dr[n*3]
  double *g2 = (double *)malloc (sizeof (double) * nc*3);
  for (i = 0; i < nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      g2[ix] = dr[ix+3] - dr[ix];
      g2[iy] = dr[iy+3] - dr[iy];
      g2[iz] = dr[iz+3] - dr[iz];
    }

  fprintf (stdout, "# i, (u' + g)^2\n");
  for (i = 0; i < nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      double d[3];

      d[0] = br->uu[ix] + g1[ix];
      d[1] = br->uu[iy] + g1[iy];
      d[2] = br->uu[iz] + g1[iz];
      double d1 = norm (d);

      d[0] = br->uu[ix] + g2[ix];
      d[1] = br->uu[iy] + g2[iy];
      d[2] = br->uu[iz] + g2[iz];
      double d2 = norm (d);

      fprintf (stdout, "%d : %e %e\n",
	       i, d1, d2);
    }


  // borrowed from BeadRod_NLEQ_for_gamma() in bead-rod.c
  double *f0 = (double *)malloc (sizeof (double) * nc*3);
  double *f1 = (double *)malloc (sizeof (double) * nc*3);
  double *f2 = (double *)malloc (sizeof (double) * nc*3);
  double *g3 = (double *)malloc (sizeof (double) * nc*3);
  for (i = 0; i < nc; i ++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // vector g[i] := gamma_j A_{ij} u_j
      g3[ix] = 0.0;
      g3[iy] = 0.0;
      g3[iz] = 0.0;

      int j;
      for (j = 0; j < nc; j ++)
	{
	  double A = BeadRod_Rouse_matrix (br, i, j);
	  if (A != 0.0)
	    {
	      int jx = j * 3;
	      g3[ix] += gamma[j] * A * br->u[jx+0];
	      g3[iy] += gamma[j] * A * br->u[jx+1];
	      g3[iz] += gamma[j] * A * br->u[jx+2];
	    }
	}

      // linear term
      double lin
	= br->uu[ix] * g3[ix]
	+ br->uu[iy] * g3[iy]
	+ br->uu[iz] * g3[iz];

      f2[i] = br->d2 * (g3[ix] * g3[ix]
			+ g3[iy] * g3[iy]
			+ g3[iz] * g3[iz]);
      f1[i] =
	- 2.0 * br->d1 * lin;
      f0[i] = 
	+ (  br->uu[ix] * br->uu[ix]
	   + br->uu[iy] * br->uu[iy]
	   + br->uu[iz] * br->uu[iz])
	- 1.0;
    }
  fprintf (stdout, "# i, f0, f1, f2, f, F\n");
  for (i = 0; i < nc; i ++)
    {
      fprintf (stdout, "%d : %e %e %e => %e <=> %e\n",
	       i, f0[i], f1[i], f2[i], f0[i] + f1[i] + f2[i], f[i]);
    }
  fprintf (stdout, "# g1, g2, g3\n");
  for (i = 0; i < nc; i ++)
    {
      fprintf (stdout, "%d (x) %e %e %e\n",
	       i, g1[i*3], g2[i*3],
	       - br->d1 * g3[i*3]);
    }

  free (g1);
  free (g2);

  free (f0);
  free (f1);
  free (f2);
  free (g3);
  */

  free (r0);
  free (rr);

  free (a);
  free (ia);
  free (ib);

  BeadRod_free (br);

  free (gamma);
  free (f);
  free (dr);
  free (r1);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


int
check_BeadRod_calc_dr (int n, double eps,
		       int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BeadRod_calc_dr (n=%d, eps=%e)\n",
	       n, eps);
    }

  int check = 0;
  double max = 0.0;


  // make straight configuration and it's deviation
  int nc = n - 1;
  double *a = (double *)malloc (sizeof (double) * nc);
  int *ia = (int *)malloc (sizeof (int) * nc);
  int *ib = (int *)malloc (sizeof (int) * nc);
  CHECK_MALLOC (a,  "check_BeadRod_calc_dr");
  CHECK_MALLOC (ia, "check_BeadRod_calc_dr");
  CHECK_MALLOC (ib, "check_BeadRod_calc_dr");

  srand48 (1);
  int i;
  for (i = 0; i < nc; i ++)
    {
      a[i] = 1.0 + 0.2 * (drand48() - 0.5);
      ia[i] = i+1;
      ib[i] = i;
    }

  double *r0 = (double *)malloc (sizeof (double) * n*3);
  double *rr = (double *)malloc (sizeof (double) * n*3);
  CHECK_MALLOC (r0,  "check_BeadRod_calc_dr");
  CHECK_MALLOC (rr,  "check_BeadRod_calc_dr");

  r0[0] = 0.0;
  r0[1] = 0.0;
  r0[2] = 0.0;

  rr[0] = 0.2 * (drand48() - 0.5);
  rr[1] = 0.2 * (drand48() - 0.5);
  rr[2] = 0.2 * (drand48() - 0.5);

  for (i = 0; i < nc; i ++)
    {
      r0[(i+1)*3+0] = r0[i*3+0] + a[i];
      r0[(i+1)*3+1] = 0.0;
      r0[(i+1)*3+2] = 0.0;

      rr[(i+1)*3+0] = r0[(i+1)*3+0] + 0.2 * (drand48() - 0.5);
      rr[(i+1)*3+1] = r0[(i+1)*3+1] + 0.2 * (drand48() - 0.5);
      rr[(i+1)*3+2] = r0[(i+1)*3+2] + 0.2 * (drand48() - 0.5);
    }

  struct stokes *sys = stokes_init ();
  stokes_set_np (sys, n, n);
  stokes_set_pos (sys, r0);
  sys->version = 0;

  struct BeadRod *br
    = BeadRod_init (sys, nc, a, ia, ib);
  BeadRod_set_scheme (br,
		      //1, // nitsol
		      0, // linear
		      eps);
  if (br->scheme == 1)
    {
      NITSOL_set_iplvl (br->nit,
			0,  // iplvl
			6); // stdout
    }

  double dt = 0.1;
  double zeta = 1.0;
  BeadRod_set_coefs (br, dt, zeta);

  BeadRod_set_u_by_r  (br, r0);
  BeadRod_set_uu_by_r (br, rr);

  double *gamma = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (gamma, "check_BeadRod_calc_dr");
  BeadRod_solve_gamma (br, gamma);

  double *g0 = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (g0, "check_BeadRod_calc_dr");
  BeadRod_update_gamma (br, gamma, g0);
  double *g1 = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (g1, "check_BeadRod_calc_dr");
  BeadRod_update_gamma_ (br, gamma, g1);
  char label[80];
  for (i = 0; i < nc; i ++)
    {
      sprintf (label, " gamma[%d]", i);
      check += compare_max (g0[i], g1[i], label,  verbose, tiny, &max);
    }

  double *f = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (f, "check_BeadRod_calc_dr");
  BeadRod_NLEQ_for_gamma (nc, gamma, f, (void *)br);

  double *f_ = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (f_, "check_BeadRod_calc_dr");
  BeadRod_NLEQ_for_gamma_ (br, gamma, f_);
  for (i = 0; i < nc; i ++)
    {
      sprintf (label, " F[%d]", i);
      check += compare_max (f[i], f_[i], label,  verbose, tiny, &max);
    }

  free (r0);
  free (rr);

  free (a);
  free (ia);
  free (ib);

  BeadRod_free (br);

  free (gamma);
  free (g0);
  free (g1);
  free (f);
  free (f_);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
