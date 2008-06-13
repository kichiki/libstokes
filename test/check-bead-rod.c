/* test code for bead-rod.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bead-rod.c,v 1.1 2008/06/13 03:12:18 kichiki Exp $
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
  int *ia = (int *)malloc (sizeof (int) * nc);
  int *ib = (int *)malloc (sizeof (int) * nc);
  CHECK_MALLOC (ia, "BeadRod_set_connections_straight");
  CHECK_MALLOC (ib, "BeadRod_set_connections_straight");
  for (i = 0; i < nc; i ++)
    {
      ia[i] = i+1;
      ib[i] = i;
    }

  struct BeadRod *br = BeadRod_init (nc,
				     NULL, // a[]
				     ia, ib);
  CHECK_MALLOC (br, "BeadRod_init_for_test");
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

  double *dr = (double *)malloc (sizeof (double) * n * 3);
  CHECK_MALLOC (dr, "check_BeadRod_constraint_displacement");
  BeadRod_constraint_displacement (br, gamma, n, dr);

  double *du = (double *)malloc (sizeof (double) * nc * 3);
  CHECK_MALLOC (du, "check_BeadRod_constraint_displacement");
  BeadRod_bead_to_connector (nc,
			     br->ia, br->ib,
			     br->a, dr, du);


  char label[80];
  double c = - 0.5 * br->c1; // = - 2 dt / zeta
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

  double *gamma = (double *)malloc (sizeof (double) * nc);
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


  double *gamma = (double *)malloc (sizeof (double) * nc);
  CHECK_MALLOC (gamma, "check_BeadRod_solve_iter_gamma");

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
