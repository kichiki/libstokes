/* test for ev-dh-grid.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ev-dh-grid.c,v 1.1 2008/10/31 05:52:43 kichiki Exp $
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
#include <stdio.h> /* printf() fprintf() */
#include <stdlib.h> /* exit() */
#include <string.h> /* strcmp() */

#include "memory-check.h" // CHECK_MALLOC
#include "check.h" // compare_max
#include "bench.h" // ptime_ms_d()

#include "ev-dh-grid.h"

/*
 * INPUT
 */
int
check_EV_DH_calc_force_grid (int np,
			     int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_EV_DH_calc_force_grid(np=%d): start\n", np);
    }

  int check = 0;
  double max = 0.0;


  // EV_DH parameters
  double length = 1.0; // [nm]
  double peclet = 1.0;
  double rd     = 3.07; // [nm]
  double T      = 298.0; // [K]
  double e      = 80.0;

  double eps = 1.0e-6;

  struct EV_DH *ev_dh = EV_DH_init (length, peclet, rd, T, e, eps, np);
  CHECK_MALLOC (ev_dh, "check_EV_DH_calc_force_grid");

  double nu = 2.43; // [e/nm]
  double l0 = 5.0; // [nm]
  int i;
  for (i = 0; i < np; i ++)
    {
      ev_dh->nu[i] = nu * l0;
    }
  double ev_dh_r = sqrt (ev_dh->r2);


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_EV_DH_calc_force_grid");
  sys->version = 0; // F version

  stokes_set_np (sys, np, np);


  double r0 = 0.5 * ev_dh_r; // (mean) lattice distance
  int nx = (int)(pow ((double)np, 1.0 / 3.0)) + 1;

  // random configuration
  double L = r0 * (double)nx;
  srand48 (0);
  for (i = 0; i < np; i ++)
    {
      int ix = i*3;
      sys->pos[ix  ] = L * drand48 ();
      sys->pos[ix+1] = L * drand48 ();
      sys->pos[ix+2] = L * drand48 ();
    }

  /*
  // array configuration
  double r0 = 0.5 * ev_dh_r; // lattice distance
  int nx = (int)(pow ((double)np, 1.0 / 3.0)) + 1;
  i = 0;
  int ix, iy, iz;
  for (ix = 0; ix < nx; ix ++)
    {
      double x = (double)ix * r0;
      for (iy = 0; iy < nx; iy ++)
	{
	  double y = (double)iy * r0;
	  for (iz = 0; iz < nx; iz ++)
	    {
	      double z = (double)iz * r0;

	      if (i >= np) break;
	      // now (i < np)

	      int ip = i*3;
	      sys->pos [ip  ] = x;
	      sys->pos [ip+1] = y;
	      sys->pos [ip+2] = z;

	      i++;
	    }
	  if (i >= np) break;
	}
      if (i >= np) break;
    }
  */


  double *f1 = (double *)calloc (sizeof (double), np * 3);
  double *f2 = (double *)calloc (sizeof (double), np * 3);
  CHECK_MALLOC (f1, "check_EV_DH_calc_force_grid");
  CHECK_MALLOC (f2, "check_EV_DH_calc_force_grid");

  double t, t0;
  t0 = ptime_ms_d ();
  EV_DH_calc_force (ev_dh, sys, f1, 0); // zero-clear and set the force
  t = ptime_ms_d ();
  double t1 = t - t0;

  t0 = ptime_ms_d ();
  struct RYUON_grid *g = GRID_init_all_by_l (sys, ev_dh_r);
  CHECK_MALLOC (g, "check_EV_DH_calc_force_grid");
  EV_DH_calc_force_grid (ev_dh, g, sys, f2, 0); // zero-clear and set the force
  GRID_free (g);
  t = ptime_ms_d ();
  double t2 = t - t0;


  char label [80];
  for (i = 0; i < np; i ++)
    {
      int ix = i * 3;
      sprintf (label, " fx [%d]", i);
      check += compare_max (f1[ix], f2[ix], label, verbose, tiny, &max);

      sprintf (label, " fy [%d]", i);
      check += compare_max (f1[ix+1], f2[ix+1], label, verbose, tiny, &max);

      sprintf (label, " fz [%d]", i);
      check += compare_max (f1[ix+2], f2[ix+2], label, verbose, tiny, &max);
    }

  EV_DH_free (ev_dh);
  stokes_free (sys);
  free (f1);
  free (f2);


  if (verbose != 0)
    {
      // benchmark
      fprintf (stdout, " CPU times normal vs grid : %f %f\n", t1, t2);
      fprintf (stdout, "\n");

      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
