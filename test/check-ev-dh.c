/* test code for ev-dh.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ev-dh.c,v 1.1 2008/04/26 17:37:08 kichiki Exp $
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
#include "memory-check.h"
#include "check.h" // compare()

#include "stokes-guile.h" // guile_load()
#include "ev-dh.h" // struct EV_DH


/* check reading SCM script
 */
int
check_EV_DH_calc_force (double r,
			int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_EV_DH_calc_force : start\n");
    }

  int check = 0;
  double max = 0.0;

  int n = 2; // two particle problem
  double *f = (double *)malloc (sizeof (double) * 3 * n);
  CHECK_MALLOC (f, "check_EV_DH_calc_force");

  // set struct stokes *sys
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_EV_DH_calc_force");
  stokes_set_np (sys, n, n);

  sys->pos[0] = 0.0;
  sys->pos[1] = 0.0;
  sys->pos[2] = 0.0;

  sys->pos[3] = r;
  sys->pos[4] = 0.0;
  sys->pos[5] = 0.0;

  // parameters
  double a  = 1.0;
  double pe = 1.0;
  double rd = 3.07; // [nm]
  double T  = 298.0; // [K]
  double e  = 80.0;

  double r2 = 2.0 * (r * r); // set it large enough

  struct EV_DH *ev_dh = EV_DH_init (a, pe, rd, T, e, r2, n);
  CHECK_MALLOC (ev_dh, "check_EV_DH_calc_force");

  double nu = 2.43; // [e/nm]
  double l0 = 5.0; // [nm]
  ev_dh->nu[0] = nu * l0;
  ev_dh->nu[1] = nu * l0;

  EV_DH_calc_force (ev_dh, sys, f, 0); // after zero-clear


  // Boltzmann constant
  double kB = 1.3806504e-23; // [N m / K]
  // elementary charge
  double ec = 1.602176487e-19; // [C]
  // permittivity of vacuum
  double e0 = 8.8541878176e-12; // [F / m] = [C^2 / N m^2]

  double a_sys = ec * ec / (4.0 * M_PI * e * e0 * pe * kB * T * a);
  double rda = rd / a;
  double nu2 = nu * nu * l0 * l0;

  check += compare_max (ev_dh->a_sys, a_sys, " a_sys",
			verbose, tiny, &max);
  check += compare_max (ev_dh->rd, rda, " rd",
			verbose, tiny, &max);
  check += compare_max (ev_dh->nu[0], nu * l0, " nu[0]",
			verbose, tiny, &max);
  check += compare_max (ev_dh->nu[1], nu * l0, " nu[1]",
			verbose, tiny, &max);

  // scalar part of the force is positive
  double kr = r / rda;
  double f_check = a_sys * nu2 / (r * r) * (1.0 + kr) * exp (- kr);

  /* note that F_j = f_check * (x_j - x_i)/|rij|
   *           F_i = - F_j
   * now rij = x_j - x_i = (r, 0, 0)
   * so that for r>0, F_j > 0 and F_i < 0 (pure repulsive).
   */
  check += compare_max (f[0], -f_check, " Fx(0,0,0)",
			verbose, tiny, &max);
  check += compare_max (f[1], 0.0, " Fy(0,0,0)",
			verbose, tiny, &max);
  check += compare_max (f[2], 0.0, " Fz(0,0,0)",
			verbose, tiny, &max);

  check += compare_max (f[3], f_check, " Fx(r,0,0)",
			verbose, tiny, &max);
  check += compare_max (f[4], 0.0, " Fy(r,0,0)",
			verbose, tiny, &max);
  check += compare_max (f[5], 0.0, " Fz(r,0,0)",
			verbose, tiny, &max);

  free (f);
  stokes_free (sys);
  EV_DH_free (ev_dh);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
