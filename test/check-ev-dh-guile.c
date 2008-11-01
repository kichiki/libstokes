/* test code for ev-dh-guile.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ev-dh-guile.c,v 1.3 2008/11/01 05:55:39 kichiki Exp $
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
#include "ev-dh-guile.h"


/* check reading SCM script
 */
int
check_EV_DH_guile_get (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_EV_DH_guile_get : start\n");
    }

  int check = 0;
  double max = 0.0;


  double length = 2.0;
  double peclet = 1.0;
  int np = 5;

  char *filename = "check-ev-dh-guile.scm";

  // read a parameter file
  guile_load (filename);

  struct EV_DH *ev_dh = EV_DH_guile_get ("ev-dh",
					 length, peclet, np);
  if (ev_dh == NULL) // FALSE
    {
      fprintf (stdout, " fail to parse ev_dh.\n");
      exit (1);
    }

  double eps = 1.0e-6;
  double T_  = 298.0;
  double e_  = 80.0;
  double rd_ = 3.07;

  double r_  = - log (eps) * rd_ / length;
  double r2_ = r_ * r_;
  check += compare_max (ev_dh->r2, r2_, " r2",
			verbose, tiny, &max);

  // Boltzmann constant
  double kB = 1.3806504e-23; // [N m / K]
  // elementary charge
  double ec = 1.602176487e-19; // [C]
  // permittivity of vacuum
  double e0 = 8.8541878176e-12; // [F / m] = [C^2 / N m^2]

  double a_sys_ = ec * ec / (4.0 * M_PI * e_ * e0 * peclet * kB * T_ * length);
  check += compare_max (ev_dh->a_sys, a_sys_, " a_sys",
			verbose, tiny, &max);
  
  rd_ /= length;
  check += compare_max (ev_dh->rd, rd_, " rd",
			verbose, tiny, &max);

  check += compare_max ((double)ev_dh->flag_grid, 1.0, " flag_grid",
			verbose, tiny, &max);

  check += compare_max ((double)ev_dh->n, (double)np, " np",
			verbose, tiny, &max);

  // check particles
  char label[80];
  int i;
  for (i = 0; i < 3; i++)
    {
      double nu = 2.43 * 5.0;
      sprintf (label, " nu[%d]", i);
      check += compare_max (ev_dh->nu[i], nu, label,
			    verbose, tiny, &max);
    }
  for (; i < 5; i++)
    {
      double nu = 2.0 * 4.0;
      sprintf (label, " nu[%d]", i);
      check += compare_max (ev_dh->nu[i], nu, label,
			    verbose, tiny, &max);
    }

  EV_DH_free (ev_dh);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
