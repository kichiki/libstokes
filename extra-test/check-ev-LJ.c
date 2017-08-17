/* test code for ev-LJ.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ev-LJ.c,v 1.1 2008/05/24 06:06:09 kichiki Exp $
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
#include "ev-LJ.h" // struct EV_DH


/* check reading SCM script
 */
int
check_EV_LJ_calc_force (double r,
			int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_EV_LJ_calc_force : start\n");
    }

  int check = 0;
  double max = 0.0;

  int n = 2; // two particle problem
  double *f = (double *)malloc (sizeof (double) * 3 * n);
  CHECK_MALLOC (f, "check_EV_LJ_calc_force");

  // set struct stokes *sys
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_EV_LJ_calc_force");
  stokes_set_np (sys, n, n);

  sys->pos[0] = 0.0;
  sys->pos[1] = 0.0;
  sys->pos[2] = 0.0;

  sys->pos[3] = r;
  sys->pos[4] = 0.0;
  sys->pos[5] = 0.0;

  struct EV_LJ *ev_LJ = EV_LJ_init (n);
  CHECK_MALLOC (ev_LJ, "check_EV_LJ_calc_force");

  // parameters
  double e  = 10.0; // in kT
  double r0 = 2.0;  // in length ==> r can be in the range (0, 2).

  ev_LJ->e[0] = e;
  ev_LJ->e[1] = e;

  ev_LJ->r0[0] = r0;
  ev_LJ->r0[1] = r0;

  // scale
  double peclet = 1.0;
  EV_LJ_scale (ev_LJ, peclet);

  check += compare_max (ev_LJ->e[0], e / peclet / r0, " e[0]",
			verbose, tiny, &max);
  check += compare_max (ev_LJ->e[1], e / peclet / r0, " e[1]",
			verbose, tiny, &max);

  EV_LJ_calc_force (ev_LJ, sys, f, 0); // after zero-clear

  check += compare_max (f[1], 0.0, " Fy(0,0,0)",
			verbose, tiny, &max);
  check += compare_max (f[2], 0.0, " Fz(0,0,0)",
			verbose, tiny, &max);

  check += compare_max (f[4], 0.0, " Fy(r,0,0)",
			verbose, tiny, &max);
  check += compare_max (f[5], 0.0, " Fz(r,0,0)",
			verbose, tiny, &max);

  // scalar part of the force is positive
  double overlap = 2.0 - r;
  double r_LJ = r0 - overlap;
  if (r_LJ < 0.0)
    {
      r_LJ = 1.0e-12;
    }
  double ri = r0 / r_LJ;
  double ri6 = pow (ri, 6.0);
  double ri7 = ri6 * ri;
  double f_check = 12.0 * e * ri7 * (ri6 - 1.0);
  if (overlap < 0.0)
    {
      check += compare_max (f[0], 0.0, " Fx(0,0,0)",
			    verbose, tiny, &max);
      check += compare_max (f[3], 0.0, " Fx(r,0,0)",
			    verbose, tiny, &max);
    }
  else
    {
      check += compare_max (f[0], -f_check, " Fx(0,0,0)",
			    verbose, tiny, &max);
      check += compare_max (f[3], f_check, " Fx(r,0,0)",
			    verbose, tiny, &max);
    }

  free (f);
  stokes_free (sys);
  EV_LJ_free (ev_LJ);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
