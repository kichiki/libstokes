/* test code for angles.c and angles-guile.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-angles.c,v 1.1 2008/04/12 19:16:31 kichiki Exp $
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
#include "angles-guile.h"
#include "angles.h"


/* check reading SCM script
 */
int
check_guile_get_angles (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_guile_get_angles : start\n");
    }

  int check = 0;
  double max = 0.0;

  char *filename = "check-angles.scm";

  // read a parameter file
  guile_load (filename);

  struct angles *ang = guile_get_angles ("angles");
  if (ang == NULL) // FALSE
    {
      fprintf (stdout, " fail to parse angles.\n");
      exit (1);
    }

  // check the result
  int n = ang->n;
  if (n != 5)
    {
      if (verbose != 0)
	{
	  fprintf (stdout, " n = %d != 5\n", n);
	}
      check ++;
    }

  int i;
  struct angle *a;

  // check triplets
  for (i = 0, a = ang->a;
       i < 5;
       i++, a++)
    {
      if (a->ia != i)
	{
	  if (verbose != 0)
	    {
	      fprintf (stdout, " ia[%d] = %d != %d\n", i, a->ia, i);
	    }
	  check ++;
	}
      if (a->ib != (i+1))
	{
	  if (verbose != 0)
	    {
	      fprintf (stdout, " ib[%d] = %d != %d\n", i, a->ib, i+1);
	    }
	  check ++;
	}
      if (a->ic != (i+2))
	{
	  if (verbose != 0)
	    {
	      fprintf (stdout, " ic[%d] = %d != %d\n", i, a->ic, i+2);
	    }
	  check ++;
	}
    }

  // check angle parameters k and t0
  // for the first angle (3 triplets)
  for (i = 0, a = ang->a;
       i < 3;
       i++, a++)
    {
      check += compare_max (a->k, 10.0, " k[0]",
			    verbose, tiny, &max);
      check += compare_max (a->t0, M_PI, " t0[0]",
			    verbose, tiny, &max);
    }
  // for the second angle (2 triplets)
  for (;
       i < 5;
       i++, a++)
    {
      check += compare_max (a->k, 20.0, " k[1]",
			    verbose, tiny, &max);
      check += compare_max (a->t0, 0.5*M_PI, " t0[1]",
			    verbose, tiny, &max);
    }

  angles_free (ang);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


/* check angle_calc_force()
 * INPUT
 *  k  : 
 *  t0 : (degree)
 */
int
check_angles_calc_force (double k, double t0,
			 int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_angles_calc_force(k=%f, t0=%f) : start\n", k, t0);
    }

  int check = 0;
  double max = 0.0;

  // convert degree to radian
  t0 = t0 * M_PI / 180.0;

  int n = 3; // three particle problem
  double *f = (double *)malloc (sizeof (double) * 3 * n);
  CHECK_MALLOC (f, "check_angles_calc_force");

  // set struct stokes *sys
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_angles_calc_force");
  stokes_set_np (sys, n, n);

  sys->pos[0] = 1.0;
  sys->pos[1] = 0.0;
  sys->pos[2] = 0.0;

  sys->pos[3] = 0.0;
  sys->pos[4] = 0.0;
  sys->pos[5] = 0.0;

  // sys->pos[6] = -cos (theta);
  // sys->pos[7] = sin (theta);
  sys->pos[8] = 0.0;

  // set struct angles *ang
  struct angles *ang = angles_init ();
  CHECK_MALLOC (ang, "check_angles_calc_force");

  angles_add (ang, 0, 1, 2, k, t0);

  int i;
  for (i = 1; i < 100; i ++)
    {
      double theta = (double)i / 100.0 * M_PI;
      sys->pos[6] = -cos (theta);
      sys->pos[7] = sin (theta);

      // calc force
      angles_calc_force (ang, sys, f, 0 /* zero-clear */);

      // particle X
      check += compare_max (f[0], 0.0, " Fx(X)",
			    verbose, tiny, &max);
      check += compare_max (f[1], -k*(theta-t0), " Fy(X)",
			    verbose, tiny, &max);
      check += compare_max (f[2], 0.0, " Fz(X)",
			    verbose, tiny, &max);
      // particle Y
      check += compare_max (f[3], +k*(theta-t0)*sin(theta), " Fx(Y)",
			    verbose, tiny, &max);
      check += compare_max (f[4], +k*(theta-t0)*(1.0+cos(theta)), " Fy(Y)",
			    verbose, tiny, &max);
      check += compare_max (f[5], 0.0, " Fz(Y)",
			    verbose, tiny, &max);
      // particle Z
      check += compare_max (f[6], -k*(theta-t0)*sin(theta), " Fx(Z)",
			    verbose, tiny, &max);
      check += compare_max (f[7], -k*(theta-t0)*cos(theta), " Fy(Z)",
			    verbose, tiny, &max);
      check += compare_max (f[8], 0.0, " Fz(Z)",
			    verbose, tiny, &max);
    }

  free (f);
  stokes_free (sys);
  angles_free (ang);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
