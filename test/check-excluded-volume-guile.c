/* test code for excluded-volume-guile.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-excluded-volume-guile.c,v 1.1 2008/07/17 21:46:46 kichiki Exp $
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
#include "excluded-volume.h" // struct EV
#include "excluded-volume-guile.h"


/* check reading SCM script
 */
int
check_EV_guile_get (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_EV_guile_get : start\n");
    }

  int check = 0;
  double max = 0.0;

  char *filename = "check-excluded-volume-guile.scm";

  // read a parameter file
  guile_load (filename);

  int np = 5;
  double length = 3.3;
  double peclet = 6.6;

  struct EV *ev = EV_guile_get ("ev", np, length, peclet);
  if (ev == NULL) // FALSE
    {
      fprintf (stdout, " fail to parse angles.\n");
      exit (1);
    }

  // check the result
  check += compare_max ((double)ev->n, (double)np, " ev->n",
			verbose, tiny, &max);

  char label[80];
  int i;
  for (i = 0; i < ev->n; i ++)
    {
      double v;
      double N;
      double b;
      if (i < 3)
	{
	  v  = 0.0012;
	  double p1 = 1.0; // Asp
	  double p2 = 2.1; // Ls

	  b = 3.0 / (peclet * p1);
	  N = p2 / b;
	}
      else
	{
	  v  = 0.002;
	  N = 19.8;  // p1
	  b = 106.0; // p2
	  b /= length;
	}

      double ls = sqrt (N * b * b / 3.0);
      double l = ls / sqrt (1.5);
      double A = pow (M_PI, -1.5) // (1 / pi)^{3/2}
	/ peclet
	* (v / (length * length * length))
	* N * N
	* pow (l, -5.0);

      sprintf (label, " l[%d]", i);
      check += compare_max (ev->l[i], l, label,
			    verbose, tiny, &max);
      sprintf (label, " A[%d]", i);
      check += compare_max (ev->A[i], A, label,
			    verbose, tiny, &max);
    }

  EV_free (ev);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
