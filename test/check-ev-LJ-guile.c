/* test code for ev-LJ-guile.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ev-LJ-guile.c,v 1.1 2008/05/24 06:06:55 kichiki Exp $
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
#include "ev-LJ.h" // struct EV_LJ
#include "ev-LJ-guile.h"


/* check reading SCM script
 */
int
check_EV_LJ_guile_get (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_EV_LJ_guile_get : start\n");
    }

  int check = 0;
  double max = 0.0;


  int np = 5;

  char *filename = "check-ev-LJ-guile.scm";

  // read a parameter file
  guile_load (filename);

  struct EV_LJ *ev_LJ = EV_LJ_guile_get ("ev-LJ", np);
  if (ev_LJ == NULL) // FALSE
    {
      fprintf (stdout, " fail to parse ev_LJ.\n");
      exit (1);
    }

  check += compare_max ((double)ev_LJ->n, (double)np, "n",
			verbose, tiny, &max);

  // check particles
  char label[80];
  int i;
  for (i = 0; i < 3; i++)
    {
      sprintf (label, " flag[%d]", i);
      check += compare_max ((double)ev_LJ->flag[i], 0.0, label,
			    verbose, tiny, &max);
      double e = 10.0;
      sprintf (label, " e[%d]", i);
      check += compare_max (ev_LJ->e[i], e, label,
			    verbose, tiny, &max);
      double r0 = 1.0;
      sprintf (label, " r0[%d]", i);
      check += compare_max (ev_LJ->r0[i], r0, label,
			    verbose, tiny, &max);
    }
  for (; i < 5; i++)
    {
      sprintf (label, " flag[%d]", i);
      check += compare_max ((double)ev_LJ->flag[i], 0.0, label,
			    verbose, tiny, &max);
      double e = 8.0;
      sprintf (label, " e[%d]", i);
      check += compare_max (ev_LJ->e[i], e, label,
			    verbose, tiny, &max);
      double r0 = 2.0;
      sprintf (label, " r0[%d]", i);
      check += compare_max (ev_LJ->r0[i], r0, label,
			    verbose, tiny, &max);
    }

  double peclet = 3.14;
  EV_LJ_scale (ev_LJ, peclet);
  for (i = 0; i < 3; i++)
    {
      sprintf (label, " flag[%d]", i);
      check += compare_max ((double)ev_LJ->flag[i], 1.0, label,
			    verbose, tiny, &max);
      double r0 = 1.0;
      double e = 10.0 / (peclet * r0);
      sprintf (label, " scaled e[%d]", i);
      check += compare_max (ev_LJ->e[i], e, label,
			    verbose, tiny, &max);
    }
  for (; i < 5; i++)
    {
      sprintf (label, " flag[%d]", i);
      check += compare_max ((double)ev_LJ->flag[i], 1.0, label,
			    verbose, tiny, &max);
      double r0 = 2.0;
      double e = 8.0 / (peclet * r0);
      sprintf (label, " scaled e[%d]", i);
      check += compare_max (ev_LJ->e[i], e, label,
			    verbose, tiny, &max);
    }

  EV_LJ_free (ev_LJ);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
