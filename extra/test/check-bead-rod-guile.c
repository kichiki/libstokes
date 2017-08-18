/* test code for bead-rod-guile.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bead-rod-guile.c,v 1.1 2008/07/17 03:04:42 kichiki Exp $
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
#include "bead-rod-guile.h"
#include "bead-rod.h"


/* check reading SCM script
 */
int
check_BeadRod_guile_get (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BeadRod_guile_get : start\n");
    }

  int check = 0;
  double max = 0.0;

  char *filename = "check-bead-rod-guile.scm";

  // read a parameter file
  guile_load (filename);

  struct BeadRod *br = BeadRod_guile_get ("constraints",
					  NULL, // struct stokes
					  1.0); // length
  if (br == NULL) // FALSE
    {
      fprintf (stdout, " fail to parse constraints.\n");
      exit (1);
    }

  // nc
  check += compare_max ((double)br->nc, 5.0, " nc", verbose, tiny, &max);
  // scheme
  check += compare_max ((double)br->scheme, 1.0,
			" scheme", verbose, tiny, &max);

  // check a[]
  char label[80];
  int i;
  for (i = 0; i < 3; i ++)
    {
      sprintf (label, " a2[%d]", i);
      check += compare_max ((double)br->a2[i], 25.0, label,
			    verbose, tiny, &max);
    }
  for (; i < 5; i ++)
    {
      sprintf (label, " a2[%d]", i);
      check += compare_max ((double)br->a2[i], 100.0, label,
			    verbose, tiny, &max);
    }

  // check ia[] and ib
  for (i = 0; i < 5; i ++)
    {
      sprintf (label, " ia[%d]", i);
      check += compare_max ((double)br->ia[i], (double)i,
			    label, verbose, tiny, &max);

      sprintf (label, " ib[%d]", i);
      check += compare_max ((double)br->ib[i], (double)(i+1),
			    label, verbose, tiny, &max);
    }

  BeadRod_free (br);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
