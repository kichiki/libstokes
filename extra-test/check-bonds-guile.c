/* test code for bonds-guile.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bonds-guile.c,v 1.1 2008/07/17 03:06:11 kichiki Exp $
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
#include "bonds.h" // struct BONDS
#include "bonds-guile.h"


/* check reading SCM script
 */
int
check_bonds_guile_get (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_bonds_guile_get : start\n");
    }

  int check = 0;
  double max = 0.0;

  char *filename = "check-bonds-guile.scm";

  // read a parameter file
  guile_load (filename);

  struct BONDS *b = BONDS_guile_get ("bonds");
  if (b == NULL) // FALSE
    {
      fprintf (stdout, " fail to parse angles.\n");
      exit (1);
    }

  // check the result
  check += compare_max ((double)b->n, 8.0, " b->n",
			verbose, tiny, &max);

  char label[80];
  int i;
  for (i = 0; i < b->n; i ++)
    {
      if (i >= 0 && i < 3)
	{
	  sprintf (label, " type[%d]", i);
	  check += compare_max ((double)b->type[i], 0.0, label,
				verbose, tiny, &max);
	  sprintf (label, " fene[%d]", i);
	  check += compare_max ((double)b->fene[i], 0.0, label,
				verbose, tiny, &max);
	  sprintf (label, " p1[%d]", i);
	  check += compare_max (b->p1[i], 1.0, label,
				verbose, tiny, &max);
	  sprintf (label, " p2[%d]", i);
	  check += compare_max (b->p2[i], 2.1, label,
				verbose, tiny, &max);
	  sprintf (label, " p3[%d]", i);
	  check += compare_max (b->p3[i], 0.0, label,
				verbose, tiny, &max);
	  sprintf (label, " ia[%d]", i);
	  check += compare_max ((double)b->ia[i], (double)i, label,
				verbose, tiny, &max);
	  sprintf (label, " ib[%d]", i);
	  check += compare_max ((double)b->ib[i], (double)(i+1), label,
				verbose, tiny, &max);
	}
      else if (i < 6)
	{
	  sprintf (label, " type[%d]", i);
	  check += compare_max ((double)b->type[i], 2.0, label,
				verbose, tiny, &max);
	  sprintf (label, " fene[%d]", i);
	  check += compare_max ((double)b->fene[i], 1.0, label,
				verbose, tiny, &max);
	  sprintf (label, " p1[%d]", i);
	  check += compare_max (b->p1[i], 19.8, label,
				verbose, tiny, &max);
	  sprintf (label, " p2[%d]", i);
	  check += compare_max (b->p2[i], 106.0, label,
				verbose, tiny, &max);
	  sprintf (label, " p3[%d]", i);
	  check += compare_max (b->p3[i], 0.0, label,
				verbose, tiny, &max);
	  sprintf (label, " ia[%d]", i);
	  check += compare_max ((double)b->ia[i], (double)(i+1), label,
				verbose, tiny, &max);
	  sprintf (label, " ib[%d]", i);
	  check += compare_max ((double)b->ib[i], (double)(i+2), label,
				verbose, tiny, &max);
	}
      else
	{
	  sprintf (label, " type[%d]", i);
	  check += compare_max ((double)b->type[i], 7.0, label,
				verbose, tiny, &max);
	  sprintf (label, " fene[%d]", i);
	  check += compare_max ((double)b->fene[i], 0.0, label,
				verbose, tiny, &max);
	  sprintf (label, " p1[%d]", i);
	  check += compare_max (b->p1[i], 1.0e6, label,
				verbose, tiny, &max);
	  sprintf (label, " p2[%d]", i);
	  check += compare_max (b->p2[i], 0.5, label,
				verbose, tiny, &max);
	  sprintf (label, " p3[%d]", i);
	  check += compare_max (b->p3[i], 0.01, label,
				verbose, tiny, &max);
	  sprintf (label, " ia[%d]", i);
	  check += compare_max ((double)b->ia[i], (double)(i+2), label,
				verbose, tiny, &max);
	  sprintf (label, " ib[%d]", i);
	  check += compare_max ((double)b->ib[i], (double)(i+3), label,
				verbose, tiny, &max);
	}
    }

  BONDS_free (b);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
