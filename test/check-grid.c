/* test for grid.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-grid.c,v 1.1 2008/10/31 05:51:54 kichiki Exp $
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

#include "grid.h"

/*
 * INPUT
 */
int
check_GRID_ixyz_to_in_to_ixyz (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_GRID_ixyz_to_in_to_ixyz: start\n");
    }

  int check = 0;
  double max = 0.0;


  double x0 = 0.0;
  double x1 = 10.0;
  double n  = 10;
  struct RYUON_grid *g = GRID_init (x0, x1,
				    x0, x1,
				    x0, x1,
				    n, n, n);
  CHECK_MALLOC (g, "GRID_ixyz_to_in_to_ixyz");

  char label [80];
  int ixx, iyy, izz;
  int ix, iy, iz;
  for (ix = 0; ix < n; ix ++)
    {
      for (iy = 0; iy < n; iy ++)
	{
	  for (iz = 0; iz < n; iz ++)
	    {
	      int in = GRID_ixyz_to_in (g, ix, iy, iz);
	      GRID_in_to_ixyz (g, in,
			       &ixx, &iyy, &izz);
	      sprintf (label, " ix for in=%d", in);
	      check += compare_max ((double)ix, (double)ixx,
				    label, verbose, tiny, &max);
	      sprintf (label, " iy for in=%d", in);
	      check += compare_max ((double)iy, (double)iyy,
				    label, verbose, tiny, &max);
	      sprintf (label, " iz for in=%d", in);
	      check += compare_max ((double)iz, (double)izz,
				    label, verbose, tiny, &max);
	    }
	}
    }

  GRID_free (g);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
