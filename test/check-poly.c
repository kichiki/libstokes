/* test code for polydisperse handling for non-periodic systems
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-poly.c,v 1.1 2007/04/03 02:36:43 kichiki Exp $
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
#include <math.h> // fabs()
#include "memory-check.h" // macro CHECK_MALLOC

#include <stokes.h> // struct stokes
#include <non-ewald.h> // scalars_nonewald_poly(), scalars_nonewald()


static int
compare (double x, double y, char *label,
	 int verbose, double tiny)
{
  int check = 0;

  double d;
  d = fabs (x - y);
  if (d > tiny)
    {
      if (verbose != 0)
	{
	  fprintf (stdout, "%s %f %f %e\n",
		   label, x, y, d);
	}
      check ++;
    }

  return (check);
}

/* 
 * INPUT
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_scalars_nonewald_poly (int verbose,
			     double tiny)
{
  int check = 0;


  double r = 2.5;

  double poly[11];
  scalars_nonewald_poly (2, /* FTS */ r, 1.0, 1.0, poly);

  double mono[11];
  scalars_nonewald (2, /* FTS */ r, mono);

  check += compare (poly [0], mono [0], "check_scalars_nonewald_poly : xa",
		    verbose, tiny);
  check += compare (poly [1], mono [1], "check_scalars_nonewald_poly : ya",
		    verbose, tiny);
  check += compare (poly [2], mono [2], "check_scalars_nonewald_poly : yb",
		    verbose, tiny);
  check += compare (poly [3], mono [3], "check_scalars_nonewald_poly : xc",
		    verbose, tiny);
  check += compare (poly [4], mono [4], "check_scalars_nonewald_poly : yc",
		    verbose, tiny);
  check += compare (poly [5], mono [5], "check_scalars_nonewald_poly : xg",
		    verbose, tiny);
  check += compare (poly [6], mono [6], "check_scalars_nonewald_poly : yg",
		    verbose, tiny);
  check += compare (poly [7], mono [7], "check_scalars_nonewald_poly : yh",
		    verbose, tiny);
  check += compare (poly [8], mono [8], "check_scalars_nonewald_poly : xm",
		    verbose, tiny);
  check += compare (poly [9], mono [9], "check_scalars_nonewald_poly : ym",
		    verbose, tiny);
  check += compare (poly[10], mono[10], "check_scalars_nonewald_poly : zm",
		    verbose, tiny);

  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_scalars_nonewald_poly : PASSED\n");
    }

  return (check);
}
