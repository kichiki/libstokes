/* test code for chebyshev.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-chebyshev.c,v 1.3 2007/12/01 18:26:27 kichiki Exp $
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
#include "check.h" // compare()
#include "memory-check.h" // macro CHECK_MALLOC

#include "chebyshev.h"


static int
check_chebyshev_sub (int ncheb, double x0, double x1, 
		     double (*func)(double),
		     char *label,
		     int verbose, double tiny, double *max)
{
  int check = 0;

  double *a_cheb = (double *)malloc (sizeof (double) * ncheb);
  CHECK_MALLOC (a_cheb, "check_chebyshev");

  chebyshev_coef (ncheb, func, x0, x1, a_cheb);

  char label_i [80];

  int i;
  for (i = 0; i < 100; i ++)
    {
      double x = x0 + (double)i/100.0*(x1 - x0);
      double y0 = func (x);
      double y1 = chebyshev_eval  (ncheb, a_cheb, x0, x1, x);
      double y2 = chebyshev_eval_ (ncheb, a_cheb, x0, x1, x);

      sprintf (label_i, "%s(%f) y0-y1 [%d]", label, x, i);
      check += compare_max (y0, y1, label_i, verbose, tiny, max);

      sprintf (label_i, "%s(%f) y1-y2 [%d]", label, x, i);
      check += compare_max (y1, y2, label_i, verbose, tiny, max);
    }

  free (a_cheb);

  return (check);
}


/*
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_chebyshev (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_chebyshev : start\n");
    }

  int check = 0;
  double max = 0.0;

  int ncheb = 100;
  check += check_chebyshev_sub (ncheb, 0.0, M_PI, cos, "cos(x)",
				verbose, tiny, &max);
  check += check_chebyshev_sub (ncheb, 1e-2, M_PI, sin, "sin(x)",
				verbose, tiny, &max);
  check += check_chebyshev_sub (ncheb, 1e-2, 1.0, sqrt, "sqrt(x)",
				verbose, tiny, &max);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
