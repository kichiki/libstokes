/* test code for BD_sqrt_by_dgeev() in brownian.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-sqrt-dgeev.c,v 1.2 2007/12/01 18:32:48 kichiki Exp $
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
#include <stdlib.h> // malloc(), free()
#include <math.h>
#include "memory-check.h"
#include "check.h" // compare()

#include <brownian.h> // BD_sqrt_by_dgeev()


int
check_BD_sqrt_by_dgeev (int n,
			int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_BD_sqrt_by_dgeev(n=%d) : start\n", n);
      fprintf (stdout, " to check sqrt(A)^2 =? A\n");
    }

  int check = 0;
  double max = 0.0;


  double *a = (double *)malloc (sizeof (double) * n * n);
  double *s = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_BD_sqrt_by_dgeev");
  CHECK_MALLOC (s, "check_BD_sqrt_by_dgeev");


  // set the matrix 'a'
  int i, j;
  for (i = 0; i < n; i ++)
    {
      a[i*n+i] = 1.0;
      for (j = i+1; j < n; j ++)
	{
	  a[i*n+j] = 1.0 / (double)(i+j+1);
	  a[j*n+i] = a[i*n+j];
	}
    }

  BD_sqrt_by_dgeev (n, a, s);

  double *a2 = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a2, "check_BD_sqrt_by_dgeev");

  int k;
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  char label [80];
	  sprintf (label, "sqrt(a)^2[%d %d]", i, j);
	  a2[i*n+j] = 0.0;
	  for (k = 0; k < n; k ++)
	    {
	      a2[i*n+j] += s[i*n+k] * s[k*n+j];
	    }
	  check += compare_max (a[i*n+j], a2[i*n+j],
				label, verbose, tiny, &max);
	}
    }
  free (a2);

  free (a);
  free (s);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
