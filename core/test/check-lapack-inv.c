/* test code for lapack_inv_() in matrix.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-lapack-inv.c,v 1.3 2007/12/01 18:25:38 kichiki Exp $
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

#include <dgetri_c.h> // lapack_inv_()
#include <matrix.h> // mul_matrices()


/*
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_lapack_inv_ (int n, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_lapack_inv_(n=%d) : start\n", n);
    }

  int check = 0;
  double max = 0.0;


  double *a  = (double *)malloc (sizeof (double) * n * n);
  double *a_ = (double *)malloc (sizeof (double) * n * n);
  double *b  = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a,  "check_lapack_inv_");
  CHECK_MALLOC (a_, "check_lapack_inv_");
  CHECK_MALLOC (b, "check_lapack_inv_");

  int i;
  srand48 (0);
  for (i = 0; i < n * n; i ++)
    {
      a [i] = drand48();
      a_[i] = a [i];
    }

  // a = a^-1
  lapack_inv_ (n, a);

  // b = a^-1 . a
  mul_matrices (a, n, n, a_, n, n, b);
  /*
  int j, k;
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  b[i*n+j] = 0.0;
	  for (k = 0; k < n; k ++)
	    {
	      b[i*n+j] += a[i*n+k] * a_[k*n+j];
	    }
	}
    }
  */

  char label[80];
  double d;
  for (i = 0; i < n; i ++)
    {
      int j;
      for (j = 0; j < n; j ++)
	{
	  if (i == j) d = fabs (b [i*n+j] - 1.0);
	  else        d = fabs (b [i*n+j]);

	  sprintf (label, "check_lapack_inv_ : [%d]", i);
	  check += compare_max (d+1.0, 1.0, label, verbose, tiny, &max);
	}
    }

  free (a);
  free (a_);
  free (b);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
