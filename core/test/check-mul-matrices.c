/* test code for mul_matrices() in matrix.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-mul-matrices.c,v 1.3 2007/12/01 18:28:04 kichiki Exp $
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

#include <matrix.h> // mul_matrices()



/* tedious version
 */
static void
mul_matrices_ (const double *A, int na1, int na2,
	       const double *B, int nb1, int nb2,
	       double *C)
{
  int i, j, k;
  for (i = 0; i < na1; i ++)
    {
      for (j = 0; j < nb2; j ++)
	{
	  C [i*na2+j] = 0.0;
	  for (k = 0; k < na2; k ++)
	    {
	      C [i*na2+j] += A [i*na2+k] * B [k*nb2+j];
	    }
	}
    }
}


int
check_mul_matrices (int n, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_mul_matrices : start\n");
    }

  int check = 0;
  double max = 0.0;


  double *a = (double *)malloc (sizeof (double) * n * n);
  double *b = (double *)malloc (sizeof (double) * n * n);
  double *c1 = (double *)malloc (sizeof (double) * n * n);
  double *c2 = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_mul_matrices");
  CHECK_MALLOC (a, "check_mul_matrices");
  CHECK_MALLOC (c1, "check_mul_matrices");
  CHECK_MALLOC (c2, "check_mul_matrices");

  int i, j;
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  a[i*n+j] = (double)(i*n+j+1);
	  b[i*n+j] = (double)(i*n+j+1 + n*n);
	}
    }

  // tedious one
  mul_matrices_ (a, n, n, b, n, n, c1);

  // mul_matrices
  mul_matrices (a, n, n, b, n, n, c2);
  
  // check for mul_matrices
  char label[80];
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j++)
	{
	  sprintf (label, "check_mul_matrices : c1-c2[%d,%d]", i, j);
	  check += compare_max (c1[i*n+j], c2[i*n+j], label, verbose, tiny,
				&max);
	}
    }

  free (a);
  free (b);
  free (c1);
  free (c2);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
