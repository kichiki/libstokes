/* test code for mul_matrices() in matrix.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-mul-matrices.c,v 1.1 2007/03/19 02:55:29 kichiki Exp $
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
check_mul_matrices (int n,
		    int verbose)
{
  int check = 0;


  double *a = NULL;
  a = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_mul_matrices");

  double *b = NULL;
  b = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_mul_matrices");

  double *c1 = NULL;
  c1 = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (c1, "check_mul_matrices");

  double *c2 = NULL;
  c2 = (double *)malloc (sizeof (double) * n * n);
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
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j++)
	{
	  double d;
	  d = fabs (c1 [i*n+j] - c2 [i*n+j]);
	  if (d > 1.0e-10)
	    {
	      if (verbose != 0)
		{
		  fprintf (stdout, "check_mul_matrices:"
			   " %d %d %f %f %e\n",
			   i, j, c1 [i*n+j], c2 [i*n+j], d);
		}
	      check ++;
	    }
	}
    }

  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_mul_matrices: PASSED\n");
    }

  free (a);
  free (b);
  free (c1);
  free (c2);

  return (check);
}
