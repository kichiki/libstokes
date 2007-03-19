/* test code for lapack_inv_() in matrix.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-lapack-inv.c,v 1.1 2007/03/19 02:51:07 kichiki Exp $
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
check_lapack_inv_ (int n,
		   int verbose)
{
  int check = 0;


  double *a = NULL;
  a = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_lapack_inv_");

  double *a_ = NULL;
  a_ = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a_, "check_lapack_inv_");

  int i;
  for (i = 0; i < n * n; i ++)
    {
      a [i] = drand48();
      a_[i] = a [i];
    }

  // a = a^-1
  lapack_inv_ (n, a);

  double *b = NULL;
  b = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (b, "check_lapack_inv_");

  // b = a^-1 . a
  mul_matrices (a, n, n, a_, n, n, b);

  double d;
  for (i = 0; i < n; i ++)
    {
      int j;
      for (j = 0; j < n; j ++)
	{
	  if (i == j) d = fabs (b [i*n+j] - 1.0);
	  else        d = fabs (b [i*n+j]);

	  if (d > 1.0e-10)
	    {
	      if (verbose != 0)
		{
		  fprintf (stdout, "check_lapack_inv_ : %d %d %f\n",
			   i, j, b[i*n+j]);
		}
	      check ++;
	    }
	}
    }

  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_lapack_inv_ : PASSED\n");
    }


  free (a);
  free (a_);
  free (b);

  return (check);
}
