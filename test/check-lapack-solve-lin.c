/* test code for lapack_solve_lin() in dgetri_c.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-lapack-solve-lin.c,v 1.1 2007/12/01 18:24:14 kichiki Exp $
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

#include "bench.h"

#include "matrix.h" // dot_prod_matrix()
#include "dgetri_c.h" // lapack_solve_lin()

void dcopy_ (int *, double *, int *, double *, int *);


/*
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_lapack_solve_lin (int n, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_lapack_solve_lin (n=%d, tiny=%e): start\n",
	       n, tiny);
    }

  int check = 0;
  double max = 0.0;


  double *a  = (double *)malloc (sizeof (double) * n * n);
  double *b  = (double *)malloc (sizeof (double) * n);
  double *ai = (double *)malloc (sizeof (double) * n * n);
  double *x1 = (double *)malloc (sizeof (double) * n);
  double *x2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (a,  "check_dgeev_c");
  CHECK_MALLOC (b,  "check_dgeev_c");
  CHECK_MALLOC (ai, "check_dgeev_c");
  CHECK_MALLOC (x1, "check_dgeev_c");
  CHECK_MALLOC (x2, "check_dgeev_c");

  int i, j;
  srand48 (0);
  for (i = 0; i < n; i ++)
    {
      a[i*n+i] = 1.0;
      for (j = 0; j < n; j ++)
	{
	  a [i*n+j] = drand48();
	}
    }
  for (i = 0; i < n; i ++)
    {
      b[i] = drand48();
    }

  int i_1 = 1;
  int nn = n * n;
  dcopy_ (&nn, a, &i_1, ai, &i_1);

  double t0 = ptime_ms_d();
  lapack_inv_ (n, ai);
  dot_prod_matrix (ai, n, n, b, x1);
  double t1 = ptime_ms_d();

  lapack_solve_lin (n, a, b, x2);
  double t2 = ptime_ms_d();


  char label[80];
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "x[%d]", i);
      check += compare_max (x1[i], x2[i], label, verbose, tiny, &max);
    }

  free (a);
  free (b);
  free (ai);
  free (x1);
  free (x2);

  if (verbose != 0)
    {
      fprintf (stdout, " CPU times for lapack_inv_() => dot_prod_matrix(): "
	       "%e [msec]\n", t1-t0);
      fprintf (stdout, " CPU times for lapack_solve_lin()                : "
	       "%e [msec]\n", t2-t1);

      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
