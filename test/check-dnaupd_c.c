/* test code for dnaupd_c.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-dnaupd_c.c,v 1.4 2007/12/12 06:34:18 kichiki Exp $
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

#include "dgeev_c.h" // dgeev_wrap()
#include "dnaupd_c.h" // dnaupd_wrap_min_max()


static void
check_dnaupd_atimes (int n, const double *x, double *b, void *user_data)
{
  double *a = (double *)user_data;
  int i, j;
  for (i = 0; i < n; i ++)
    {
      b[i] = 0.0;
      for (j = 0; j < n; j ++)
	{
	  b[i] += a[i*n+j]*x[j];
	}
    }
}

/*
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_dnaupd_c (int n, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_dnaupd_c(n=%d) : start\n", n);
    }

  int check = 0;
  double max = 0.0;


  double *a  = (double *)malloc (sizeof (double) * n * n);
  double *wr = (double *)malloc (sizeof (double) * n);
  double *wi = (double *)malloc (sizeof (double) * n);
  double *v  = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a,  "check_dnaupd_c");
  CHECK_MALLOC (wr, "check_dnaupd_c");
  CHECK_MALLOC (wi, "check_dnaupd_c");
  CHECK_MALLOC (v,  "check_dnaupd_c");

  int i, j;
  for (i = 0; i < n; i ++)
    {
      a[i*n+i] = 1.0;
      for (j = i+1; j < n; j ++)
	{
	  a[i*n+j] = (double)j;
	  a[j*n+i] = a[i*n+j];
	}
    }

  // set matrix
  // a[n*n]
  dgeev_wrap (n, a, wr, wi, v);

  double lmin = wr[0];
  double lmax = wr[0];
  for (i = 1; i < n; i ++)
    {
      if (lmin > wr[i]) lmin = wr[i];
      if (lmax < wr[i]) lmax = wr[i];
    }

  double l[2];
  dnaupd_wrap_min_max (n, l, check_dnaupd_atimes, (void *)a, tiny);

  check += compare_max (lmin, l[0], "lmin", verbose, tiny, &max);
  check += compare_max (lmax, l[1], "lmax", verbose, tiny, &max);

  free (a);
  free (wr);
  free (wi);
  free (v);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
