/* test code for dpotrf_c.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-dpotrf_c.c,v 1.2 2007/12/01 18:32:03 kichiki Exp $
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
#include <math.h>
#include "memory-check.h"
#include "check.h" // compare()

#include <dpotrf_c.h>


int
check_dpotrf_c (int n,
		int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_dpotrf_c(n=%d) : start\n", n);
    }

  int check = 0;
  double max = 0.0;


  double *a = (double *)malloc (sizeof (double) * n * n);
  double *l = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_dpotrf_c");
  CHECK_MALLOC (l, "check_dpotrf_c");


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

  //dpotf2_wrap (n, a, l);
  dpotrf_wrap (n, a, l);

  char label[80];
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  double b = 0.0;
	  int k;
	  for (k = 0; k < n; k ++)
	    {
	      b += l[i*n+k] * l[j*n+k];
	    }
	  sprintf (label, "a[%d,%d]", i, j);
	  check += compare_max (a[i*n+j], b, label, verbose, tiny, &max);
	}
    }


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

int
check_dpotf2_c (int n,
		int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_dpotf2_c(n=%d) : start\n", n);
    }

  int check = 0;
  double max = 0.0;


  double *a = (double *)malloc (sizeof (double) * n * n);
  double *l = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a, "check_dpotf2_c");
  CHECK_MALLOC (l, "check_dpotf2_c");


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

  dpotf2_wrap (n, a, l);

  char label[80];
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  double b = 0.0;
	  int k;
	  for (k = 0; k < n; k ++)
	    {
	      b += l[i*n+k] * l[j*n+k];
	    }
	  sprintf (label, "a[%d,%d]", i, j);
	  check += compare_max (a[i*n+j], b, label, verbose, tiny, &max);
	}
    }


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
