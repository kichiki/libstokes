/* test code for dgeev_c.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-dgeev_c.c,v 1.1 2007/09/30 03:50:57 kichiki Exp $
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


/*
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_dgeev_c (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_dgeev_c : start\n");
    }

  int check = 0;


  int n = 3;
  double *a  = (double *)malloc (sizeof (double) * n * n);
  double *wr = (double *)malloc (sizeof (double) * n);
  double *wi = (double *)malloc (sizeof (double) * n);
  double *v  = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a,  "check_dgeev_c");
  CHECK_MALLOC (wr, "check_dgeev_c");
  CHECK_MALLOC (wi, "check_dgeev_c");
  CHECK_MALLOC (v,  "check_dgeev_c");

  a[0] = +6.0;
  a[1] = -1.0;
  a[2] = +5.0;

  a[3] = -3.0;
  a[4] = +2.0;
  a[5] = -3.0;

  a[6] = -7.0;
  a[7] = +1.0;
  a[8] = -6.0;

  // answeres
  double *wr0 = (double *)malloc (sizeof (double) * n);
  double *wi0 = (double *)malloc (sizeof (double) * n);
  double *v0  = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (wr0, "check_dgeev_c");
  CHECK_MALLOC (wi0, "check_dgeev_c");
  CHECK_MALLOC (v0,  "check_dgeev_c");
  wr0[0] = -9.999999999999980e-01;
  wr0[1] = 9.999999999999993e-01;
  wr0[2] = 2.000000000000000e+00;

  wi0[0] = 0.000000000000000e+00;
  wi0[1] = 0.000000000000000e+00;
  wi0[2] = 0.000000000000000e+00;

  v0[0] = -5.345224838248490e-01;
  v0[1] = 2.672612419124243e-01;
  v0[2] = 8.017837257372731e-01;

  v0[3] = 7.071067811865475e-01;
  v0[4] = -1.897943516187247e-16;
  v0[5] = -7.071067811865475e-01;

  v0[6] = 5.773502691896260e-01;
  v0[7] = -5.773502691896257e-01;
  v0[8] = -5.773502691896257e-01;


  // set matrix
  // a[n*n]
  dgeev_wrap (n, a, wr, wi, v);


  int i;
  char label[80];
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "wr[%d]", i);
      check += compare (wr[i], wr0[i], label, verbose, tiny);

      sprintf (label, "wi[%d]", i);
      check += compare (wi[i]+1.0, wi0[i]+1.0, label, verbose, tiny);
    }

  int j;
  for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
	{
	  sprintf (label, "v%d[%d]", i, j);
	  check += compare (v[i*n+j], v0[i*n+j], label, verbose, tiny);
	}
    }

  free (a);
  free (wr);
  free (wi);
  free (v);
  free (wr0);
  free (wi0);
  free (v0);

  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
