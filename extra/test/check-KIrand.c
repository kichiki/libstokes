/* test code for KIrand.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-KIrand.c,v 1.2 2008/06/05 03:27:49 kichiki Exp $
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
#include "memory-check.h" // macro CHECK_MALLOC
#include "check.h" // compare()

#include "KIrand.h"


/* check KIrand_Gaussian()
 */
int
check_KIrand_Gaussian (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_KIrand_Gaussian : start\n");
    }

  int check = 0;
  double max = 0.0;

  int n = 1000000;
  double *x0 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x0, "check_KIrand_Gaussian");

  int i;
  unsigned long seed = 0;
  struct KIrand *rng = KIrand_init ();
  CHECK_MALLOC (rng, "check_KIrand_Gaussian");
  KIrand_init_genrand (rng, seed);
  for (i = 0; i < n; i ++)
    {
      x0[i] = KIrand_Gaussian (rng);
    }
  KIrand_free (rng);

  // check for Gaussian properties
  double mean = 0.0;
  for (i = 0; i < n; i ++)
    {
      mean += x0[i];
    }
  mean /= (double)n;
  double vari = 0.0;
  for (i = 0; i < n; i ++)
    {
      vari += (x0[i] - mean) * (x0[i] - mean);
    }
  vari /= (double)n;

  check += compare_max (mean+1.0, 1.0, " mean", verbose, tiny, &max);
  check += compare_max (vari,     1.0, " vari", verbose, tiny, &max);

  free (x0);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

/* check the continuation of KIrand_Gaussian()
 */
int
check_KIrand_Gaussian_cont (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_KIrand_Gaussian_cont : start\n");
    }

  int check = 0;
  double max = 0.0;

  int n = 1000;
  double *x0 = (double *)malloc (sizeof (double) * n);
  double *x1 = (double *)malloc (sizeof (double) * n);
  double *x2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x0, "check_KIrand_Gaussian_cont");
  CHECK_MALLOC (x1, "check_KIrand_Gaussian_cont");
  CHECK_MALLOC (x2, "check_KIrand_Gaussian_cont");

  int i;
  unsigned long seed = 0;
  struct KIrand *rng = NULL;


  // reference sequence
  rng = KIrand_init ();
  CHECK_MALLOC (rng, "check_KIrand_Gaussian_cont");
  KIrand_init_genrand (rng, seed);
  for (i = 0; i < n; i ++)
    {
      x0[i] = KIrand_Gaussian (rng);
    }
  KIrand_free (rng);

  // even stop
  int n1 = 500;
  rng = KIrand_init ();
  CHECK_MALLOC (rng, "check_KIrand_Gaussian_cont");
  KIrand_init_genrand (rng, seed);
  for (i = 0; i < n1; i ++)
    {
      x1[i] = KIrand_Gaussian (rng);
    }
  unsigned long mt1[MTRNG_N];
  for (i = 0; i < MTRNG_N; i ++)
    {
      mt1[i] = rng->mt[i];
    }
  int mti1 = rng->mti;
  int saved1 = rng->Gaussian_has_saved;
  double x_saved1 = rng->Gaussian_saved;
  KIrand_free (rng);

  // odd stop
  int n2 = 501;
  rng = KIrand_init ();
  CHECK_MALLOC (rng, "check_KIrand_Gaussian_cont");
  KIrand_init_genrand (rng, seed);
  for (i = 0; i < n2; i ++)
    {
      x2[i] = KIrand_Gaussian (rng);
    }
  unsigned long mt2[MTRNG_N];
  for (i = 0; i < MTRNG_N; i ++)
    {
      mt2[i] = rng->mt[i];
    }
  int mti2 = rng->mti;
  int saved2 = rng->Gaussian_has_saved;
  double x_saved2 = rng->Gaussian_saved;
  KIrand_free (rng);

  // continuation for even stop
  rng = KIrand_init ();
  CHECK_MALLOC (rng, "check_KIrand_Gaussian_cont");
  KIrand_load_Gaussian (rng, mt1, mti1, saved1, x_saved1);
  for (i = n1; i < n; i ++)
    {
      x1[i] = KIrand_Gaussian (rng);
    }
  KIrand_free (rng);

  // continuation for odd stop
  rng = KIrand_init ();
  CHECK_MALLOC (rng, "check_KIrand_Gaussian_cont");
  KIrand_load_Gaussian (rng, mt2, mti2, saved2, x_saved2);
  for (i = n2; i < n; i ++)
    {
      x2[i] = KIrand_Gaussian (rng);
    }
  KIrand_free (rng);

  // check for continuation
  char label[80];
  for (i = 0; i < n; i ++)
    {
      sprintf (label, " even-cont x[%d]", i);
      check += compare_max (x0[i], x1[i], label, verbose, tiny, &max);

      sprintf (label, " odd-cont x[%d]", i);
      check += compare_max (x0[i], x2[i], label, verbose, tiny, &max);
    }

  free (x0);
  free (x1);
  free (x2);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
