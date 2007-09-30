/* test code for Brownian dynamics scheme
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bd.c,v 1.2 2007/09/30 03:54:47 kichiki Exp $
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
#include <math.h> // sqrt()
#include <libstokes.h> // struct stokes
#include <ewald.h> // atimes_3all()
#include <KIrand.h> // struct KIrand

#include "check.h" // compare()
#include "memory-check.h" // macro CHECK_MALLOC

#include "chebyshev.h" // chebyshev_coef(), etc.
#include "dsaupd_c.h" // dsaupd_wrap_min_max()


static double
check_bd_my_inv_sqrt (double x)
{
  return (1.0 / sqrt (x));
}


/*
 * INPUT
 *  verbose : if non-zero, print results
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_bd (int ncheb, int np, double r,
	  int verbose, double tiny)
{
  char label[80];
  sprintf (label, "check-bd (ncheb=%d, eps=%e)", ncheb, tiny);

  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "%s : start\n", label);
    }

  int check = 0;


  struct KIrand *rng = KIrand_init ();

  struct stokes *sys = stokes_init ();
  stokes_set_np (sys, np, np); // all mobile particles
  sys->lubmin2 = 4.0000000001;
  sys->lubmax  = 4.0;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 1, stdout);

  double *pos = (double *)malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (pos, "test-bd");

  double x0 = (double)np * r / 2.0;
  int i;
  for (i = 0; i < np; i ++)
    {
      pos [i*3  ] = (double)i * r - x0;
      pos [i*3+1] = 0.0;
      pos [i*3+2] = 0.0;
    }
  stokes_set_pos (sys, pos);


  // F version
  sys->version = 0;
  int n = sys->np * 3;

  double eig[2];
  dsaupd_wrap_min_max (n, eig, atimes_3all, (void *)sys, 1.0e-12);
  //dnaupd_wrap_min_max (n, eig, atimes_3all, (void *)sys, 1.0e-12);


  // prepare chebyshev coefficients for f(x) = 1/sqrt(x)
  double *a_cheb = (double *)malloc (sizeof (double) * ncheb);
  CHECK_MALLOC (a_cheb, "test_bd");
  chebyshev_coef (ncheb, check_bd_my_inv_sqrt, eig[0], eig[1], a_cheb);

  double *y = (double *)malloc (sizeof (double) * n);
  double *z = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y, "test_bd");
  CHECK_MALLOC (z, "test_bd");
  for (i = 0; i < n; i ++)
    {
      y[i] = KIrand_Gaussian (rng);
    }

  // calc vector z = (M)^{-1/2}.y
  chebyshev_eval_atimes (ncheb, a_cheb,
			 n, y, z, 
			 eig[0], eig[1],
			 atimes_3all, (void *)sys);
  double err_cheb
    = chebyshev_error_minvsqrt (n, y, z, 
				atimes_3all, (void *)sys);

  free (y);
  free (z);
  free (a_cheb);

  free (pos);
  stokes_free (sys);

  KIrand_free (rng);


  check += compare (err_cheb+1.0, 1.0, label, verbose, tiny);

  if (verbose != 0)
    {
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
