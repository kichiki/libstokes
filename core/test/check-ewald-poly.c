/* test code for polydisperse code in ewald.c
 * Copyright (C) 2007-2017 Kengo Ichiki <kengoichiki@gmail.com>
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

#include <stokes.h> // struct stokes
#include <ewald.h> // atimes_ewald_3all()
#include <bench.h> // ptime_ms_d()

#include "check.h" // compare()


/** check routines **/

/* check atimes_ewald_3all() for mono and poly(a=1) with SC config with N=1
 * INPUT
 *  version   : 0 (F), 1 (FT), 2 (FTS)
 *  phi       : volume fraction, that is, phi = (4/3)pi a^3/l^3
 *  ewald_tr  :
 *  ewald_eps :
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_atimes_ewald_3all_poly_SC_1 (int version,
				   double phi,
				   double ewald_tr, double ewald_eps,
				   int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_atimes_ewald_3all_poly_SC_1\n"
	       "(%d,phi=%f,tr=%f,eps=%f)"
	       ": start\n",
	       version, phi, ewald_tr, ewald_eps);
    }

  // initialize struct stokes *sys
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_atimes_ewald_3all_poly_SC_1");
  sys->version = version;

  int np = 1;
  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;


  sys->lubmin2 = 4.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  // periodic
  sys->periodic = 1;
  double l;
  l = pow (M_PI / 0.75 / phi, 1.0 / 3.0);
  stokes_set_l (sys, l, l, l);
  
  double xi = xi_by_tratio (sys, ewald_tr);
  stokes_set_xi (sys, xi, ewald_eps);

  int n;
  if (version == 0)      n =  3 * np;
  else if (version == 1) n =  6 * np;
  else                   n = 11 * np;

  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_atimes_ewald_3all_poly_SC_1");
  int i;
  for (i = 0; i < n; i ++)
    {
      x[i] = 1.0;
    }

  // case 1) -- mono code
  double *y1 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y1, "check_atimes_ewald_3all_poly_SC_1");

  double t0, t;
  t0 = ptime_ms_d ();
  atimes_ewald_3all (n, x, y1, (void *)sys);
  t = ptime_ms_d ();
  double ptime_mono = t - t0;


  // set sys->a to test poly version
  double a [1] = {1.0};
  stokes_set_radius (sys, a);

  // case 2) -- poly code
  double *y2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y2, "check_atimes_ewald_3all_poly_SC_1");

  t0 = ptime_ms_d ();
  atimes_ewald_3all (n, x, y2, (void *)sys);
  t = ptime_ms_d ();
  double ptime_poly = t - t0;


  // compare
  int check = 0;
  double max = 0.0;
  char label [80];
  for (i = 0; i < n; i ++)
    {
      sprintf (label, " (%d) ", i);
      check += compare_max (y1[i], y2[i], label, verbose, tiny, &max);
    }

  free (x);
  free (y1);
  free (y2);
  stokes_free (sys);


  if (verbose != 0)
    {
      fprintf (stdout, " ptime mono, poly = %.3f %.3f, poly/mono = %f\n",
	       ptime_mono, ptime_poly, ptime_poly / ptime_mono);

      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

/* check atimes_ewald_3all() for mono and poly(a=1) with SC config with N=2
 * INPUT
 *  version   : 0 (F), 1 (FT), 2 (FTS)
 *  dir       : direction of the config, 0 (x), 1 (y), 2(z).
 *  phi       : volume fraction, that is, phi = (4/3)pi a^3/l^3
 *  ewald_tr  :
 *  ewald_eps :
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_atimes_ewald_3all_poly_SC_2 (int version,
				   int dir,
				   double phi,
				   double ewald_tr, double ewald_eps,
				   int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_atimes_ewald_3all_poly_SC_2\n"
	       "(%d,dir=%d,phi=%f,tr=%f,eps=%f): start\n",
	       version, dir, phi, ewald_tr, ewald_eps);
    }

  // initialize struct stokes *sys
  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_atimes_ewald_3all_poly_SC_2");
  sys->version = version;

  int np = 2;
  stokes_set_np (sys, np, np);

  //sys->lubmin = 2.0000000001;
  sys->lubmin2 = 4.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 1, stderr);

  // periodic
  sys->periodic = 1;
  double l;
  l = pow (M_PI / 0.75 / phi, 1.0 / 3.0);

  if (dir == 1)
    {
      // x direction
      stokes_set_l (sys, 2.0*l, l, l);
  
      sys->pos [0] = 0.0;
      sys->pos [1] = 0.0;
      sys->pos [2] = 0.0;
      sys->pos [3] = l;
      sys->pos [4] = 0.0;
      sys->pos [5] = 0.0;
    }
  else if (dir == 2)
    {
      // y direction
      stokes_set_l (sys, l, 2.0*l, l);
  
      sys->pos [0] = 0.0;
      sys->pos [1] = 0.0;
      sys->pos [2] = 0.0;
      sys->pos [3] = 0.0;
      sys->pos [4] = l;
      sys->pos [5] = 0.0;
    }
  else
    {
      // z direction
      stokes_set_l (sys, l, l, 2.0*l);
  
      sys->pos [0] = 0.0;
      sys->pos [1] = 0.0;
      sys->pos [2] = 0.0;
      sys->pos [3] = 0.0;
      sys->pos [4] = 0.0;
      sys->pos [5] = l;
    }


  double xi = xi_by_tratio (sys, ewald_tr);
  stokes_set_xi (sys, xi, ewald_eps);

  int n;
  if (version == 0)      n =  3 * np;
  else if (version == 1) n =  6 * np;
  else                   n = 11 * np;

  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_atimes_ewald_3all_poly_SC_2");
  int i;
  for (i = 0; i < n; i ++)
    {
      x[i] = 1.0;
    }

  // case 1) -- mono code
  double *y1 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y1, "check_atimes_ewald_3all_poly_SC_2");

  double t0, t;
  t0 = ptime_ms_d ();
  atimes_ewald_3all (n, x, y1, (void *)sys);
  t = ptime_ms_d ();
  double ptime_mono = t - t0;


  // set sys->a to test poly version
  double a [2] = {1.0, 1.0};
  stokes_set_radius (sys, a);

  // case 2) -- poly code
  double *y2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y2, "check_atimes_ewald_3all_poly_SC_2");

  t0 = ptime_ms_d ();
  atimes_ewald_3all (n, x, y2, (void *)sys);
  t = ptime_ms_d ();
  double ptime_poly = t - t0;


  // compare
  int check = 0;
  double max = 0.0;
  char label [80];
  for (i = 0; i < n; i ++)
    {
      sprintf (label, " (%d) ", i);
      check += compare_max (y1[i], y2[i], label, verbose, tiny, &max);
    }

  free (x);
  free (y1);
  free (y2);
  stokes_free (sys);


  if (verbose != 0)
    {
      fprintf (stdout, " ptime mono, poly = %.3f %.3f, poly/mono = %f\n",
	       ptime_mono, ptime_poly, ptime_poly / ptime_mono);

      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
