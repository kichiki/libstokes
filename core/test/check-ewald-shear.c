/* test code for ewald.c in the new shear mode
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ewald-shear.c,v 1.1 2007/12/22 18:24:23 kichiki Exp $
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
#include <matrix.h> // multiply_extmat_with_extvec_3fts()

#include "check.h" // compare()


/** check routines **/

/* compare atimes_3all in periodic B.C. with shear-mode 1 or 2
 * based on the fact that for the simple-cubic lattice with 
 * the size (2L, L, L) containing 2 particles, the shift with 
 * (+/-)L in shear B.C. is the same without the shift.
 * INPUT
 *  version   : 0 (F), 1 (FT), 2 (FTS)
 *  shear_mode: 1 for (x = flow dir, y = grad dir)
 *              2 for (x = flow dir, z = grad dir)
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
check_atimes_3all_ewald_shear (int version,
			       int shear_mode,
			       double phi,
			       double ewald_tr, double ewald_eps,
			       int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_atimes_3all_ewald_shear\n"
	       "(%d,shear-mode=%d,phi=%f,tr=%f,eps=%e)"
	       ": start\n",
	       version, shear_mode, phi, ewald_tr, ewald_eps);
    }

  int check = 0;
  double max = 0.0;


  double l = pow (M_PI / 0.75 / phi, 1.0 / 3.0);

  int n0;
  if      (version == 0) n0 =  3;
  else if (version == 1) n0 =  6;
  else                   n0 = 11;


  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_atimes_3all_ewald_shear");
  sys->version = version;

  // N=2
  int np = 2;
  int n = n0 * np;
  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;

  sys->pos [3] = l;
  sys->pos [4] = 0.0;
  sys->pos [5] = 0.0;

  sys->periodic = 1;
  stokes_set_l (sys, 2.0*l, l, l);
  
  double xi = xi_by_tratio (sys, ewald_tr);
  stokes_set_xi (sys, xi, ewald_eps);

  sys->lubmin2 = 4.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-12, 1, stderr);


  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_atimes_3all_ewald_shear");
  int i;
  for (i = 0; i < n; i ++)
    {
      x[i] = 1.0;
    }

  // shift = 0
  double *y1 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y1, "check_atimes_3all_ewald_shear");

  // set shear-mode
  stokes_set_shear (sys, shear_mode, 1.0); // shear_rate = 1.0
  sys->shear_shift = 0.0;

  atimes_ewald_3all (n, x, y1, (void *)sys);


  // shift = + L
  double *y2 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y2,  "check_atimes_3all_ewald_shear");

  sys->shear_shift = +l;

  atimes_ewald_3all (n, x, y2, (void *)sys);

  char label [80];
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "+L (%d) ", i);
      check += compare_max (y1[i], y2[i], label, verbose, tiny, &max);
    }


  // shift = - L
  sys->shear_shift = -l;

  atimes_ewald_3all (n, x, y2, (void *)sys);

  // between y1 and y2
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "-L (%d) ", i);
      check += compare_max (y1[i], y2[i], label, verbose, tiny, &max);
    }

  free (y1);
  free (y2);
  free (x);
  stokes_free (sys);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
