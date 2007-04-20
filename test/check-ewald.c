/* test code for ewald.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ewald.c,v 1.1 2007/04/20 01:59:39 kichiki Exp $
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

/* compare atimes and matrix processes for ewald_3all
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
check_ewald_3all_atimes_matrix_SC (int version,
				   double phi,
				   double ewald_tr, double ewald_eps,
				   int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout, "check_ewald_3all_atimes_matrix_SC"
	       " Tr=%f: start\n", ewald_tr);
    }

  double l;
  l = pow (M_PI / 0.75 / phi, 1.0 / 3.0);

  int n0;
  if (version == 0)      n0 =  3;
  else if (version == 1) n0 =  6;
  else                   n0 = 11;


  struct stokes * sys = NULL;
  sys = stokes_init ();
  sys->version = version;

  /*
  // N=1
  int np = 1;
  int n = n0 * np;
  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;


  sys->periodic = 1;
  stokes_set_l (sys, l, l, l);
  */
  // N=8
  int np = 8;
  int n = n0 * np;
  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;

  sys->pos [3] = l;
  sys->pos [4] = 0.0;
  sys->pos [5] = 0.0;

  sys->pos [6] = 0.0;
  sys->pos [7] = l;
  sys->pos [8] = 0.0;

  sys->pos [9] = 0.0;
  sys->pos[10] = 0.0;
  sys->pos[11] = l;

  sys->pos[12] = l;
  sys->pos[13] = l;
  sys->pos[14] = 0.0;

  sys->pos[15] = l;
  sys->pos[16] = 0.0;
  sys->pos[17] = l;

  sys->pos[18] = 0.0;
  sys->pos[19] = l;
  sys->pos[20] = l;

  sys->pos[21] = l;
  sys->pos[22] = l;
  sys->pos[23] = l;


  sys->periodic = 1;
  stokes_set_l (sys, 2.0*l, 2.0*l, 2.0*l);
  
  double xi = xi_by_tratio (sys, ewald_tr);
  stokes_set_xi (sys, xi, ewald_eps);

  sys->lubmin = 2.0000000001;
  stokes_set_iter (sys, "gmres", 2000, 20, 1.0e-6, 1, stderr);



  double *x = NULL;
  x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (x, "check_ewald_3all_atimes_matrix_SC");
  int i;
  for (i = 0; i < n; i ++)
    {
      x[i] = 1.0;
    }

  // atimes version
  double *y1 = NULL;
  y1 = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (y1, "check_ewald_3all_atimes_matrix_SC");

  atimes_ewald_3all (n, x, y1, (void *)sys);


  // matrix version
  double *y2 = NULL;
  double *mat = NULL;
  y2 = (double *)malloc (sizeof (double) * n);
  mat = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (y2,  "check_ewald_3all_atimes_matrix_SC");
  CHECK_MALLOC (mat, "check_ewald_3all_atimes_matrix_SC");

  make_matrix_mob_ewald_3all (sys, mat);
  if (version == 2)
    {
      // FTS
      multiply_extmat_with_extvec_3fts (np, mat, x, y2);
    }
  else
    {
      dot_prod_matrix (mat, n, n, x, y2);
    }


  // compare
  int check = 0;

  char label [80];
  // between y1 and y2
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "check_ewald_3all_atimes_matrix_SC"
	       " Tr=%f: (%d) ", ewald_tr, i);
      check += compare (y1[i], y2[i], label, verbose, tiny);
    }


  if (verbose != 0)
    {
      if (check == 0)
	fprintf (stdout, "check_ewald_3all_atimes_matrix_SC"
		 " Tr=%f: PASSED\n\n", ewald_tr);
      else
	fprintf (stdout, "check_ewald_3all_atimes_matrix_SC"
		 " Tr=%f: FAILED\n\n", ewald_tr);
    }

  free (y1);
  free (y2);
  free (mat);
  free (x);
  stokes_free (sys);

  return (check);
}
