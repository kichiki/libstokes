/* test code for lubrication for polydisperse systems
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-lub-poly.c,v 1.1 2007/04/14 00:36:45 kichiki Exp $
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
#include <stokes.h> // struct stokes
#include <fts.h>
#include <matrix.h> // multiply_extmat_with_extvec_3fts()


#include "check.h" // compare()


/* check calc_lub_fts_2b_poly() with a1=a2=a
 * comparing with calc_lub_fts_2b()
 */
int
check_lub_fts_2b_poly (double r, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout, "check_lub_fts_2b_poly : start\n");
    }

  struct stokes *sys = NULL;
  sys = stokes_init ();
  sys->version = 2; // FTS

  sys->lubcut = 2.0 + 1.0e-15;
  sys->twobody_nmax = 100;
  sys->twobody_lub = 1; // lub form

  int np = 2;
  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;
  sys->pos [3] = 0.0;
  sys->pos [4] = 0.0;
  sys->pos [5] = r;
  // so that r12 = r2 - r1 = (0, 0, r)

  double uoe[22];
  double fts_mono[22];
  double fts_poly[22];

  int i;
  for (i = 0; i < 22; i ++)
    {
      uoe [i] = 1.0;
      fts_mono [i] = 0.0;
      fts_poly [i] = 0.0;
    }

  calc_lub_fts_2b (sys,
		   &uoe[0], &uoe[11],
		   &(sys->pos[0]), &(sys->pos[3]),
		   &fts_mono[0], &fts_mono[11]);

  calc_lub_fts_2b_poly (sys,
			&uoe[0], &uoe[11],
			&(sys->pos[0]), &(sys->pos[3]),
			1.0, 1.0,
			&fts_poly[0], &fts_poly[11]);

  int check = 0;
  check += compare (fts_mono [0], fts_poly [0],  "check_lub_fts_2b_poly: Fx1", verbose, tiny);
  check += compare (fts_mono [1], fts_poly [1],  "check_lub_fts_2b_poly: Fy1", verbose, tiny);
  check += compare (fts_mono [2], fts_poly [2],  "check_lub_fts_2b_poly: Fz1", verbose, tiny);
  check += compare (fts_mono [3], fts_poly [3],  "check_lub_fts_2b_poly: Tx1", verbose, tiny);
  check += compare (fts_mono [4], fts_poly [4],  "check_lub_fts_2b_poly: Ty1", verbose, tiny);
  check += compare (fts_mono [5], fts_poly [5],  "check_lub_fts_2b_poly: Tz1", verbose, tiny);
  check += compare (fts_mono [6], fts_poly [6],  "check_lub_fts_2b_poly: Sxx1", verbose, tiny);
  check += compare (fts_mono [7], fts_poly [7],  "check_lub_fts_2b_poly: Sxy1", verbose, tiny);
  check += compare (fts_mono [8], fts_poly [8],  "check_lub_fts_2b_poly: Sxz1", verbose, tiny);
  check += compare (fts_mono [9], fts_poly [9],  "check_lub_fts_2b_poly: Syz1", verbose, tiny);
  check += compare (fts_mono[10], fts_poly[10],  "check_lub_fts_2b_poly: Syy1", verbose, tiny);
  check += compare (fts_mono[11], fts_poly[11],  "check_lub_fts_2b_poly: Fx2", verbose, tiny);
  check += compare (fts_mono[12], fts_poly[12],  "check_lub_fts_2b_poly: Fy2", verbose, tiny);
  check += compare (fts_mono[13], fts_poly[13],  "check_lub_fts_2b_poly: Fz2", verbose, tiny);
  check += compare (fts_mono[14], fts_poly[14],  "check_lub_fts_2b_poly: Tx2", verbose, tiny);
  check += compare (fts_mono[15], fts_poly[15],  "check_lub_fts_2b_poly: Ty2", verbose, tiny);
  check += compare (fts_mono[16], fts_poly[16],  "check_lub_fts_2b_poly: Tz2", verbose, tiny);
  check += compare (fts_mono[17], fts_poly[17],  "check_lub_fts_2b_poly: Sxx2", verbose, tiny);
  check += compare (fts_mono[18], fts_poly[18],  "check_lub_fts_2b_poly: Sxy2", verbose, tiny);
  check += compare (fts_mono[19], fts_poly[19],  "check_lub_fts_2b_poly: Sxz2", verbose, tiny);
  check += compare (fts_mono[20], fts_poly[20],  "check_lub_fts_2b_poly: Syz2", verbose, tiny);
  check += compare (fts_mono[21], fts_poly[21],  "check_lub_fts_2b_poly: Syy2", verbose, tiny);


  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_lub_fts_2b_poly : PASSED\n\n");
    }

  stokes_free (sys);

  return (check);
}

/* check matrix_lub_fts_2b_poly() with a1=a2=a
 * comparing with matrix_lub_fts_2b()
 */
int
check_matrix_lub_fts_2b_poly (double r, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout, "check_matrix_lub_fts_2b_poly : start\n");
    }

  struct stokes *sys = NULL;
  sys = stokes_init ();
  sys->version = 2; // FTS

  sys->lubcut = 2.0 + 1.0e-15;
  sys->twobody_nmax = 100;
  sys->twobody_lub = 1; // lub form

  int np = 2;
  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;
  sys->pos [3] = 0.0;
  sys->pos [4] = 0.0;
  sys->pos [5] = r;
  // so that r12 = r2 - r1 = (0, 0, r)

  double uoe[22];
  double mat_mono[22*22];
  double mat_poly[22*22];

  int i;
  for (i = 0; i < 22; i ++)
    {
      uoe [i] = 1.0;
    }
  for (i = 0; i < 22*22; i ++)
    {
      mat_mono [i] = 0.0;
      mat_poly [i] = 0.0;
    }

  matrix_lub_fts_2b (sys,
		     0, 1,
		     &(sys->pos[0]), &(sys->pos[3]),
		     22, mat_mono);

  matrix_lub_fts_2b_poly (sys,
			  0, 1,
			  &(sys->pos[0]), &(sys->pos[3]),
			  1.0, 1.0,
			  22, mat_poly);

  int check = 0;
  char label [80];
  int j;
  for (i = 0; i < 22; i ++)
    {
      for (j = 0; j < 22; j ++)
	{
	  sprintf (label, "check_matrix_lub_fts_2b_poly: (%d, %d) ", i, j);
	  check += compare (mat_mono [i*22+j], mat_poly [i*22+j], label, verbose, tiny);
	}
    }


  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_matrix_lub_fts_2b_poly : PASSED\n\n");
    }

  stokes_free (sys);

  return (check);
}

/* compare calc_lub_fts_2b_poly() and matrix_lub_fts_2b_poly()
 */
int
check_atimes_matrix_lub_fts_2b_poly (double r, double a1, double a2,
				     int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout, "check_atimes_matrix_lub_fts_2b_poly : start\n");
    }

  struct stokes *sys = NULL;
  sys = stokes_init ();
  sys->version = 2; // FTS

  sys->lubcut = 2.0 + 1.0e-15;
  sys->twobody_nmax = 100;
  sys->twobody_lub = 1; // lub form

  int np = 2;
  stokes_set_np (sys, np, np);

  sys->pos [0] = 0.0;
  sys->pos [1] = 0.0;
  sys->pos [2] = 0.0;
  sys->pos [3] = 0.0;
  sys->pos [4] = 0.0;
  sys->pos [5] = r;
  // so that r12 = r2 - r1 = (0, 0, r)

  double uoe[22];
  double fts_atimes[22];
  double fts_matrix[22];
  double mat[22*22];

  int i;
  for (i = 0; i < 22; i ++)
    {
      uoe [i] = 1.0;
      fts_atimes [i] = 0.0;
      fts_matrix [i] = 0.0;
    }
  for (i = 0; i < 22*22; i ++)
    {
      mat [i] = 0.0;
    }

  // atimes version
  calc_lub_fts_2b_poly (sys,
			&uoe[0], &uoe[11],
			&(sys->pos[0]), &(sys->pos[3]),
			a1, a2,
			&fts_atimes[0], &fts_atimes[11]);

  // matrix version
  matrix_lub_fts_2b_poly (sys,
			  0, 1,
			  &(sys->pos[0]), &(sys->pos[3]),
			  a1, a2,
			  22, mat);
  multiply_extmat_with_extvec_3fts (np, mat, uoe,
				    fts_matrix);

  int check = 0;
  check += compare (fts_atimes [0], fts_matrix [0],  "check_atimes_matrix_lub_fts_2b_poly: Fx1", verbose, tiny);
  check += compare (fts_atimes [1], fts_matrix [1],  "check_atimes_matrix_lub_fts_2b_poly: Fy1", verbose, tiny);
  check += compare (fts_atimes [2], fts_matrix [2],  "check_atimes_matrix_lub_fts_2b_poly: Fz1", verbose, tiny);
  check += compare (fts_atimes [3], fts_matrix [3],  "check_atimes_matrix_lub_fts_2b_poly: Tx1", verbose, tiny);
  check += compare (fts_atimes [4], fts_matrix [4],  "check_atimes_matrix_lub_fts_2b_poly: Ty1", verbose, tiny);
  check += compare (fts_atimes [5], fts_matrix [5],  "check_atimes_matrix_lub_fts_2b_poly: Tz1", verbose, tiny);
  check += compare (fts_atimes [6], fts_matrix [6],  "check_atimes_matrix_lub_fts_2b_poly: Sxx1", verbose, tiny);
  check += compare (fts_atimes [7], fts_matrix [7],  "check_atimes_matrix_lub_fts_2b_poly: Sxy1", verbose, tiny);
  check += compare (fts_atimes [8], fts_matrix [8],  "check_atimes_matrix_lub_fts_2b_poly: Sxz1", verbose, tiny);
  check += compare (fts_atimes [9], fts_matrix [9],  "check_atimes_matrix_lub_fts_2b_poly: Syz1", verbose, tiny);
  check += compare (fts_atimes[10], fts_matrix[10],  "check_atimes_matrix_lub_fts_2b_poly: Syy1", verbose, tiny);
  check += compare (fts_atimes[11], fts_matrix[11],  "check_atimes_matrix_lub_fts_2b_poly: Fx2", verbose, tiny);
  check += compare (fts_atimes[12], fts_matrix[12],  "check_atimes_matrix_lub_fts_2b_poly: Fy2", verbose, tiny);
  check += compare (fts_atimes[13], fts_matrix[13],  "check_atimes_matrix_lub_fts_2b_poly: Fz2", verbose, tiny);
  check += compare (fts_atimes[14], fts_matrix[14],  "check_atimes_matrix_lub_fts_2b_poly: Tx2", verbose, tiny);
  check += compare (fts_atimes[15], fts_matrix[15],  "check_atimes_matrix_lub_fts_2b_poly: Ty2", verbose, tiny);
  check += compare (fts_atimes[16], fts_matrix[16],  "check_atimes_matrix_lub_fts_2b_poly: Tz2", verbose, tiny);
  check += compare (fts_atimes[17], fts_matrix[17],  "check_atimes_matrix_lub_fts_2b_poly: Sxx2", verbose, tiny);
  check += compare (fts_atimes[18], fts_matrix[18],  "check_atimes_matrix_lub_fts_2b_poly: Sxy2", verbose, tiny);
  check += compare (fts_atimes[19], fts_matrix[19],  "check_atimes_matrix_lub_fts_2b_poly: Sxz2", verbose, tiny);
  check += compare (fts_atimes[20], fts_matrix[20],  "check_atimes_matrix_lub_fts_2b_poly: Syz2", verbose, tiny);
  check += compare (fts_atimes[21], fts_matrix[21],  "check_atimes_matrix_lub_fts_2b_poly: Syy2", verbose, tiny);


  if (check == 0 && verbose != 0)
    {
      fprintf (stdout, "check_atimes_matrix_lub_fts_2b_poly : PASSED\n\n");
    }

  stokes_free (sys);

  return (check);
}

