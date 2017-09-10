/* test code for lubrication for polydisperse systems
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
#include <stokes.h> // struct stokes
#include <fts.h>
#include <matrix.h> // multiply_extmat_with_extvec_3fts()
#include <two-body-res.h> // scalar_two_body_res ()
#include <minv-poly.h> // scalars_lub_poly_full() scalars_res_poly_scale_SD()
#include <bench.h> // ptime_ms_d()


#include "check.h" // compare()


/* compare lub scalar functions
 *  scalarwo-body_res() + scalar_minv_fts() for mono
 *  and scalars_lub_poly_full() for poly
 * for equal-sphere case
 * INPUT
 *  r       := x_2 - x_1
 *  nmax    : for scalars_lub_poly_full() in twobody.c
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_lub_scalars_poly (double r, int nmax,
			int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_lub_scalars_poly(r=%f,nmax=%d) : start\n",
	       r, nmax);
    }

  int check = 0;
  double max = 0.0;


  double poly[44];
  double t0, t;
  t0 = ptime_ms_d ();
  scalars_lub_poly_full (2, // FTS
			 r, 1.0, 1.0, NULL, NULL,
			 nmax, 1, /* lub form */
			 poly);
  scalars_res_poly_scale_SD (2, // FTS version
			     1.0, 1.0, poly);
  t = ptime_ms_d ();
  double ptime_poly = t - t0;

  double res[22];
  double minv[22];
  t0 = ptime_ms_d ();
  scalar_two_body_res (r, res);
  scalar_minv_fts (r, minv);
  t = ptime_ms_d ();
  double ptime_mono = t - t0;

  check += compare_max (poly [0], res [0]-minv [0], " xa11", verbose, tiny, &max);
  check += compare_max (poly [1], res [1]-minv [1], " xa12", verbose, tiny, &max);
  check += compare_max (poly [4], res [2]-minv [2], " ya11", verbose, tiny, &max);
  check += compare_max (poly [5], res [3]-minv [3], " ya12", verbose, tiny, &max);
  check += compare_max (poly [8], res [4]-minv [4], " yb11", verbose, tiny, &max);
  check += compare_max (poly [9], res [5]-minv [5], " yb12", verbose, tiny, &max);
  check += compare_max (poly[12], res [6]-minv [6], " xc11", verbose, tiny, &max);
  check += compare_max (poly[13], res [7]-minv [7], " xc12", verbose, tiny, &max);
  check += compare_max (poly[16], res [8]-minv [8], " yc11", verbose, tiny, &max);
  check += compare_max (poly[17], res [9]-minv [9], " yc12", verbose, tiny, &max);
  check += compare_max (poly[20], res[10]-minv[10], " xg11", verbose, tiny, &max);
  check += compare_max (poly[21], res[11]-minv[11], " xg12", verbose, tiny, &max);
  check += compare_max (poly[24], res[12]-minv[12], " yg11", verbose, tiny, &max);
  check += compare_max (poly[25], res[13]-minv[13], " yg12", verbose, tiny, &max);
  check += compare_max (poly[28], res[14]-minv[14], " yh11", verbose, tiny, &max);
  check += compare_max (poly[29], res[15]-minv[15], " yh12", verbose, tiny, &max);
  check += compare_max (poly[32], res[16]-minv[16], " xm11", verbose, tiny, &max);
  check += compare_max (poly[33], res[17]-minv[17], " xm12", verbose, tiny, &max);
  check += compare_max (poly[36], res[18]-minv[18], " ym11", verbose, tiny, &max);
  check += compare_max (poly[37], res[19]-minv[19], " ym12", verbose, tiny, &max);
  check += compare_max (poly[40], res[20]-minv[20], " zm11", verbose, tiny, &max);
  check += compare_max (poly[41], res[21]-minv[21], " zm12", verbose, tiny, &max);

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


/* check matrix_lub_fts_2b_poly() with a1=a2=a
 * comparing with matrix_lub_fts_2b()
 */
int
check_matrix_lub_fts_2b_poly (double r, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_matrix_lub_fts_2b_poly(r=%f) : start\n",
	       r);
    }

  struct stokes *sys = NULL;
  sys = stokes_init ();
  sys->version = 2; // FTS

  //sys->lubmin = 2.0 + 1.0e-15;
  sys->lubmin2 = 2.0 + 1.0e-15;
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

  double t0, t;
  t0 = ptime_ms_d ();
  matrix_lub_fts_2b (sys,
		     0, 1,
		     &(sys->pos[0]), &(sys->pos[3]),
		     22, mat_mono);
  t = ptime_ms_d ();
  double ptime_mono = t - t0;


  // set sys->a to test poly version
  double a [2] = {1.0, 1.0};
  stokes_set_radius (sys, a);

  t0 = ptime_ms_d ();
  matrix_lub_fts_2b_poly (sys,
			  0, 1,
			  &(sys->pos[0]), &(sys->pos[3]),
			  //1.0, 1.0,
			  0, 1,
			  22, mat_poly);
  t = ptime_ms_d ();
  double ptime_poly = t - t0;

  stokes_free (sys);


  int check = 0;
  double max = 0.0;
  char label [80];
  int j;
  for (i = 0; i < 22; i ++)
    {
      for (j = 0; j < 22; j ++)
	{
	  sprintf (label, " (%d, %d) ", i, j);
	  check += compare_max (mat_mono [i*22+j], mat_poly [i*22+j],
				label, verbose, tiny, &max);
	}
    }


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


/* compare calc_lub_fts_2b_poly() and matrix_lub_fts_2b_poly()
 */
int
check_atimes_matrix_lub_fts_2b_poly (double r, double a1, double a2,
				     int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_atimes_matrix_lub_fts_2b_poly\n"
	       "(r=%f,a1=%f,a2=%f) : start\n",
	       r, a1, a2);
    }

  struct stokes *sys = NULL;
  sys = stokes_init ();
  sys->version = 2; // FTS

  //sys->lubmin = 2.0 + 1.0e-15;
  sys->lubmin2 = 2.0 + 1.0e-15;
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

  // set sys->a to test poly version
  double a [2] = {1.0, 1.0};
  stokes_set_radius (sys, a);

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

  stokes_free (sys);


  int check = 0;
  double max = 0.0;
  check += compare_max (fts_atimes [0], fts_matrix [0], " Fx1", verbose, tiny, &max);
  check += compare_max (fts_atimes [1], fts_matrix [1], " Fy1", verbose, tiny, &max);
  check += compare_max (fts_atimes [2], fts_matrix [2], " Fz1", verbose, tiny, &max);
  check += compare_max (fts_atimes [3], fts_matrix [3], " Tx1", verbose, tiny, &max);
  check += compare_max (fts_atimes [4], fts_matrix [4], " Ty1", verbose, tiny, &max);
  check += compare_max (fts_atimes [5], fts_matrix [5], " Tz1", verbose, tiny, &max);
  check += compare_max (fts_atimes [6], fts_matrix [6], " Sxx1", verbose, tiny, &max);
  check += compare_max (fts_atimes [7], fts_matrix [7], " Sxy1", verbose, tiny, &max);
  check += compare_max (fts_atimes [8], fts_matrix [8], " Sxz1", verbose, tiny, &max);
  check += compare_max (fts_atimes [9], fts_matrix [9], " Syz1", verbose, tiny, &max);
  check += compare_max (fts_atimes[10], fts_matrix[10], " Syy1", verbose, tiny, &max);
  check += compare_max (fts_atimes[11], fts_matrix[11], " Fx2", verbose, tiny, &max);
  check += compare_max (fts_atimes[12], fts_matrix[12], " Fy2", verbose, tiny, &max);
  check += compare_max (fts_atimes[13], fts_matrix[13], " Fz2", verbose, tiny, &max);
  check += compare_max (fts_atimes[14], fts_matrix[14], " Tx2", verbose, tiny, &max);
  check += compare_max (fts_atimes[15], fts_matrix[15], " Ty2", verbose, tiny, &max);
  check += compare_max (fts_atimes[16], fts_matrix[16], " Tz2", verbose, tiny, &max);
  check += compare_max (fts_atimes[17], fts_matrix[17], " Sxx2", verbose, tiny, &max);
  check += compare_max (fts_atimes[18], fts_matrix[18], " Sxy2", verbose, tiny, &max);
  check += compare_max (fts_atimes[19], fts_matrix[19], " Sxz2", verbose, tiny, &max);
  check += compare_max (fts_atimes[20], fts_matrix[20], " Syz2", verbose, tiny, &max);
  check += compare_max (fts_atimes[21], fts_matrix[21], " Syy2", verbose, tiny, &max);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}

