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
#include <bench.h> // ptime_ms_d()


#include "check.h" // compare()


/* check calc_lub_fts_2b_poly() with a1=a2=a
 * comparing with calc_lub_fts_2b()
 */
int
check_lub_fts_2b_poly (double r, int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_lub_fts_2b_poly(r=%f) : start\n",
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
  double fts_mono[22];
  double fts_poly[22];

  int i;
  for (i = 0; i < 22; i ++)
    {
      uoe [i] = 1.0;
      fts_mono [i] = 0.0;
      fts_poly [i] = 0.0;
    }

  double t0, t;
  t0 = ptime_ms_d ();
  calc_lub_fts_2b (sys,
		   &uoe[0], &uoe[11],
		   &(sys->pos[0]), &(sys->pos[3]),
		   &fts_mono[0], &fts_mono[11]);
  t = ptime_ms_d ();
  double ptime_mono = t - t0;


  // set sys->a to test poly version
  double a [2] = {1.0, 1.0};
  stokes_set_radius (sys, a);

  t0 = ptime_ms_d ();
  calc_lub_fts_2b_poly (sys,
			&uoe[0], &uoe[11],
			&(sys->pos[0]), &(sys->pos[3]),
			//1.0, 1.0,
			0, 1,
			&fts_poly[0], &fts_poly[11]);
  t = ptime_ms_d ();
  double ptime_poly = t - t0;

  stokes_free (sys);


  int check = 0;
  double max = 0.0;

  check += compare_max (fts_mono [0], fts_poly [0], " Fx1", verbose, tiny, &max);
  check += compare_max (fts_mono [1], fts_poly [1], " Fy1", verbose, tiny, &max);
  check += compare_max (fts_mono [2], fts_poly [2], " Fz1", verbose, tiny, &max);
  check += compare_max (fts_mono [3], fts_poly [3], " Tx1", verbose, tiny, &max);
  check += compare_max (fts_mono [4], fts_poly [4], " Ty1", verbose, tiny, &max);
  check += compare_max (fts_mono [5], fts_poly [5], " Tz1", verbose, tiny, &max);
  check += compare_max (fts_mono [6], fts_poly [6], " Sxx1", verbose, tiny, &max);
  check += compare_max (fts_mono [7], fts_poly [7], " Sxy1", verbose, tiny, &max);
  check += compare_max (fts_mono [8], fts_poly [8], " Sxz1", verbose, tiny, &max);
  check += compare_max (fts_mono [9], fts_poly [9], " Syz1", verbose, tiny, &max);
  check += compare_max (fts_mono[10], fts_poly[10], " Syy1", verbose, tiny, &max);
  check += compare_max (fts_mono[11], fts_poly[11], " Fx2", verbose, tiny, &max);
  check += compare_max (fts_mono[12], fts_poly[12], " Fy2", verbose, tiny, &max);
  check += compare_max (fts_mono[13], fts_poly[13], " Fz2", verbose, tiny, &max);
  check += compare_max (fts_mono[14], fts_poly[14], " Tx2", verbose, tiny, &max);
  check += compare_max (fts_mono[15], fts_poly[15], " Ty2", verbose, tiny, &max);
  check += compare_max (fts_mono[16], fts_poly[16], " Tz2", verbose, tiny, &max);
  check += compare_max (fts_mono[17], fts_poly[17], " Sxx2", verbose, tiny, &max);
  check += compare_max (fts_mono[18], fts_poly[18], " Sxy2", verbose, tiny, &max);
  check += compare_max (fts_mono[19], fts_poly[19], " Sxz2", verbose, tiny, &max);
  check += compare_max (fts_mono[20], fts_poly[20], " Syz2", verbose, tiny, &max);
  check += compare_max (fts_mono[21], fts_poly[21], " Syy2", verbose, tiny, &max);

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
