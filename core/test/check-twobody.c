/* test code for twobody.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-twobody.c,v 1.3 2007/12/01 18:29:09 kichiki Exp $
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

#include <twobody.h>
#include <two-body-res.h>
#include <bench.h> // ptime_ms_d()

#include "check.h" // compare()


/* check twobody_lub()
 */
int
check_twobody_lub (int version, double r, double a1, double a2, int nmax,
		   int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_twobody_lub : start\n");
    }

  int check = 0;
  double max = 0.0;


  double l = a2 / a1;
  double s = 2.0 * r / (a1 + a2);

  double lub1[22];
  double t0, t;
  t0 = ptime_ms_d ();
  twobody_lub (version, nmax, l, s, lub1);
  t = ptime_ms_d ();
  double ptime_lub1 = t - t0;

  double lub2[22];
  struct twobody_f *f = twobody_f_init (nmax, l);
  t0 = ptime_ms_d ();
  twobody_lub_with_f (version,
		      f, // struct twobody_f
		      nmax, l, s, lub2);
  t = ptime_ms_d ();
  double ptime_lub2 = t - t0;
  twobody_f_free (f);
  // poly gives the results in dimensional form but for a1=1, it's the same

  check += compare_max (lub1 [0], lub2 [0], " XA11", verbose, tiny, &max);
  check += compare_max (lub1 [1], lub2 [1], " XA12", verbose, tiny, &max);
  check += compare_max (lub1 [2], lub2 [2], " YA11", verbose, tiny, &max);
  check += compare_max (lub1 [3], lub2 [3], " YA12", verbose, tiny, &max);
  check += compare_max (lub1 [4], lub2 [4], " YB11", verbose, tiny, &max);
  check += compare_max (lub1 [5], lub2 [5], " YB12", verbose, tiny, &max);
  check += compare_max (lub1 [6], lub2 [6], " XC11", verbose, tiny, &max);
  check += compare_max (lub1 [7], lub2 [7], " XC12", verbose, tiny, &max);
  check += compare_max (lub1 [8], lub2 [8], " YC11", verbose, tiny, &max);
  check += compare_max (lub1 [9], lub2 [9], " YC12", verbose, tiny, &max);
  check += compare_max (lub1[10], lub2[10], " XG11", verbose, tiny, &max);
  check += compare_max (lub1[11], lub2[11], " XG12", verbose, tiny, &max);
  check += compare_max (lub1[12], lub2[12], " YG11", verbose, tiny, &max);
  check += compare_max (lub1[13], lub2[13], " YG12", verbose, tiny, &max);
  check += compare_max (lub1[14], lub2[14], " YH11", verbose, tiny, &max);
  check += compare_max (lub1[15], lub2[15], " YH12", verbose, tiny, &max);
  check += compare_max (lub1[16], lub2[16], " XM11", verbose, tiny, &max);
  check += compare_max (lub1[17], lub2[17], " XM12", verbose, tiny, &max);
  check += compare_max (lub1[18], lub2[18], " YM11", verbose, tiny, &max);
  check += compare_max (lub1[19], lub2[19], " YM12", verbose, tiny, &max);
  check += compare_max (lub1[20], lub2[20], " ZM11", verbose, tiny, &max);
  check += compare_max (lub1[21], lub2[21], " ZM12", verbose, tiny, &max);

  if (verbose != 0)
    {
      fprintf (stdout, " ptime lub1, lub2 = %.3f %.3f, lub1/lub2 = %f\n",
	       ptime_lub1, ptime_lub2, ptime_lub1 / ptime_lub2);

      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}



/* check scalars_minv_f_poly() with scalar_minv_f() for equal sphere
 */
int
check_twobody_scalars_res_with_equal (double r, int nmax,
				      int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_twobody_scalars_res_with_equal: start\n");
    }

  int check = 0;
  double max = 0.0;


  double mono[22];
  double t0, t;
  t0 = ptime_ms_d ();
  scalar_two_body_res (r, mono);
  t = ptime_ms_d ();
  double ptime_mono = t - t0;

  double poly[22];
  t0 = ptime_ms_d ();
  twobody_scalars_res (2, // FTS
		       r, 1.0, 1.0, NULL, nmax,
		       1, // lub form
		       1, // dimensional form
		       poly);
  t = ptime_ms_d ();
  double ptime_poly = t - t0;
  // poly gives the results in dimensional form but for a1=1, it's the same

  check += compare_max (mono [0], poly [0], " XA11", verbose, tiny, &max);
  check += compare_max (mono [1], poly [1], " XA12", verbose, tiny, &max);
  check += compare_max (mono [2], poly [2], " YA11", verbose, tiny, &max);
  check += compare_max (mono [3], poly [3], " YA12", verbose, tiny, &max);
  check += compare_max (mono [4], poly [4], " YB11", verbose, tiny, &max);
  check += compare_max (mono [5], poly [5], " YB12", verbose, tiny, &max);
  check += compare_max (mono [6], poly [6], " XC11", verbose, tiny, &max);
  check += compare_max (mono [7], poly [7], " XC12", verbose, tiny, &max);
  check += compare_max (mono [8], poly [8], " YC11", verbose, tiny, &max);
  check += compare_max (mono [9], poly [9], " YC12", verbose, tiny, &max);
  check += compare_max (mono[10], poly[10], " XG11", verbose, tiny, &max);
  check += compare_max (mono[11], poly[11], " XG12", verbose, tiny, &max);
  check += compare_max (mono[12], poly[12], " YG11", verbose, tiny, &max);
  check += compare_max (mono[13], poly[13], " YG12", verbose, tiny, &max);
  check += compare_max (mono[14], poly[14], " YH11", verbose, tiny, &max);
  check += compare_max (mono[15], poly[15], " YH12", verbose, tiny, &max);
  check += compare_max (mono[16], poly[16], " XM11", verbose, tiny, &max);
  check += compare_max (mono[17], poly[17], " XM12", verbose, tiny, &max);
  check += compare_max (mono[18], poly[18], " YM11", verbose, tiny, &max);
  check += compare_max (mono[19], poly[19], " YM12", verbose, tiny, &max);
  check += compare_max (mono[20], poly[20], " ZM11", verbose, tiny, &max);
  check += compare_max (mono[21], poly[21], " ZM12", verbose, tiny, &max);

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
