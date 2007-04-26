/* test code for twobody.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-twobody.c,v 1.2 2007/04/26 05:19:13 kichiki Exp $
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

  check += compare (lub1 [0], lub2 [0], "check_twobody_lub: XA11", verbose, tiny);
  check += compare (lub1 [1], lub2 [1], "check_twobody_lub: XA12", verbose, tiny);
  check += compare (lub1 [2], lub2 [2], "check_twobody_lub: YA11", verbose, tiny);
  check += compare (lub1 [3], lub2 [3], "check_twobody_lub: YA12", verbose, tiny);
  check += compare (lub1 [4], lub2 [4], "check_twobody_lub: YB11", verbose, tiny);
  check += compare (lub1 [5], lub2 [5], "check_twobody_lub: YB12", verbose, tiny);
  check += compare (lub1 [6], lub2 [6], "check_twobody_lub: XC11", verbose, tiny);
  check += compare (lub1 [7], lub2 [7], "check_twobody_lub: XC12", verbose, tiny);
  check += compare (lub1 [8], lub2 [8], "check_twobody_lub: YC11", verbose, tiny);
  check += compare (lub1 [9], lub2 [9], "check_twobody_lub: YC12", verbose, tiny);
  check += compare (lub1[10], lub2[10], "check_twobody_lub: XG11", verbose, tiny);
  check += compare (lub1[11], lub2[11], "check_twobody_lub: XG12", verbose, tiny);
  check += compare (lub1[12], lub2[12], "check_twobody_lub: YG11", verbose, tiny);
  check += compare (lub1[13], lub2[13], "check_twobody_lub: YG12", verbose, tiny);
  check += compare (lub1[14], lub2[14], "check_twobody_lub: YH11", verbose, tiny);
  check += compare (lub1[15], lub2[15], "check_twobody_lub: YH12", verbose, tiny);
  check += compare (lub1[16], lub2[16], "check_twobody_lub: XM11", verbose, tiny);
  check += compare (lub1[17], lub2[17], "check_twobody_lub: XM12", verbose, tiny);
  check += compare (lub1[18], lub2[18], "check_twobody_lub: YM11", verbose, tiny);
  check += compare (lub1[19], lub2[19], "check_twobody_lub: YM12", verbose, tiny);
  check += compare (lub1[20], lub2[20], "check_twobody_lub: ZM11", verbose, tiny);
  check += compare (lub1[21], lub2[21], "check_twobody_lub: ZM12", verbose, tiny);

  if (verbose != 0)
    {
      fprintf (stdout, "check_twobody_lub :"
	       " ptime lub1, lub2 = %.3f %.3f, lub1/lub2 = %f\n",
	       ptime_lub1, ptime_lub2, ptime_lub1 / ptime_lub2);

      if (check == 0)
	fprintf (stdout, "check_twobody_lub :"
		 " PASSED\n\n");
      else
	fprintf (stdout, "check_twobody_lub :"
		 " FAILED\n\n");
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

  check += compare (mono [0], poly [0], "check_twobody_scalars_res_with_equal: XA11", verbose, tiny);
  check += compare (mono [1], poly [1], "check_twobody_scalars_res_with_equal: XA12", verbose, tiny);
  check += compare (mono [2], poly [2], "check_twobody_scalars_res_with_equal: YA11", verbose, tiny);
  check += compare (mono [3], poly [3], "check_twobody_scalars_res_with_equal: YA12", verbose, tiny);
  check += compare (mono [4], poly [4], "check_twobody_scalars_res_with_equal: YB11", verbose, tiny);
  check += compare (mono [5], poly [5], "check_twobody_scalars_res_with_equal: YB12", verbose, tiny);
  check += compare (mono [6], poly [6], "check_twobody_scalars_res_with_equal: XC11", verbose, tiny);
  check += compare (mono [7], poly [7], "check_twobody_scalars_res_with_equal: XC12", verbose, tiny);
  check += compare (mono [8], poly [8], "check_twobody_scalars_res_with_equal: YC11", verbose, tiny);
  check += compare (mono [9], poly [9], "check_twobody_scalars_res_with_equal: YC12", verbose, tiny);
  check += compare (mono[10], poly[10], "check_twobody_scalars_res_with_equal: XG11", verbose, tiny);
  check += compare (mono[11], poly[11], "check_twobody_scalars_res_with_equal: XG12", verbose, tiny);
  check += compare (mono[12], poly[12], "check_twobody_scalars_res_with_equal: YG11", verbose, tiny);
  check += compare (mono[13], poly[13], "check_twobody_scalars_res_with_equal: YG12", verbose, tiny);
  check += compare (mono[14], poly[14], "check_twobody_scalars_res_with_equal: YH11", verbose, tiny);
  check += compare (mono[15], poly[15], "check_twobody_scalars_res_with_equal: YH12", verbose, tiny);
  check += compare (mono[16], poly[16], "check_twobody_scalars_res_with_equal: XM11", verbose, tiny);
  check += compare (mono[17], poly[17], "check_twobody_scalars_res_with_equal: XM12", verbose, tiny);
  check += compare (mono[18], poly[18], "check_twobody_scalars_res_with_equal: YM11", verbose, tiny);
  check += compare (mono[19], poly[19], "check_twobody_scalars_res_with_equal: YM12", verbose, tiny);
  check += compare (mono[20], poly[20], "check_twobody_scalars_res_with_equal: ZM11", verbose, tiny);
  check += compare (mono[21], poly[21], "check_twobody_scalars_res_with_equal: ZM12", verbose, tiny);

  if (verbose != 0)
    {
      fprintf (stdout, "check_twobody_scalars_res_with_equal:"
	       " ptime mono, poly = %.3f %.3f, poly/mono = %f\n",
	       ptime_mono, ptime_poly, ptime_poly / ptime_mono);

      if (check == 0)
	fprintf (stdout, "check_twobody_scalars_res_with_equal:"
		 " PASSED\n\n");
      else
	fprintf (stdout, "check_twobody_scalars_res_with_equal:"
		 " FAILED\n\n");
    }

  return (check);
}
