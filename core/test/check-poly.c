/* test code for polydisperse handling for non-periodic systems
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-poly.c,v 1.5 2007/12/01 18:29:32 kichiki Exp $
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
#include <math.h> // fabs()
#include "memory-check.h" // macro CHECK_MALLOC

#include <stokes.h> // struct stokes
#include <non-ewald.h> // scalars_nonewald_poly(), scalars_nonewald()
#include <ewald.h> // scalars_ewald_real_poly(), scalars_ewald_real()
#include <bench.h> // ptime_ms_d()

#include "check.h" // compare()


/* compare scalar functions for two routines
 *  scalars_nonewald() and scalars_nonewald_poly()
 * for equal-sphere case
 * INPUT
 *  r       := x_2 - x_1
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_scalars_nonewald_poly (double r,
			     int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_scalars_nonewald_poly : start\n");
    }

  int check = 0;
  double max = 0.0;


  double poly[11];
  double t0, t;
  t0 = ptime_ms_d ();
  scalars_nonewald_poly (2, /* FTS */ r, 1.0, 1.0, poly);
  t = ptime_ms_d ();
  double ptime_poly = t - t0;

  double mono[11];
  t0 = ptime_ms_d ();
  scalars_nonewald (2, /* FTS */ r, mono);
  t = ptime_ms_d ();
  double ptime_mono = t - t0;

  check += compare_max (poly [0], mono [0], " xa", verbose, tiny, &max);
  check += compare_max (poly [1], mono [1], " ya", verbose, tiny, &max);
  check += compare_max (poly [2], mono [2], " yb", verbose, tiny, &max);
  check += compare_max (poly [3], mono [3], " xc", verbose, tiny, &max);
  check += compare_max (poly [4], mono [4], " yc", verbose, tiny, &max);
  check += compare_max (poly [5], mono [5], " xg", verbose, tiny, &max);
  check += compare_max (poly [6], mono [6], " yg", verbose, tiny, &max);
  check += compare_max (poly [7], mono [7], " yh", verbose, tiny, &max);
  check += compare_max (poly [8], mono [8], " xm", verbose, tiny, &max);
  check += compare_max (poly [9], mono [9], " ym", verbose, tiny, &max);
  check += compare_max (poly[10], mono[10], " zm", verbose, tiny, &max);

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


/* check the symmetry between (ab) and (ba) for M^inf
 * INPUT
 *  r       := x_2 - x_1
 *  a1, a2  : radius of particle 1 and 2
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_scalars_nonewald_poly_symmetry (double r, double a1, double a2,
				      int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_scalars_nonewald_poly_symmetry : start\n");
    }

  int check = 0;
  double max = 0.0;


  // for (12)-interaction
  // scalar[] is in the dimensional form (not the SD scaling)
  double scalar[11];
  scalars_nonewald_poly (2, // FTS version
			 r, a1, a2, scalar);

  double xa12, ya12;
  double yb12;
  double xc12, yc12;
  double xg12, yg12;
  double yh12;
  double xm12, ym12, zm12;
  xa12 = scalar [0];
  ya12 = scalar [1];
  yb12 = scalar [2];
  xc12 = scalar [3];
  yc12 = scalar [4];
  xg12 = scalar [5];
  yg12 = scalar [6];
  yh12 = scalar [7];
  xm12 = scalar [8];
  ym12 = scalar [9];
  zm12 = scalar[10];


  // for (21)-interaction
  // scalar[] is in the dimensional form (not the SD scaling)
  scalars_nonewald_poly (2, // FTS version
			 r, a2, a1, scalar);

  double xa21, ya21;
  double yb21;
  double xc21, yc21;
  double xg21, yg21;
  double yh21;
  double xm21, ym21, zm21;
  xa21 = scalar [0];
  ya21 = scalar [1];
  yb21 = scalar [2];
  xc21 = scalar [3];
  yc21 = scalar [4];
  xg21 = scalar [5];
  yg21 = scalar [6];
  yh21 = scalar [7];
  xm21 = scalar [8];
  ym21 = scalar [9];
  zm21 = scalar[10];


  // compare (12) and (21)
  check += compare_max (xa12, xa21, " xa12", verbose, tiny, &max);
  check += compare_max (ya12, ya21, " ya12", verbose, tiny, &max);
  check += compare_max (xc12, xc21, " xc12", verbose, tiny, &max);
  check += compare_max (yc12, yc21, " yc12", verbose, tiny, &max);
  check += compare_max (xm12, xm21, " xm12", verbose, tiny, &max);
  check += compare_max (ym12, ym21, " ym12", verbose, tiny, &max);
  check += compare_max (zm12, zm21, " zm12", verbose, tiny, &max);

  /* because of the finite size effects in G part,
   * [xy]g12 != [xy]g21 for unequal spheres in general.
   * note that for other parts such as b and h,
   * mobility functions for 12 is equal to 21 even for unequal spheres
   * because there is no finite size effects in M^inf.
   */
  //check += compare_max (xg12, xg21, " xg12", verbose, tiny, &max);
  //check += compare_max (yg12, yg21, " yg12", verbose, tiny, &max);
  check += compare_max (yb12, yb21, " yb12", verbose, tiny, &max);
  check += compare_max (yh12, yh21, " yh12", verbose, tiny, &max);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}



/* compare scalar functions for two routines
 *  scalars_ewald_real() and scalars_ewald_real_poly()
 * for equal-sphere case
 * INPUT
 *  r       := x_2 - x_1
 *  xi      : splitting parameter for ewald summation
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_scalars_ewald_real_poly (double r, double xi,
			       int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_scalars_ewald_real_poly : start\n");
    }

  int check = 0;
  double max = 0.0;


  double poly[11];
  double t0, t;
  t0 = ptime_ms_d ();
  scalars_ewald_real_poly (2, // FTS
			   xi, r, 1.0, 1.0,
			   poly + 0,
			   poly + 1,
			   poly + 2,
			   poly + 3,
			   poly + 4,
			   poly + 5,
			   poly + 6,
			   poly + 7,
			   poly + 8,
			   poly + 9,
			   poly +10);
  t = ptime_ms_d ();
  double ptime_poly = t - t0;

  double mono[11];
  t0 = ptime_ms_d ();
  scalars_ewald_real (2, // FTS
		      xi, r,
		      mono + 0,
		      mono + 1,
		      mono + 2,
		      mono + 3,
		      mono + 4,
		      mono + 5,
		      mono + 6,
		      mono + 7,
		      mono + 8,
		      mono + 9,
		      mono +10);
  t = ptime_ms_d ();
  double ptime_mono = t - t0;

  check += compare_max (poly [0], mono [0], " xa", verbose, tiny, &max);
  check += compare_max (poly [1], mono [1], " ya", verbose, tiny, &max);
  check += compare_max (poly [2], mono [2], " yb", verbose, tiny, &max);
  check += compare_max (poly [3], mono [3], " xc", verbose, tiny, &max);
  check += compare_max (poly [4], mono [4], " yc", verbose, tiny, &max);
  check += compare_max (poly [5], mono [5], " xg", verbose, tiny, &max);
  check += compare_max (poly [6], mono [6], " yg", verbose, tiny, &max);
  check += compare_max (poly [7], mono [7], " yh", verbose, tiny, &max);
  check += compare_max (poly [8], mono [8], " xm", verbose, tiny, &max);
  check += compare_max (poly [9], mono [9], " ym", verbose, tiny, &max);
  check += compare_max (poly[10], mono[10], " zm", verbose, tiny, &max);

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
