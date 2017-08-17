/* test code for twobody.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: $
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
#include <minv-poly.h> // scalars_minv_fts_poly()
#include "twobody-slip.h"
#include "twobody.h"

#include "memory-check.h"


void
check_twobody_slip (double r,
		    struct twobody_slip_f *f,
		    double a1, double a2)
{
  double *minv = (double *)malloc (sizeof (double) * 44);
  CHECK_MALLOC (minv, "check_twobody_slip");

  if (f->hat_g1 > 0.0 && f->hat_g2 > 0.0)
    {
      scalars_minv_fts_poly (r, f->slip_a1, f->slip_a2, minv);

      fprintf (stdout, "XA11-minv %f %d %.15e\n", r, 0, minv [0]);
      fprintf (stdout, "XA12-minv %f %d %.15e\n", r, 0, minv [1]);
      fprintf (stdout, "YA11-minv %f %d %.15e\n", r, 0, minv [4]);
      fprintf (stdout, "YA12-minv %f %d %.15e\n", r, 0, minv [5]);
      fprintf (stdout, "YB11-minv %f %d %.15e\n", r, 0, minv [8]);
      fprintf (stdout, "YB12-minv %f %d %.15e\n", r, 0, minv [9]);
      fprintf (stdout, "XC11-minv %f %d %.15e\n", r, 0, minv [12]);
      fprintf (stdout, "XC12-minv %f %d %.15e\n", r, 0, minv [13]);
      fprintf (stdout, "YC11-minv %f %d %.15e\n", r, 0, minv [16]);
      fprintf (stdout, "YC12-minv %f %d %.15e\n", r, 0, minv [17]);
      fprintf (stdout, "XG11-minv %f %d %.15e\n", r, 0, minv [20]);
      fprintf (stdout, "XG12-minv %f %d %.15e\n", r, 0, minv [21]);
      fprintf (stdout, "YG11-minv %f %d %.15e\n", r, 0, minv [24]);
      fprintf (stdout, "YG12-minv %f %d %.15e\n", r, 0, minv [25]);
      fprintf (stdout, "YH11-minv %f %d %.15e\n", r, 0, minv [28]);
      fprintf (stdout, "YH12-minv %f %d %.15e\n", r, 0, minv [29]);
      fprintf (stdout, "XM11-minv %f %d %.15e\n", r, 0, minv [32]);
      fprintf (stdout, "XM12-minv %f %d %.15e\n", r, 0, minv [33]);
      fprintf (stdout, "YM11-minv %f %d %.15e\n", r, 0, minv [36]);
      fprintf (stdout, "YM12-minv %f %d %.15e\n", r, 0, minv [37]);
      fprintf (stdout, "ZM11-minv %f %d %.15e\n", r, 0, minv [40]);
      fprintf (stdout, "ZM12-minv %f %d %.15e\n", r, 0, minv [41]);
    }


  double *res12 = (double *)malloc (sizeof (double) * 22);
  CHECK_MALLOC (res12, "check_twobody_slip");

  // currently slip is implemented only up to k=15. (2007/08/17)
  int i;
  for (i = 2; i <= 15; i += 2)
    {
      twobody_slip_scalars_res (2, /* FTS */
				r,
				a1, a2, f,
				i, 0, /* far form */ 1, /* dimensional */
				res12);
      fprintf (stdout, "XA11-res %f %d %.15e\n", r, i, res12 [0]);
      fprintf (stdout, "XA12-res %f %d %.15e\n", r, i, res12 [1]);
      fprintf (stdout, "YA11-res %f %d %.15e\n", r, i, res12 [2]);
      fprintf (stdout, "YA12-res %f %d %.15e\n", r, i, res12 [3]);
      fprintf (stdout, "YB11-res %f %d %.15e\n", r, i, res12 [4]);
      fprintf (stdout, "YB12-res %f %d %.15e\n", r, i, res12 [5]);
      fprintf (stdout, "XC11-res %f %d %.15e\n", r, i, res12 [6]);
      fprintf (stdout, "XC12-res %f %d %.15e\n", r, i, res12 [7]);
      fprintf (stdout, "YC11-res %f %d %.15e\n", r, i, res12 [8]);
      fprintf (stdout, "YC12-res %f %d %.15e\n", r, i, res12 [9]);
      fprintf (stdout, "XG11-res %f %d %.15e\n", r, i, res12 [10]);
      fprintf (stdout, "XG12-res %f %d %.15e\n", r, i, res12 [11]);
      fprintf (stdout, "YG11-res %f %d %.15e\n", r, i, res12 [12]);
      fprintf (stdout, "YG12-res %f %d %.15e\n", r, i, res12 [13]);
      fprintf (stdout, "YH11-res %f %d %.15e\n", r, i, res12 [14]);
      fprintf (stdout, "YH12-res %f %d %.15e\n", r, i, res12 [15]);
      fprintf (stdout, "XM11-res %f %d %.15e\n", r, i, res12 [16]);
      fprintf (stdout, "XM12-res %f %d %.15e\n", r, i, res12 [17]);
      fprintf (stdout, "YM11-res %f %d %.15e\n", r, i, res12 [18]);
      fprintf (stdout, "YM12-res %f %d %.15e\n", r, i, res12 [19]);
      fprintf (stdout, "ZM11-res %f %d %.15e\n", r, i, res12 [20]);
      fprintf (stdout, "ZM12-res %f %d %.15e\n", r, i, res12 [21]);
    }

  free (minv);
  free (res12);
}

void
check_twobody_noslip (double r,
		      double a1, double a2)
{
  double *minv = (double *)malloc (sizeof (double) * 44);
  CHECK_MALLOC (minv, "check_twobody_slip");

  scalars_minv_fts_poly (r, a1, a2, minv);

  fprintf (stdout, "XA11-minv %f %d %.15e\n", r, 0, minv [0]);
  fprintf (stdout, "XA12-minv %f %d %.15e\n", r, 0, minv [1]);
  fprintf (stdout, "YA11-minv %f %d %.15e\n", r, 0, minv [4]);
  fprintf (stdout, "YA12-minv %f %d %.15e\n", r, 0, minv [5]);
  fprintf (stdout, "YB11-minv %f %d %.15e\n", r, 0, minv [8]);
  fprintf (stdout, "YB12-minv %f %d %.15e\n", r, 0, minv [9]);
  fprintf (stdout, "XC11-minv %f %d %.15e\n", r, 0, minv [12]);
  fprintf (stdout, "XC12-minv %f %d %.15e\n", r, 0, minv [13]);
  fprintf (stdout, "YC11-minv %f %d %.15e\n", r, 0, minv [16]);
  fprintf (stdout, "YC12-minv %f %d %.15e\n", r, 0, minv [17]);
  fprintf (stdout, "XG11-minv %f %d %.15e\n", r, 0, minv [20]);
  fprintf (stdout, "XG12-minv %f %d %.15e\n", r, 0, minv [21]);
  fprintf (stdout, "YG11-minv %f %d %.15e\n", r, 0, minv [24]);
  fprintf (stdout, "YG12-minv %f %d %.15e\n", r, 0, minv [25]);
  fprintf (stdout, "YH11-minv %f %d %.15e\n", r, 0, minv [28]);
  fprintf (stdout, "YH12-minv %f %d %.15e\n", r, 0, minv [29]);
  fprintf (stdout, "XM11-minv %f %d %.15e\n", r, 0, minv [32]);
  fprintf (stdout, "XM12-minv %f %d %.15e\n", r, 0, minv [33]);
  fprintf (stdout, "YM11-minv %f %d %.15e\n", r, 0, minv [36]);
  fprintf (stdout, "YM12-minv %f %d %.15e\n", r, 0, minv [37]);
  fprintf (stdout, "ZM11-minv %f %d %.15e\n", r, 0, minv [40]);
  fprintf (stdout, "ZM12-minv %f %d %.15e\n", r, 0, minv [41]);


  double *res12 = (double *)malloc (sizeof (double) * 22);
  CHECK_MALLOC (res12, "check_twobody_slip");

  int i;
  //for (i = 2; i <= 150; i += 2)
  for (i = 2; i <= 15; i += 2)
    {
      twobody_scalars_res (2, /* FTS */
			   r,
			   a1, a2, NULL,
			   i, 0, /* far form */ 1, /* dimensional */
			   res12);
      fprintf (stdout, "XA11-res %f %d %.15e\n", r, i, res12 [0]);
      fprintf (stdout, "XA12-res %f %d %.15e\n", r, i, res12 [1]);
      fprintf (stdout, "YA11-res %f %d %.15e\n", r, i, res12 [2]);
      fprintf (stdout, "YA12-res %f %d %.15e\n", r, i, res12 [3]);
      fprintf (stdout, "YB11-res %f %d %.15e\n", r, i, res12 [4]);
      fprintf (stdout, "YB12-res %f %d %.15e\n", r, i, res12 [5]);
      fprintf (stdout, "XC11-res %f %d %.15e\n", r, i, res12 [6]);
      fprintf (stdout, "XC12-res %f %d %.15e\n", r, i, res12 [7]);
      fprintf (stdout, "YC11-res %f %d %.15e\n", r, i, res12 [8]);
      fprintf (stdout, "YC12-res %f %d %.15e\n", r, i, res12 [9]);
      fprintf (stdout, "XG11-res %f %d %.15e\n", r, i, res12 [10]);
      fprintf (stdout, "XG12-res %f %d %.15e\n", r, i, res12 [11]);
      fprintf (stdout, "YG11-res %f %d %.15e\n", r, i, res12 [12]);
      fprintf (stdout, "YG12-res %f %d %.15e\n", r, i, res12 [13]);
      fprintf (stdout, "YH11-res %f %d %.15e\n", r, i, res12 [14]);
      fprintf (stdout, "YH12-res %f %d %.15e\n", r, i, res12 [15]);
      fprintf (stdout, "XM11-res %f %d %.15e\n", r, i, res12 [16]);
      fprintf (stdout, "XM12-res %f %d %.15e\n", r, i, res12 [17]);
      fprintf (stdout, "YM11-res %f %d %.15e\n", r, i, res12 [18]);
      fprintf (stdout, "YM12-res %f %d %.15e\n", r, i, res12 [19]);
      fprintf (stdout, "ZM11-res %f %d %.15e\n", r, i, res12 [20]);
      fprintf (stdout, "ZM12-res %f %d %.15e\n", r, i, res12 [21]);
    }

  free (minv);
  free (res12);
}


/* main program */
int
main (int argc, char** argv)
{
  double a1 = 1.0;
  double a2 = 1.0;
  // noslip
  /*
  double gamma1 = 0.0;
  double gamma2 = 0.0;
  */
  // perfect slip
  double gamma1 = -1.0;
  double gamma2 = -1.0;

  double hat_g1 = gamma1 / a1;
  double hat_g2 = gamma2 / a2;

  double slip_a1 = a1 / sqrt (1.0 + 2.0 * hat_g1);
  double slip_a2 = a2 / sqrt (1.0 + 2.0 * hat_g2);

  double l = a2 / a1;

  struct twobody_slip_f *f
    = twobody_slip_f_init (15, l, hat_g1, hat_g2,
			   slip_a1, slip_a2);

  int i;
  for (i = 0; i < 10; i ++)
    {
      double r = 2.0 + (double)i / 10.0;
      check_twobody_slip (r, f, a1, a2);
      //check_twobody_noslip (r, a1, a2);
    }

  return 0;
}
