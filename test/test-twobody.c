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
#include "twobody.h"
#include "two-body-res.h"

#include "memory-check.h"


void
check_twobody (double r)
{
  double l = 1.0;
  double s = 2.0 * r / (1.0 + l);

  double *new = (double *)malloc (sizeof (double) * 22);
  CHECK_MALLOC (new, "check_twobody");

  double *lub = (double *)malloc (sizeof (double) * 22);
  CHECK_MALLOC (lub, "check_twobody");

  double *old = (double *)malloc (sizeof (double) * 22);
  CHECK_MALLOC (old, "check_twobody");

  scalar_two_body_res (r, old);
  fprintf (stdout, "XA11-old %f %d %.15e\n", s, 0, old [0]);
  fprintf (stdout, "XA12-old %f %d %.15e\n", s, 0, old [1]);
  fprintf (stdout, "YA11-old %f %d %.15e\n", s, 0, old [2]);
  fprintf (stdout, "YA12-old %f %d %.15e\n", s, 0, old [3]);
  fprintf (stdout, "YB11-old %f %d %.15e\n", s, 0, old [4]);
  fprintf (stdout, "YB12-old %f %d %.15e\n", s, 0, old [5]);
  fprintf (stdout, "XC11-old %f %d %.15e\n", s, 0, old [6]);
  fprintf (stdout, "XC12-old %f %d %.15e\n", s, 0, old [7]);
  fprintf (stdout, "YC11-old %f %d %.15e\n", s, 0, old [8]);
  fprintf (stdout, "YC12-old %f %d %.15e\n", s, 0, old [9]);
  fprintf (stdout, "XG11-old %f %d %.15e\n", s, 0, old [10]);
  fprintf (stdout, "XG12-old %f %d %.15e\n", s, 0, old [11]);
  fprintf (stdout, "YG11-old %f %d %.15e\n", s, 0, old [12]);
  fprintf (stdout, "YG12-old %f %d %.15e\n", s, 0, old [13]);
  fprintf (stdout, "YH11-old %f %d %.15e\n", s, 0, old [14]);
  fprintf (stdout, "YH12-old %f %d %.15e\n", s, 0, old [15]);
  fprintf (stdout, "XM11-old %f %d %.15e\n", s, 0, old [16]);
  fprintf (stdout, "XM12-old %f %d %.15e\n", s, 0, old [17]);
  fprintf (stdout, "YM11-old %f %d %.15e\n", s, 0, old [18]);
  fprintf (stdout, "YM12-old %f %d %.15e\n", s, 0, old [19]);
  fprintf (stdout, "ZM11-old %f %d %.15e\n", s, 0, old [20]);
  fprintf (stdout, "ZM12-old %f %d %.15e\n", s, 0, old [21]);



  int i;
  for (i = 2; i <= 150; i += 2)
    {
      twobody_far (2, /* FTS */ i, l, s, new);
      // scale from Jeffrey-Onishi to Stokesian dynamics
      //twobody_scale_SD (2, new, l);
      twobody_scale (2, new, 1.0, l); // dimensional

      fprintf (stdout, "XA11 %f %d %.15e\n", s, i, new [0]);
      fprintf (stdout, "XA12 %f %d %.15e\n", s, i, new [1]);
      fprintf (stdout, "YA11 %f %d %.15e\n", s, i, new [2]);
      fprintf (stdout, "YA12 %f %d %.15e\n", s, i, new [3]);
      fprintf (stdout, "YB11 %f %d %.15e\n", s, i, new [4]);
      fprintf (stdout, "YB12 %f %d %.15e\n", s, i, new [5]);
      fprintf (stdout, "XC11 %f %d %.15e\n", s, i, new [6]);
      fprintf (stdout, "XC12 %f %d %.15e\n", s, i, new [7]);
      fprintf (stdout, "YC11 %f %d %.15e\n", s, i, new [8]);
      fprintf (stdout, "YC12 %f %d %.15e\n", s, i, new [9]);
      fprintf (stdout, "XG11 %f %d %.15e\n", s, i, new [10]);
      fprintf (stdout, "XG12 %f %d %.15e\n", s, i, new [11]);
      fprintf (stdout, "YG11 %f %d %.15e\n", s, i, new [12]);
      fprintf (stdout, "YG12 %f %d %.15e\n", s, i, new [13]);
      fprintf (stdout, "YH11 %f %d %.15e\n", s, i, new [14]);
      fprintf (stdout, "YH12 %f %d %.15e\n", s, i, new [15]);
      fprintf (stdout, "XM11 %f %d %.15e\n", s, i, new [16]);
      fprintf (stdout, "XM12 %f %d %.15e\n", s, i, new [17]);
      fprintf (stdout, "YM11 %f %d %.15e\n", s, i, new [18]);
      fprintf (stdout, "YM12 %f %d %.15e\n", s, i, new [19]);
      fprintf (stdout, "ZM11 %f %d %.15e\n", s, i, new [20]);
      fprintf (stdout, "ZM12 %f %d %.15e\n", s, i, new [21]);


      twobody_lub (2, /* FTS */ i, l, s, lub);
      // scale from Jeffrey-Onishi to Stokesian dynamics
      //twobody_scale_SD (2, new, l);
      twobody_scale (2, new, 1.0, l); // dimensional

      fprintf (stdout, "XA11-lub %f %d %.15e\n", s, i, lub [0]);
      fprintf (stdout, "XA12-lub %f %d %.15e\n", s, i, lub [1]);
      fprintf (stdout, "YA11-lub %f %d %.15e\n", s, i, lub [2]);
      fprintf (stdout, "YA12-lub %f %d %.15e\n", s, i, lub [3]);
      fprintf (stdout, "YB11-lub %f %d %.15e\n", s, i, lub [4]);
      fprintf (stdout, "YB12-lub %f %d %.15e\n", s, i, lub [5]);
      fprintf (stdout, "XC11-lub %f %d %.15e\n", s, i, lub [6]);
      fprintf (stdout, "XC12-lub %f %d %.15e\n", s, i, lub [7]);
      fprintf (stdout, "YC11-lub %f %d %.15e\n", s, i, lub [8]);
      fprintf (stdout, "YC12-lub %f %d %.15e\n", s, i, lub [9]);
      fprintf (stdout, "XG11-lub %f %d %.15e\n", s, i, lub [10]);
      fprintf (stdout, "XG12-lub %f %d %.15e\n", s, i, lub [11]);
      fprintf (stdout, "YG11-lub %f %d %.15e\n", s, i, lub [12]);
      fprintf (stdout, "YG12-lub %f %d %.15e\n", s, i, lub [13]);
      fprintf (stdout, "YH11-lub %f %d %.15e\n", s, i, lub [14]);
      fprintf (stdout, "YH12-lub %f %d %.15e\n", s, i, lub [15]);
      fprintf (stdout, "XM11-lub %f %d %.15e\n", s, i, lub [16]);
      fprintf (stdout, "XM12-lub %f %d %.15e\n", s, i, lub [17]);
      fprintf (stdout, "YM11-lub %f %d %.15e\n", s, i, lub [18]);
      fprintf (stdout, "YM12-lub %f %d %.15e\n", s, i, lub [19]);
      fprintf (stdout, "ZM11-lub %f %d %.15e\n", s, i, lub [20]);
      fprintf (stdout, "ZM12-lub %f %d %.15e\n", s, i, lub [21]);
    }


  free (new);
  free (lub);
  free (old);
}


/* main program */
int
main (int argc, char** argv)
{
  check_twobody (2.1);

  return 0;
}
