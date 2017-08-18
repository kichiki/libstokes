/* check for twobody-slip.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-twobody-slip.c,v 1.2 2007/12/01 18:31:41 kichiki Exp $
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
#include <twobody.h>
#include <twobody-slip.h>

#include "check.h" // compare()
#include "memory-check.h"


/* compare twobody_slip_far() with twobody_far() with zero slip length.
 */
int
check_twobody_slip_with_noslip (double s, double l,
				int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_twobody_slip_with_noslip (s=%f,l=%f): start\n",
	       s, l);
    }

  int check = 0;
  double max = 0.0;


  double *slip = (double *)malloc (sizeof (double) * 22);
  CHECK_MALLOC (slip, "check_twobody");

  double *noslip = (double *)malloc (sizeof (double) * 22);
  CHECK_MALLOC (noslip, "check_twobody");

  // currently nmax = 15 for slip (2007/08/17)
  twobody_slip_far (2, /* FTS */ 15, l, s, 0.0, 0.0, slip);
  twobody_far      (2, /* FTS */ 15, l, s,           noslip);

  check += compare_max (slip [0], noslip [0], " XA11", verbose, tiny, &max);
  check += compare_max (slip [1], noslip [1], " XA12", verbose, tiny, &max);
  check += compare_max (slip [2], noslip [2], " YA11", verbose, tiny, &max);
  check += compare_max (slip [3], noslip [3], " YA12", verbose, tiny, &max);
  check += compare_max (slip [4], noslip [4], " YB11", verbose, tiny, &max);
  check += compare_max (slip [5], noslip [5], " YB12", verbose, tiny, &max);
  check += compare_max (slip [6], noslip [6], " XC11", verbose, tiny, &max);
  check += compare_max (slip [7], noslip [7], " XC12", verbose, tiny, &max);
  check += compare_max (slip [8], noslip [8], " YC11", verbose, tiny, &max);
  check += compare_max (slip [9], noslip [9], " YC12", verbose, tiny, &max);
  check += compare_max (slip[10], noslip[10], " XG11", verbose, tiny, &max);
  check += compare_max (slip[11], noslip[11], " XG12", verbose, tiny, &max);
  check += compare_max (slip[12], noslip[12], " YG11", verbose, tiny, &max);
  check += compare_max (slip[13], noslip[13], " YG12", verbose, tiny, &max);
  check += compare_max (slip[14], noslip[14], " YH11", verbose, tiny, &max);
  check += compare_max (slip[15], noslip[15], " YH12", verbose, tiny, &max);
  check += compare_max (slip[16], noslip[16], " XM11", verbose, tiny, &max);
  check += compare_max (slip[17], noslip[17], " XM12", verbose, tiny, &max);
  check += compare_max (slip[18], noslip[18], " YM11", verbose, tiny, &max);
  check += compare_max (slip[19], noslip[19], " YM12", verbose, tiny, &max);
  check += compare_max (slip[20], noslip[20], " ZM11", verbose, tiny, &max);
  check += compare_max (slip[21], noslip[21], " ZM12", verbose, tiny, &max);

  free (slip);
  free (noslip);


  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
