/* test code for ode-quaternion.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-ode-quaternion.c,v 1.2 2007/12/01 18:31:18 kichiki Exp $
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
#include <ode-quaternion.h>


#include "check.h" // compare()


/* check quaternion_Winv() with quaternion_Winv_lapack()
 * INPUT
 *  verbose : if non-zero, print results
 *  tiny    : small number for check
 * OUTPUT
 *  (returned value) : 0 => passed
 *                     otherwise => failed
 */
int
check_quaternion_Winv (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_quaternion_Winv : start\n");
    }

  int check = 0;
  double max = 0.0;


  /* this is not proper quaternion, but here
   * it is to test the inversion
   */
  double Q[4];
  Q[0] = 1.0;
  Q[1] = 2.0;
  Q[2] = 3.0;
  Q[3] = 4.0;

  double W1[16];
  quaternion_Winv_lapack (Q, W1);


  double W2[16];
  quaternion_Winv (Q, W2);

  int i;
  for (i = 0; i < 16; i ++)
    {
      char label[80];
      sprintf (label, " Winv[%d] :", i);
      check += compare_max (W1[i], W2[i], label, verbose, tiny, &max);
    }

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
