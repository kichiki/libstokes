/* utility routines for check
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check.c,v 1.3 2007/12/12 06:30:57 kichiki Exp $
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


/* compare x and y
 */
int
compare (double x, double y, char *label,
	 int verbose, double tiny)
{
  int check = 0;

  double d;
  if (fabs (x) > tiny)
    {
      // see the relative error
      d = fabs ((x - y) / x);
    }
  else if (fabs (y) > tiny)
    {
      // see the relative error
      d = fabs ((x - y) / y);
    }
  else
    {
      // see the absolute error
      d = fabs (x - y);
    }

  if (d > tiny)
    {
      if (verbose != 0)
	{
	  fprintf (stdout, "%s %e %e %e\n",
		   label, x, y, d);
	}
      check ++;
    }

  return (check);
}

/* compare x and y and keep the max error
 */
int
compare_max (double x, double y, char *label,
	     int verbose, double tiny,
	     double *max)
{
  int check = 0;

  double d;
  if (fabs (x) > tiny)
    {
      // see the relative error
      d = fabs ((x - y) / x);
    }
  else if (fabs (y) > tiny)
    {
      // see the relative error
      d = fabs ((x - y) / y);
    }
  else
    {
      // see the absolute error
      d = fabs (x - y);
    }

  if (d > *max) *max = d;

  if (d > tiny)
    {
      if (verbose != 0)
	{
	  fprintf (stdout, "%s %e %e %e\n",
		   label, x, y, d);
	}
      check ++;
    }

  return (check);
}
