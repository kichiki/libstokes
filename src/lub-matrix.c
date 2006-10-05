/* lubrication routines -- MATRIX procedure
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: lub-matrix.c,v 1.2 2006/10/05 21:26:46 ichiki Exp $
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
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */

#include "stokes.h"
#include "bench.h"
#include "f.h"
#include "ft.h"
#include "fts.h"

#include "matrix.h"
#include "lub-matrix.h"


/* condition for lubrication
 * INPUT
 *  x1 [3], x2 [3] : position
 * OUTPUT (return value)
 *  0 : r != 0 and r < 3.0
 *  1 : otherwise
 */
static int
cond_lub (const double * x1, const double * x2)
{
  double x, y, z;
  double r2;


  x = x1 [0] - x2 [0];
  y = x1 [1] - x2 [1];
  z = x1 [2] - x2 [2];

  r2 = x * x
    + y * y
    + z * z;

  if (r2 != 0.0
      && r2 < 9.0) // r = 3.0 is the critical separation for lubrication now.
    {
      return 0;
    }
  else
    {
      return 1;
    }
}

/* make lubrication matrix for F version for all particles
 * under the periodic boundary condition
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 * OUTPUT
 *  mat [np * 3 * np * 3] :
 */
void
make_matrix_lub_ewald_3f (struct stokes * sys,
			  double * mat)
{
  int np;
  int i, j, k;
  int i3;
  int j3;
  int n;

  double tmp_pos [3];


  np = sys->np;
  n = np * 3;

  /* clear result */
  for (i = 0; i < n * n; ++i)
    {
      mat [i] = 0.0;
    }

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  /* all image cells */
	  for (k = 0; k < 27; ++k)
	    {
	      tmp_pos [0] = sys->pos [j3 + 0] + sys->llx [k];
	      tmp_pos [1] = sys->pos [j3 + 1] + sys->lly [k];
	      tmp_pos [2] = sys->pos [j3 + 2] + sys->llz [k];
	      if (cond_lub (sys->pos + i3, tmp_pos) == 0)
		{
		  matrix_lub_f_2b (sys,
				   i, j,
				   sys->pos + i3, tmp_pos,
				   n, mat);
		}
	    }
	}
    }

  free (tmp_pos);
}

/* make lubrication matrix for FT version for all particles
 * under the periodic boundary condition
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 * OUTPUT
 *  mat [np * 6 * np * 6] :
 */
void
make_matrix_lub_ewald_3ft (struct stokes * sys,
			   double * mat)
{
  int np;
  int i, j, k;
  int i3, j3;
  int n;

  double tmp_pos [3];


  np = sys->np;
  n = np * 6;

  /* clear result */
  for (i = 0; i < n * n; ++i)
    {
      mat [i] = 0.0;
    }

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  /* all image cells */
	  for (k = 0; k < 27; ++k)
	    {
	      tmp_pos [0] = sys->pos [j3 + 0] + sys->llx [k];
	      tmp_pos [1] = sys->pos [j3 + 1] + sys->lly [k];
	      tmp_pos [2] = sys->pos [j3 + 2] + sys->llz [k];
	      if (cond_lub (sys->pos + i3, tmp_pos) == 0)
		{
		  matrix_lub_ft_2b (sys,
				    i, j,
				    sys->pos + i3, tmp_pos,
				    n, mat);
		}
	    }
	}
    }

  free (tmp_pos);
}

/* make lubrication matrix for FTS version for all particles
 * under the periodic boundary condition
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 * OUTPUT
 *  mat [np * 11 * np * 11] :
 */
void
make_matrix_lub_ewald_3fts (struct stokes * sys,
			    double * mat)
{
  int np;
  int i, j, k;
  int i3, j3;
  int n;

  double * tmp_pos;


  np = sys->np;

  tmp_pos = malloc (sizeof (double) * 3);
  if (tmp_pos == NULL)
    {
      fprintf (stderr, "allocation error in calc_lub_ewald_3fts().\n");
      exit (1);
    }

  n = np * 11;

  /* clear result */
  for (i = 0; i < n * n; ++i)
    {
      mat [i] = 0.0;
    }

  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;
	  /* all image cells */
	  for (k = 0; k < 27; ++k)
	    {
	      tmp_pos [0] = sys->pos [j3 + 0] + sys->llx [k];
	      tmp_pos [1] = sys->pos [j3 + 1] + sys->lly [k];
	      tmp_pos [2] = sys->pos [j3 + 2] + sys->llz [k];
	      if (cond_lub (sys->pos + i3, tmp_pos) == 0)
		{
		  matrix_lub_fts_2b (sys,
				     i, j,
				     sys->pos + i3, tmp_pos,
				     n, mat);
		}
	    }
	}
    }

  free (tmp_pos);
}

