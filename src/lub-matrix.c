/* lubrication routines -- MATRIX procedure
 * Copyright (C) 1993-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: lub-matrix.c,v 1.10 2007/04/20 01:55:24 kichiki Exp $
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
#include "f.h"
#include "ft.h"
#include "fts.h"

#include "matrix.h"
#include "lub-matrix.h"


/* condition for lubrication
 * INPUT
 *  x1 [3], x2 [3] : position
 *  lubmax2 : square of the max distance (0 means no limit)
 * OUTPUT (return value)
 *  0 : r != 0 and r < 3.0
 *  1 : otherwise
 */
static int
cond_lub (const double *x1, const double *x2, double lubmax2)
{
  double x, y, z;
  x = x1 [0] - x2 [0];
  y = x1 [1] - x2 [1];
  z = x1 [2] - x2 [2];

  double r2 = x*x + y*y + z*z;

  if (r2 == 0.0)
    {
      // self part
      return 1; // false
    }
  else if (lubmax2 <= 0.0)
    {
      // no limit
      return 0; // true
    }
  else if (r2 < lubmax2)
    {
      // within the limit
      return 0; // true
    }
  else
    {
      // out of the limit
      return 1; // false
    }
}

/* condition for lubrication for polydisperse system
 * INPUT
 *  x1 [3], x2 [3] : position
 *  a1, a2         : radii for particles 1 and 2
 *  lubmax2        : square of the max distance (0 means no limit)
 * OUTPUT (return value)
 *  0 : r != 0 and r < 3.0
 *  1 : otherwise
 */
static int
cond_lub_poly (const double *x1, const double *x2,
	       double a1, double a2,
	       double lubmax2)
{
  double x, y, z;
  x = x1 [0] - x2 [0];
  y = x1 [1] - x2 [1];
  z = x1 [2] - x2 [2];

  double r2 = x*x + y*y + z*z;

  if (r2 == 0.0)
    {
      // self part
      return 1; // false
    }
  else if (lubmax2 <= 0.0)
    {
      // no limit
      return 0; // true
    }
  else
    {
      // max2 should be compared with r2 = (2 r / (a1 + a2))^2
      double fac = 2.0 / (a1 + a2);
      double s2 = fac * fac * r2;
      if (s2 < lubmax2)
	{
	  // within the limit
	  return 0; // true
	}
      else
	{
	  // out of the limit
	  return 1; // false
	}
    }
}

/*
 * INPUT
 *  sys : struct stokes *sys.
 *        lx, ly, lz, and lubmax are used.
 * OUTPUT
 *  *imax[xyz] : range for the image-cell loop
 */
static void
set_imax_lub_periodic (struct stokes *sys,
		       int *imaxx, int *imaxy, int *imaxz)
{
  double lubmax;
  if (sys->lubmax <= 0.0)
    {
      // if lubmax2 is not define, at least, take symmetric in x, y, z.
      // set lubmax by max (lx, ly, lz)
      lubmax = sys->lx;
      if (sys->ly > lubmax) lubmax = sys->ly;
      if (sys->lz > lubmax) lubmax = sys->lz;
    }
  else
    {
      lubmax = sys->lubmax;
    }

  *imaxx = (int)(lubmax / sys->lx) + 1;
  *imaxy = (int)(lubmax / sys->ly) + 1;
  *imaxz = (int)(lubmax / sys->lz) + 1;
}


/* make lubrication matrix for F version for all particles
 * for both periodic and non-periodic boundary conditions
 * polydisperse system can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 * OUTPUT
 *  mat [np * 3 * np * 3] :
 */
void
make_matrix_lub_3f (struct stokes *sys,
		    double *mat)
{
  int i, j;
  int i3;
  int j3;
  int n;

  double tmp_pos [3];

  int np = sys->np;
  n = np * 3;
  /* clear result */
  for (i = 0; i < n * n; ++i)
    {
      mat [i] = 0.0;
    }

  // set lubmax2 for cond_lub()
  double lubmax2 = sys->lubmax * sys->lubmax;

  // set imax[xyz] for periodic systems
  // which covers enough images for sys->lubmax
  int imaxx = 0;
  int imaxy = 0;
  int imaxz = 0;
  if (sys->periodic != 0)
    {
      set_imax_lub_periodic (sys, &imaxx, &imaxy, &imaxz);
    }


  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;

	  if (sys->periodic == 0)
	    {
	      // non-periodic
	      if (sys->a == NULL)
		{
		  // monodisperse
		  if (cond_lub (sys->pos + i3, sys->pos + j3,
				lubmax2) == 0)
		    {
		      matrix_lub_f_2b (sys,
				       i, j,
				       sys->pos + i3, sys->pos + j3,
				       n, mat);
		    }
		}
	      else
		{
		  // polydisperse
		  if (cond_lub_poly (sys->pos + i3, sys->pos + j3,
				     sys->a[i], sys->a[j],
				     lubmax2) == 0)
		    {
		      matrix_lub_f_2b_poly (sys,
					    i, j,
					    sys->pos + i3, sys->pos + j3,
					    sys->a[i], sys->a[j],
					    n, mat);
		    }
		}
	    }
	  else
	    {
	      int ix, iy, iz;
	      for (ix = -imaxx; ix <= imaxx; ix++)
		{
		  tmp_pos[0] = sys->pos[j3 + 0] + sys->lx * (double)ix;
		  for (iy = -imaxy; iy <= imaxy; iy++)
		    {
		      tmp_pos[1] = sys->pos[j3 + 1] + sys->ly * (double)iy;
		      for (iz = -imaxz; iz <= imaxz; iz++)
			{
			  tmp_pos[2] = sys->pos[j3 + 2] + sys->lz * (double)iz;

			  if (sys->a == NULL)
			    {
			      // monodisperse
			      if (cond_lub (sys->pos + i3, tmp_pos,
					    lubmax2) == 0)
				{
				  matrix_lub_f_2b
				    (sys,
				     i, j,
				     sys->pos + i3, tmp_pos,
				     n, mat);
				}
			    }
			  else
			    {
			      // polydisperse
			      if (cond_lub_poly (sys->pos + i3, sys->pos + j3,
						 sys->a[i], sys->a[j],
						 lubmax2) == 0)
				{
				  matrix_lub_f_2b_poly
				    (sys,
				     i, j,
				     sys->pos + i3, tmp_pos,
				     sys->a[i], sys->a[j],
				     n, mat);
				}
			    }
			}
		    }
		}
	      // endif for periodic
	    }
	}
    }
}

/* make lubrication matrix for FT version for all particles
 * for both periodic and non-periodic boundary conditions
 * polydisperse system can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 * OUTPUT
 *  mat [np * 6 * np * 6] :
 */
void
make_matrix_lub_3ft (struct stokes *sys,
		     double *mat)
{
  int i, j;
  int i3, j3;
  int n;

  double tmp_pos [3];

  int np = sys->np;
  n = np * 6;
  /* clear result */
  for (i = 0; i < n * n; ++i)
    {
      mat [i] = 0.0;
    }

  // set lubmax2 for cond_lub()
  double lubmax2 = sys->lubmax * sys->lubmax;

  // set imax[xyz] for periodic systems
  // which covers enough images for sys->lubmax
  int imaxx = 0;
  int imaxy = 0;
  int imaxz = 0;
  if (sys->periodic != 0)
    {
      set_imax_lub_periodic (sys, &imaxx, &imaxy, &imaxz);
    }


  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;

	  if (sys->periodic == 0)
	    {
	      // non-periodic
	      if (sys->a == NULL)
		{
		  // monodisperse
		  if (cond_lub (sys->pos + i3, sys->pos + j3,
				lubmax2) == 0)
		    {
		      matrix_lub_ft_2b (sys,
					i, j,
					sys->pos + i3, sys->pos + j3,
					n, mat);
		    }
		}
	      else
		{
		  // polydisperse
		  if (cond_lub_poly (sys->pos + i3, sys->pos + j3,
				     sys->a[i], sys->a[j],
				     lubmax2) == 0)
		    {
		      matrix_lub_ft_2b_poly (sys,
					     i, j,
					     sys->pos + i3, sys->pos + j3,
					     sys->a[i], sys->a[j],
					     n, mat);
		    }
		}
	    }
	  else
	    {
	      int ix, iy, iz;
	      for (ix = -imaxx; ix <= imaxx; ix++)
		{
		  tmp_pos[0] = sys->pos[j3 + 0] + sys->lx * (double)ix;
		  for (iy = -imaxy; iy <= imaxy; iy++)
		    {
		      tmp_pos[1] = sys->pos[j3 + 1] + sys->ly * (double)iy;
		      for (iz = -imaxz; iz <= imaxz; iz++)
			{
			  tmp_pos[2] = sys->pos[j3 + 2] + sys->lz * (double)iz;

			  if (sys->a == NULL)
			    {
			      // monodisperse
			      if (cond_lub (sys->pos + i3, tmp_pos,
					    lubmax2) == 0)
				{
				  matrix_lub_ft_2b
				    (sys,
				     i, j,
				     sys->pos + i3, tmp_pos,
				     n, mat);
				}
			    }
			  else
			    {
			      // polydisperse
			      if (cond_lub_poly (sys->pos + i3, sys->pos + j3,
						 sys->a[i], sys->a[j],
						 lubmax2) == 0)
				{
				  matrix_lub_ft_2b_poly
				    (sys,
				     i, j,
				     sys->pos + i3, tmp_pos,
				     sys->a[i], sys->a[j],
				     n, mat);
				}
			    }
			}
		    }
		}
	      // endif for periodic
	    }
	}
    }
}

/* make lubrication matrix for FTS version for all particles
 * for both periodic and non-periodic boundary conditions
 * polydisperse system can be handled.
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 * OUTPUT
 *  mat [np * 11 * np * 11] :
 */
void
make_matrix_lub_3fts (struct stokes *sys,
		      double *mat)
{
  int i, j;
  int i3, j3;
  int n;

  double tmp_pos [3];

  int np = sys->np;
  n = np * 11;
  /* clear result */
  for (i = 0; i < n * n; ++i)
    {
      mat [i] = 0.0;
    }

  // set lubmax2 for cond_lub()
  double lubmax2 = sys->lubmax * sys->lubmax;

  // set imax[xyz] for periodic systems
  // which covers enough images for sys->lubmax
  int imaxx = 0;
  int imaxy = 0;
  int imaxz = 0;
  if (sys->periodic != 0)
    {
      set_imax_lub_periodic (sys, &imaxx, &imaxy, &imaxz);
    }


  for (i = 0; i < np; ++i)
    {
      i3 = i * 3;
      for (j = i; j < np; ++j)
	{
	  j3 = j * 3;

	  if (sys->periodic == 0)
	    {
	      // non-periodic
	      if (sys->a == NULL)
		{
		  // monodisperse
		  if (cond_lub (sys->pos + i3, sys->pos + j3,
				lubmax2) == 0)
		    {
		      matrix_lub_fts_2b (sys,
					 i, j,
					 sys->pos + i3, sys->pos + j3,
					 n, mat);
		    }
		}
	      else
		{
		  // polydisperse
		  if (cond_lub_poly (sys->pos + i3, sys->pos + j3,
				     sys->a[i], sys->a[j],
				     lubmax2) == 0)
		    {
		      matrix_lub_fts_2b_poly (sys,
					      i, j,
					      sys->pos + i3, sys->pos + j3,
					      sys->a[i], sys->a[j],
					      n, mat);
		    }
		}
	    }
	  else
	    {
	      int ix, iy, iz;
	      for (ix = -imaxx; ix <= imaxx; ix++)
		{
		  tmp_pos[0] = sys->pos[j3 + 0] + sys->lx * (double)ix;
		  for (iy = -imaxy; iy <= imaxy; iy++)
		    {
		      tmp_pos[1] = sys->pos[j3 + 1] + sys->ly * (double)iy;
		      for (iz = -imaxz; iz <= imaxz; iz++)
			{
			  tmp_pos[2] = sys->pos[j3 + 2] + sys->lz * (double)iz;

			  if (sys->a == NULL)
			    {
			      // monodisperse
			      if (cond_lub (sys->pos + i3, tmp_pos,
					    lubmax2) == 0)
				{
				  matrix_lub_fts_2b
				    (sys,
				     i, j,
				     sys->pos + i3, tmp_pos,
				     n, mat);
				}
			    }
			  else
			    {
			      // polydisperse
			      if (cond_lub_poly (sys->pos + i3, sys->pos + j3,
						 sys->a[i], sys->a[j],
						 lubmax2) == 0)
				{
				  matrix_lub_fts_2b_poly
				    (sys,
				     i, j,
				     sys->pos + i3, tmp_pos,
				     sys->a[i], sys->a[j],
				     n, mat);
				}
			    }
			}
		    }
		}
	    }
	}
    }
}
