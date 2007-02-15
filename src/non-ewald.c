/* utility for non-Ewald routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: non-ewald.c,v 1.1 2007/02/15 03:25:55 kichiki Exp $
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
#include "f.h"      // matrix_f_ij() matrix_f_atimes()
#include "ft.h"     // matrix_ft_ij() matrix_ft_atimes()
#include "fts.h"    // matrix_fts_ij() matrix_fts_atimes()
#include "matrix.h" // dot_prod_matrix() multiply_extmat_with_extvec_3fts()

#include "non-ewald.h"


/* calculate scalar functions under no periodic boundary condition
 * all functions are explicitly shown in Durlofsky-Brady (1987).
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 * OUTPUT
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 */
void
scalars_nonewald (int version,
		  double r,
		  double *xa, double *ya,
		  double *yb,
		  double *xc, double *yc,
		  double *xg, double *yg,
		  double *yh,
		  double *xm, double *ym, double *zm)
{
  double a2;
  double c2;

  double s, s2;

  // zero clear for FT and higher
  *yb = 0.0;
  *xc = 0.0;
  *yc = 0.0;
  *xg = 0.0;
  *yg = 0.0;
  *yh = 0.0;
  *xm = 0.0;
  *ym = 0.0;
  *zm = 0.0;

  s  = r;
  s2 = r * r;

  *ya = (0.75 + 0.5 / s2) / s;
  a2 = (0.75 - 1.5 / s2) / s;
  *xa = a2 + *ya;

  // check point for F version
  if (version == 0) return;


  *yb = - 0.75 / s2;

  *yc = - 3.0 / 8.0 / s2 / s;
  c2 = 9.0 / 8.0 / s2 / s;
  *xc = c2 + *yc;

  // check point for FT version
  if (version == 1) return;

  *xg = (2.25 - 3.6 / s2) / s2;
  *yg = 1.2 / s2 / s2;

  *yh = - 9.0 / 8.0 / s2 / s;

  *xm = (- 4.5 + 10.8 / s2) / s / s2;
  *ym = (+ 2.25 - 7.2 / s2) / s / s2;
  *zm = + 1.8 / s2 / s / s2;
}


/* ATIMES of calc plain mobility for F/FT/FTS versions
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_3all (int n, const double *x, double *y, void * user_data)
{
  struct stokes * sys;

  double xa, ya; 
  double yb;
  double xc, yc;
  double xg, yg;
  double yh;
  double xm, ym, zm;

  double ex, ey, ez;

  double xx, yy, zz, rr;
  double r;

  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;


  sys = (struct stokes *) user_data;

  /* clear result */
  for (i = 0; i < n; i ++)
    {
      y [i] = 0.0;
    }

  /* diagonal part ( self part ) */
  for (i = 0; i < sys->np; i++)
    {
      if (sys->version == 0) // F version
	{
	  matrix_f_atimes (x + i*3, y + i*3,
			   0.0, 0.0, 0.0,
			   1.0, 1.0); // a part
	}
      else if (sys->version == 1) // FT version
	{
	  matrix_ft_atimes (x + i*6, y + i*6,
			    0.0, 0.0, 0.0,
			    1.0, 1.0, // a part
			    0.0,
			    0.75, 0.75); // c part
	}
      else // FTS version
	{
	  matrix_fts_atimes (x + i*11, y + i*11,
			     0.0, 0.0, 0.0,
			     1.0, 1.0, // a part
			     0.0,
			     0.75, 0.75, // c part
			     0.0, 0.0,
			     0.0,
			     0.9, 0.9, 0.9); // m part
	}
    }

  /* loop for all particle-particle pairs without the self */
  for (i = 0; i < sys->np; i++)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = 0; j < sys->np; j++)
	{
	  // exclude the self part
	  if (i == j) continue;

	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;

	  xx = sys->pos [jx] - sys->pos [ix];
	  yy = sys->pos [jy] - sys->pos [iy];
	  zz = sys->pos [jz] - sys->pos [iz];
	  rr = xx * xx + yy * yy + zz * zz;

	  r = sqrt (rr);
	  ex = xx / r;
	  ey = yy / r;
	  ez = zz / r;

	  scalars_nonewald (sys->version,
			    r,
			    &xa, &ya,
			    &yb,
			    &xc, &yc,
			    &xg, &yg,
			    &yh,
			    &xm, &ym, &zm);

	  if (sys->version == 0) // F version
	    {
	      matrix_f_atimes (x + i*3, y + j*3,
			       ex, ey, ez,
			       xa, ya);
	    }
	  else if (sys->version == 1) // FT version
	    {
	      matrix_ft_atimes (x + i*6, y + j*6,
				ex, ey, ez,
				xa, ya,
				yb,
				xc, yc);
	    }
	  else // FTS version
	    {
	      matrix_fts_atimes (x + i*11, y + j*11,
				 ex, ey, ez,
				 xa, ya,
				 yb,
				 xc, yc,
				 xg, yg,
				 yh,
				 xm, ym, zm);
	    }
	}
    }
}

/* make plain mobility matrix for F/FT/FTS versions
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_3all (struct stokes * sys, double * mat)
{
  double xa, ya; 
  double yb;
  double xc, yc;
  double xg, yg;
  double yh;
  double xm, ym, zm;

  double ex, ey, ez;

  double xx, yy, zz, rr;
  double r;

  int n;
  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;


  if (sys->version == 0) // F version
    {
      n = sys->np * 3;
    }
  else if (sys->version == 1) // FT version
    {
      n = sys->np * 6;
    }
  else // FTS version
    {
      n = sys->np * 11;
    }

  /* clear result */
  for (i = 0; i < n * n; ++i)
    {
      mat [i] = 0.0;
    }

  /* diagonal part ( self part ) */
  for (i = 0; i < sys->np; i++)
    {
      if (sys->version == 0) // F version
	{
	  matrix_f_ij (i, i,
		       0.0, 0.0, 0.0,
		       1.0, 1.0, // a part
		       n, mat);
	}
      else if (sys->version == 1) // FT version
	{
	  matrix_ft_ij (i, i,
			0.0, 0.0, 0.0,
			1.0, 1.0, // a part
			0.0,
			0.75, 0.75, // c part
			n, mat);
	}
      else // FTS version
	{
	  matrix_fts_ij (i, i,
			 0.0, 0.0, 0.0,
			 1.0, 1.0, // a part
			 0.0,
			 0.75, 0.75, // c part
			 0.0, 0.0,
			 0.0,
			 0.9, 0.9, 0.9, // m part
			 n, mat);
	}
    }

  /* loop for all particle-particle pairs without the self */
  for (i = 0; i < sys->np; i++)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = 0; j < sys->np; j++)
	{
	  // exclude the self part
	  if (i == j) continue;

	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;

	  xx = sys->pos [jx] - sys->pos [ix];
	  yy = sys->pos [jy] - sys->pos [iy];
	  zz = sys->pos [jz] - sys->pos [iz];
	  rr = xx * xx + yy * yy + zz * zz;

	  r = sqrt (rr);
	  ex = xx / r;
	  ey = yy / r;
	  ez = zz / r;

	  scalars_nonewald (sys->version,
			    r,
			    &xa, &ya,
			    &yb,
			    &xc, &yc,
			    &xg, &yg,
			    &yh,
			    &xm, &ym, &zm);

	  if (sys->version == 0) // F version
	    {
	      matrix_f_ij (i, j,
			   ex, ey, ez,
			   xa, ya,
			   n, mat);
	    }
	  else if (sys->version == 1) // FT version
	    {
	      matrix_ft_ij (i, j,
			    ex, ey, ez,
			    xa, ya,
			    yb,
			    xc, yc,
			    n, mat);
	    }
	  else // FTS version
	    {
	      matrix_fts_ij (i, j,
			     ex, ey, ez,
			     xa, ya,
			     yb,
			     xc, yc,
			     xg, yg,
			     yh,
			     xm, ym, zm,
			     n, mat);
	    }
	}
    }
}

/* ATIMES of calc plain mobility for F/FT/FTS versions
 * through matrix with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_3all_matrix (int n, const double *x,
		    double *y, void * user_data)
{
  struct stokes * sys;
  int np;
  double * mat;


  sys = (struct stokes *) user_data;

  np = sys->np;
  mat = (double *) malloc (sizeof (double) * n * n);
  if (mat == NULL)
    {
      fprintf (stderr, "allocation error in atimes_ewald_3fts_matrix ().\n");
      exit (1);
    }

  make_matrix_mob_3all (sys, mat);
  if (sys->version == 0 // F version
      || sys->version == 1) // FT version
    {
      dot_prod_matrix (mat, n, n, x, y);
    }
  else // FTS version
    {
      multiply_extmat_with_extvec_3fts (sys->np, mat, x, y);
    }

  free (mat);
}
