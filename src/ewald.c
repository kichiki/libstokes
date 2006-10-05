/* utility for Ewald summation calculation
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald.c,v 1.1 2006/10/05 00:36:19 ichiki Exp $
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
#include <libiter.h> /* solve_iter() */

#include <libstokes.h> /* struct stokeks */
#include "/home/ichiki/WORK/SF/ryuon/libstokes/bench.h"
#include "/home/ichiki/WORK/SF/ryuon/libstokes/f.h"
#include "/home/ichiki/WORK/SF/ryuon/libstokes/ft.h"
#include "/home/ichiki/WORK/SF/ryuon/libstokes/fts.h"
#include "matrix.h"

#include "ewald.h"


/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_table (int n, const double *x, double *y, void * user_data)
{
  struct stokes * sys;

  double cpu0, cpu; /* for ptime_ms_d() */

  double xa, ya; 
  double yb;
  double xc, yc;
  double xg, yg;
  double yh;
  double xm, ym, zm;

  double ex, ey, ez;

  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;

  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int m;

  double k1, k2, k3;
  double cf, sf;

  double erfczr;
  double expzr2;

  double a2, c2;


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
			   sys->self_a, sys->self_a);
	}
      else if (sys->version == 1) // FT version
	{
	  matrix_ft_atimes (x + i*6, y + i*6,
			    0.0, 0.0, 0.0,
			    sys->self_a, sys->self_a,
			    0.0,
			    sys->self_c, sys->self_c);
	}
      else // FTS version
	{
	  matrix_fts_atimes (x + i*11, y + i*11,
			     0.0, 0.0, 0.0,
			     sys->self_a, sys->self_a,
			     0.0,
			     sys->self_c, sys->self_c,
			     0.0, 0.0,
			     0.0,
			     sys->self_m, sys->self_m, sys->self_m);
	}
    }

  /* for zeta code to measure CPU times */
  cpu0 = ptime_ms_d ();

  /* first Ewald part ( real space ) */
  for (i = 0; i < sys->np; i++)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = 0; j < sys->np; j++)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;

	  for (m = 0; m < sys->nr; m ++)
	    {
	      xx = sys->pos [jx] - sys->pos [ix] + sys->rlx [m];
	      yy = sys->pos [jy] - sys->pos [iy] + sys->rly [m];
	      zz = sys->pos [jz] - sys->pos [iz] + sys->rlz [m];
	      rr = sqrt (xx * xx + yy * yy + zz * zz);

	      if (rr == 0.0) continue; // to exclude the self part

	      if (sys->rmax != 0.0 && rr > sys->rmax) continue;

	      zr = sys->zeta * rr;
	      zr2 = zr * zr;
	      s  = rr;
	      s2 = s * s;

	      erfczr = erfc (zr);
	      expzr2 = sys->zaspi * exp (- zr2);

	      ex = xx / rr;
	      ey = yy / rr;
	      ez = zz / rr;

	      ya = (0.75 + 0.5 / s2) / s * erfczr
		+ ((1.0 + zr2 *
		    (14.0 + 4.0 * zr2 *
		     (- 5.0 + zr2))) / s2
		   - 4.5 + 3.0 * zr2)
		* expzr2;
	      a2 = (0.75 - 1.5 / s2) / s * erfczr
		+ ((- 3.0 + zr2 *
		    (- 2.0 + 4.0 * zr2 *
		     (4.0 - zr2))) / s2
		   + 1.5 - 3.0 * zr2)
		* expzr2;
	      xa = a2 + ya;

	      // check point for F version
	      if (sys->version == 0)
		{
		  matrix_f_atimes (x + i*3, y + j*3,
				   ex, ey, ez,
				   xa, ya);
		  continue;
		}

	      yb = - 0.75 / s2 * erfczr
		- 1.5 * (+ 1.0 + zr2 *
			 (- 6.0 + zr2 *
			  (+ 2.0)))
		/ s * expzr2;

	      yc = - 3.0 / 8.0 / s2 / s * erfczr
		- 0.75 * (+ 1.0 + zr2 *
			  (+ 14.0 + zr2 *
			   (-20.0 + zr2 *
			    ( + 4.0))))
		/ s2 * expzr2;
	      c2 = 9.0 / 8.0 / s2 / s * erfczr
		- 0.75 * (- 3.0 + zr2 *
			  (- 2.0 + zr2 *
			   (+ 16.0 + zr2 *
			    (- 4.0))))
		/ s2 * expzr2;
	      xc = c2 + yc;

	      // check point for FT version
	      if (sys->version == 1)
		{
		  matrix_ft_atimes (x + i*6, y + j*6,
				    ex, ey, ez,
				    xa, ya,
				    yb,
				    xc, yc);
		  continue;
		}

	      xg = (2.25 - 3.6 / s2) / s2 * erfczr
		+ (- 1.5 * (- 3.0 + zr2 *
			    (+ 6.0))
		   - 0.8 * (+ 9.0 + zr2 *
			    (+ 6.0 + zr2 *
			     (- 48.0 + zr2 *
			      (+ 12.0)))) / s2)
		/ s * expzr2;
	      yg = 1.2 / s2 / s2 * erfczr
		+ (- 3.0 * ( zr2 *
			     (2.0 + zr2 *
			      (- 1.0)))
		   - 0.8 * (- 3.0 + zr2 *
			    (- 2.0 + zr2 *
			     (- 26.0 + zr2 *
			      (+ 26.0 + zr2 *
			       (- 4.0))))) / s2)
		/ s * expzr2;

	      yh = - 9.0 / 8.0 / s2 / s * erfczr
		+ 1.5 * (- 1.5 + zr2 *
			 (- 1.0 + zr2 *
			  (+ 8.0 + zr2 *
			   (- 2.0))))
		/ s2 * expzr2;

	      xm = (- 4.5 + 10.8 / s2) / s / s2 * erfczr
		+ (+ 1.5 * (- 6.0 +  zr2 *
			    (- 12.0 + zr2 *
			     (+ 12.0)))
		   + 1.2 * (+ 18.0 + zr2 *
			    (+ 12.0 + zr2 *
			     (+ 30.0 + zr2 *
			      (- 66.0 + zr2 *
			       (+ 12.0))))) / s2)
		/ s2 * expzr2;
	      ym = (+ 2.25 - 7.2 / s2) / s / s2 * erfczr
		+ (- 1.5 * (- 3.0 +  zr2 *
			    (+ 6.0 + zr2 *
			     (- 12.0 + zr2 *
			      (+ 4.0))))
		   - 1.2 * (+ 12.0 + zr2 *
			    (+ 8.0 + zr2 *
			     (- 22.0 + zr2 *
			      (+ 58.0 + zr2 *
			       (- 34.0 + zr2 *
				(+ 4.0)))))) / s2)
		/ s2 * expzr2;
	      zm = + 1.8 / s2 / s / s2 * erfczr
		+ (- 1.5 * (+ 0.0 +  zr2 *
			    (+ 8.0 + zr2 *
			     (- 4.0)))
		   - 1.2 * (- 3.0 + zr2 *
			    (- 2.0 + zr2 *
			     (- 26.0 + zr2 *
			      (+ 26.0 + zr2 *
			       (- 4.0))))) / s2)
		/ s2 * expzr2;

	      // for FTS version
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

  /* for zeta code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu2 = cpu - cpu0;
  cpu0 = cpu;

  /* Second Ewald part ( reciprocal space ) */
  for (m = 0; m < sys->nk; m ++)
    {
      ex = sys->ex[m];
      ey = sys->ey[m];
      ez = sys->ez[m];

      k1 = sys->k1[m];
      k2 = sys->k2[m];
      k3 = sys->k3[m];

      ya = sys->ya[m];
      yb = sys->yb[m];
      yc = sys->yc[m];
      yg = sys->yg[m];
      yh = sys->yh[m];
      ym = sys->ym[m];

      //if (sys->kmax != 0.0 && sys->k[m] > sys->kmax) continue;
      // this is already satisfied.

      for (i = 0; i < sys->np; i++)
	{
	  ix = i * 3;
	  iy = ix + 1;
	  iz = ix + 2;
	  for (j = 0; j < sys->np; j++)
	    {
	      jx = j * 3;
	      jy = jx + 1;
	      jz = jx + 2;

	      xx = sys->pos [jx] - sys->pos [ix];
	      yy = sys->pos [jy] - sys->pos [iy];
	      zz = sys->pos [jz] - sys->pos [iz];

	      cf = cos (+ k1 * xx
			+ k2 * yy
			+ k3 * zz);

	      if (sys->version == 0) // F version
		{
		  matrix_f_atimes (x + i*3, y + j*3,
				   ex, ey, ez,
				   0.0, cf * ya);
		}
	      else if (sys->version == 1) // FT version
		{
		  sf = - sin (+ k1 * xx
			      + k2 * yy
			      + k3 * zz);

		  matrix_ft_atimes (x + i*6, y + j*6,
				    ex, ey, ez,
				    0.0, cf * ya,
				    sf * yb,
				    0.0, cf * yc);
		}
	      else // FTS version
		{
		  sf = - sin (+ k1 * xx
			      + k2 * yy
			      + k3 * zz);

		  matrix_fts_atimes (x + i*11, y + j*11,
				     ex, ey, ez,
				     0.0, cf * ya,
				     sf * yb,
				     0.0, cf * yc,
				     0.0, sf * yg,
				     cf * yh,
				     0.0, cf * ym, 0.0);
		}
	    }
	}
    }

  /* for zeta code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu3 = cpu - cpu0;
  sys->cpu1 = sys->cpu2 + sys->cpu3;
}

/* make ewald-summed mobility matrix for F/FT/FTS versions
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_ewald_3all_table (struct stokes * sys, double * mat)
{
  double cpu0, cpu; /* for ptime_ms_d() */

  double xa, ya; 
  double yb;
  double xc, yc;
  double xg, yg;
  double yh;
  double xm, ym, zm;

  double ex, ey, ez;

  double xx, yy, zz, rr;
  double zr, zr2;
  double s, s2;

  int n;
  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int m;

  double k1, k2, k3;
  double cf, sf;

  double erfczr;
  double expzr2;

  double a2, c2;


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
		       sys->self_a, sys->self_a,
		       n, mat);
	}
      else if (sys->version == 1) // FT version
	{
	  matrix_ft_ij (i, i,
			0.0, 0.0, 0.0,
			sys->self_a, sys->self_a,
			0.0,
			sys->self_c, sys->self_c,
			n, mat);
	}
      else // FTS version
	{
	  matrix_fts_ij (i, i,
			 0.0, 0.0, 0.0,
			 sys->self_a, sys->self_a,
			 0.0,
			 sys->self_c, sys->self_c,
			 0.0, 0.0,
			 0.0,
			 sys->self_m, sys->self_m, sys->self_m,
			 n, mat);
	}
    }

  /* for zeta code to measure CPU times */
  cpu0 = ptime_ms_d ();

  /* first Ewald part ( real space ) */
  for (i = 0; i < sys->np; i++)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = 0; j < sys->np; j++)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;

	  for (m = 0; m < sys->nr; m ++)
	    {
	      xx = sys->pos [jx] - sys->pos [ix] + sys->rlx [m];
	      yy = sys->pos [jy] - sys->pos [iy] + sys->rly [m];
	      zz = sys->pos [jz] - sys->pos [iz] + sys->rlz [m];
	      rr = sqrt (xx * xx + yy * yy + zz * zz);

	      if (rr == 0.0) continue; // to exclude the self part

	      if (sys->rmax != 0.0 && rr > sys->rmax) continue;

	      zr = sys->zeta * rr;
	      zr2 = zr * zr;
	      s  = rr;
	      s2 = s * s;

	      erfczr = erfc (zr);
	      expzr2 = sys->zaspi * exp (- zr2);

	      ex = xx / rr;
	      ey = yy / rr;
	      ez = zz / rr;

	      ya = (0.75 + 0.5 / s2) / s * erfczr
		+ ((1.0 + zr2 *
		    (14.0 + 4.0 * zr2 *
		     (- 5.0 + zr2))) / s2
		   - 4.5 + 3.0 * zr2)
		* expzr2;
	      a2 = (0.75 - 1.5 / s2) / s * erfczr
		+ ((- 3.0 + zr2 *
		    (- 2.0 + 4.0 * zr2 *
		     (4.0 - zr2))) / s2
		   + 1.5 - 3.0 * zr2)
		* expzr2;
	      xa = a2 + ya;

	      // check point for F version
	      if (sys->version == 0)
		{
		  matrix_f_ij (i, j,
			       ex, ey, ez,
			       xa, ya,
			       n, mat);
		  continue;
		}

	      yb = - 0.75 / s2 * erfczr
		- 1.5 * (+ 1.0 + zr2 *
			 (- 6.0 + zr2 *
			  (+ 2.0)))
		/ s * expzr2;

	      yc = - 3.0 / 8.0 / s2 / s * erfczr
		- 0.75 * (+ 1.0 + zr2 *
			  (+ 14.0 + zr2 *
			   (-20.0 + zr2 *
			    ( + 4.0))))
		/ s2 * expzr2;
	      c2 = 9.0 / 8.0 / s2 / s * erfczr
		- 0.75 * (- 3.0 + zr2 *
			  (- 2.0 + zr2 *
			   (+ 16.0 + zr2 *
			    (- 4.0))))
		/ s2 * expzr2;
	      xc = c2 + yc;

	      // check point for FT version
	      if (sys->version == 1)
		{
		  matrix_ft_ij (i, j,
				ex, ey, ez,
				xa, ya,
				yb,
				xc, yc,
				n, mat);
		  continue;
		}

	      xg = (2.25 - 3.6 / s2) / s2 * erfczr
		+ (- 1.5 * (- 3.0 + zr2 *
			    (+ 6.0))
		   - 0.8 * (+ 9.0 + zr2 *
			    (+ 6.0 + zr2 *
			     (- 48.0 + zr2 *
			      (+ 12.0)))) / s2)
		/ s * expzr2;
	      yg = 1.2 / s2 / s2 * erfczr
		+ (- 3.0 * ( zr2 *
			     (2.0 + zr2 *
			      (- 1.0)))
		   - 0.8 * (- 3.0 + zr2 *
			    (- 2.0 + zr2 *
			     (- 26.0 + zr2 *
			      (+ 26.0 + zr2 *
			       (- 4.0))))) / s2)
		/ s * expzr2;

	      yh = - 9.0 / 8.0 / s2 / s * erfczr
		+ 1.5 * (- 1.5 + zr2 *
			 (- 1.0 + zr2 *
			  (+ 8.0 + zr2 *
			   (- 2.0))))
		/ s2 * expzr2;

	      xm = (- 4.5 + 10.8 / s2) / s / s2 * erfczr
		+ (+ 1.5 * (- 6.0 +  zr2 *
			    (- 12.0 + zr2 *
			     (+ 12.0)))
		   + 1.2 * (+ 18.0 + zr2 *
			    (+ 12.0 + zr2 *
			     (+ 30.0 + zr2 *
			      (- 66.0 + zr2 *
			       (+ 12.0))))) / s2)
		/ s2 * expzr2;
	      ym = (+ 2.25 - 7.2 / s2) / s / s2 * erfczr
		+ (- 1.5 * (- 3.0 +  zr2 *
			    (+ 6.0 + zr2 *
			     (- 12.0 + zr2 *
			      (+ 4.0))))
		   - 1.2 * (+ 12.0 + zr2 *
			    (+ 8.0 + zr2 *
			     (- 22.0 + zr2 *
			      (+ 58.0 + zr2 *
			       (- 34.0 + zr2 *
				(+ 4.0)))))) / s2)
		/ s2 * expzr2;
	      zm = + 1.8 / s2 / s / s2 * erfczr
		+ (- 1.5 * (+ 0.0 +  zr2 *
			    (+ 8.0 + zr2 *
			     (- 4.0)))
		   - 1.2 * (- 3.0 + zr2 *
			    (- 2.0 + zr2 *
			     (- 26.0 + zr2 *
			      (+ 26.0 + zr2 *
			       (- 4.0))))) / s2)
		/ s2 * expzr2;

	      // for FTS version
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

  /* for zeta code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu2 = cpu - cpu0;
  cpu0 = cpu;

  /* Second Ewald part ( reciprocal space ) */
  for (m = 0; m < sys->nk; m ++)
    {
      ex = sys->ex[m];
      ey = sys->ey[m];
      ez = sys->ez[m];

      k1 = sys->k1[m];
      k2 = sys->k2[m];
      k3 = sys->k3[m];

      ya = sys->ya[m];
      yb = sys->yb[m];
      yc = sys->yc[m];
      yg = sys->yg[m];
      yh = sys->yh[m];
      ym = sys->ym[m];

      for (i = 0; i < sys->np; i++)
	{
	  ix = i * 3;
	  iy = ix + 1;
	  iz = ix + 2;
	  for (j = 0; j < sys->np; j++)
	    {
	      jx = j * 3;
	      jy = jx + 1;
	      jz = jx + 2;

	      xx = sys->pos [jx] - sys->pos [ix];
	      yy = sys->pos [jy] - sys->pos [iy];
	      zz = sys->pos [jz] - sys->pos [iz];

	      cf = cos (+ k1 * xx
			+ k2 * yy
			+ k3 * zz);

	      if (sys->version == 0) // F version
		{
		  matrix_f_ij (i, j,
			       ex, ey, ez,
			       0.0, cf * ya,
			       n, mat);
		}
	      else if (sys->version == 1) // FT version
		{
		  sf = - sin (+ k1 * xx
			      + k2 * yy
			      + k3 * zz);

		  matrix_ft_ij (i, j,
				ex, ey, ez,
				0.0, cf * ya,
				sf * yb,
				0.0, cf * yc,
				n, mat);
		}
	      else // FTS version
		{
		  sf = - sin (+ k1 * xx
			      + k2 * yy
			      + k3 * zz);

		  matrix_fts_ij (i, j,
				 ex, ey, ez,
				 0.0, cf * ya,
				 sf * yb,
				 0.0, cf * yc,
				 0.0, sf * yg,
				 cf * yh,
				 0.0, cf * ym, 0.0,
				 n, mat);
		}
	    }
	}
    }

  /* for zeta code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu3 = cpu - cpu0;
  sys->cpu1 = sys->cpu2 + sys->cpu3;
}


/** copy from test-fts-atimes.c Rev 1.5 **/
/* utility routine for matrix in the extracted form
 * INPUT
 *  np : # particles (not # elements!)
 *  m [np *11 * np *11] : matrix in the extracted form
 *  x [np *11] : vector in the extracted form
 * INPUT
 *  y [np *11] : output vector in the extracted form (:= m.x)
 */
static void
multiply_extmat_with_extvec_3fts (int np, const double * m, const double * x,
				  double * y)
{
  int n11;
  int i;
  int j;
  int j11;
  int jj;


  n11 = np * 11;

  for (i = 0; i < n11; ++i)
    {
      y [i] = 0.0;
      for (j = 0; j < np; ++j)
	{
	  j11 = j * 11;
	  for (jj = 0; jj < 6; ++jj)
	    {
	      y [i] += m [i * n11 + j11 + jj] * x [j11 + jj];
	    }
	  y [i] += m [i * n11 + j11 + 6]
	    * (2.0 * x [j11 + 6] + x [j11 + 10]);
	  y [i] += m [i * n11 + j11 + 7] * 2.0 * x [j11 + 7];
	  y [i] += m [i * n11 + j11 + 8] * 2.0 * x [j11 + 8];
	  y [i] += m [i * n11 + j11 + 9] * 2.0 * x [j11 + 9];
	  y [i] += m [i * n11 + j11 + 10]
	    * (2.0 * x [j11 + 10] + x [j11 + 6]);
	}
    }
}
/* ATIMES version (for O(N^2) scheme) of
 * calc ewald-summed mobility for F/FT/FTS versions
 * INPUT
 *  n := np * 11
 *  x [n * 11] : FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n * 11] : UOE
 */
void
atimes_ewald_3all_table_matrix (int n, const double *x,
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

  make_matrix_mob_ewald_3all_table (sys, mat);
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
