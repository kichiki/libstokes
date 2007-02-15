/* utility for Ewald summation calculation
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald.c,v 1.4 2007/02/15 03:27:07 kichiki Exp $
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
#include "f.h"      // matrix_f_ij() matrix_f_atimes()
#include "ft.h"     // matrix_ft_ij() matrix_ft_atimes()
#include "fts.h"    // matrix_fts_ij() matrix_fts_atimes()
#include "matrix.h" // dot_prod_matrix() multiply_extmat_with_extvec_3fts()
#include "non-ewald.h" // atimes_3all() make_matrix_mob_3all()

#include "ewald.h"


/*
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
scalars_ewald_real (int version,
		    double xi, double r,
		    double *xa, double *ya,
		    double *yb,
		    double *xc, double *yc,
		    double *xg, double *yg,
		    double *yh,
		    double *xm, double *ym, double *zm)
{
  double a2;
  double c2;

  double xir, xir2;
  double s, s2;
  double erfcxir;
  double expxir2;


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

  xir = xi * r;
  xir2 = xir * xir;
  s  = r;
  s2 = r * r;

  erfcxir = erfc (xir);
  expxir2 = xi / sqrt (M_PI) * exp (- xir2);


  *ya = (0.75 + 0.5 / s2) / s * erfcxir
    + ((1.0 + xir2 *
	(14.0 + 4.0 * xir2 *
	 (- 5.0 + xir2))) / s2
       - 4.5 + 3.0 * xir2)
    * expxir2;
  a2 = (0.75 - 1.5 / s2) / s * erfcxir
    + ((- 3.0 + xir2 *
	(- 2.0 + 4.0 * xir2 *
	 (4.0 - xir2))) / s2
       + 1.5 - 3.0 * xir2)
    * expxir2;
  *xa = a2 + *ya;

  // check point for F version
  if (version == 0) return;


  *yb = - 0.75 / s2 * erfcxir
    - 1.5 * (+ 1.0 + xir2 *
	     (- 6.0 + xir2 *
	      (+ 2.0)))
    / s * expxir2;

  *yc = - 3.0 / 8.0 / s2 / s * erfcxir
    - 0.75 * (+ 1.0 + xir2 *
	      (+ 14.0 + xir2 *
	       (-20.0 + xir2 *
		( + 4.0))))
    / s2 * expxir2;
  c2 = 9.0 / 8.0 / s2 / s * erfcxir
    - 0.75 * (- 3.0 + xir2 *
	      (- 2.0 + xir2 *
	       (+ 16.0 + xir2 *
		(- 4.0))))
    / s2 * expxir2;
  *xc = c2 + *yc;

  // check point for FT version
  if (version == 1) return;

  *xg = (2.25 - 3.6 / s2) / s2 * erfcxir
    + (- 1.5 * (- 3.0 + xir2 *
		(+ 6.0))
       - 0.8 * (+ 9.0 + xir2 *
		(+ 6.0 + xir2 *
		 (- 48.0 + xir2 *
		  (+ 12.0)))) / s2)
    / s * expxir2;
  *yg = 1.2 / s2 / s2 * erfcxir
    + (- 3.0 * ( xir2 *
		 (2.0 + xir2 *
		  (- 1.0)))
       - 0.8 * (- 3.0 + xir2 *
		(- 2.0 + xir2 *
		 (- 26.0 + xir2 *
		  (+ 26.0 + xir2 *
		   (- 4.0))))) / s2)
    / s * expxir2;

  *yh = - 9.0 / 8.0 / s2 / s * erfcxir
    + 1.5 * (- 1.5 + xir2 *
	     (- 1.0 + xir2 *
	      (+ 8.0 + xir2 *
	       (- 2.0))))
    / s2 * expxir2;

  *xm = (- 4.5 + 10.8 / s2) / s / s2 * erfcxir
    + (+ 1.5 * (- 6.0 +  xir2 *
		(- 12.0 + xir2 *
		 (+ 12.0)))
       + 1.2 * (+ 18.0 + xir2 *
		(+ 12.0 + xir2 *
		 (+ 30.0 + xir2 *
		  (- 66.0 + xir2 *
		   (+ 12.0))))) / s2)
    / s2 * expxir2;
  *ym = (+ 2.25 - 7.2 / s2) / s / s2 * erfcxir
    + (- 1.5 * (- 3.0 +  xir2 *
		(+ 6.0 + xir2 *
		 (- 12.0 + xir2 *
		  (+ 4.0))))
       - 1.2 * (+ 12.0 + xir2 *
		(+ 8.0 + xir2 *
		 (- 22.0 + xir2 *
		  (+ 58.0 + xir2 *
		   (- 34.0 + xir2 *
		    (+ 4.0)))))) / s2)
    / s2 * expxir2;
  *zm = + 1.8 / s2 / s / s2 * erfcxir
    + (- 1.5 * (+ 0.0 +  xir2 *
		(+ 8.0 + xir2 *
		 (- 4.0)))
       - 1.2 * (- 3.0 + xir2 *
		(- 2.0 + xir2 *
		 (- 26.0 + xir2 *
		  (+ 26.0 + xir2 *
		   (- 4.0))))) / s2)
    / s2 * expxir2;
}


/** table version **/

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
 * this routine also can handle non-periodic case seamlessly
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all (int n, const double *x, double *y, void * user_data)
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
  double r;

  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int m;

  double k1, k2, k3;
  double cf, sf;


  sys = (struct stokes *) user_data;

  if (sys->periodic == 0)
    {
      // non-periodic
      atimes_3all (n, x, y, user_data);
      return;
    }


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

  /* for xi code to measure CPU times */
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
	      rr = xx * xx + yy * yy + zz * zz;

	      if (rr == 0.0) continue; // to exclude the self part
	      if (rr > sys->rmax2) continue;

	      r = sqrt (rr);
	      ex = xx / r;
	      ey = yy / r;
	      ez = zz / r;

	      scalars_ewald_real (sys->version,
				  sys->xi, r,
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

  /* for xi code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu2 = cpu - cpu0;
  cpu0 = cpu;

  /* Second Ewald part ( reciprocal space ) */
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

  /* for xi code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu3 = cpu - cpu0;
  sys->cpu1 = sys->cpu2 + sys->cpu3;
}

/* make ewald-summed mobility matrix for F/FT/FTS versions
 * with the ewald table
 * this routine also can handle non-periodic case seamlessly
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_ewald_3all (struct stokes * sys, double * mat)
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
  double r;

  int n;
  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int m;

  double k1, k2, k3;
  double cf, sf;


  if (sys->periodic == 0)
    {
      // non-periodic
      make_matrix_mob_3all (sys, mat);
      return;
    }


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

  /* for xi code to measure CPU times */
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
	      rr = xx * xx + yy * yy + zz * zz;

	      if (rr == 0.0) continue; // to exclude the self part
	      if (rr > sys->rmax2) continue;

	      r = sqrt (rr);
	      ex = xx / r;
	      ey = yy / r;
	      ez = zz / r;

	      scalars_ewald_real (sys->version,
				  sys->xi, r,
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

  /* for xi code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu2 = cpu - cpu0;
  cpu0 = cpu;

  /* Second Ewald part ( reciprocal space ) */
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

  /* for xi code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu3 = cpu - cpu0;
  sys->cpu1 = sys->cpu2 + sys->cpu3;
}


/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * through matrix with the ewald table
 * this routine also can handle non-periodic case seamlessly
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_matrix (int n, const double *x,
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

  make_matrix_mob_ewald_3all (sys, mat);
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


/** non-table version **/

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_notbl (int n, const double *x,
			 double *y, void * user_data)
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
  double r;
  double rlx, rly, rlz;

  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double k;
  double cf, sf;
  double kexp;


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

  /* for xi code to measure CPU times */
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

	  for (m1 = - sys->rmaxx; m1 <= sys->rmaxx; m1++)
	    {
	      rlx = sys->lx * (double) m1;
	      for (m2 = - sys->rmaxy; m2 <= sys->rmaxy; m2++)
		{
		  rly = sys->ly * (double) m2;
		  for (m3 = - sys->rmaxz; m3 <= sys->rmaxz; m3++)
		    {
		      rlz = sys->lz * (double) m3;
  
		      xx = sys->pos [jx] - sys->pos [ix] + rlx;
		      yy = sys->pos [jy] - sys->pos [iy] + rly;
		      zz = sys->pos [jz] - sys->pos [iz] + rlz;
		      rr = xx * xx + yy * yy + zz * zz;

		      if (rr == 0.0) continue; // to exclude the self part
		      if (rr > sys->rmax2) continue;

		      r = sqrt (rr);
		      ex = xx / r;
		      ey = yy / r;
		      ez = zz / r;

		      scalars_ewald_real (sys->version,
					  sys->xi, r,
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
	}
    }

  /* for xi code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu2 = cpu - cpu0;
  cpu0 = cpu;

  /* Second Ewald part ( reciprocal space ) */
  for (m1 = - sys->kmaxx; m1 <= sys->kmaxx; m1++)
    {
      k1 = 2.0 * M_PI * (double) m1 / sys->lx;
      for (m2 = - sys->kmaxy; m2 <= sys->kmaxy; m2++)
	{
	  k2 = 2.0 * M_PI * (double) m2 / sys->ly;
	  for (m3 = - sys->kmaxz; m3 <= sys->kmaxz; m3++)
	    {
	      k3 = 2.0 * M_PI * (double) m3 / sys->lz;
	      if (m1 != 0 || m2 != 0 || m3 != 0)
		{
		  kk = k1 * k1 + k2 * k2 + k3 * k3;
		  k = sqrt (kk);

		  if (sys->kmax != 0.0 && k > sys->kmax) continue;

		  k4z = kk / 4.0 / sys->xi2;
		  kexp = sys->pivol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  ex = k1 / k;
		  ey = k2 / k;
		  ez = k3 / k;

		  ya = 6.0 * (1.0 - kk / 3.0) * kexp;
		  yb = 3.0 * k * kexp;
		  yc = 3.0 / 2.0 * kk * kexp;
		  yg = 3.0 * (1.0 - 4.0 / 15.0 * kk) * k * kexp;
		  yh = 3.0 / 2.0 * kk * kexp;
		  ym = 3.0 * (1.0 - kk / 5.0) * kk * kexp;
      
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
	    }
	}
    }

  /* for xi code to measure CPU times */
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
make_matrix_mob_ewald_3all_notbl (struct stokes * sys, double * mat)
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
  double r;
  double rlx, rly, rlz;

  int n;
  int i, j;
  int ix, iy, iz;
  int jx, jy, jz;
  int m1, m2, m3;

  double k1, k2, k3, kk, k4z;
  double k;
  double cf, sf;
  double kexp;


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

  /* for xi code to measure CPU times */
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

	  for (m1 = - sys->rmaxx; m1 <= sys->rmaxx; m1++)
	    {
	      rlx = sys->lx * (double) m1;
	      for (m2 = - sys->rmaxy; m2 <= sys->rmaxy; m2++)
		{
		  rly = sys->ly * (double) m2;
		  for (m3 = - sys->rmaxz; m3 <= sys->rmaxz; m3++)
		    {
		      rlz = sys->lz * (double) m3;
  
		      xx = sys->pos [jx] - sys->pos [ix] + rlx;
		      yy = sys->pos [jy] - sys->pos [iy] + rly;
		      zz = sys->pos [jz] - sys->pos [iz] + rlz;
		      rr = xx * xx + yy * yy + zz * zz;

		      if (rr == 0.0) continue; // to exclude the self part
		      if (rr > sys->rmax2) continue;

		      r = sqrt (rr);
		      ex = xx / r;
		      ey = yy / r;
		      ez = zz / r;

		      scalars_ewald_real (sys->version,
					  sys->xi, r,
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
	}
    }

  /* for xi code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu2 = cpu - cpu0;
  cpu0 = cpu;

  /* Second Ewald part ( reciprocal space ) */
  for (m1 = - sys->kmaxx; m1 <= sys->kmaxx; m1++)
    {
      k1 = 2.0 * M_PI * (double) m1 / sys->lx;
      for (m2 = - sys->kmaxy; m2 <= sys->kmaxy; m2++)
	{
	  k2 = 2.0 * M_PI * (double) m2 / sys->ly;
	  for (m3 = - sys->kmaxz; m3 <= sys->kmaxz; m3++)
	    {
	      k3 = 2.0 * M_PI * (double) m3 / sys->lz;
	      if (m1 != 0 || m2 != 0 || m3 != 0)
		{
		  kk = k1 * k1 + k2 * k2 + k3 * k3;
		  k = sqrt (kk);

		  if (sys->kmax != 0.0 && k > sys->kmax) continue;

		  k4z = kk / 4.0 / sys->xi2;
		  kexp = sys->pivol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  ex = k1 / k;
		  ey = k2 / k;
		  ez = k3 / k;

		  ya = 6.0 * (1.0 - kk / 3.0) * kexp;
		  yb = 3.0 * k * kexp;
		  yc = 3.0 / 2.0 * kk * kexp;
		  yg = 3.0 * (1.0 - 4.0 / 15.0 * kk) * k * kexp;
		  yh = 3.0 / 2.0 * kk * kexp;
		  ym = 3.0 * (1.0 - kk / 5.0) * kk * kexp;
      
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
	    }
	}
    }

  /* for xi code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu3 = cpu - cpu0;
  sys->cpu1 = sys->cpu2 + sys->cpu3;
}
/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * through matrix
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_matrix_notbl (int n, const double *x,
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

  // non-table version!
  make_matrix_mob_ewald_3all_notbl (sys, mat);
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
