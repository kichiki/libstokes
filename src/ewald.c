/* utility for Ewald summation calculation
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ewald.c,v 1.8 2007/08/12 18:09:49 kichiki Exp $
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
#include "non-ewald.h" // atimes_nonewald_3all() make_matrix_mob_3all()

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

/* calculate scalar functions of (12)-interaction for unequal spheres
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  r      := x_2 - x_1
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 */
void
scalars_ewald_real_poly (int version,
			 double xi, double r,
			 double aa, double ab,
			 double *xa, double *ya,
			 double *yb,
			 double *xc, double *yc,
			 double *xg, double *yg,
			 double *yh,
			 double *xm, double *ym, double *zm)
{
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


  double xir = xi * r;
  double xir2 = xir * xir;
  double r2 = r * r;
  double r3 = r2 * r;
  double r4 = r2 * r2;
  double r5 = r3 * r2;

  double aa2 = aa * aa;
  double ab2 = ab * ab;

  double erfcxir = erfc (xir);
  double expxir2 = xi / sqrt (M_PI) * exp (- xir2);

  double a0 = 
    2.0 * (-3.0 +xir2* 2.0) * expxir2
    + erfcxir / r;
  double ad0 = 
    -2.0 / r * (1.0 +xir2* (-10.0 +xir2 * 4.0)) * expxir2
    - erfcxir / r2;
  double add0 = 
    4.0 / r2 * (1.0 +xir2* (6.0 +xir2* (-16.0 +xir2 * 4.0))) * expxir2
    + 2.0 / r3 * erfcxir;

  double b0 = 
    2.0 * (1.0 +xir2* (-2.0)) * expxir2
    + erfcxir / r;
  double bd0 = 
    -2.0 / r * (1.0 +xir2* (6.0 +xir2* (-4.0))) * expxir2
    - erfcxir / r2;
  /* not used
  double bdd0 = 
    4.0 / r2 * (1.0 +xir2* (-2.0 +xir2* (12.0 +xir2* (-4.0)))) * expxir2
    + 2.0 / r3 * erfcxir;
  */

  double a2 = 
    4.0 / r2 * (1.0 +xir2* (14.0 +xir2* (-20.0 +xir2 * 4.0))) * expxir2
    +2.0 / r3 * erfcxir;
  double ad2 = 
    -4.0 / r3 * (3.0 +xir2* (2.0 +xir2* (68.0 +xir2 * (-56.0 +xir2* 8.0)))) * expxir2
    -6.0 / r4 * erfcxir;
  double add2 = 
    16.0 / r4 * (3.0 +xir2* (2.0 +xir2* (-16.0 +xir2 * (76.0 +xir2* (-38.0 +xir2 * 4.0))))) * expxir2
    +24.0 / r5 * erfcxir;

  double b2 = 
    -4.0 / r2 * (3.0 +xir2* (2.0 +xir2* (-16.0 +xir2 * 4.0))) * expxir2
    -6.0 / r3 * erfcxir;
  double bd2 = 
    4.0 / r3 * (9.0 +xir2* (6.0 +xir2* (36.0 +xir2 * (-48.0 +xir2* 8.0)))) * expxir2
    +18.0 / r4 * erfcxir;
  /* not used
  double bdd2 = 
    -16.0 / r4 * (9.0 +xir2* (6.0 +xir2* (-6.0 +xir2 * (54.0 +xir2* (-34.0 +xir2 * 4.0))))) * expxir2
    -72.0 / r5 * erfcxir;
  */

  double a10 = a0;
  double a12 = a2 * (aa2+ab2)/6.0;
  double a20 = b0;
  double a22 = b2 * (aa2+ab2)/6.0;
  *ya = 0.75 * (a10 + a12);
  *xa = 0.75 * (a20 + a22) + (*ya);

  // check point for F version
  if (version == 0) return;

  double b10 = -0.375 * (-ad0 + b0/r);
  *yb = b10;

  double c1 = 0.1875 * (-add0 + (-ad0 + bd0)/r);
  double c2 = 0.1875 * ( add0 + (-ad0 - bd0 + 2.0*b0/r) /r);
  // note 3/16 = 0.1875
  *yc = c1;
  *xc = c2 + (*yc);

  // check point for FT version
  if (version == 1) return;

  double g10 = -0.75/r * b0;
  double g12 = -0.75/r * b2 * (aa2/10.0 + ab2/6.0);
  double g20 = -0.375 * (ad0 + b0/r);
  double g22 = -0.375 * (ad2 + b2/r) * (aa2/10.0 + ab2/6.0);
  *xg = -3.0 * (g10 + g12);
  *yg =         g20 + g22;

  /*
  // for consistency check
  double d;
  double g30 = -0.75 * (bd0 - 2.0*b0/r);
  double g32 = -0.75 * (bd2 - 2.0*b2/r)* (aa2/10.0 + ab2/6.0);
  // compare (g30 + g32) with (xg - 2 yg) !!
  d = fabs (g30 + g32 - ((*xg) - 2.0 * (*yg)));
  if (d > 1.0e-15)
    {
      fprintf (stderr, "inconsistency on g part: %e\n", d);
    }
  */

  double h10 = -0.1875 * (add0 + (-ad0 - bd0 + 2.0*b0/r) /r);
  // note 3/16 = 0.1875
  *yh = h10;

  double m50 = -0.375 * (ad0 + b0/r) /r;
  double m52 = -0.375 * (ad2 + b2/r) /r * (aa2+ab2)/10.0;
  *zm = 2.0 * (m50 + m52);

  double m30 = -0.75 * b0 / r2;
  double m32 = -0.75 * b2 / r2 * (aa2+ab2)/10.0;
  *xm = 6.0 * (m30 + m32) + 3.0 * (*zm);

  double m40 = -0.1875 * (add0 + (- ad0 + 3.0*bd0 - 6.0*b0/r) /r);
  double m42 = -0.1875 * (add2 + (- ad2 + 3.0*bd2 - 6.0*b2/r) /r) * (aa2 + ab2)/10.0;
  *ym = 2.0 * (m40 + m42) + *zm;

  /*
  // for consistency check
  double m10 = -0.75 * (bdd0 + (-5.0*bd0 + 8.0*b0/r) /r);
  double m12 = -0.75 * (bdd2 + (-5.0*bd2 + 8.0*b2/r) /r) * (aa2 + ab2)/10.0;
  d = fabs (m10 + m12 - (1.5*(*xm)-2.0*(*ym)+0.5*(*zm)));
  if (d > 1.0e-15)
    {
      fprintf (stderr, "inconsistency on m part (M1): %e\n", d);
    }

  double m20 = -0.75 * (bd0 - 2.0*b0/r) /r;
  double m22 = -0.75 * (bd2 - 2.0*b2/r) /r * (aa2+ab2)/10.0;
  d = fabs (-2.0 * (m20 + m22) + *zm - *xm);
  if (d > 1.0e-15)
    {
      fprintf (stderr, "inconsistency on m part(M2): %e\n", d);
    }
  */
}

/* convert scalar functions for mobility from dimensional to SD form
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  a1      : radius for the particle 1
 *            Note that the scalar functions are for (12)-interaction.
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 * OUTPUT
 *  xa, ya : for F version
 *  yb,
 *  xc, yc : for FT version
 *  xg, yg,
 *  yh,
 *  xm, ym, zm : for FTS version
 */
void
scalars_mob_poly_scale_SD_ (int version,
			   double a1,
			    double *xa, double *ya,
			    double *yb,
			    double *xc, double *yc,
			    double *xg, double *yg,
			    double *yh,
			    double *xm, double *ym, double *zm)
{
  // xa12, ya12
  (*xa) *= a1;
  (*ya) *= a1;
  // check point for F version
  if (version == 0) return;

  double a12 = a1*a1;
  double a13 = a12*a1;
  // yb12
  (*yb) *= a12;

  // xc12, yc12
  (*xc) *= a13;
  (*yc) *= a13;
  // check point for FT version
  if (version == 1) return;

  // xg12, yg12
  (*xg) *= a12;
  (*yg) *= a12;

  // yh12
  (*yh) *= a13;

  // xm12, ym12, zm12
  (*xm) *= a13;
  (*ym) *= a13;
  (*zm) *= a13;
}

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * this is a wrapper for non-periodic and periodic cases
 * also polydisperse systems for non-periodic
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

  sys = (struct stokes *) user_data;

  if (sys->periodic == 0)
    {
      // non-periodic
      atimes_nonewald_3all (n, x, y, user_data);
    }
  else
    {
      // periodic
      atimes_ewald_3all (n, x, y, user_data);
    }
}

/* make mobility matrix for F/FT/FTS versions
 * this is a wrapper for non-periodic and periodic cases
 * also polydisperse systems for non-periodic
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
  if (sys->periodic == 0)
    {
      // non-periodic
      make_matrix_mob_nonewald_3all (sys, mat);
    }
  else
    {
      // periodic
      make_matrix_mob_ewald_3all (sys, mat);
    }
}

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions through matrix
 * for both periodic and non-periodic boundary conditions
 * (now polydisperse can be handled for non-periodic case)
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
      fprintf (stderr, "allocation error in atimes_3all_matrix ().\n");
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


/** table version **/

/* ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
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

  /* clear result */
  for (i = 0; i < n; i ++)
    {
      y [i] = 0.0;
    }

  /* diagonal part ( self part ) */
  if (sys->a == NULL)
    {
      // monodisperse
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
    }
  else
    {
      // polydisperse
      double self_a;
      double self_c;
      double self_m;

      for (i = 0; i < sys->np; i++)
	{
	  double a = sys->a [i];
	  double a2 = a*a;
	  double a3 = a2*a;
	  self_a = 1.0
	    - a * sys->xiaspi * (6.0 - 40.0 / 3.0 * sys->xia2 * a2);
	  self_c = 0.75
	    - a3 * sys->xiaspi * sys->xia2 * 10.0;
	  self_m = 0.9
	    - a3 * sys->xiaspi * sys->xia2 * (12.0 - 30.24 * sys->xia2 * a2);
	  // in SD scaling

	  if (sys->version == 0) // F version
	    {
	      matrix_f_atimes (x + i*3, y + i*3,
			       0.0, 0.0, 0.0,
			       self_a, self_a);
	    }
	  else if (sys->version == 1) // FT version
	    {
	      matrix_ft_atimes (x + i*6, y + i*6,
				0.0, 0.0, 0.0,
				self_a, self_a,
				0.0,
				self_c, self_c);
	    }
	  else // FTS version
	    {
	      matrix_fts_atimes (x + i*11, y + i*11,
				 0.0, 0.0, 0.0,
				 self_a, self_a,
				 0.0,
				 self_c, self_c,
				 0.0, 0.0,
				 0.0,
				 self_m, self_m, self_m);
	    }
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

	      if (sys->a == NULL)
		{
		  // monodisperse
		  scalars_ewald_real (sys->version,
				      sys->xi, r,
				      &xa, &ya,
				      &yb,
				      &xc, &yc,
				      &xg, &yg,
				      &yh,
				      &xm, &ym, &zm);
		}
	      else
		{
		  // polydisperse
		  scalars_ewald_real_poly (sys->version,
					   sys->xi, r,
					   sys->a[i], sys->a[j],
					   &xa, &ya,
					   &yb,
					   &xc, &yc,
					   &xg, &yg,
					   &yh,
					   &xm, &ym, &zm);
		  scalars_mob_poly_scale_SD_ (sys->version,
					      sys->a[i],
					      &xa, &ya,
					      &yb,
					      &xc, &yc,
					      &xg, &yg,
					      &yh,
					      &xm, &ym, &zm);
		  // now scalars are in the SD form
		}

	      // note that interaction (i,j) should be for (U[i], F[j])
	      if (sys->version == 0) // F version
		{
		  matrix_f_atimes (x + j*3, y + i*3,
				   ex, ey, ez,
				   xa, ya);
		}
	      else if (sys->version == 1) // FT version
		{
		  matrix_ft_atimes (x + j*6, y + i*6,
				    ex, ey, ez,
				    xa, ya,
				    yb,
				    xc, yc);
		}
	      else // FTS version
		{
		  matrix_fts_atimes (x + j*11, y + i*11,
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

	      if (sys->a == NULL)
		{
		  // monodisperse
		  ya = sys->ya[m];
		  yb = sys->yb[m];
		  yc = sys->yc[m];
		  yg = sys->yg[m];
		  yh = sys->yh[m];
		  ym = sys->ym[m];
		}
	      else
		{
		  // polydisperse
		  double aa;
		  double aa3;
		  double aa2;
		  aa = sys->a [i];
		  aa2 = aa * aa;
		  aa3 = aa2 * aa;
		  double ab2;
		  ab2 = sys->a [j];
		  ab2 *= ab2;

		  double k = sys->k [m];
		  double kk = k * k;
		  double k4z = kk / 4.0 / sys->xi2;
		  double kexp = sys->pivol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  ya = aa * 6.0 * (1.0 - kk * (aa2+ab2) / 6.0) * kexp;
		  yb = aa2 * 3.0 * k * kexp;
		  yc = aa3 * 3.0 / 2.0 * kk * kexp;
		  yg = aa2 * 3.0 * (1.0 - kk * (aa2/10.0 + ab2/6.0)) * k * kexp;
		  yh = aa3 * 3.0 / 2.0 * kk * kexp;
		  ym = aa3 * 3.0 * (1.0 - kk * (aa2+ab2) / 10.0) * kk * kexp;
		  // in SD scaling
		}

	      cf = cos (+ k1 * xx
			+ k2 * yy
			+ k3 * zz);

	      // note that interaction (i,j) should be for (U[i], F[j])
	      if (sys->version == 0) // F version
		{
		  matrix_f_atimes (x + j*3, y + i*3,
				   ex, ey, ez,
				   0.0, cf * ya);
		}
	      else if (sys->version == 1) // FT version
		{
		  sf = - sin (+ k1 * xx
			      + k2 * yy
			      + k3 * zz);

		  matrix_ft_atimes (x + j*6, y + i*6,
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

		  matrix_fts_atimes (x + j*11, y + i*11,
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
  if (sys->a == NULL)
    {
      // monodisperse
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
    }
  else
    {
      // polydisperse
      double self_a;
      double self_c;
      double self_m;

      for (i = 0; i < sys->np; i++)
	{
	  double a = sys->a [i];
	  double a2 = a*a;
	  double a3 = a2*a;
	  self_a = 1.0
	    - a * sys->xiaspi * (6.0 - 40.0 / 3.0 * sys->xia2 * a2);
	  self_c = 0.75
	    - a3 * sys->xiaspi * sys->xia2 * 10.0;
	  self_m = 0.9
	    - a3 * sys->xiaspi * sys->xia2 * (12.0 - 30.24 * sys->xia2 * a2);
	  // in SD scaling

	  if (sys->version == 0) // F version
	    {
	      matrix_f_ij (i, i,
			   0.0, 0.0, 0.0,
			   self_a, self_a,
			   n, mat);
	    }
	  else if (sys->version == 1) // FT version
	    {
	      matrix_ft_ij (i, i,
			    0.0, 0.0, 0.0,
			    self_a, self_a,
			    0.0,
			    self_c, self_c,
			    n, mat);
	    }
	  else // FTS version
	    {
	      matrix_fts_ij (i, i,
			     0.0, 0.0, 0.0,
			     self_a, self_a,
			     0.0,
			     self_c, self_c,
			     0.0, 0.0,
			     0.0,
			     self_m, self_m, self_m,
			     n, mat);
	    }
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

	      if (sys->a == NULL)
		{
		  // monodisperse
		  scalars_ewald_real (sys->version,
				      sys->xi, r,
				      &xa, &ya,
				      &yb,
				      &xc, &yc,
				      &xg, &yg,
				      &yh,
				      &xm, &ym, &zm);
		}
	      else
		{
		  // polydisperse
		  scalars_ewald_real_poly (sys->version,
					   sys->xi, r,
					   sys->a[i], sys->a[j],
					   &xa, &ya,
					   &yb,
					   &xc, &yc,
					   &xg, &yg,
					   &yh,
					   &xm, &ym, &zm);
		  scalars_mob_poly_scale_SD_ (sys->version,
					      sys->a[i],
					      &xa, &ya,
					      &yb,
					      &xc, &yc,
					      &xg, &yg,
					      &yh,
					      &xm, &ym, &zm);
		  // now scalars are in the SD form
		}

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

	      if (sys->a == NULL)
		{
		  // monodisperse
		  ya = sys->ya[m];
		  yb = sys->yb[m];
		  yc = sys->yc[m];
		  yg = sys->yg[m];
		  yh = sys->yh[m];
		  ym = sys->ym[m];
		}
	      else
		{
		  // polydisperse
		  double aa;
		  double aa3;
		  double aa2;
		  aa = sys->a [i];
		  aa2 = aa * aa;
		  aa3 = aa2 * aa;
		  double ab2;
		  ab2 = sys->a [j];
		  ab2 *= ab2;

		  double k = sys->k [m];
		  double kk = k * k;
		  double k4z = kk / 4.0 / sys->xi2;
		  double kexp = sys->pivol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  ya = aa * 6.0 * (1.0 - kk * (aa2+ab2) / 6.0) * kexp;
		  yb = aa2 * 3.0 * k * kexp;
		  yc = aa3 * 3.0 / 2.0 * kk * kexp;
		  yg = aa2 * 3.0 * (1.0 - kk * (aa2/10.0 + ab2/6.0)) * k * kexp;
		  yh = aa3 * 3.0 / 2.0 * kk * kexp;
		  ym = aa3 * 3.0 * (1.0 - kk * (aa2+ab2) / 10.0) * kk * kexp;
		  // in SD scaling
		}

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
