/* utility for non-Ewald routines
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: non-ewald.c,v 1.7 2007/08/13 00:31:11 kichiki Exp $
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
 *  scalar [11]:
 *   0, 1,    : (xa12, ya12) for F version
 *   2,       : (yb12)
 *   3, 4,    : (xc12, yc12) for FT version
 *   5, 6,    : (xg12, yg12)
 *   7,       : (yh12)
 *   8, 9, 10 : (xm12, ym12, zm12) for FTS version
 */
void
scalars_nonewald (int version,
		  double r,
		  double *scalar)
{
  // zero clear for FT and higher
  int i;
  for (i = 0; i < 11; i ++) scalar [i] = 0.0;

  double r2;
  r2 = r * r;

  // xa12, ya12
  scalar [0] = (1.5  - 1.0 / r2) / r;
  scalar [1] = (0.75 + 0.5 / r2) / r;

  // check point for F version
  if (version == 0) return;


  // yb12
  scalar [2] = - 0.75 / r2;

  // xc, yc
  scalar [3] = + 0.75  / r2 / r;
  scalar [4] = - 0.375 / r2 / r;

  // check point for FT version
  if (version == 1) return;

  // xg, yg
  scalar [5] = (2.25 - 3.6 / r2) / r2;
  scalar [6] = 1.2 / r2 / r2;

  // yh
  scalar [7] = - 1.125 / r2 / r;

  // xm, ym, zm
  scalar [8] = (- 4.5 + 10.8 / r2) / r2 / r;
  scalar [9] = (+ 2.25 - 7.2 / r2) / r2 / r;
  scalar[10] = + 1.8 / r2 / r2 / r;
}

/* calculate scalar functions for unequal spheres
 * under no periodic boundary condition in dimensional form
 * to convert them in the SD form, use scalars_mob_poly_scale_SD ().
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar [11]:
 *   0, 1,    : (xa12, ya12) for F version
 *   2,       : (yb12)
 *   3, 4,    : (xc12, yc12) for FT version
 *   5, 6,    : (xg12, yg12)
 *   7,       : (yh12)
 *   8, 9, 10 : (xm12, ym12, zm12) for FTS version
 */
void
scalars_nonewald_poly (int version,
		       double r,
		       double aa, double ab,
		       double *scalar)
{
  // zero clear for FT and higher
  int i;
  for (i = 0; i < 11; i ++) scalar [i] = 0.0;

  double r2 = r * r;

  double aa2 = aa * aa;
  double ab2 = ab * ab;

  // xa12, ya12
  scalar [0] = (1.5 - (aa2 + ab2) * 0.5 / r2) / r;
  scalar [1] = (0.75 + (aa2 + ab2) * 0.25 / r2) / r;

  // check point for F version
  if (version == 0) return;

  // yb12
  scalar [2] = - 0.75 / r2;

  // xc12, yc12
  scalar [3] = 0.75 / r2 / r;
  scalar [4] = - 3.0 / 8.0 / r2 / r;

  // check point for FT version
  if (version == 1) return;

  // xg12, yg12
  scalar [5] = (2.25 - 13.5 * (aa2 / 10.0 + ab2 / 6.0) / r2) / r2;
  scalar [6] = 4.5 * (aa2 / 10.0 + ab2 / 6.0) / r2 / r2;

  // yh12
  scalar [7] = - 9.0 / 8.0 / r2 / r;

  // xm12, ym12, zm12
  scalar [8] = (- 4.5 + 5.4 * (aa2 + ab2) / r2) / r2 / r;
  scalar [9] = (2.25 - 3.6 * (aa2 + ab2) / r2) / r2 / r;
  scalar[10] = 0.9 * (aa2 + ab2) / r2 / r2 / r;
}

/* convert scalar functions for mobility from dimensional to SD form
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  a1      : radius for the particle 1
 *            Note that the scalar functions are for (12)-interaction.
 *  scalar [11]:
 *    0, 1,    : (xa12, ya12) for F version
 *    2,       : (yb12)
 *    3, 4,    : (xc12, yc12) for FT version
 *    5, 6,    : (xg12, yg12)
 *    7,       : (yh12)
 *    8, 9, 10 : (xm12, ym12, zm12) for FTS version
 * OUTPUT
 *  scalar [11]: scaled
 */
void
scalars_mob_poly_scale_SD (int version,
			   double a1,
			   double *scalar)
{
  // xa12, ya12
  scalar [0] *= a1;
  scalar [1] *= a1;
  // check point for F version
  if (version == 0) return;

  double a12 = a1*a1;
  double a13 = a12*a1;
  // yb12
  scalar [2] *= a12;

  // xc12, yc12
  scalar [3] *= a13;
  scalar [4] *= a13;
  // check point for FT version
  if (version == 1) return;

  // xg12, yg12
  scalar [5] *= a12;
  scalar [6] *= a12;

  // yh12
  scalar [7] *= a13;

  // xm12, ym12, zm12
  scalar [8] *= a13;
  scalar [9] *= a13;
  scalar[10] *= a13;
}

/* calculate scalar functions for unequal spheres
 * under no periodic boundary condition in dimensional form
 * to convert them in the SD form, use scalars_mob_poly_scale_SD ().
 * INPUT
 *  version : 0 = F version
 *            1 = FT version
 *            2 = FTS version
 *  r      := x_b - x_a
 *  aa, ab : radius of particle a and b
 * OUTPUT
 *  scalar [44] : scalar functions in dimensional form!
 *   0, 1, 2, 3 : (xa11, xa12, xa21, xa22)
 *   4, 5, 6, 7 : (ya11, ya12, ya21, ya22)
 *   8, 9,10,11 : (yb11, yb12, yb21, yb22)
 *  12,13,14,15 : (xc11, xc12, xc21, xc22)
 *  16,17,18,19 : (yc11, yc12, yc21, yc22)
 *  20,21,22,23 : (xg11, xg12, xg21, xg22)
 *  24,25,26,27 : (yg11, yg12, yg21, yg22)
 *  28,29,30,31 : (yh11, yh12, yh21, yh22)
 *  32,33,34,35 : (xm11, xm12, xm21, xm22)
 *  36,37,38,39 : (ym11, ym12, ym21, ym22)
 *  40,41,42,43 : (zm11, zm12, zm21, zm22)
 */
void
scalars_nonewald_poly_full (int version,
			    double r,
			    double aa, double ab,
			    double *scalar)
{
  // zero clear for FT and higher
  int i;
  for (i = 0; i < 44; i ++) scalar [i] = 0.0;

  double r2 = r * r;

  double aa2 = aa * aa;
  double ab2 = ab * ab;

  double aa3 = aa2 * aa;
  double ab3 = ab2 * ab;

  // xa part
  scalar [0] = 1.0 / aa;
  scalar [1] = (1.5 - (aa2 + ab2) * 0.5 / r2) / r;
  scalar [2] = scalar [1];
  scalar [3] = 1.0 / ab;

  // ya part
  scalar [4] = 1.0 / aa;
  scalar [5] = (0.75 + (aa2 + ab2) * 0.25 / r2) / r;
  scalar [6] = scalar [5];
  scalar [7] = 1.0 / ab;

  // check point for F version
  if (version == 0) return;

  // yb part
  //scalar [8] = 0.0;
  scalar [9] = - 0.75 / r2;
  scalar[10] = scalar [9];
  //scalar[11] = 0.0;

  // xc part
  scalar[12] = 0.75 / aa3;
  scalar[13] = 0.75 / r2 / r;
  scalar[14] = scalar[13];
  scalar[15] = 0.75 / ab3;

  // yc part
  scalar[16] = 0.75 / aa3;
  scalar[17] = - 3.0 / 8.0 / r2 / r;
  scalar[18] = scalar[17];
  scalar[19] = 0.75 / ab3;

  // check point for FT version
  if (version == 1) return;

  // xg part
  //scalar[20] = 0.0;
  scalar[21] = (2.25 - 13.5 * (aa2 / 10.0 + ab2 / 6.0) / r2) / r2;
  scalar[22] = (2.25 - 13.5 * (ab2 / 10.0 + aa2 / 6.0) / r2) / r2;
  //scalar[23] = 0.0;

  // yg part
  //scalar[24] = 0.0;
  scalar[25] = 4.5 * (aa2 / 10.0 + ab2 / 6.0) / r2 / r2;
  scalar[26] = 4.5 * (ab2 / 10.0 + aa2 / 6.0) / r2 / r2;
  //scalar[27] = 0.0;

  // yh part
  //scalar[28] = 0.0;
  scalar[29] = - 9.0 / 8.0 / r2 / r;
  scalar[30] = scalar[29];
  //scalar[31] = 0.0;

  // xm part
  scalar[32] = 0.9 / aa3;
  scalar[33] = (- 4.5 + 5.4 * (aa2 + ab2) / r2) / r2 / r;
  scalar[34] = scalar[33];
  scalar[35] = 0.9 / ab3;

  // ym part
  scalar[36] = 0.9 / aa3;
  scalar[37] = (2.25 - 3.6 * (aa2 + ab2) / r2) / r2 / r;
  scalar[38] = scalar[37];
  scalar[39] = 0.9 / ab3;

  // zm part
  scalar[40] = 0.9 / aa3;
  scalar[41] = 0.9 * (aa2 + ab2) / r2 / r2 / r;
  scalar[42] = scalar[41];
  scalar[43] = 0.9 / ab3;
}

/* ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_nonewald_3all (int n, const double *x, double *y, void *user_data)
{
  struct stokes * sys;

  double mob [22];

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
  // the default values are for no-slip case
  double self_a = 1.0;
  double self_c = 0.75;
  double self_m = 0.9;

  for (i = 0; i < sys->np; i++)
    {
      if (sys->slip != NULL) // slip
	{
	  self_a = 1.0  * sys->slip_G32 [i]; // Lambda(3,2) = 1/Lambda(2,3)
	  self_c = 0.75 * sys->slip_G30 [i]; // Lambda(3,0) = 1/Lambda(0,3)
	  self_m = 0.9  * sys->slip_G52 [i]; // Lambda(5,2) = 1/Lambda(2,5)
	}

      if (sys->version == 0) // F version
	{
	  matrix_f_atimes (x + i*3, y + i*3,
			   0.0, 0.0, 0.0,
			   self_a, self_a); // a part
	}
      else if (sys->version == 1) // FT version
	{
	  matrix_ft_atimes (x + i*6, y + i*6,
			    0.0, 0.0, 0.0,
			    self_a, self_a, // a part
			    0.0,
			    self_c, self_c); // c part
	}
      else // FTS version
	{
	  matrix_fts_atimes (x + i*11, y + i*11,
			     0.0, 0.0, 0.0,
			     self_a, self_a, // a part
			     0.0,
			     self_c, self_c, // c part
			     0.0, 0.0,
			     0.0,
			     self_m, self_m, self_m); // m part
	}
    }

  /* loop for all particle-particle pairs without the self */
  for (i = 0; i < sys->np; i++)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;

      // used in the slip case
      double ai;
      if (sys->a == NULL) ai = 1.0;
      else                ai = sys->a [i];

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

	  if (sys->slip == NULL // no-slip
	      && sys->a == NULL) // monodisperse
	    {
	      scalars_nonewald (sys->version, r, mob);
	    }
	  else if (sys->slip == NULL) // no-slip polydisperse
	    {
	      scalars_nonewald_poly (sys->version, r,
				     sys->a[i], sys->a[j],
				     mob);
	      scalars_mob_poly_scale_SD (sys->version,
					 sys->a[i],
					 mob);
	      // now mob is in the SD form
	    }
	  else // slip (for both mono- and poly-disperse systems)
	    {
	      // use the effective radius in slip->a[] defined by
	      // slip->a[i] = sys->a[i] * sqrt (Lambda^(i)(0,2)).
	      // for both monodisperse and polydisperse.
	      // (slip->a[] is defined for both cases properly.)
	      scalars_nonewald_poly (sys->version, r,
				     sys->slip_a[i], sys->slip_a[j],
				     mob);
	      // scale is done by the real radius
	      scalars_mob_poly_scale_SD (sys->version,
					 //sys->a[i],
					 ai,
					 mob);
	      // now scalars are in the SD form
	    }

	  // note that interaction (i,j) should be for (U[i], F[j])
	  if (sys->version == 0) // F version
	    {
	      matrix_f_atimes (x + j*3, y + i*3,
			       ex, ey, ez,
			       mob[0], mob[1]); // xa, ya
	    }
	  else if (sys->version == 1) // FT version
	    {
	      matrix_ft_atimes (x + j*6, y + i*6,
				ex, ey, ez,
				mob[0], mob[1],  // xa, ya
				mob[2],          // yb
				mob[3], mob[4]); // xc, yc
	    }
	  else // FTS version
	    {
	      matrix_fts_atimes (x + j*11, y + i*11,
				 ex, ey, ez,
				 mob[0], mob[1], // xa, ya
				 mob[2],         // yb
				 mob[3], mob[4], // xc, yc
				 mob[5], mob[6], // xg, yg
				 mob[7],         // yh
				 mob[8], mob[9], mob[10]); // xm, ym, zm
	    }
	}
    }
}

/* make plain mobility matrix for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_nonewald_3all (struct stokes *sys, double *mat)
{
  double mob [22];

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
  // the default values are for no-slip case
  double self_a = 1.0;
  double self_c = 0.75;
  double self_m = 0.9;

  for (i = 0; i < sys->np; i++)
    {
      if (sys->slip != NULL) // slip
	{
	  self_a = 1.0  * sys->slip_G32 [i]; // Lambda(3,2) = 1/Lambda(2,3)
	  self_c = 0.75 * sys->slip_G30 [i]; // Lambda(3,0) = 1/Lambda(0,3)
	  self_m = 0.9  * sys->slip_G52 [i]; // Lambda(5,2) = 1/Lambda(2,5)
	}

      if (sys->version == 0) // F version
	{
	  matrix_f_ij (i, i,
		       0.0, 0.0, 0.0,
		       self_a, self_a, // a part
		       n, mat);
	}
      else if (sys->version == 1) // FT version
	{
	  matrix_ft_ij (i, i,
			0.0, 0.0, 0.0,
			self_a, self_a, // a part
			0.0,
			self_c, self_c, // c part
			n, mat);
	}
      else // FTS version
	{
	  matrix_fts_ij (i, i,
			 0.0, 0.0, 0.0,
			 self_a, self_a, // a part
			 0.0,
			 self_c, self_c, // c part
			 0.0, 0.0,
			 0.0,
			 self_m, self_m, self_m, // m part
			 n, mat);
	}
    }

  /* loop for all particle-particle pairs without the self */
  for (i = 0; i < sys->np; i++)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;

      // used in the slip case
      double ai;
      if (sys->a == NULL) ai = 1.0;
      else                ai = sys->a [i];

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

	  if (sys->slip == NULL // no-slip
	      && sys->a == NULL) // monodisperse
	    {
	      scalars_nonewald (sys->version, r, mob);
	    }
	  else if (sys->slip == NULL) // no-slip polydisperse
	    {
	      scalars_nonewald_poly (sys->version, r,
				     sys->a[i], sys->a[j],
				     mob);
	      scalars_mob_poly_scale_SD (sys->version,
					 sys->a[i],
					 mob);
	      // now mob is in the SD form
	    }
	  else // slip (for both mono- and poly-disperse systems)
	    {
	      // use the effective radius in slip->a[] defined by
	      // slip->a[i] = sys->a[i] * sqrt (Lambda^(i)(0,2)).
	      // for both monodisperse and polydisperse.
	      // (slip->a[] is defined for both cases properly.)
	      scalars_nonewald_poly (sys->version, r,
				     sys->slip_a[i], sys->slip_a[j],
				     mob);
	      // scale is done by the real radius
	      scalars_mob_poly_scale_SD (sys->version,
					 //sys->a[i],
					 ai,
					 mob);
	      // now scalars are in the SD form
	    }

	  if (sys->version == 0) // F version
	    {
	      matrix_f_ij (i, j,
			   ex, ey, ez,
			   mob[0], mob[1], // xa, ya
			   n, mat);
	    }
	  else if (sys->version == 1) // FT version
	    {
	      matrix_ft_ij (i, j,
			    ex, ey, ez,
			    mob[0], mob[1], // xa, ya
			    mob[2],         // yb
			    mob[3], mob[4], // xc, yc
			    n, mat);
	    }
	  else // FTS version
	    {
	      matrix_fts_ij (i, j,
			     ex, ey, ez,
			     mob[0], mob[1], // xa, ya
			     mob[2],         // yb
			     mob[3], mob[4], // xc, yc
			     mob[5], mob[6], // xg, yg
			     mob[7],         // yh
			     mob[8], mob[9], mob[10], // xm, ym, zm
			     n, mat);
	    }
	}
    }
}

/* ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * through matrix procedure
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_nonewald_3all_matrix (int n, const double *x,
			     double *y, void *user_data)
{
  struct stokes * sys;
  int np;
  double * mat;


  sys = (struct stokes *) user_data;

  np = sys->np;
  mat = (double *) malloc (sizeof (double) * n * n);
  if (mat == NULL)
    {
      fprintf (stderr, "allocation error in atimes_nonewald_3fts_matrix ().\n");
      exit (1);
    }

  make_matrix_mob_nonewald_3all (sys, mat);
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
