/* backup of bug fixing for polydisperse systems
 * utility for non-Ewald routines
 * Copyright (C) 2007-2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include <math.h> // sqrt()
//#include <stdio.h> /* for printf() */
//#include <stdlib.h> /* for exit() */
//#include "memory-check.h" // CHECK_MALLOC

#include "stokes.h"
#include "f.h"      // matrix_f_ij() matrix_f_atimes()
#include "ft.h"     // matrix_ft_ij() matrix_ft_atimes()
#include "fts.h"    // matrix_fts_ij() matrix_fts_atimes()
#include "matrix.h" // dot_prod_matrix() multiply_extmat_with_extvec_3fts()

#include "non-ewald.h"


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


/* backup of bug fixing for polydisperse systems of
 * ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_nonewald_3all_old
(int n, const double *x, double *y, void *user_data)
{
  struct stokes *sys = (struct stokes *)user_data;

  double mob [22];

  int i, j;
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
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // used in the slip case
      double ai;
      if (sys->a == NULL) ai = 1.0;
      else                ai = sys->a [i];

      for (j = 0; j < sys->np; j++)
	{
	  // exclude the self part
	  if (i == j) continue;

	  double aj;
	  if (sys->a == NULL) aj = 1.0;
	  else                aj = sys->a [j];

	  int jx = j * 3;
	  int jy = jx + 1;
	  int jz = jx + 2;

	  double xx = sys->pos [jx] - sys->pos [ix];
	  double yy = sys->pos [jy] - sys->pos [iy];
	  double zz = sys->pos [jz] - sys->pos [iz];
	  double rr = xx * xx + yy * yy + zz * zz;

	  double r = sqrt (rr);
	  double rmin = (ai + aj) * sys->rmin;
	  if (r < rmin) r = rmin;

	  double ex = xx / r;
	  double ey = yy / r;
	  double ez = zz / r;

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

/* backup of bug fixing for polydisperse systems of
 * make plain mobility matrix for F/FT/FTS versions for non-periodic case
 * for both monodisplerse and polydisperse systems (given by sys->a)
 * INPUT
 * sys : system parameters
 * OUTPUT
 *  mat [np * 3  * np * 3 ] : for F version
 *  mat [np * 6  * np * 6 ] : for FT version
 *  mat [np * 11 * np * 11] : for FTS version
 */
void
make_matrix_mob_nonewald_3all_old
(struct stokes *sys, double *mat)
{
  double mob [22];

  int n;
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

  int i, j;
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
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // used in the slip case
      double ai;
      if (sys->a == NULL) ai = 1.0;
      else                ai = sys->a [i];

      for (j = 0; j < sys->np; j++)
	{
	  // exclude the self part
	  if (i == j) continue;

	  double aj;
	  if (sys->a == NULL) aj = 1.0;
	  else                aj = sys->a [j];

	  int jx = j * 3;
	  int jy = jx + 1;
	  int jz = jx + 2;

	  double xx = sys->pos [jx] - sys->pos [ix];
	  double yy = sys->pos [jy] - sys->pos [iy];
	  double zz = sys->pos [jz] - sys->pos [iz];
	  double rr = xx * xx + yy * yy + zz * zz;

	  double r = sqrt (rr);
	  double rmin = (ai + aj) * sys->rmin;
	  if (r < rmin) r = rmin;

	  double ex = xx / r;
	  double ey = yy / r;
	  double ez = zz / r;

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

/* backup of bug fixing for polydisperse systems of
 * ATIMES of calc plain mobility for F/FT/FTS versions for non-periodic case
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
atimes_nonewald_3all_matrix_old
(int n, const double *x,
 double *y, void *user_data)
{
  struct stokes *sys = (struct stokes *)user_data;

  double *mat = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (mat, "atimes_nonewald_3all_matrix");

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
