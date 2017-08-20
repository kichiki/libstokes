/* bug fixing for polydisperse systems
 * utility for non-Ewald routines
 * Copyright (C) 2017 Kengo Ichiki <kengoichiki@gmail.com>
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

#include "stokes.h" // struct stokes
#include "f.h"      // matrix_f_ij() matrix_f_atimes()
#include "ft.h"     // matrix_ft_ij() matrix_ft_atimes()
#include "fts.h"    // matrix_fts_ij() matrix_fts_atimes()
#include "non-ewald.h" // scalars_nonewald() scalars_nonewald_poly()

#include "non-ewald-new.h"


/* fixed version for polydisperse systems of
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
atimes_nonewald_3all_new (int n, const double *x, double *y, void *user_data)
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
  double self_a_0 = 1.0;
  double self_c_0 = 0.75;
  double self_m_0 = 0.9;
  for (i = 0; i < sys->np; i++)
    {
      double self_a = self_a_0;
      double self_c = self_c_0;
      double self_m = self_m_0;

      if (sys->slip != NULL) // slip
	{
	  self_a = 1.0  * sys->slip_G32 [i]; // Lambda(3,2) = 1/Lambda(2,3)
	  self_c = 0.75 * sys->slip_G30 [i]; // Lambda(3,0) = 1/Lambda(0,3)
	  self_m = 0.9  * sys->slip_G52 [i]; // Lambda(5,2) = 1/Lambda(2,5)
	}

      // adjust scaling for polydisperse systems
      if (sys->a != NULL)
	{
	  double ai = sys->a[i];
	  double ai3 = ai * ai * ai;
	  self_a /= ai;
	  self_c /= ai3;
	  self_m /= ai3;
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
	      /*
	      scalars_mob_poly_scale_SD (sys->version,
					 sys->a[i],
					 mob);
	      // now mob is in the SD form
	      */
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
	      /*
	      // scale is done by the real radius
	      scalars_mob_poly_scale_SD (sys->version,
					 //sys->a[i],
					 ai,
					 mob);
	      // now scalars are in the SD form
	      */
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
