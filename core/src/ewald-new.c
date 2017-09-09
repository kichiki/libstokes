/* bug fixing for polydisperse systems
 * utility for Ewald summation calculation
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
//#include <stdio.h> /* for printf() */
//#include <stdlib.h> /* for exit() */
//#include "memory-check.h" // CHECK_MALLOC

#include "stokes.h" // struct stokes
#include "bench.h"
#include "f.h"      // matrix_f_ij() matrix_f_atimes()
#include "ft.h"     // matrix_ft_ij() matrix_ft_atimes()
#include "fts.h"    // matrix_fts_ij() matrix_fts_atimes()
#include "non-ewald-new.h" // atimes_nonewald_3all_new()

#include "ewald.h"  // scalars_ewald_real() scalars_ewald_real_poly()

#include "ewald-new.h"




/* fixed version for polydisperse systems of
 * ATIMES of calc ewald-summed mobility for F/FT/FTS versions
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
atimes_3all_new (int n, const double *x, double *y, void * user_data)
{
  struct stokes *sys = (struct stokes *)user_data;

  if (sys->periodic == 0)
    {
      // non-periodic
      atimes_nonewald_3all_new (n, x, y, user_data);
    }
  else
    {
      // periodic
      atimes_ewald_3all_new (n, x, y, user_data);
    }
}




/* fixed version for polydisperse systems of
 * ATIMES of calc ewald-summed mobility for F/FT/FTS versions
 * with the ewald table
 * INPUT
 *  n := np*3 (F), np*6 (FT), or np*11 (FTS)
 *  x [n] : F, FT, or FTS
 *  user_data = (struct stokes *) sys
 * OUTPUT
 *  y [n] : U, UO, or UOE
 */
void
atimes_ewald_3all_new (int n, const double *x, double *y, void * user_data)
{
  struct stokes *sys = (struct stokes *)user_data;

  double cpu0, cpu; /* for ptime_ms_d() */

  double xa, ya; 
  double yb;
  double xc, yc;
  double xg, yg;
  double yh;
  double xm, ym, zm;

  int i, j;
  /* clear result */
  for (i = 0; i < n; i ++)
    {
      y [i] = 0.0;
    }

  /* diagonal part ( self part ) */
  if (sys->slip == NULL // no-slip
      && sys->a == NULL) // monodisperse
    {
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
  else // slip or no-slip polydisperse
    {
      double self_a;
      double self_c;
      double self_m;

      for (i = 0; i < sys->np; i++)
	{
	  double a;
	  double a2;
	  double a3;

	  if (sys->slip == NULL) // no-slip
	    {
	      a = sys->a [i];
	      a2 = a*a;
	      a3 = a2*a;

	      /*
	      self_a = 1.0
		- a * sys->xiaspi * (6.0 - 40.0 / 3.0 * sys->xia2 * a2);
	      self_c = 0.75
		- a3 * sys->xiaspi * sys->xia2 * 10.0;
	      self_m = 0.9
		- a3 * sys->xiaspi * sys->xia2
		* (12.0 - 30.24 * sys->xia2 * a2);
	      // in SD scaling
	      */

	      // *** DIDN'T CHECK YET ***
	      self_a = 1.0 / a
		- sys->xiaspi * (6.0 - 40.0 / 3.0 * sys->xia2 * a2);
	      self_c = 0.75 / a3
		- sys->xiaspi * sys->xia2 * 10.0;
	      self_m = 0.9 / a3
		- sys->xiaspi * sys->xia2
		* (12.0 - 30.24 * sys->xia2 * a2);
	    }
	  else // slip (for both mono- and poly-disperse systems)
	    {
	      if (sys->a == NULL) a = 1.0;
	      else                a = sys->a [i];
	      a2 = a*a;
	      a3 = a2*a;

	      /*
	      self_a = 1.0  * sys->slip_G32 [i] // Lambda(3,2) = 1/Lambda(2,3)
		- a * sys->xiaspi * (6.0 - 40.0 / 3.0 * sys->xia2 * a2);
	      self_c = 0.75 * sys->slip_G30 [i] // Lambda(3,0) = 1/Lambda(0,3)
		- a3 * sys->xiaspi * sys->xia2 * 10.0;
	      self_m = 0.9  * sys->slip_G52 [i] // Lambda(5,2) = 1/Lambda(2,5)
		- a3 * sys->xiaspi * sys->xia2
		* (12.0 - 30.24 * sys->xia2 * a2);
	      // in SD scaling
	      */

	      // *** DIDN'T CHECK YET ***
	      self_a = 1.0 / a
		* sys->slip_G32 [i] // Lambda(3,2) = 1/Lambda(2,3)
		- sys->xiaspi * (6.0 - 40.0 / 3.0 * sys->xia2 * a2);
	      self_c = 0.75 / a3
		* sys->slip_G30 [i] // Lambda(3,0) = 1/Lambda(0,3)
		- sys->xiaspi * sys->xia2 * 10.0;
	      self_m = 0.9 / a3
		* sys->slip_G52 [i] // Lambda(5,2) = 1/Lambda(2,5)
		- sys->xiaspi * sys->xia2
		* (12.0 - 30.24 * sys->xia2 * a2);
	    }

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
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;

      // used in the slip case
      double ai;
      if (sys->a == NULL) ai = 1.0;
      else                ai = sys->a [i];

      for (j = 0; j < sys->np; j++)
	{
	  double aj;
	  if (sys->a == NULL) aj = 1.0;
	  else                aj = sys->a [j];

	  int jx = j * 3;
	  int jy = jx + 1;
	  int jz = jx + 2;

	  int m;
	  for (m = 0; m < sys->nr; m ++)
	    {
	      double xx = sys->pos [jx] - sys->pos [ix] + sys->rlx [m];
	      double yy = sys->pos [jy] - sys->pos [iy] + sys->rly [m];
	      double zz = sys->pos [jz] - sys->pos [iz] + sys->rlz [m];

	      // shift for shear
	      if (sys->shear_mode == 1)
		{
		  xx += (double)sys->rmy[m] * sys->shear_shift;
		}
	      else if (sys->shear_mode == 2)
		{
		  xx += (double)sys->rmz[m] * sys->shear_shift;
		}

	      double rr = xx * xx + yy * yy + zz * zz;

	      if (rr == 0.0) continue; // to exclude the self part
	      if (rr > sys->rmax2) continue;

	      double r = sqrt (rr);
	      double rmin = (ai + aj) * sys->rmin;
	      if (r < rmin) r = rmin;

	      double ex = xx / r;
	      double ey = yy / r;
	      double ez = zz / r;

	      if (sys->slip == NULL // no-slip
		  && sys->a == NULL) // monodisperse
		{
		  scalars_ewald_real (sys->version,
				      sys->xi, r,
				      &xa, &ya,
				      &yb,
				      &xc, &yc,
				      &xg, &yg,
				      &yh,
				      &xm, &ym, &zm);
		}
	      else if (sys->slip == NULL) // no-slip polydisperse
		{
		  scalars_ewald_real_poly (sys->version,
					   sys->xi, r,
					   sys->a[i], sys->a[j],
					   &xa, &ya,
					   &yb,
					   &xc, &yc,
					   &xg, &yg,
					   &yh,
					   &xm, &ym, &zm);
		  /*
		  scalars_mob_poly_scale_SD_ (sys->version,
					      sys->a[i],
					      &xa, &ya,
					      &yb,
					      &xc, &yc,
					      &xg, &yg,
					      &yh,
					      &xm, &ym, &zm);
		  // now scalars are in the SD form
		  */
		}
	      else // slip (for both mono- and poly-disperse systems)
		{
		  // use the effective radius in sys->slip_a[] defined by
		  // sys->slip_a[i] := sys->a[i] * sqrt (Lambda^(i)(0,2))
		  // for both monodisperse and polydisperse.
		  // (sys->slip_a[] is defined for both cases properly.)
		  scalars_ewald_real_poly (sys->version,
					   sys->xi, r,
					   sys->slip_a[i], sys->slip_a[j],
					   &xa, &ya,
					   &yb,
					   &xc, &yc,
					   &xg, &yg,
					   &yh,
					   &xm, &ym, &zm);
		  /*
		  scalars_mob_poly_scale_SD_ (sys->version,
					      //sys->a[i],
					      ai,
					      &xa, &ya,
					      &yb,
					      &xc, &yc,
					      &xg, &yg,
					      &yh,
					      &xm, &ym, &zm);
		  // now scalars are in the SD form
		  */
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

  /* Second Ewald part ( reciprocal space )
   * -- table version
  for (i = 0; i < sys->np; i++)
    {
      int ix = i * 3;
      int iy = ix + 1;
      int iz = ix + 2;
      for (j = 0; j < sys->np; j++)
	{
	  int jx = j * 3;
	  int jy = jx + 1;
	  int jz = jx + 2;

	  double xx = sys->pos [jx] - sys->pos [ix];
	  double yy = sys->pos [jy] - sys->pos [iy];
	  double zz = sys->pos [jz] - sys->pos [iz];

	  int m;
	  for (m = 0; m < sys->nk; m ++)
	    {
	      double ex = sys->ex[m];
	      double ey = sys->ey[m];
	      double ez = sys->ez[m];

	      double k1 = sys->k1[m];
	      double k2 = sys->k2[m];
	      double k3 = sys->k3[m];

	      if (sys->slip == NULL // no-slip
		  && sys->a == NULL) // monodisperse
		{
		  ya = sys->ya[m];
		  yb = sys->yb[m];
		  yc = sys->yc[m];
		  yg = sys->yg[m];
		  yh = sys->yh[m];
		  ym = sys->ym[m];
		}
	      else // slip or no-slip polydisperse
		{
		  double aa;
		  double aa3;
		  double aa2;
		  if (sys->a == NULL) aa = 1.0;
		  else                aa = sys->a [i];
		  aa2 = aa * aa;
		  aa3 = aa2 * aa;

		  // a^2 in nabla^2 term
		  double a2a;
		  double a2b;
		  if (sys->slip == NULL) // no-slip
		    {
		      a2a = aa2;

		      if (sys->a == NULL) a2b = 1.0;
		      else                a2b = sys->a [j];
		      a2b *= a2b;
		    }
		  else // slip
		    {
		      a2a = sys->slip_a [i];
		      a2a *= a2a;

		      a2b = sys->slip_a [j];
		      a2b *= a2b;
		    }

		  double k = sys->k [m];
		  double kk = k * k;
		  double k4z = kk / 4.0 / sys->xi2;
		  double kexp = sys->pivol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  ya = aa * 6.0 * (1.0 - kk * (a2a+a2b) / 6.0) * kexp;
		  yb = aa2 * 3.0 * k * kexp;
		  yc = aa3 * 3.0 / 2.0 * kk * kexp;
		  yg = aa2 * 3.0 * (1.0 - kk * (a2a/10.0 + a2b/6.0)) * k * kexp;
		  yh = aa3 * 3.0 / 2.0 * kk * kexp;
		  ym = aa3 * 3.0 * (1.0 - kk * (a2a+a2b) / 10.0) * kk * kexp;
		  // in SD scaling
		}

	      double cf = cos (+ k1 * xx
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
		  double sf = - sin (+ k1 * xx
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
		  double sf = - sin (+ k1 * xx
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
   */

  /* Second Ewald part ( reciprocal space ) */
  int m1, m2, m3;
  double k0x = 2.0 * M_PI / sys->lx;
  double k0y = 2.0 * M_PI / sys->ly;
  double k0z = 2.0 * M_PI / sys->lz;
  double sl = sys->shear_shift / sys->lx;
  for (m1 = - sys->kmaxx; m1 <= sys->kmaxx; m1++)
    {
      double k1 = k0x * (double) m1;
      for (m2 = - sys->kmaxy; m2 <= sys->kmaxy; m2++)
	{
	  double k2 = k0y * (double) m2;
	  for (m3 = - sys->kmaxz; m3 <= sys->kmaxz; m3++)
	    {
	      double k3 = k0z * (double) m3;
	      if (m1 != 0 || m2 != 0 || m3 != 0)
		{
		  // shift for shear
		  if (sys->shear_mode == 1)
		    {
		      k2 = k0y * ((double)m2 - sl * (double)m1);
		    }
		  else if (sys->shear_mode == 2)
		    {
		      k3 = k0z * ((double)m3 - sl * (double)m1);
		    }

		  double kk = k1 * k1 + k2 * k2 + k3 * k3;
		  double k = sqrt (kk);

		  if (sys->kmax != 0.0 && k > sys->kmax) continue;

		  double k4z = kk / 4.0 / sys->xi2;
		  double kexp = sys->pivol
		    * (1.0 + k4z * (1.0 + 2.0 * k4z))
		    / kk * exp (- k4z);

		  double ex = k1 / k;
		  double ey = k2 / k;
		  double ez = k3 / k;

		  for (i = 0; i < sys->np; i++)
		    {
		      int ix = i * 3;
		      int iy = ix + 1;
		      int iz = ix + 2;
		      for (j = 0; j < sys->np; j++)
			{
			  int jx = j * 3;
			  int jy = jx + 1;
			  int jz = jx + 2;

			  if (sys->slip == NULL // no-slip
			      && sys->a == NULL) // monodisperse
			    {
			      ya = 6.0 * (1.0 - kk / 3.0) * kexp;
			      yb = 3.0 * k * kexp;
			      yc = 3.0 / 2.0 * kk * kexp;
			      yg = 3.0 * (1.0 - 4.0 / 15.0 * kk) * k * kexp;
			      yh = 3.0 / 2.0 * kk * kexp;
			      ym = 3.0 * (1.0 - kk / 5.0) * kk * kexp;
			    }
			  else // slip or no-slip polydisperse
			    {
			      double aa;
			      double aa3;
			      double aa2;
			      if (sys->a == NULL) aa = 1.0;
			      else                aa = sys->a [i];
			      aa2 = aa * aa;
			      aa3 = aa2 * aa;

			      // a^2 in nabla^2 term
			      double a2a;
			      double a2b;
			      if (sys->slip == NULL) // no-slip
				{
				  a2a = aa2;

				  if (sys->a == NULL) a2b = 1.0;
				  else                a2b = sys->a [j];
				  a2b *= a2b;
				}
			      else // slip
				{
				  a2a = sys->slip_a [i];
				  a2a *= a2a;

				  a2b = sys->slip_a [j];
				  a2b *= a2b;
				}

			      /*
			      ya = aa * 6.0 * (1.0 - kk * (a2a+a2b) / 6.0)
				* kexp;
			      yb = aa2 * 3.0 * k * kexp;
			      yc = aa3 * 3.0 / 2.0 * kk * kexp;
			      yg = aa2 * 3.0 * (1.0 - kk * (a2a/10.0 + a2b/6.0))
				* k * kexp;
			      yh = aa3 * 3.0 / 2.0 * kk * kexp;
			      ym = aa3 * 3.0 * (1.0 - kk * (a2a+a2b) / 10.0)
				* kk * kexp;
			      // in SD scaling
			      */
			      ya = 6.0 * (1.0 - kk * (a2a+a2b) / 6.0)
				* kexp;
			      yb = 3.0 * k * kexp;
			      yc = 3.0 / 2.0 * kk * kexp;
			      yg = 3.0 * (1.0 - kk * (a2a/10.0 + a2b/6.0))
				* k * kexp;
			      yh = 3.0 / 2.0 * kk * kexp;
			      ym = 3.0 * (1.0 - kk * (a2a+a2b) / 10.0)
				* kk * kexp;
			    }
      
			  double xx = sys->pos [jx] - sys->pos [ix];
			  double yy = sys->pos [jy] - sys->pos [iy];
			  double zz = sys->pos [jz] - sys->pos [iz];

			  double cf = cos (+ k1 * xx
					   + k2 * yy
					   + k3 * zz);

			  // note that interaction (i,j) should be
			  // for (U[i], F[j])
			  if (sys->version == 0) // F version
			    {
			      matrix_f_atimes (x + j*3, y + i*3,
					       ex, ey, ez,
					       0.0, cf * ya);
			    }
			  else if (sys->version == 1) // FT version
			    {
			      double sf = - sin (+ k1 * xx
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
			      double sf = - sin (+ k1 * xx
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
	    }
	}
    }

  /* for xi code to measure CPU times */
  cpu = ptime_ms_d ();
  sys->cpu3 = cpu - cpu0;
  sys->cpu1 = sys->cpu2 + sys->cpu3;
}
