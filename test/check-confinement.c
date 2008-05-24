/* test code for confinement.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-confinement.c,v 1.1 2008/05/24 05:52:57 kichiki Exp $
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory-check.h"
#include "check.h" // compare()

#include "confinement.h"


int
check_CF_init (int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_CF_init : start\n");
    }

  int check = 0;
  double max = 0.0;

  struct confinement *cf;
  int type;
  double R;
  double r;
  double x, y, z;
  double R2;
  double L;
  int flag_LJ;
  double e;
  double r0;

  double theta;
  double theta2;
  double center1; // x component of the center of left cavity
  double center2; // x component of the center of right cavity
  double Lx;
  double Ly;
  double Lz;

  // type == 0 : sphere
  type = 0;
  R  = 10.0;
  r  = 0.0;
  x  = 0.0;
  y  = 0.0;
  z  = 0.0;
  R2 = 0.0;
  L  = 0.0;
  flag_LJ = 0;
  e  = 10.0;
  r0 = 0.8;
  cf = CF_init (type, R, r, x, y, z, R2, L, flag_LJ, e, r0);
  CHECK_MALLOC (cf, "check_CF_init");

  check += compare_max (cf->R, R, " sphere R", verbose, tiny, &max);
  check += compare_max (cf->r, r, " sphere r", verbose, tiny, &max);
  check += compare_max (cf->x, x, " sphere x", verbose, tiny, &max);
  check += compare_max (cf->y, y, " sphere y", verbose, tiny, &max);
  check += compare_max (cf->z, z, " sphere z", verbose, tiny, &max);
  check += compare_max (cf->R2, R2, " sphere R2", verbose, tiny, &max);
  check += compare_max (cf->L, L, " sphere L", verbose, tiny, &max);
  check += compare_max ((double)cf->flag_LJ, (double)flag_LJ, " sphere LJ", verbose, tiny, &max);
  check += compare_max (cf->e, e, " sphere LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, r0, " sphere LJ-r0", verbose, tiny, &max);

  theta = asin (r / R);
  check += compare_max (cf->theta, theta, " sphere theta", verbose, tiny, &max);
  check += compare_max (cf->theta2, 0.0, " sphere theta2", verbose, tiny, &max);
  check += compare_max (cf->center1, 0.0, " sphere center1", verbose, tiny, &max);
  check += compare_max (cf->center2, 0.0, " sphere center2", verbose, tiny, &max);
  check += compare_max (cf->Lx, 0.0, " sphere Lx", verbose, tiny, &max);
  check += compare_max (cf->Ly, 0.0, " sphere Ly", verbose, tiny, &max);
  check += compare_max (cf->Lz, 0.0, " sphere Lz", verbose, tiny, &max);

  CF_free (cf);


  // type == 1 : sphere + hole
  type = 1;
  R  = 5.0;
  r  = 1.25;
  x  = 0.0;
  y  = 0.0;
  z  = 0.0;
  R2 = 0.0;
  L  = 0.0;
  flag_LJ = 1;
  e  = 9.0;
  r0 = 0.7;
  cf = CF_init (type, R, r, x, y, z, R2, L, flag_LJ, e, r0);
  CHECK_MALLOC (cf, "check_CF_init");

  check += compare_max (cf->R, R, " sph+hole R", verbose, tiny, &max);
  check += compare_max (cf->r, r, " sph+hole r", verbose, tiny, &max);
  check += compare_max (cf->x, x, " sph+hole x", verbose, tiny, &max);
  check += compare_max (cf->y, y, " sph+hole y", verbose, tiny, &max);
  check += compare_max (cf->z, z, " sph+hole z", verbose, tiny, &max);
  check += compare_max (cf->R2, R2, " sph+hole R2", verbose, tiny, &max);
  check += compare_max (cf->L, L, " sph+hole L", verbose, tiny, &max);
  check += compare_max ((double)cf->flag_LJ, (double)flag_LJ, " sph+hole LJ", verbose, tiny, &max);
  check += compare_max (cf->e, e, " sph+hole LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, r0, " sph+hole LJ-r0", verbose, tiny, &max);

  theta = asin (r / R);
  check += compare_max (cf->theta, theta, " sph+hole theta", verbose, tiny, &max);
  check += compare_max (cf->theta2, 0.0, " sph+hole theta2", verbose, tiny, &max);
  check += compare_max (cf->center1, 0.0, " sph+hole center1", verbose, tiny, &max);
  check += compare_max (cf->center2, 0.0, " sph+hole center2", verbose, tiny, &max);
  check += compare_max (cf->Lx, 0.0, " sph+hole Lx", verbose, tiny, &max);
  check += compare_max (cf->Ly, 0.0, " sph+hole Ly", verbose, tiny, &max);
  check += compare_max (cf->Lz, 0.0, " sph+hole Lz", verbose, tiny, &max);

  CF_free (cf);


  // type == 2 : cylinder
  type = 2;
  R  = 0.0;
  r  = 2.5;
  x  = 0.5;
  y  = 0.6;
  z  = 0.7;
  R2 = 0.0;
  L  = 0.0;
  flag_LJ = 0;
  e  = 11.0;
  r0 = 1.1;
  cf = CF_init (type, R, r, x, y, z, R2, L, flag_LJ, e, r0);
  CHECK_MALLOC (cf, "check_CF_init");

  check += compare_max (cf->R, R, " cylinder R", verbose, tiny, &max);
  check += compare_max (cf->r, r, " cylinder r", verbose, tiny, &max);
  check += compare_max (cf->x, x, " cylinder x", verbose, tiny, &max);
  check += compare_max (cf->y, y, " cylinder y", verbose, tiny, &max);
  check += compare_max (cf->z, z, " cylinder z", verbose, tiny, &max);
  check += compare_max (cf->R2, R2, " cylinder R2", verbose, tiny, &max);
  check += compare_max (cf->L, L, " cylinder L", verbose, tiny, &max);
  check += compare_max ((double)cf->flag_LJ, (double)flag_LJ, " cylinder LJ", verbose, tiny, &max);
  check += compare_max (cf->e, e, " cylinder LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, r0, " cylinder LJ-r0", verbose, tiny, &max);

  theta = asin (r / R);
  check += compare_max (cf->theta, theta, " cylinder theta", verbose, tiny, &max);
  check += compare_max (cf->theta2, 0.0, " cylinder theta2", verbose, tiny, &max);
  check += compare_max (cf->center1, 0.0, " cylinder center1", verbose, tiny, &max);
  check += compare_max (cf->center2, 0.0, " cylinder center2", verbose, tiny, &max);
  check += compare_max (cf->Lx, 0.0, " cylinder Lx", verbose, tiny, &max);
  check += compare_max (cf->Ly, 0.0, " cylinder Ly", verbose, tiny, &max);
  check += compare_max (cf->Lz, 0.0, " cylinder Lz", verbose, tiny, &max);

  CF_free (cf);


  // type == 3 : dumbbell
  type = 3;
  R  = 5.0;
  r  = 1.3;
  x  = 0.0;
  y  = 0.0;
  z  = 0.0;
  R2 = 10.0;
  L  = 2.0;
  flag_LJ = 1;
  e  = 15.0;
  r0 = 0.5;
  cf = CF_init (type, R, r, x, y, z, R2, L, flag_LJ, e, r0);
  CHECK_MALLOC (cf, "check_CF_init");

  check += compare_max (cf->R, R, " dumbbell R", verbose, tiny, &max);
  check += compare_max (cf->r, r, " dumbbell r", verbose, tiny, &max);
  check += compare_max (cf->x, x, " dumbbell x", verbose, tiny, &max);
  check += compare_max (cf->y, y, " dumbbell y", verbose, tiny, &max);
  check += compare_max (cf->z, z, " dumbbell z", verbose, tiny, &max);
  check += compare_max (cf->R2, R2, " dumbbell R2", verbose, tiny, &max);
  check += compare_max (cf->L, L, " dumbbell L", verbose, tiny, &max);
  check += compare_max ((double)cf->flag_LJ, (double)flag_LJ, " dumbbell LJ", verbose, tiny, &max);
  check += compare_max (cf->e, e, " dumbbell LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, r0, " dumbbell LJ-r0", verbose, tiny, &max);

  theta = asin (r / R);
  check += compare_max (cf->theta, theta, " dumbbell theta", verbose, tiny, &max);
  theta2 = asin (r / R2);
  check += compare_max (cf->theta2, theta2, " dumbbell theta2", verbose, tiny, &max);
  center1 = -(0.5 * L + R * cos (theta));
  center2 = +(0.5 * L + R2* cos (theta2));
  check += compare_max (cf->center1, center1, " dumbbell center1", verbose, tiny, &max);
  check += compare_max (cf->center2, center2, " dumbbell center2", verbose, tiny, &max);
  check += compare_max (cf->Lx, 0.0, " dumbbell Lx", verbose, tiny, &max);
  check += compare_max (cf->Ly, 0.0, " dumbbell Ly", verbose, tiny, &max);
  check += compare_max (cf->Lz, 0.0, " dumbbell Lz", verbose, tiny, &max);

  CF_free (cf);


  // type == 4 : hex2d
  type = 4;
  R  = 5.0;
  r  = 1.3;
  x  = 0.0;
  y  = 0.0;
  z  = 0.0;
  R2 = 0.0;
  L  = 5.5;
  flag_LJ = 0;
  e  = 11.1;
  r0 = 0.81;
  cf = CF_init (type, R, r, x, y, z, R2, L, flag_LJ, e, r0);
  CHECK_MALLOC (cf, "check_CF_init");

  check += compare_max (cf->R, R, " hex2d R", verbose, tiny, &max);
  check += compare_max (cf->r, r, " hex2d r", verbose, tiny, &max);
  check += compare_max (cf->x, x, " hex2d x", verbose, tiny, &max);
  check += compare_max (cf->y, y, " hex2d y", verbose, tiny, &max);
  check += compare_max (cf->z, z, " hex2d z", verbose, tiny, &max);
  check += compare_max (cf->R2, R2, " hex2d R2", verbose, tiny, &max);
  check += compare_max (cf->L, L, " hex2d L", verbose, tiny, &max);
  check += compare_max ((double)cf->flag_LJ, (double)flag_LJ, " hex2d LJ", verbose, tiny, &max);
  check += compare_max (cf->e, e, " hex2d LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, r0, " hex2d LJ-r0", verbose, tiny, &max);

  theta = asin (r / R);
  check += compare_max (cf->theta, theta, " hex2d theta", verbose, tiny, &max);
  check += compare_max (cf->theta2, 0.0, " hex2d theta2", verbose, tiny, &max);
  check += compare_max (cf->center1, 0.0, " hex2d center1", verbose, tiny, &max);
  check += compare_max (cf->center2, 0.0, " hex2d center2", verbose, tiny, &max);
  Lx = 0.5 * L;
  Ly = sqrt (3.0) * 0.5 * L;
  check += compare_max (cf->Lx, Lx, " hex2d Lx", verbose, tiny, &max);
  check += compare_max (cf->Ly, Ly, " hex2d Ly", verbose, tiny, &max);
  check += compare_max (cf->Lz, 0.0, " hex2d Lz", verbose, tiny, &max);

  CF_free (cf);

  // type == 5 : porous
  type = 5;
  R  = 5.0;
  r  = 0.0;
  x  = 0.0;
  y  = 0.0;
  z  = 0.0;
  R2 = 0.0;
  L  = 10.0;
  flag_LJ = 0;
  e  = 12.3;
  r0 = 0.76;
  cf = CF_init (type, R, r, x, y, z, R2, L, flag_LJ, e, r0);
  CHECK_MALLOC (cf, "check_CF_init");

  check += compare_max (cf->R, R, " porous R", verbose, tiny, &max);
  check += compare_max (cf->r, r, " porous r", verbose, tiny, &max);
  check += compare_max (cf->x, x, " porous x", verbose, tiny, &max);
  check += compare_max (cf->y, y, " porous y", verbose, tiny, &max);
  check += compare_max (cf->z, z, " porous z", verbose, tiny, &max);
  check += compare_max (cf->R2, R2, " porous R2", verbose, tiny, &max);
  check += compare_max (cf->L, L, " porous L", verbose, tiny, &max);
  check += compare_max ((double)cf->flag_LJ, (double)flag_LJ, " porous LJ", verbose, tiny, &max);
  check += compare_max (cf->e, e, " porous LJ-e", verbose, tiny, &max);
  check += compare_max (cf->r0, r0, " porous LJ-r0", verbose, tiny, &max);

  check += compare_max (cf->theta, 0.0, " porous theta", verbose, tiny, &max);
  check += compare_max (cf->theta2, 0.0, " porous theta2", verbose, tiny, &max);
  check += compare_max (cf->center1, 0.0, " porous center1", verbose, tiny, &max);
  check += compare_max (cf->center2, 0.0, " porous center2", verbose, tiny, &max);
  Lx = L;
  Ly = L * sqrt (3.0);
  Lz = L * sqrt (6.0) * 2.0 / 3.0;
  check += compare_max (cf->Lx, Lx, " porous Lx", verbose, tiny, &max);
  check += compare_max (cf->Ly, Ly, " porous Ly", verbose, tiny, &max);
  check += compare_max (cf->Lz, Lz, " porous Lz", verbose, tiny, &max);

  CF_free (cf);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}


// place the proto-type here because it is not in "confinement.h"
void
CF_sphere_calc_force (struct confinement *cf,
		      struct stokes *sys,
		      double *f,
		      int flag_add);

int
check_CF_sphere_calc_force (double R,
			    int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_CF_sphere_calc_force (R=%f) : start\n", R);
    }

  int check = 0;
  double max = 0.0;

  struct confinement *cf
    = CF_init (0, // sphere
	       R,
	       0.0, // r,
	       0.0, // x,
	       0.0, // y,
	       0.0, // z,
	       0.0, // R2,
	       0.0, // L,
	       1, // flag_LJ
	       10.0, // e
	       R + 2.0); // r0
  // note now r is in [0, 2R) and a = 1, so it is to prevent overflow.
  CHECK_MALLOC (cf, "check_CF_sphere_calc_force");

  // set struct stokes *sys
  int np = 1;
  double *f = (double *)malloc (sizeof (double) * np * 3);
  CHECK_MALLOC (f, "check_CF_sphere_calc_force");

  struct stokes *sys = stokes_init ();
  CHECK_MALLOC (sys, "check_CF_sphere_calc_force");
  stokes_set_np (sys, np, np);

  char label[80];
  int i;
  for (i = 0; i < 100; i ++)
    {
      double r = (double)i / 100.0 * 2.0 * R;
      int j;
      for (j = 0; j < 10; j ++)
	{
	  double theta = M_PI * (double)j / 10.0;
	  int k;
	  for (k = 0; k < 20; k ++)
	    {
	      double phi = 2.0 * M_PI * (double) k / 20.0;
	      sys->pos[0] = r * sin(theta) * cos(phi);
	      sys->pos[1] = r * sin(theta) * sin(phi);
	      sys->pos[2] = r * cos(theta);

	      CF_sphere_calc_force (cf, sys, f, 0 /* zero-clear */);

	      if (r <= (R-1.0))
		{
		  sprintf (label, " sphere (%f,%5.1f,%5.1f) fx",
			   r, 180.0/M_PI*theta, 180.0/M_PI*phi);
		  check += compare_max (f[0], 0.0, label, verbose, tiny, &max);

		  sprintf (label, " sphere (%f,%5.1f,%5.1f) fy",
			   r, 180.0/M_PI*theta, 180.0/M_PI*phi);
		  check += compare_max (f[1], 0.0, label, verbose, tiny, &max);

		  sprintf (label, " sphere (%f,%5.1f,%5.1f) fz",
			   r, 180.0/M_PI*theta, 180.0/M_PI*phi);
		  check += compare_max (f[2], 0.0, label, verbose, tiny, &max);
		}
	      else
		{
		  double x = r - (R-1.0); // depth of overlap
		  double LJ_r = cf->r0 - x;
		  double fr = 12.0 * cf->e * pow (cf->r0 / LJ_r, 7.0)
		    * (pow (cf->r0 / LJ_r, 6.0) - 1.0);
		  double fx = -fr * sin(theta) * cos(phi);
		  double fy = -fr * sin(theta) * sin(phi);
		  double fz = -fr * cos(theta);

		  sprintf (label, " sphere (%f,%5.1f,%5.1f) fx",
			   r, 180.0/M_PI*theta, 180.0/M_PI*phi);
		  check += compare_max (f[0], fx, label, verbose, tiny, &max);

		  sprintf (label, " sphere (%f,%5.1f,%5.1f) fy",
			   r, 180.0/M_PI*theta, 180.0/M_PI*phi);
		  check += compare_max (f[1], fy, label, verbose, tiny, &max);

		  sprintf (label, " sphere (%f,%5.1f,%5.1f) fz",
			   r, 180.0/M_PI*theta, 180.0/M_PI*phi);
		  check += compare_max (f[2], fz, label, verbose, tiny, &max);
		}
	    }
	}
    }

  CF_free (cf);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
