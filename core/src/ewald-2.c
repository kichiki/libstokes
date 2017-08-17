/* utility routines for Ewald-summation code in 2D
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: ewald-2.c,v 1.1 2001/02/05 09:38:47 ichiki Exp $
 */
#include <math.h>
#include <stdlib.h> // srand48(), drand48()

#include "ewald-2.h"


static int
check_overlap (int np, double lx, double ly, double lz, double * pos);


/* initialize random configuration and the primary cell in monolayer
 * INPUT
 *  seed : seed for random
 *  phi : volume fraction of particles
 *  np  : # particles
 * OUTPUT
 *  pos [np * 3] : positions of particles
 *  lx, ly, lz   : geometry of the primary cell, where lz = 2.0 is fixed.
 */
void
init_config_random_2d (long seed,
		       double phi, int np,
		       double *pos, double *lx, double *ly, double *lz)
{
  int i;
  int ix, iy, iz;


  srand48 (seed);

  (*lz) = 2.0;
  (*lx) = (*ly)
    = pow (4.0 * M_PI * (double) np / phi / 3.0 / (*lz),
	   1.0 / 2.0);

  pos [0] = (*lx) * drand48 ();
  pos [1] = (*ly) * drand48 ();
  pos [2] = 0.0;
  for (i = 1; i < np; i ++)
    {
    retry_init_config_random_2d:
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;

      pos [ix] = (*lx) * drand48 ();
      pos [iy] = (*ly) * drand48 ();
      pos [iz] = 0.0;

      if (check_overlap (i, (*lx), (*ly), (*lz), pos) != 0)
	{
	  if (i > 1) --i;
	  goto retry_init_config_random_2d;
	}
    }
}

/* check overlap
 * INPUT
 * OUTPUT (return value)
 *  0 : no-overlap
 *  1 : overlap
 */
static int
check_overlap (int np, double lx, double ly, double lz, double * pos)
{
  int i;
  int ix, iy, iz;
  int j;
  int jx, jy, jz;
  int kx, ky, kz;

  double rr;
  double x, y, z;


  for (i = 0; i < np; ++i)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;
      for (j = i + 1; j < np; ++j)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;
	  for (kx = -1; kx <= 1; ++kx)
	    {
	      for (ky = -1; ky <= 1; ++ky)
		{
		  for (kz = -1; kz <= 1; ++kz)
		    {
		      x = pos [ix] - pos [jx] + (double) kx * lx;
		      y = pos [iy] - pos [jy] + (double) ky * ly;
		      z = pos [iz] - pos [jz] + (double) kz * lz;
		      rr = x * x
			+ y * y
			+ z * z;
		      if (rr <= 4.0)
			{
			  return 1;
			}
		    }
		}
	    }
	}
    }
  return 0;
}

