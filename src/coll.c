/* collision handling routines
 * Copyright (C) 1995-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: coll.c,v 1.4 2007/05/04 01:19:53 kichiki Exp $
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
#include "stokes.h" // struct stokes
#include "coll.h"

/*
 * INPUT
 *  sys : system parameters
 *  x [np * 3] : position of particles for BOTH mobile and fixed particles
 *  v [nm * 3] : velocity of particles before collisions
 *               only mobile particles required.
 *               assumed that the velocity for fixed particles are zero.
 *  en : elastic constant
 * OUTPUT
 *  v [nm * 3] : velocity of particles after collisions
 */
void
collide_particles (struct stokes *sys,
		   const double *x, double *v, double en)
{
  int np, nm;
  
  double R = 4.0; /* diameter^2 */

  double e;
  double xx, yy, zz;
  double r;
  double vx, vy, vz;
  double B;
  int i, j, k;
  int ix, iy, iz;
  int jx, jy, jz;


  np = sys->np;
  nm = sys->nm;

  e = (1.0 + en) / 2.0;

  for (i = 0; i < nm; ++i)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;

      // loop for mobile particles
      for(j = 0; j < nm; ++j)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;

	  for (k = 0; k < 27; ++k)
	    {
	      xx = x [ix] - x [jx] - sys->llx [k];
	      yy = x [iy] - x [jy] - sys->lly [k];
	      zz = x [iz] - x [jz] - sys->llz [k];

	      vx = v [ix] - v [jx];
	      vy = v [iy] - v [jy];
	      vz = v [iz] - v [jz];

	      B = xx * vx + yy * vy + zz * vz;
	      r = xx * xx + yy * yy + zz * zz;

	      if (B < 0.0
		  && r < R)
		{
		  B = e * B / r;
		  v [ix] += - xx * B;
		  v [iy] += - yy * B;
		  v [iz] += - zz * B;
		  v [jx] -= - xx * B;
		  v [jy] -= - yy * B;
		  v [jz] -= - zz * B;
		}
	    }
	}

      // loop for fixed particles
      for(j = nm; j < np; ++j)
	{
	  jx = j * 3;
	  jy = jx + 1;
	  jz = jx + 2;

	  for (k = 0; k < 27; ++k)
	    {
	      xx = x [ix] - x [jx] - sys->llx [k];
	      yy = x [iy] - x [jy] - sys->lly [k];
	      zz = x [iz] - x [jz] - sys->llz [k];

	      vx = v [ix];
	      vy = v [iy];
	      vz = v [iz];

	      B = xx * vx + yy * vy + zz * vz;
	      r = xx * xx + yy * yy + zz * zz;

	      if (B < 0.0
		  && r < R)
		{
		  B = e * B / r;
		  v [ix] += - 2.0 * xx * B;
		  v [iy] += - 2.0 * yy * B;
		  v [iz] += - 2.0 * zz * B;
		}
	    }
	}
    }
}

/*
 * INPUT
 *  xx : relative position
 *  vx : relative velocity
 *  e := (1 + en) / 2,
 *      where 'en' is the restitution coefficient for normal component
 * OUTPUT
 *  *v : velocity to collide
 */
static void
collide_1 (double xx, double vx, double *v, double e)
{
  double R = 1.0; /* radius^2 */
  double B;
  double r;


  B = xx * vx;
  r = xx * xx;

  if (B < 0.0
      && r < R)
    {
      *v += - xx * 2.0 * e * B / r;
    }
}
/*
 * INPUT
 *  sys : system parameters
 *  x [np * 3] : position of particles for BOTH mobile and fixed particles
 *  v [nm * 3] : velocity of particles before collisions
 *               only mobile particles required.
 *               assumed that the velocity for fixed particles are zero.
 *  en : elastic constant
 *  x_wall : position of the wall
 *  v_wall : 
 * OUTPUT
 *  v [nm * 3] : velocity of particles after collisions
 */
void
collide_wall_x (struct stokes *sys,
		const double *x, double *v, double en,
		double x_wall, double v_wall)
{
  int nm;

  double e;
  double xx;
  double vx;
  int i;
  int ix;


  nm = sys->nm;

  e = (1.0 + en) / 2.0;

  for (i = 0; i < nm; ++i)
    {
      ix = i * 3;

      vx = v [ix] - v_wall;

      xx = x [ix] - x_wall;
      collide_1 (xx, vx, & v [ix], e);

      xx = x [ix] - x_wall - sys->lx;
      collide_1 (xx, vx, & v [ix], e);

      xx = x [ix] - x_wall + sys->lx;
      collide_1 (xx, vx, & v [ix], e);
    }
}

void
collide_wall_y (struct stokes *sys,
		const double *x, double *v, double en,
		double y_wall, double v_wall)
{
  int nm;

  double e;
  double yy;
  double vy;
  int i;
  int iy;


  nm = sys->nm;

  e = (1.0 + en) / 2.0;

  for (i = 0; i < nm; ++i)
    {
      iy = i * 3 + 1;

      vy = v [iy] - v_wall;

      yy = x [iy] - y_wall;
      collide_1 (yy, vy, & v [iy], e);

      yy = x [iy] - y_wall - sys->ly;
      collide_1 (yy, vy, & v [iy], e);

      yy = x [iy] - y_wall + sys->ly;
      collide_1 (yy, vy, & v [iy], e);
    }
}

void
collide_wall_z (struct stokes *sys,
		const double *x, double *v, double en,
		double z_wall, double v_wall)
{
  int nm;

  double e;
  double zz;
  double vz;
  int i;
  int iz;


  nm = sys->nm;

  e = (1.0 + en) / 2.0;

  for (i = 0; i < nm; ++i)
    {
      iz = i * 3 + 2;

      vz = v [iz] - v_wall;

      zz = x [iz] - z_wall;
      collide_1 (zz, vz, & v [iz], e);

      zz = x [iz] - z_wall - sys->lz;
      collide_1 (zz, vz, & v [iz], e);

      zz = x [iz] - z_wall + sys->lz;
      collide_1 (zz, vz, & v [iz], e);
    }
}

/*
 * INPUT
 *  sys : system parameters
 *  x [np * 2] : position of particles for BOTH mobile and fixed particles
 *  v [nm * 2] : velocity of particles before collisions
 *               only mobile particles required.
 *               assumed that the velocity for fixed particles are zero.
 *  en : elastic constant
 * OUTPUT
 *  v [nm * 2] : velocity of particles after collisions
 */
void
collide_particles_2d (struct stokes *sys,
		      const double *x, double *v, double en)
{
  int np, nm;

  double R = 4.0; /* diameter^2 */

  double e;
  double xx, yy;
  double r;
  double vx, vy;
  double B;
  int i, j, k;
  int ix, iy;
  int jx, jy;


  np = sys->np;
  nm = sys->nm;

  e = (1.0 + en) / 2.0;

  for (i = 0; i < nm; ++i)
    {
      ix = i * 2;
      iy = ix + 1;

      // loop for mobile particles
      for(j = 0; j < nm; ++j)
	{
	  jx = j * 2;
	  jy = jx + 1;

	  for (k = 0; k < 9; ++k)
	    {
	      xx = x [ix] - x [jx] - sys->llx [k];
	      yy = x [iy] - x [jy] - sys->lly [k];

	      vx = v [ix] - v [jx];
	      vy = v [iy] - v [jy];

	      B = xx * vx + yy * vy;
	      r = xx * xx + yy * yy;

	      if (B < 0.0
		  && r < R)
		{
		  B = e * B / r;
		  v [ix] += - xx * B;
		  v [iy] += - yy * B;
		  v [jx] -= - xx * B;
		  v [jy] -= - yy * B;
		}
	    }
	}

      // loop for fixed particles
      for(j = nm; j < np; ++j)
	{
	  jx = j * 2;
	  jy = jx + 1;

	  for (k = 0; k < 9; ++k)
	    {
	      xx = x [ix] - x [jx] - sys->llx [k];
	      yy = x [iy] - x [jy] - sys->lly [k];

	      vx = v [ix];
	      vy = v [iy];

	      B = xx * vx + yy * vy;
	      r = xx * xx + yy * yy;

	      if (B < 0.0
		  && r < R)
		{
		  B = e * B / r;
		  v [ix] += - 2.0 * xx * B;
		  v [iy] += - 2.0 * yy * B;
		}
	    }
	}
    }
}
