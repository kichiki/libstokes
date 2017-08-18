/* periodicity handling routines
 * Copyright (C) 1995-2007,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#include <math.h> /* M_PI */
#include <libstokes-core.h> // struct stokes

#include "periodicity.h"

/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  x [np * 3] : position to check; the range is [0, l[xyz])
 */
void
check_periodic (struct stokes *sys,
		double *x)
{
  int i;
  int ix, iy, iz;
  int nm;
  double lx, ly, lz;


  nm = sys->nm;
  lx = sys->lx;
  ly = sys->ly;
  lz = sys->lz;
  for (i = 0; i < nm; ++i)
    {
      ix = i * 3;
      iy = ix + 1;
      iz = ix + 2;

      /* periodic condition in horizontal */
      for (;x [ix] < 0.0;) x [ix] += lx;
      for (;x [ix] >= lx;) x [ix] -= lx;

      /* periodic condition in horizontal */
      if (sys->shear_mode == 1)
	{
	  for (;x [iy] < 0.0;)
	    {
	      x [iy] += ly;
	      x [ix] += sys->shear_shift;
	    }
	  for (;x [iy] >= ly;)
	    {
	      x [iy] -= ly;
	      x [ix] -= sys->shear_shift;
	    }
	}
      else
	{
	  for (;x [iy] < 0.0;) x [iy] += ly;
	  for (;x [iy] >= ly;) x [iy] -= ly;
	}

      /* periodic condition in vertical */
      if (sys->shear_mode == 2)
	{
	  for (;x [iz] < 0.0;)
	    {
	      x [iz] += lz;
	      x [ix] += sys->shear_shift;
	    }
	  for (;x [iz] >= lz;)
	    {
	      x [iz] -= lz;
	      x [ix] -= sys->shear_shift;
	    }
	}
      else
	{
	  for (;x [iz] < 0.0;) x [iz] += lz;
	  for (;x [iz] >= lz;) x [iz] -= lz;
	}
    }
}

/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  angle [np * 3] : angle to check; the range is [0, 2\pi)
 */
void
check_angle (struct stokes *sys, double *angle)
{
  int nm;
  int i;

  nm = sys->nm;
  for (i = 0; i < nm * 3; ++i)
    {
      for (;angle [i] < 0.0;)         angle [i] += M_PI * 2.0;
      for (;angle [i] >= M_PI * 2.0;) angle [i] -= M_PI * 2.0;
    }
}


/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  x [np * 2] : position to check; the range is [0, l[xy])
 */
void
check_periodic_2d (struct stokes *sys,
		   double *x)
{
  int i;
  int ix, iy;

  int nm;
  double lx, ly;

  nm = sys->nm;
  lx = sys->lx;
  ly = sys->ly;

  for (i = 0; i < nm; ++i)
    {
      ix = i * 2;
      iy = ix + 1;

      /* periodic condition in horizontal */
      for (;x [ix] < 0.0;) x [ix] += lx;
      for (;x [ix] >= lx;) x [ix] -= lx;

      /* periodic condition in vertical */
      if (sys->shear_mode == 1)
	{
	  for (;x [iy] < 0.0;)
	    {
	      x [iy] += ly;
	      x [ix] += sys->shear_shift;
	    }
	  for (;x [iy] >= ly;)
	    {
	      x [iy] -= ly;
	      x [ix] -= sys->shear_shift;
	    }
	}
      else
	{
	  for (;x [iy] < 0.0;) x [iy] += ly;
	  for (;x [iy] >= ly;) x [iy] -= ly;
	}
    }
}

/*
 * INTPUT
 *  sys : system parameters
 * OUTPUT
 *  angle [np] : angle to check; the range is [0, 2\pi)
 */
void
check_angle_2d (struct stokes *sys, double *angle)
{
  int nm;
  int i;

  nm = sys->nm;
  for (i = 0; i < nm; ++i)
    {
      for (;angle [i] < 0.0;)         angle [i] += M_PI * 2.0;
      for (;angle [i] >= M_PI * 2.0;) angle [i] -= M_PI * 2.0;
    }
}
