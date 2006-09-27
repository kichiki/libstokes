/* structure for system parameters of stokes library.
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes.c,v 2.2 2006/09/27 00:02:03 ichiki Exp $
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
#include <stdlib.h>
#include <math.h> /* log() */

#include "stokes.h"

/* all elements are zero-cleared
 */
struct stokes *
stokes_init (void)
{
  struct stokes * sys = NULL;

  sys = (struct stokes *) malloc (sizeof (struct stokes));

  sys->np = 0;
  sys->nm = 0;
  sys->pos = NULL;

  /* # of cell in real space */
  sys->pcellx = 0;
  sys->pcelly = 0;
  sys->pcellz = 0;

  /* # of cell in reciprocal space */
  sys->kmaxx = 0;
  sys->kmaxy = 0;
  sys->kmaxz = 0;

  sys->zeta = 0.0;
  sys->zeta2 = 0.0;
  sys->zaspi = 0.0;
  sys->za2 = 0.0;

  sys->pi2 = 0.0;
  sys->pivol = 0.0;

  sys->lx = 0.0;
  sys->ly = 0.0;
  sys->lz = 0.0;

  int i;
  for (i = 0; i < 27; i ++)
    {
      sys->llx [i] = 0.0;
      sys->lly [i] = 0.0;
      sys->llz [i] = 0.0;
    }

  /* for zeta program */
  sys->cpu1 = 0.0;
  sys->cpu2 = 0.0;
  sys->cpu3 = 0.0;

  return (sys);
}

void
stokes_free (struct stokes * sys)
{
  if (sys != NULL)
    {
      if (sys->pos != NULL) free (sys->pos);
      free (sys);
    }
}


void
stokes_set_np (struct stokes * sys,
	       int np, int nm)
{
  sys->np = np;
  sys->nm = nm;

  if (sys->pos != NULL) free (sys->pos);
  sys->pos = (double *) malloc (sizeof (double) * np * 3);
}


void
stokes_set_ll (struct stokes * sys,
	       double lx, double ly, double lz)
{
  sys->pi2 = M_PI * 2.0;
  sys->pivol = M_PI / lx / ly / lz;

  sys->lx = lx;
  sys->ly = ly;
  sys->lz = lz;

  /* 0 is always the primary cell */
  sys->llx [ 0] = 0.0; sys->lly [ 0] = 0.0; sys->llz [ 0] = 0.0;
  sys->llx [ 1] = -lx; sys->lly [ 1] = 0.0; sys->llz [ 1] = 0.0;
  sys->llx [ 2] = 0.0; sys->lly [ 2] = -ly; sys->llz [ 2] = 0.0;
  sys->llx [ 3] = +lx; sys->lly [ 3] = 0.0; sys->llz [ 3] = 0.0;
  sys->llx [ 4] = 0.0; sys->lly [ 4] = +ly; sys->llz [ 4] = 0.0;
  sys->llx [ 5] = -lx; sys->lly [ 5] = -ly; sys->llz [ 5] = 0.0;
  sys->llx [ 6] = -lx; sys->lly [ 6] = +ly; sys->llz [ 6] = 0.0;
  sys->llx [ 7] = +lx; sys->lly [ 7] = -ly; sys->llz [ 7] = 0.0;
  sys->llx [ 8] = +lx; sys->lly [ 8] = +ly; sys->llz [ 8] = 0.0;
  /* up to 8, the monolayer (xy plain)
   * note: y is the vertical direction for 2D monolayer case! */

  sys->llx [ 9] = 0.0; sys->lly [ 9] = 0.0; sys->llz [ 9] = - lz;
  sys->llx [10] = -lx; sys->lly [10] = 0.0; sys->llz [10] = - lz;
  sys->llx [11] = 0.0; sys->lly [11] = -ly; sys->llz [11] = - lz;
  sys->llx [12] = +lx; sys->lly [12] = 0.0; sys->llz [12] = - lz;
  sys->llx [13] = 0.0; sys->lly [13] = +ly; sys->llz [13] = - lz;
  sys->llx [14] = -lx; sys->lly [14] = -ly; sys->llz [14] = - lz;
  sys->llx [15] = -lx; sys->lly [15] = +ly; sys->llz [15] = - lz;
  sys->llx [16] = +lx; sys->lly [16] = -ly; sys->llz [16] = - lz;
  sys->llx [17] = +lx; sys->lly [17] = +ly; sys->llz [17] = - lz;

  sys->llx [18] = 0.0; sys->lly [18] = 0.0; sys->llz [18] = + lz;
  sys->llx [19] = -lx; sys->lly [19] = 0.0; sys->llz [19] = + lz;
  sys->llx [20] = 0.0; sys->lly [20] = -ly; sys->llz [20] = + lz;
  sys->llx [21] = +lx; sys->lly [21] = 0.0; sys->llz [21] = + lz;
  sys->llx [22] = 0.0; sys->lly [22] = +ly; sys->llz [22] = + lz;
  sys->llx [23] = -lx; sys->lly [23] = -ly; sys->llz [23] = + lz;
  sys->llx [24] = -lx; sys->lly [24] = +ly; sys->llz [24] = + lz;
  sys->llx [25] = +lx; sys->lly [25] = -ly; sys->llz [25] = + lz;
  sys->llx [26] = +lx; sys->lly [26] = +ly; sys->llz [26] = + lz;
}

void
stokes_set_zeta (struct stokes * sys,
		 double zeta, double cutlim)
{
  sys->zeta  = zeta;
  sys->zeta2 = zeta * zeta;
  sys->za2   = sys->zeta2;
  sys->zaspi = zeta / sqrt (M_PI);

  /* define # of cells  */
  /* in real space */
  sys->pcellx = (int)(sqrt (-log(cutlim)) / zeta / sys->lx) + 1;
  sys->pcelly = (int)(sqrt (-log(cutlim)) / zeta / sys->ly) + 1;
  sys->pcellz = (int)(sqrt (-log(cutlim)) / zeta / sys->lz) + 1;
  /* in reciprocal space */
  sys->kmaxx = (int)(sqrt (-log(cutlim)) * zeta * sys->lx / M_PI);
  sys->kmaxy = (int)(sqrt (-log(cutlim)) * zeta * sys->ly / M_PI);
  sys->kmaxz = (int)(sqrt (-log(cutlim)) * zeta * sys->lz / M_PI);
}


double
zeta_by_tratio (struct stokes * sys,
		double tratio)
{
  double zeta;

  zeta = pow (tratio, 1.0 / 6.0)
    * sqrt (M_PI)
    / pow (sys->lx * sys->ly * sys->lz, 1.0 / 3.0);

  return (zeta);
}
