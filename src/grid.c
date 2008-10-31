/* RYUON_grid routines
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: grid.c,v 1.1 2008/10/31 05:38:43 kichiki Exp $
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
#include <stdio.h>  // fprintf()
#include <stdlib.h> // malloc()
#include <math.h>   // floor()
#include "memory-check.h" // CHECK_MALLOC

#include "grid.h"


/* initialize RYUON_grid
 * by bounding box (x0,x1, y0,y1, z0,z1)
 * and cell number (nx, ny, nz)
 * NOTE: particles are not assigned yet.
 * use GRID_add_particle() or GRID_assign_particles().
 * alternatively, use GRID_init_all_by_l() in the first place.
 * INPUT
 *  x0,x1, y0,y1, z0,z1 : bounding box
 *  nx, ny, nz : cell number
 * OUTPUT
 *  (returned value) : struct RYUON_grid
 */
struct RYUON_grid *
GRID_init (double x0, double x1,
	   double y0, double y1,
	   double z0, double z1,
	   int nx, int ny, int nz)
{
  struct RYUON_grid *g
    = (struct RYUON_grid *)malloc (sizeof (struct RYUON_grid));
  CHECK_MALLOC (g, "GRID_init");

  g->x0 = x0;
  g->x1 = x1;
  g->y0 = y0;
  g->y1 = y1;
  g->z0 = z0;
  g->z1 = z1;

  g->nx = nx;
  g->ny = ny;
  g->nz = nz;
  g->n  = nx * ny * nz;

  //double safe_factor = 1.0 + 1.0e-12;
  double safe_factor = 1.0 + 1.0e-15;
  g->lx = (x1 - x0) / (double)nx * safe_factor;
  g->ly = (y1 - y0) / (double)ny * safe_factor;
  g->lz = (z1 - z0) / (double)nz * safe_factor;

  g->np = (int *)calloc (sizeof (int), g->n);
  CHECK_MALLOC (g->np, "GRID_init");

  g->ip = (int **)malloc (sizeof (int *) * g->n);
  int i;
  for (i = 0; i < g->n; i ++)
    {
      g->ip[i] = NULL;
    }

  return (g);
}

/* initialize RYUON_grid
 * by bounding box (x0,x1, y0,y1, z0,z1)
 * and cell size l
 * NOTE: particles are not assigned yet.
 * use GRID_add_particle() or GRID_assign_particles().
 * alternatively, use GRID_init_all_by_l() in the first place.
 * INPUT
 *  x0,x1, y0,y1, z0,z1 : bounding box
 *  l : cell size
 * OUTPUT
 *  (returned value) : struct RYUON_grid
 */
struct RYUON_grid *
GRID_init_by_l (double x0, double x1,
		double y0, double y1,
		double z0, double z1,
		double l)
{
  struct RYUON_grid *g
    = (struct RYUON_grid *)malloc (sizeof (struct RYUON_grid));
  CHECK_MALLOC (g, "GRID_init_by_l");

  g->x0 = x0;
  g->x1 = x1;
  g->y0 = y0;
  g->y1 = y1;
  g->z0 = z0;
  g->z1 = z1;

  g->lx = l;
  g->ly = l;
  g->lz = l;

  double safe_factor = 1.0 + 1.0e-15;
  g->nx = (int)(safe_factor * (x1 - x0) / l) + 1;
  g->ny = (int)(safe_factor * (y1 - y0) / l) + 1;
  g->nz = (int)(safe_factor * (z1 - z0) / l) + 1;
  g->n  = g->nx * g->ny * g->nz;

  g->np = (int *)calloc (sizeof (int), g->n);
  CHECK_MALLOC (g->np, "GRID_init_by_l");

  g->ip = (int **)malloc (sizeof (int *) * g->n);
  int i;
  for (i = 0; i < g->n; i ++)
    {
      g->ip[i] = NULL;
    }

  return (g);
}

void
GRID_free (struct RYUON_grid *g)
{
  if (g == NULL) return;
  if (g->np != NULL) free (g->np);
  if (g->ip != NULL)
    {
      int i;
      for (i = 0; i < g->n; i ++)
	{
	  if (g->ip[i] != NULL) free (g->ip[i]);
	}
      free (g->ip);
    }
  free (g);
}


/* clear all particles in RYUON_grid
 */
void
GRID_clear (struct RYUON_grid *g)
{
  int i;
  for (i = 0; i < g->n; i ++)
    {
      if (g->np[i] > 0)
	{
	  g->np[i] = 0;
	  free (g->ip[i]);
	  g->ip[i] = NULL;
	}
    }
}

int
GRID_ixyz_to_in (struct RYUON_grid *g,
		 int ix, int iy, int iz)
{
  int in = ix + g->nx * (iy + g->ny * iz);
  return (in);
}

void
GRID_in_to_ixyz (struct RYUON_grid *g,
		 int in,
		 int *ix, int *iy, int *iz)
{
  *ix = in % g->nx;
  int iyy = in / g->nx;
  *iy = iyy % g->ny;
  *iz = iyy / g->ny;
}


/* assign one particle into RYUON_grid
 * INPUT
 *  g : struct RYUON_grid
 *  x[3] : position of the particle
 * OUTPUT
 *  (returned value) : 0 (false) == the particle is out of range
 *                     1 (true)  == the particle is assigned successfully
 *  g :
 */
int
GRID_add_particle (struct RYUON_grid *g,
		   const double *x,
		   int ip)
{
  int ix = (int)floor ((x[0] - g->x0) / g->lx);
  int iy = (int)floor ((x[1] - g->y0) / g->ly);
  int iz = (int)floor ((x[2] - g->z0) / g->lz);
  if (ix < 0 || ix >= g->nx ||
      iy < 0 || iy >= g->ny ||
      iz < 0 || iz >= g->nz)
    {
      return (0); // false
    }

  int in = GRID_ixyz_to_in (g, ix, iy, iz);

  g->np[in] ++;
  g->ip[in] = (int *)realloc (g->ip[in], sizeof (int) * g->np[in]);
  CHECK_MALLOC (g->ip[in], "GRID_add_particle");

  g->ip[in][g->np[in] - 1] = ip;

  return (1); // true
}


/* set particles to RYUON_grid.
 * note that g is first clear and assigned
 * INPUT
 *  np : number of particles
 *  x[np*3] : configuration
 */
void
GRID_assign_particles (struct RYUON_grid *g,
		       int np,
		       const double *x)
{
  // first, clear all particles in RYUON_grid
  GRID_clear (g);

  int i;
  for (i = 0; i < np; i ++)
    {
      int ix = i * 3;
      if (GRID_add_particle (g, x + ix, i) == 0) // false
	{
	  fprintf (stderr, "GRID_assign_particles():"
		   "fail to assign particle %d\n", i);
	}
    }
}


/* initialize RYUON_grid
 * by struct stokes *sys and cell size l
 * NOTE: particles are assigned too.
 * INPUT
 *  sys : struct stokes
 *  l : cell size
 * OUTPUT
 *  (returned value) : struct RYUON_grid
 */
struct RYUON_grid *
GRID_init_all_by_l (struct stokes *sys,
		    double l)
{
  struct RYUON_grid *g
    = (struct RYUON_grid *)malloc (sizeof (struct RYUON_grid));
  CHECK_MALLOC (g, "GRID_init_by_l");

  int i;
  if (sys->periodic == 0)
    {
      // non-periodic system
      g->x0 = g->x1 = sys->pos[0];
      g->y0 = g->y1 = sys->pos[1];
      g->z0 = g->z1 = sys->pos[2];
      for (i = 1; i < sys->np; i ++)
	{
	  int ix = i * 3;
	  if (g->x0 > sys->pos[ix  ]) g->x0 = sys->pos[ix  ];
	  if (g->x1 < sys->pos[ix  ]) g->x1 = sys->pos[ix  ];

	  if (g->y0 > sys->pos[ix+1]) g->y0 = sys->pos[ix+1];
	  if (g->y1 < sys->pos[ix+1]) g->y1 = sys->pos[ix+1];

	  if (g->z0 > sys->pos[ix+2]) g->z0 = sys->pos[ix+2];
	  if (g->z1 < sys->pos[ix+2]) g->z1 = sys->pos[ix+2];
	}
    }
  else
    {
      // periodic system
      g->x0 = 0.0;
      g->y0 = 0.0;
      g->z0 = 0.0;
      g->x1 = sys->lx;
      g->y1 = sys->ly;
      g->z1 = sys->lz;
    }

  // these are the CELL size (not the whole system size as sys->lx,ly,lz)
  g->lx = l;
  g->ly = l;
  g->lz = l;

  double safe_factor = 1.0 + 1.0e-15;
  g->nx = (int)(safe_factor * (g->x1 - g->x0) / l) + 1;
  g->ny = (int)(safe_factor * (g->y1 - g->y0) / l) + 1;
  g->nz = (int)(safe_factor * (g->z1 - g->z0) / l) + 1;
  g->n  = g->nx * g->ny * g->nz;

  g->np = (int *)calloc (sizeof (int), g->n);
  CHECK_MALLOC (g->np, "GRID_init_by_l");

  g->ip = (int **)malloc (sizeof (int *) * g->n);
  for (i = 0; i < g->n; i ++)
    {
      g->ip[i] = NULL;
    }


  GRID_assign_particles (g, sys->np, sys->pos);


  return (g);
}

