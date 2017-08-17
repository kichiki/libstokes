/* RYUON_grid routines
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: grid.c,v 1.4 2008/11/02 06:14:17 kichiki Exp $
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
  g->y0 = y0;
  g->z0 = z0;

  g->nx = nx;
  g->ny = ny;
  g->nz = nz;
  g->n  = nx * ny * nz;

  //double safe_factor = 1.0 + 1.0e-12;
  double safe_factor = 1.0 + 1.0e-15;
  g->lx = (x1 - x0) / (double)nx * safe_factor;
  g->ly = (y1 - y0) / (double)ny * safe_factor;
  g->lz = (z1 - z0) / (double)nz * safe_factor;


  // particle list for each cell
  g->np = (int *)calloc (sizeof (int), g->n);
  CHECK_MALLOC (g->np, "GRID_init");

  g->ip = (int **)malloc (sizeof (int *) * g->n);
  int i;
  for (i = 0; i < g->n; i ++)
    {
      g->ip[i] = NULL;
    }


  // nearest neighbor list
  g->nofp = 0;
  g->nnn = NULL;
  g->nnp = NULL;


  return (g);
}

/* initialize RYUON_grid
 * by bounding box (x0,x1, y0,y1, z0,z1)
 * and cell size l
 * NOTE: particles are not assigned yet.
 * use GRID_add_particle() or GRID_assign_particles().
 * alternatively, use GRID_init_all_by_cutoff() in the first place.
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
  g->y0 = y0;
  g->z0 = z0;

  g->lx = l;
  g->ly = l;
  g->lz = l;

  double safe_factor = 1.0 + 1.0e-15;
  g->nx = (int)(safe_factor * (x1 - x0) / l) + 1;
  g->ny = (int)(safe_factor * (y1 - y0) / l) + 1;
  g->nz = (int)(safe_factor * (z1 - z0) / l) + 1;
  g->n  = g->nx * g->ny * g->nz;


  // particle list for each cell
  g->np = (int *)calloc (sizeof (int), g->n);
  CHECK_MALLOC (g->np, "GRID_init_by_l");

  g->ip = (int **)malloc (sizeof (int *) * g->n);
  int i;
  for (i = 0; i < g->n; i ++)
    {
      g->ip[i] = NULL;
    }


  // nearest neighbor list
  g->nofp = 0;
  g->nnn = NULL;
  g->nnp = NULL;


  return (g);
}

void
GRID_free (struct RYUON_grid *g)
{
  if (g == NULL) return;


  // particle list for each cell
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


  // nearest neighbor list
  GRID_clear_NN (g);


  free (g);
}


/* clear all particles in RYUON_grid
 */
void
GRID_clear (struct RYUON_grid *g)
{
  // particle list for each cell
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

  // nearest neighbor list
  GRID_clear_NN (g);
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
      fprintf (stderr, "# (%f, %f, %f) => (%d %d %d)\n",
	       x[0], x[1], x[2], ix, iy, iz);
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


void
GRID_clear_NN (struct RYUON_grid *g)
{
  if (g->nnn != NULL)
    {
      int i;
      for (i = 0; i < g->nofp; i ++)
	{
	  if (g->nnp[i] != NULL) free (g->nnp[i]);
	}
      free (g->nnp);
      free (g->nnn);
      g->nnp = NULL;
      g->nnn = NULL;
    }
  g->nofp = 0;
}
	     
void
GRID_set_NN (struct RYUON_grid *g,
	     int np, const double *pos,
	     double r_cutoff2)
{
  // first clearn NN list
  GRID_clear_NN (g);


  g->nofp = np;
  g->nnn = (int *)calloc (sizeof (int), np);
  g->nnp = (int **)malloc (sizeof (int *) * np);
  int i;
  for (i = 0; i < np; i ++)
    {
      g->nnp[i] = NULL;
    }

  int ic;
  for (ic = 0; ic < g->n; ic ++)
    {
      int ipp;
      for (ipp = 0; ipp < g->np[ic]; ipp++)
	{
	  int ip = g->ip[ic][ipp];

	  // loop for the nearest neighbor cells
	  int ix, iy, iz;
	  GRID_in_to_ixyz (g, ic, &ix, &iy, &iz);
	  int jx;
	  for (jx = ix - 1; jx <= ix + 1; jx++)
	    {
	      if (jx < 0 || jx >= g->nx) continue; // out of range
	      int jy;
	      for (jy = iy - 1; jy <= iy + 1; jy++)
		{
		  if (jy < 0 || jy >= g->ny) continue; // out of range
		  int jz;
		  for (jz = iz - 1; jz <= iz + 1; jz++)
		    {
		      if (jz < 0 || jz >= g->nz) continue; // out of range

		      int jc = GRID_ixyz_to_in (g, jx, jy, jz);
		      int jpp;
		      for (jpp = 0; jpp < g->np[jc]; jpp ++)
			{
			  int jp = g->ip[jc][jpp];
			  // skip the self part and pair of (ip > jp)
			  if (ip >= jp) continue;

			  // j is in the nearest neighbor cell
			  int i3 = ip * 3;
			  int j3 = jp * 3;
			  double x = pos [i3  ] - pos [j3  ];
			  double y = pos [i3+1] - pos [j3+1];
			  double z = pos [i3+2] - pos [j3+2];

			  double r2 = x*x + y*y + z*z;
			  if (r2 < r_cutoff2)
			    {
			      g->nnn[ip] ++;
			      g->nnp[ip] = (int *)realloc
				(g->nnp[ip], sizeof (int) * g->nnn[ip]);
			      CHECK_MALLOC
				(g->nnp[ip], "GRID_set_NN");
			      g->nnp[ip][g->nnn[ip] - 1] = jp;
			    }
			  // enf of if (r2 < r_cutoff2)
			}
		      // end of for (jpp = 0; jpp < g->np[jc]; jpp ++)
		    }
		  // end of for (jz)
		}
	      // end of for (jy)
	    }
	  // end of for (jz)
	  // end of nearest neighbor loop
	}
      // end of for (ipp)
    }
  // end of for (ic)
}
	     
void
GRID_set_NN_periodic (struct RYUON_grid *g,
		      int np, const double *pos,
		      double lx, double ly, double lz,
		      double r_cutoff2)
{
  /**
   * some remarks (2008/10/31)
   * this is not completed yet
   * -- shift vector (lx,ly,lz) is also necessary to form loop.
   */

  // first clearn NN list
  GRID_clear_NN (g);


  g->nofp = np;
  g->nnn = (int *)calloc (sizeof (int), np);
  g->nnp = (int **)malloc (sizeof (int *) * np);
  int i;
  for (i = 0; i < np; i ++)
    {
      g->nnp[i] = NULL;
    }

  int ic;
  for (ic = 0; ic < g->n; ic ++)
    {
      int ipp;
      for (ipp = 0; ipp < g->np[ic]; ipp++)
	{
	  int ip = g->ip[ic][ipp];

	  // loop for the nearest neighbor cells
	  int ix, iy, iz;
	  double llx = 0.0;
	  double lly = 0.0;
	  double llz = 0.0;
	  GRID_in_to_ixyz (g, ic, &ix, &iy, &iz);
	  int jx;
	  for (jx = ix - 1; jx <= ix + 1; jx++)
	    {
	      if (jx < 0)
		{
		  jx += g->nx;
		  llx -= lx;
		}
	      if (jx >= g->nx)
		{
		  jx -= g->nx;
		  llx += lx;
		}
	      int jy;
	      for (jy = iy - 1; jy <= iy + 1; jy++)
		{
		  if (jy < 0)
		    {
		      jy += g->ny;
		      lly -= ly;
		    }
		  if (jy >= g->ny)
		    {
		      jy -= g->ny;
		      lly += ly;
		    }
		  int jz;
		  for (jz = iz - 1; jz <= iz + 1; jz++)
		    {
		      if (jz < 0)
			{
			  jz += g->nz;
			  llz -= lz;
			}
		      if (jz >= g->nz)
			{
			  jz -= g->nz;
			  llz += lz;
			}

		      int jc = GRID_ixyz_to_in (g, jx, jy, jz);
		      int jpp;
		      for (jpp = 0; jpp < g->np[jc]; jpp ++)
			{
			  int jp = g->ip[jc][jpp];
			  // skip the self part and pair of (ip > jp)
			  if (ip >= jp) continue;

			  // j is in the nearest neighbor cell
			  int i3 = ip * 3;
			  int j3 = jp * 3;
			  double x = pos [i3  ] - (pos [j3  ] + llx);
			  double y = pos [i3+1] - (pos [j3+1] + lly);
			  double z = pos [i3+2] - (pos [j3+2] + llz);

			  double r2 = x*x + y*y + z*z;
			  if (r2 < r_cutoff2)
			    {
			      g->nnn[ip] ++;
			      g->nnp[ip] = (int *)realloc
				(g->nnp[ip], sizeof (int) * g->nnn[ip]);
			      CHECK_MALLOC
				(g->nnp[ip], "GRID_set_NN_periodic");
			      g->nnp[ip][g->nnn[ip] - 1] = jp;
			    }
			  // enf of if (r2 < r_cutoff2)
			}
		      // end of for (jpp = 0; jpp < g->np[jc]; jpp ++)
		    }
		  // end of for (jz)
		}
	      // end of for (jy)
	    }
	  // end of for (jz)
	  // end of nearest neighbor loop
	}
      // end of for (ipp)
    }
  // end of for (ic)
}
	     

/* initialize RYUON_grid
 *   by struct stokes *sys
 *   and the cut-off distance r_cutoff
 * NOTE: both particle list for each cell
 *   as well as nearest neighbor list are formed.
 *   the nearest neighbor list contains pairs of particles
 *   (i,j) of (i < j) only, that is, no self part.
 * INPUT
 *  sys      : struct stokes
 *  r_cutoff : cut-off length
 * OUTPUT
 *  (returned value) : struct RYUON_grid
 */
struct RYUON_grid *
GRID_init_all_by_cutoff (struct stokes *sys,
			 double r_cutoff)
{
  struct RYUON_grid *g
    = (struct RYUON_grid *)malloc (sizeof (struct RYUON_grid));
  CHECK_MALLOC (g, "GRID_init_all_by_cutoff");

  double x1, y1, z1; // max limits
  int i;
  if (sys->periodic == 0)
    {
      // non-periodic system
      g->x0 = x1 = sys->pos[0];
      g->y0 = y1 = sys->pos[1];
      g->z0 = z1 = sys->pos[2];
      for (i = 1; i < sys->np; i ++)
	{
	  int ix = i * 3;
	  if (g->x0 > sys->pos[ix  ]) g->x0 = sys->pos[ix  ];
	  if (x1    < sys->pos[ix  ]) x1    = sys->pos[ix  ];

	  if (g->y0 > sys->pos[ix+1]) g->y0 = sys->pos[ix+1];
	  if (y1    < sys->pos[ix+1]) y1    = sys->pos[ix+1];

	  if (g->z0 > sys->pos[ix+2]) g->z0 = sys->pos[ix+2];
	  if (z1    < sys->pos[ix+2]) z1    = sys->pos[ix+2];
	}
    }
  else
    {
      // periodic system
      g->x0 = 0.0;
      g->y0 = 0.0;
      g->z0 = 0.0;

      x1 = sys->lx;
      y1 = sys->ly;
      z1 = sys->lz;
    }

  // these are the CELL size (not the whole system size as sys->lx,ly,lz)
  g->lx = r_cutoff;
  g->ly = r_cutoff;
  g->lz = r_cutoff;

  double safe_factor = 1.0 + 1.0e-15;
  g->nx = (int)(safe_factor * (x1 - g->x0) / r_cutoff) + 1;
  g->ny = (int)(safe_factor * (y1 - g->y0) / r_cutoff) + 1;
  g->nz = (int)(safe_factor * (z1 - g->z0) / r_cutoff) + 1;
  g->n  = g->nx * g->ny * g->nz;

  // shift half step for the extra 1 cell
  if (sys->periodic != 0)
    {
      g->x0 -= 0.5 * r_cutoff;
      g->y0 -= 0.5 * r_cutoff;
      g->z0 -= 0.5 * r_cutoff;
    }


  // particle list for each cell
  g->np = (int *)calloc (sizeof (int), g->n);
  CHECK_MALLOC (g->np, "GRID_init_all_by_cutoff");

  g->ip = (int **)malloc (sizeof (int *) * g->n);
  for (i = 0; i < g->n; i ++)
    {
      g->ip[i] = NULL;
    }

  // nearest neighbor list
  g->nofp = 0;
  g->nnn = NULL;
  g->nnp = NULL;


  GRID_assign_particles (g, sys->np, sys->pos);

  if (sys->periodic == 0)
    {
      GRID_set_NN (g, sys->np, sys->pos, r_cutoff * r_cutoff);
    }
  else
    {
      fprintf (stderr, "# GRID_init_all_by_cutoff() :"
	       "implementation is not finished for the periodic system.\n");
      GRID_set_NN_periodic (g,
			    sys->np, sys->pos,
			    sys->lx, sys->ly, sys->lz,
			    r_cutoff * r_cutoff);
    }


  return (g);
}

