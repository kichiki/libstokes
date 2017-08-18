/* header file for grid.c --
 * RYUON_grid routines
 * Copyright (C) 2008,2017 Kengo Ichiki <kengoichiki@gmail.com>
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
#ifndef	_GRID_H_
#define	_GRID_H_


#include <libstokes-core.h> // struct stokes


struct RYUON_grid {
  // the origin for the grid generation
  double x0;
  double y0;
  double z0;

  // number of grids
  int nx, ny, nz;
  int n; // := nx * ny * nz
  // sizes of the cell
  double lx, ly, lz;

  // particle list for each cell
  int *np;  // np[nx*ny*nz]       : number of particles in each cell
  int **ip; // ip[nx*ny*nz][np[]] : particle indices for each cell

  // nearest neighbor list
  int nofp;  // number of particles
  int *nnn;  // nnn[nofp]        : number of NN particles for each particle
  int **nnp; // nnp[nofp][nnn[]] : NN particle indices for each particle
};


/* initialize RYUON_grid
 * by bounding box (x0,x1, y0,y1, z0,z1)
 * and cell number (nx, ny, nz)
 * NOTE: particles are not assigned yet.
 * use GRID_add_particle() or GRID_assign_particles().
 * alternatively, use GRID_init_all_by_cutoff() in the first place.
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
	   int nx, int ny, int nz);

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
GRID_init_by_l (double x0, double x1,
		double y0, double y1,
		double z0, double z1,
		double l);

void
GRID_free (struct RYUON_grid *g);


/* clear all particles in RYUON_grid
 */
void
GRID_clear (struct RYUON_grid *g);

int
GRID_ixyz_to_in (struct RYUON_grid *g,
		 int ix, int iy, int iz);

void
GRID_in_to_ixyz (struct RYUON_grid *g,
		 int in,
		 int *ix, int *iy, int *iz);


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
		   int ip);

/* set particles to RYUON_grid.
 * note that g is first clear and assigned
 * INPUT
 *  np : number of particles
 *  x[np*3] : configuration
 */
void
GRID_assign_particles (struct RYUON_grid *g,
		       int np,
		       const double *x);


void
GRID_clear_NN (struct RYUON_grid *g);


void
GRID_set_NN (struct RYUON_grid *g,
	     int np, const double *pos,
	     double r_cutoff2);

void
GRID_set_NN_periodic (struct RYUON_grid *g,
		      int np, const double *pos,
		      double lx, double ly, double lz,
		      double r_cutoff2);


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
			 double r_cutoff);


#endif /* !_GRID_H_ */
