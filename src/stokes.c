/* utility routines for Ewald-summation code in 3D
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: stokes.c,v 1.5 2001/02/02 08:38:41 ichiki Exp $
 */
#include <math.h>

#include "ewald-3.h"

/** global variables **/
int pcellx, pcelly, pcellz; /* # of cell in real space */
int kmaxx, kmaxy, kmaxz; /* # of cell in reciprocal space */
double zeta, zeta2, zaspi, za2;
double pi2;
double pivol;
double lx, ly, lz; /* cell size */
double llx [27], lly [27], llz [27]; /* for regist and lub */



/* initialize parameters used in Ewald-summation code in 3D
 * INPUT
 *  lx, ly, lz : geometry of the primary cell
 *  tratio : a parameter by zeta-code to determine pcell? and kmax?
 *  cutlim : cut-off limit to determine pcell? and kmax?
 * OUTPUT
 *  (global) pcellx, pcelly, pcellz;
 *  (global) kmaxx, kmaxy, kmaxz;
 *  (global) zeta, zeta2, zaspi, za2;
 *  (global) pi2;
 *  (global) pivol;
 *  (global) lx, ly, lz;
 */
void
init_ewald_3d (double lx, double ly, double lz,
	       double tratio, double cutlim)
{
  extern int pcellx, pcelly, pcellz;
  extern int kmaxx, kmaxy, kmaxz;

  extern double zeta, zeta2, zaspi, za2;
  extern double pi2;
  extern double pivol;
  //extern double lx, ly, lz; /* cell size */
  extern double llx [27], lly [27], llz [27];


  /* initialization */
  pi2 = M_PI * 2.0;
  pivol = M_PI / lx / ly / lz;

  llx [ 0] = 0.0; lly [ 0] = 0.0; llz [ 0] = 0.0;
  llx [ 1] = -lx; lly [ 1] = 0.0; llz [ 1] = 0.0;
  llx [ 2] = 0.0; lly [ 2] = -ly; llz [ 2] = 0.0;
  llx [ 3] = +lx; lly [ 3] = 0.0; llz [ 3] = 0.0;
  llx [ 4] = 0.0; lly [ 4] = +ly; llz [ 4] = 0.0;
  llx [ 5] = -lx; lly [ 5] = -ly; llz [ 5] = 0.0;
  llx [ 6] = -lx; lly [ 6] = +ly; llz [ 6] = 0.0;
  llx [ 7] = +lx; lly [ 7] = -ly; llz [ 7] = 0.0;
  llx [ 8] = +lx; lly [ 8] = +ly; llz [ 8] = 0.0;

  llx [ 9] = 0.0; lly [ 9] = 0.0; llz [ 9] = - lz;
  llx [10] = -lx; lly [10] = 0.0; llz [10] = - lz;
  llx [11] = 0.0; lly [11] = -ly; llz [11] = - lz;
  llx [12] = +lx; lly [12] = 0.0; llz [12] = - lz;
  llx [13] = 0.0; lly [13] = +ly; llz [13] = - lz;
  llx [14] = -lx; lly [14] = -ly; llz [14] = - lz;
  llx [15] = -lx; lly [15] = +ly; llz [15] = - lz;
  llx [16] = +lx; lly [16] = -ly; llz [16] = - lz;
  llx [17] = +lx; lly [17] = +ly; llz [17] = - lz;

  llx [18] = 0.0; lly [18] = 0.0; llz [18] = + lz;
  llx [19] = -lx; lly [19] = 0.0; llz [19] = + lz;
  llx [20] = 0.0; lly [20] = -ly; llz [20] = + lz;
  llx [21] = +lx; lly [21] = 0.0; llz [21] = + lz;
  llx [22] = 0.0; lly [22] = +ly; llz [22] = + lz;
  llx [23] = -lx; lly [23] = -ly; llz [23] = + lz;
  llx [24] = -lx; lly [24] = +ly; llz [24] = + lz;
  llx [25] = +lx; lly [25] = -ly; llz [25] = + lz;
  llx [26] = +lx; lly [26] = +ly; llz [26] = + lz;

  /* define zeta */
  /*zeta = 0.01;*/
  zeta = pow (tratio, 1.0 / 6.0)
    * sqrt (M_PI)
    / pow (lx * ly * lz, 1.0 / 3.0);
  zeta2 = zeta * zeta;
  za2 = zeta2;
  zaspi = zeta / sqrt (M_PI);

  /* define # of cells  */
  /* in real space */
  pcellx = (int) (sqrt (- log (cutlim)) / zeta / lx) + 1;
  pcelly = (int) (sqrt (- log (cutlim)) / zeta / ly) + 1;
  pcellz = (int) (sqrt (- log (cutlim)) / zeta / lz) + 1;
  /* in reciprocal space */
  kmaxx = (int) (sqrt (- log (cutlim)) * zeta * lx / M_PI);
  kmaxy = (int) (sqrt (- log (cutlim)) * zeta * ly / M_PI);
  kmaxz = (int) (sqrt (- log (cutlim)) * zeta * lz / M_PI);
}


/* initialize configuration and the primary cell of simple cubic lattic
 * INPUT
 *  phi        : volume fraction of particles
 *  nx, ny, nz : # particles in each direction
 * OUTPUT
 *  pos [(nx * ny * nz) * 3] : positions of particles
 *  lx, ly, lz : geometry of the primary cell
 */
void
init_config_SC (double phi, int nx, int ny, int nz,
		double *pos, double *lx, double *ly, double *lz)
{
  int j;
  int ix, iy, iz;


  (*lx) = (*ly) = (*lz) = pow (4.0 * M_PI / phi / 3.0, 1.0 / 3.0);

  /* extend the primary cell with 1 particle to that with "n" particles*/
  j = 0;
  for (ix = 0; ix < nx; ix ++)
    {
      for (iy = 0; iy < ny; iy ++)
	{
	  for (iz = 0; iz < nz; iz ++)
	    {
	      pos [j * 3 + 0] = (double) ix * (*lx);
	      pos [j * 3 + 1] = (double) iy * (*ly);
	      pos [j * 3 + 2] = (double) iz * (*lz);
	      j ++;
	    }
	}
    }
  (*lx) *= (double) nx;
  (*ly) *= (double) ny;
  (*lz) *= (double) nz;
}

/* initialize configuration and the primary cell of BCC
 * INPUT
 *  phi        : volume fraction of particles
 *  nx, ny, nz : # single cell in each direction
 *             : so that total # particles is 2*nx*ny*nz
 * OUTPUT
 *  pos [(nx * ny * nz) 2 * * 3] : positions of particles
 *  lx, ly, lz : geometry of the primary cell
 */
void
init_config_BCC (double phi, int nx, int ny, int nz,
		 double *pos, double *lx, double *ly, double *lz)
{
  int j;
  int ix, iy, iz;


  (*lx) = (*ly) = (*lz) = pow (8.0 * M_PI / phi / 3.0, 1.0 / 3.0);

  /* extend the primary cell with 1 particle to that with "n" particles*/
  j = 0;
  for (ix = 0; ix < nx; ix ++)
    {
      for (iy = 0; iy < ny; iy ++)
	{
	  for (iz = 0; iz < nz; iz ++)
	    {
	      pos [j * 3 + 0] = (double) ix * (*lx);
	      pos [j * 3 + 1] = (double) iy * (*ly);
	      pos [j * 3 + 2] = (double) iz * (*lz);
	      j ++;
	      pos [j * 3 + 0] = ((double) ix + 0.5) * (*lx);
	      pos [j * 3 + 1] = ((double) iy + 0.5) * (*ly);
	      pos [j * 3 + 2] = ((double) iz + 0.5) * (*lz);
	      j ++;
	    }
	}
    }
  (*lx) *= (double) nx;
  (*ly) *= (double) ny;
  (*lz) *= (double) nz;
}

/* initialize configuration and the primary cell of FCC
 * INPUT
 *  phi        : volume fraction of particles
 *  nx, ny, nz : # single cell in each direction
 *             : so that total # particles is 4*nx*ny*nz
 * OUTPUT
 *  pos [(nx * ny * nz) 4 * * 3] : positions of particles
 *  lx, ly, lz : geometry of the primary cell
 */
void
init_config_FCC (double phi, int nx, int ny, int nz,
		 double *pos, double *lx, double *ly, double *lz)
{
  int j;
  int ix, iy, iz;


  (*lx) = (*ly) = (*lz) = pow (16.0 * M_PI / phi / 3.0, 1.0 / 3.0);

  /* extend the primary cell with 1 particle to that with "n" particles*/
  j = 0;
  for (ix = 0; ix < nx; ix ++)
    {
      for (iy = 0; iy < ny; iy ++)
	{
	  for (iz = 0; iz < nz; iz ++)
	    {
	      pos [j * 3 + 0] = (double) ix * (*lx);
	      pos [j * 3 + 1] = (double) iy * (*ly);
	      pos [j * 3 + 2] = (double) iz * (*lz);
	      j ++;
	      pos [j * 3 + 0] = ((double) ix + 0.0) * (*lx);
	      pos [j * 3 + 1] = ((double) iy + 0.5) * (*ly);
	      pos [j * 3 + 2] = ((double) iz + 0.5) * (*lz);
	      j ++;
	      pos [j * 3 + 0] = ((double) ix + 0.5) * (*lx);
	      pos [j * 3 + 1] = ((double) iy + 0.0) * (*ly);
	      pos [j * 3 + 2] = ((double) iz + 0.5) * (*lz);
	      j ++;
	      pos [j * 3 + 0] = ((double) ix + 0.5) * (*lx);
	      pos [j * 3 + 1] = ((double) iy + 0.5) * (*ly);
	      pos [j * 3 + 2] = ((double) iz + 0.0) * (*lz);
	      j ++;
	    }
	}
    }
  (*lx) *= (double) nx;
  (*ly) *= (double) ny;
  (*lz) *= (double) nz;
}
