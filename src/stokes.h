/* header file for stokes.c --
 * structure for system parameters of stokes library.
 * Copyright (C) 2001-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes.h,v 1.13 2007/03/29 02:26:32 kichiki Exp $
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
#ifndef	_STOKES_H_
#define	_STOKES_H_

#include <stdio.h> // FILE


struct stokes {
  int np;      /* number of all particles  */
  int nm;      /* number of mobile particles  */
  double *pos; /* position of particles  */
  double *a;   /* radius of particles
		* Note : NULL (default) is for monodisperse system  */

  int version; /* 0 = F, 1 = FT, 2 = FTS  */

  /* imposed flow */
  double Ui[3];
  double Oi[3];
  double Ei[5];

  int periodic; // 0 = non periodic, 1 = periodic

  /* for ewald codes */
  double rmax2;
  int rmaxx, rmaxy, rmaxz; /* # of cell in real space */
  double kmax;
  int kmaxx, kmaxy, kmaxz; /* # of cell in reciprocal space */

  double xi, xi2, xiaspi, xia2;
  double pivol;

  double lx, ly, lz;
  double llx[27], lly[27], llz[27]; /* for regist and lub */

  // self part
  double self_a;
  double self_c;
  double self_m;

  // table for lattice summation
  int flag_table; // 0 = inactive, 1 = active
  // real space
  int nr; // number of lattice points
  double * rlx;
  double * rly;
  double * rlz;
  // reciprocal space
  int nk; // number of lattice points
  double * ex;
  double * ey;
  double * ez;
  double * k;
  double * k1;
  double * k2;
  double * k3;
  double * ya;
  double * yb;
  double * yc;
  double * yg;
  double * yh;
  double * ym;

  /* for lubrication */
  double lubcut;
  double lubmax2; // square of the max distance (0 means no limit)

  /* for zeta program */
  double cpu1, cpu2, cpu3;

  /* for iterative solvers */
  struct iter * it;
};


/* all elements are zero-cleared
 */
struct stokes *
stokes_init (void);

void
stokes_free (struct stokes * sys);

/* set np and nm and allocate the memory for pos[np*3]
 */
void
stokes_set_np (struct stokes * sys,
	       int np, int nm);

void
stokes_set_Ui (struct stokes * sys,
	       double uix, double uiy, double uiz);
void
stokes_set_Oi (struct stokes * sys,
	       double oix, double oiy, double oiz);
void
stokes_set_Ei (struct stokes * sys,
	       double eixx, double eixy, double eixz,
	       double eiyz, double eiyy);

void
stokes_set_l (struct stokes * sys,
	      double lx, double ly, double lz);

void
stokes_set_xi (struct stokes * sys,
	       double xi, double ewald_eps);

double
xi_by_tratio (struct stokes * sys,
	      double tratio);


/* set iter param
 * INPUT
 *   solver : string indicating the solver
 *            sta, sta2, gpb, otmk, or gmres (default)
 *   eps and log10_eps
 *   max (and restart)
 *   debug = 0 : no debug info
 *         = 1 : iteration numbs and residue
 *   out   : FILE * to output debug info.
 */
void
stokes_set_iter (struct stokes * sys,
		 const char * solver,
		 int max,
		 int restart,
		 double eps,
		 int debug,
		 FILE * out);

/* set pos for all particles safely by another array
 * INPUT
 *  pos[np*3] :
 * OUTPUT
 *  sys->pos[i] for (i = 0; i < np*3)
 */
void
stokes_set_pos (struct stokes *sys,
		const double *pos);

/* set pos for mobile particles safely by another array
 * INPUT
 *  pos[nm*3] :
 * OUTPUT
 *  sys->pos[i] for (i = 0; i < nm*3)
 */
void
stokes_set_pos_mobile (struct stokes *sys,
		       const double *pos);

/* set pos for fixed particles safely by another array
 * INPUT
 *  pos[nf*3] : only fixed particles are set, where nf = np - nm
 * OUTPUT
 *  sys->pos[i] for (i = nm*3; i < np*3)
 */
void
stokes_set_pos_fixed (struct stokes *sys,
		      const double *pos);

/* set radius (sys->a[]).
 * Note that the default setting (sys->a == NULL) is for monodisperse system
 * where a=1 for all particles
 * INPUT
 *  a[np] :
 */
void
stokes_set_radius (struct stokes *sys,
		   const double *a);


#endif /* !_STOKES_H_ */
