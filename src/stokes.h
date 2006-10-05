/* header file for stokes.c --
 * structure for system parameters of stokes library.
 * Copyright (C) 2001-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes.h,v 1.4 2006/10/05 00:25:54 kichiki Exp $
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

struct stokes {
  int np; /* number of all particles */
  int nm; /* number of mobile particles */
  double * pos; /* position of particles */

  int version; /* 0 = F, 1 = FT, 2 = FTS */

  /* for ewald codes */
  double rmax;
  int rmaxx, rmaxy, rmaxz; /* # of cell in real space */
  double kmax;
  int kmaxx, kmaxy, kmaxz; /* # of cell in reciprocal space */

  double zeta,zeta2,zaspi,za2;
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

void
stokes_set_np (struct stokes * sys,
	       int np, int nm);

void
stokes_set_ll (struct stokes * sys,
	       double lx, double ly, double lz);

void
stokes_set_zeta (struct stokes * sys,
		 double zeta, double cutlim);

double
zeta_by_tratio (struct stokes * sys,
		double tratio);

/** table for lattice summation **/
/* free ewald-table related memory from sys
 */
void
stokes_ewald_table_free (struct stokes * sys);

/* make ewald-summation table
 * INPUT
 *  (struct stokes *) sys :
 *  cutlim   :
 * OUTPUT
 *  (struct stokes *) sys : 
 */
void
stokes_ewald_table_make (struct stokes * sys, double cutlim);


#endif /* !_STOKES_H_ */
